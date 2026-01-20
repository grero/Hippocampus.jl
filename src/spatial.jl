using Makie
import Makie.SpecApi as S
using Clustering
using MultivariateStats
using LinearAlgebra
using LinearRegressionUtils

abstract type AbstractMap  end
abstract type AbstractSpatialMap <: AbstractMap end

"""
Contains information about the total time spent in each spatial bin.
"""
struct SpatialOccupancy{T<:Real} <: AbstractSpatialMap
    xbins::AbstractVector{T}
    ybins::AbstractVector{T}
    weight::Array{T,3}
    bincount::Array{T,3}
end

function SpatialOccupancy(udata::UnityData, xbins::AbstractVector{T}, ybins::AbstractVector{T};trial_start=1,min_speed=0.0, min_duration=0.0, min_num_observations=0,gidx::Union{Vector{Vector{Bool}},Nothing}=nothing) where T <: Real
    nt = numtrials(udata)
    weight = zeros(T, length(xbins)-1, length(ybins)-1, nt)
    w = fill!(similar(weight), zero(T))
    bincount = zeros(T, size(weight)...)
    for i in 1:nt
        # keep track of observations per trial
        fill!(w, zero(T))
        tu, posx, posy, _ = get_trial(udata, i;trial_start=trial_start)
        if gidx !== nothing
            _gidx = gidx[i]
        else
            _gidx = fill(true, length(tu))
        end
        for j in 2:length(tu)
            if !_gidx[j]
                continue
            end
            Δt = tu[j]-tu[j-1]
            ds = sqrt((posx[j] - posx[j-1])^2 + (posy[j]-posy[j-1])^2)
            ds /= Δt
            if ds > min_speed
                xidx = searchsortedlast(xbins, posx[j-1])
                yidx = searchsortedlast(ybins, posy[j-1])
                if 0 < xidx <= size(weight,1) && 0 < yidx <= size(weight,2)
                    w[xidx,yidx,i] += Δt
                    bincount[xidx,yidx,i] += 1.0
                end
            end
        end
        # check number of observations per bin
        weight .+= w
    end
    SpatialOccupancy(xbins, ybins, weight, bincount)
end

function SpatialOccupancy(xbins::AbstractVector{T},ybins=xbins;kwargs...) where T <: Real
    udata = UnityData()
    SpatialOccupancy(udata, xbins, ybins;kwargs...)
end

struct SpatialOccupancyNew{T<:Real} <: AbstractSpatialMap
    mm::SimpleMesh
    weight::Matrix{T}
    bincount::Matrix{T}
end

function SpatialOccupancyNew(udata::UnityData,mm::SimpleMesh;trial_start=1,min_speed=0.0)
    T = eltype(udata.position)
    nt = numtrials(udata)
    ne = nelements(mm)
    weight = zeros(T, ne, nt)
    w = zeros(T,ne) 
    bincount = zeros(T, size(weight)...)
    kn = KNearestSearch(mm, 1)
    for i in 1:nt
        # keep track of observations per trial
        fill!(w, zero(T))
        tu, posx, posy, _ = get_trial(udata, i;trial_start=trial_start)
        for j in 2:length(tu)
            Δt = tu[j]-tu[j-1]
            ds = sqrt((posx[j] - posx[j-1])^2 + (posy[j]-posy[j-1])^2)
            ds /= Δt
            if ds > min_speed
                _idx,_dd = searchdists(Meshes.Point(posx[j-1], posy[j-1],0.0), kn)
                fidx = first(_idx)
                w[fidx] += Δt
                bincount[fidx,i] += 1.0
            end
        end
        # check number of observations per bin
        weight[:,i] .= w
    end
    SpatialOccupancyNew{T}(mm, weight, bincount)
end

function Makie.convert_arguments(::Type{<:AbstractPlot}, spm::SpatialOccupancyNew)
end

function SpatialOccupancyNew(;kwargs...)
    m_floor = Hippocampus.floor_topology3();
    udata = UnityData()
    SpatialOccupancyNew(udata,m_floor;kwargs...)
end

function explore(spm::SpatialOccupancyNew{<:Real};min_duration=0.05, min_num_obs=5)
    goodbinidx = findall(dropdims(sum(spm.weight .> min_duration,dims=2),dims=2) .>= min_num_obs)
    tcolor = zeros(size(spm.weight,1))
    tcolor[goodbinidx] = dropdims(sum(spm.weight[goodbinidx,:],dims=2),dims=2)
    
    mms = Shadow("xy")(spm.mm)
    fig,ax = viz(mms;color=tcolor, showsegments=true)
    Colorbar(fig[1,2], colorrange=extrema(tcolor), label="Duration [s]")
    display(fig)
    fig,ax
end

abstract type AbstractRepresentation{T1<:Real,T2<:Real} end
"""
A spatial representation of events
"""
struct SpatialRepresentation{T1<:Real,T2<:Real} <: AbstractRepresentation{T1,T2}
    position::Vector{Vector{Point{2, T1}}}
    timestamp::Vector{Vector{T2}}
    time_window::Vector{Vector{T2}}
    event::Vector{Vector{T2}}
end

DPHT.filename(::Type{SpatialRepresentation}) = "spatial_representation.jld2"
DPHT.filename(::Type{SpatialRepresentation{T1, T2}}) where T1 <: Real where T2 <: Real = "spatial_representation.jld2"
DPHT.filename(X::SpatialRepresentation{<:Real, <:Real}) = "spatial_representation.jld2"
DPHT.level(::Type{SpatialRepresentation}) = "cell"

"""
Return a matrix of all positions
"""
function get_positions(spr::SpatialRepresentation{T1,T2}) where T1 <: Real where T2 <: Real
    nspikes = sum(length.(spr.event))
    Y = zeros(T1, 2, nspikes)
    offset = 0
    for pps in spr.position
        for (j,pp)  in enumerate(pps)
            Y[:, offset+j] .= pp
        end
        offset += length(pps)
    end
    Y
end

get_rep(sr::SpatialRepresentation) = sr.position

function SpatialRepresentation(spikes::Spiketrain, rp::RippleData, udata::UnityData;kwargs...)
    sp = spikes.timestamps/1000.0 #convert to seconds
    SpatialRepresentation(sp, rp, udata;kwargs...)
end

function SpatialRepresentation(sp::AbstractVector{T}, rp::RippleData, udata::UnityData;min_speed=0.0,trial_start=1, gidx::Union{Vector{Vector{Bool}}, Nothing}=nothing) where T <: Real
    nt = numtrials(udata)
    position = Vector{Vector{Point2f}}(undef, nt)
    timestamp = Vector{Vector{Float64}}(undef, nt)
    time_window = Vector{Vector{Float64}}(undef, nt)
    events = Vector{Vector{Float64}}(undef, nt)
    for i in 1:nt
        tp,posx,posy,_ = get_trial(udata,i;trial_start=trial_start)
        tp .-= tp[1]
        use_bin = fill(true,length(tp))
        if min_speed > 0
            # mark bins as invalid if the speed is too low
            vv = sqrt.(diff(posx).^2 + diff(posy).^2)./diff(tp)
            fidx = findall(vv .< min_speed)
            use_bin[fidx] .= false
        end
        timestamps = rp.timestamps[i,:]
        idx0 = searchsortedfirst(sp, timestamps[trial_start])
        idx1 = searchsortedlast(sp, timestamps[3])
        # align to trial start
        sp_trial = sp[idx0:idx1] .- timestamps[trial_start]
        nspikes = idx1-idx0+1
        events[i] = Float64[] 
        timestamp[i] = Float64[]
        time_window[i] = Float64[]
        position[i] = Point2f[]
        for j in 1:nspikes
            k = searchsortedlast(tp,sp_trial[j])
            if 0 < k < length(posx)
                if use_bin[k]
                    push!(position[i], Point2f(posx[k],posy[k]))
                    push!(events[i], sp_trial[j])
                    push!(timestamp[i], tp[k])
                    push!(time_window[i], tp[k+1]-tp[k])
                end
            end
        end
    end
    SpatialRepresentation(position,timestamp, time_window, events)
end

function SpatialRepresentation(;redo=false, do_save=true, kwargs...)
    fname = DPHT.filename(SpatialRepresentation)
    if !redo && isfile(fname)
        spr = load_jld2(SpatialRepresentation)
    else
        rdata = cd(DPHT.process_level(level(RippleData))) do
            RippleData()
        end
        udata = cd(DPHT.process_level(level(UnityData))) do
            UnityData()
        end
        spr = SpatialRepresentation(rdata, udata;kwargs...)
        if do_save
            save_jld2(spr)
        end
    end
    spr
end

function SpatialRepresentation(rp::RippleData, udata::UnityData;kwargs...)
    sptrain = Spiketrain()
    SpatialRepresentation(sptrain, rp, udata;kwargs...)
end

each(X::AbstractVector{<:Any}) = X
each(X::AbstractMatrix{<:Any}) = eachcol(X)
each(X::Tuple{AbstractVector{<:Any}, AbstractVector{<:Any}}) = zip(X[1], X[2])

function get_population_representation(spr::Vector{T3}) where T3 <: AbstractRepresentation{T1,T2} where T1 <: Real where T2 <: Real
    # save in a dictionary for now
    nr = length(spr)
    Q = Dict()
    for (i,_spr) in enumerate(spr)
        _position = get_rep(_spr)
        _time = _spr.time_window
        for (kk,(pos,tt)) in enumerate(zip(_position,_time))
            for (j,p) in enumerate(each(pos))
                if !(p in keys(Q))
                    Q[p] = Dict("count"=> fill(0, nr), "time" => 0.0)
                end
                Q[p]["count"][i] += 1
                Q[p]["time"] +=tt[j] 
                 # TODO: We should also keep track of time windows here
            end
        end
    end
    nq = length(Q) 
    pq = first(keys(Q))
    tq = false
    if typeof(pq) == Point{2,T1}
        d = 2
    elseif typeof(pq) == Point{3,T1}
        d = 3
    elseif typeof(pq) == Tuple{Point{3,T1}, Point{2,T1}}
        d = 5
        tq = true
    else
        d = length(pq)
    end
    position = zeros(T1, d, nq)
    time_window = zeros(T1, nq)
    X = zeros(T1, nr, nq)
    for (i,(q,v)) in enumerate(Q)
        if tq 
            position[1:3,i] = q[1]
            position[4:5,i] = q[2]
        else
            position[:,i] = q
        end
        X[:,i] .= v["count"]
        time_window[i] = v["time"]
    end
    X,position,time_window
end

numtrials(spr::SpatialRepresentation) = length(spr.position)

function Makie.convert_arguments(::Type{<:AbstractPlot}, spr::SpatialRepresentation)
    nt = numtrials(spr)
    points = Point2f[]
    for pos in spr.position
        append!(points, Point2f.(pos))
    end
    S.Scatter(points)
end

function visualize!(lscene, spr::SpatialRepresentation;trial::Observable{Trial}=Observable(Trial(1)), current_time::Observable{Float64}=Observable(0.0))
    nt = numtrials(spr)
    trial_events = lift(trial) do _trial
        if 0 < _trial.i <= nt
            return spr.position[_trial.i]
        else
            return [Point2f(NaN)]
        end
    end
    scatter!(lscene, trial_events)
end

"""
Contains information about the total time spent in each spatial bin.
"""

struct SpatialMap{T<:Real} <: AbstractSpatialMap
    xbins::AbstractVector{T}
    ybins::AbstractVector{T}
    weight::Matrix{T}
    occupancy::Matrix{T}
end

struct SpatialMapNew{T<:Real} <: AbstractSpatialMap
    mm::SimpleMesh
    weight::Vector{T}
    occupancy::Vector{T}
end

DPHT.filename(::Type{SpatialMapNew}) = "spatial_map.jld2"
DPHT.filename(X::SpatialMapNew{T}) where T <: Real = "spatial_map.jld2"
DPHT.level(::Type{SpatialMapNew}) = "cell"
DPHT.level(X::SpatialMapNew{T}) where T <: Real = "cell"

struct SmoothedSpatialMap{T<:Real} <: AbstractSpatialMap
    xbins::AbstractVector{T}
    ybins::AbstractVector{T}
    weight::Matrix{T}
    occupancy::Matrix{T}
    unvisited::Vector{CartesianIndex{2}}
    α::T
end

struct SmoothedMap{T<:Real} <: AbstractMap
    mm::SimpleMesh
    weight::Vector{T}
    occupancy::Vector{T}
    unvisited::Vector{Int64}
    smooth_params::NamedTuple
end

function SmoothedMap(spm::AbstractMap;method=:gaussian, σ=5, m=4, edge_correct=false, α=1000.0^2,kwargs...)
    if method == :gaussian
        Zg, Xg, Yg = gaussian_smoothing(spm.weight, spm.occupancy, spm.mm, σ;m=m,kwargs...)
        smooth_params = (method=method, σ=σ, m=m, edge_correct=edge_correct)
    elseif method == :adaptive
        Zg, Xg, Yg = adaptive_smoothing(spm.weight, spm.occupancy, spm.mm, α;kwargs...)
        smooth_params = (method=method,α=α) 
    else
        error("Unkonwn smoothing method $method")
    end
    unvisited = findall(spm.occupancy .== 0)
    SmoothedMap(spm.mm, Xg, Yg, unvisited, smooth_params)
end

DPHT.level(::Type{<:AbstractSpatialMap}) = "cell"

function SpatialMap(spr::SpatialRepresentation{<:Real,<:Real}, spoc::SpatialOccupancy{T2};kwargs...) where T2 <: Real
    xbins = spoc.xbins
    ybins = spoc.ybins
    spatial_count = zeros(T2, length(xbins)-1, length(ybins)-1)
    nt = numtrials(spr)
    goodbinidx = findall(dropdims(sum(spoc.weight .> 0.05,dims=3),dims=3) .>= 5)
    for i in 1:nt
        position = spr.position[i]
        xpos = [pos[1] for pos in position]
        ypos = [pos[2] for pos in position]
        h = StatsBase.fit(Histogram, (xpos,ypos), (xbins, ybins))
        # remove non-valid bin counts
        spatial_count[goodbinidx] .+= h.weights[goodbinidx]
    end
    spoc_weight = zeros(T2, size(spoc.weight,1), size(spoc.weight,2))
    spoc_weight[goodbinidx] .= dropdims(sum(spoc.weight[goodbinidx, :],dims=2),dims=2)
    SpatialMap(xbins,ybins, spatial_count, spoc_weight)
end

function process_kwargs(::Type{SpatialMapNew};min_duration=0.05, min_n_obs=5, kwargs...)
    h = UInt32(0)
    if min_duration != 0.05
        h = crc32c(string(min_duration=>min_duration),h)
    end
    if min_n_obs != 5
        h = crc32c(string(min_n_obs=>min_n_obs))
    end
    h
end

function SpatialMapNew(spr::SpatialRepresentation{T,<:Real}, spoc::SpatialOccupancyNew{T2};min_duration=0.05, min_n_obs=5, kwargs...) where T2 <: Real where T <: Real
    h = process_kwargs(SpatialMapNew, min_duration=min_duration, min_n_obs=min_n_obs)
    mm = spoc.mm
    spatial_count = zeros(T2, size(spoc.weight,1))
    nt = numtrials(spr)
    # TODO: Do we need to worry about view bins here?
    goodbinidx = findall(dropdims(sum(spoc.weight .> min_duration,dims=2),dims=2) .>= min_n_obs)
    f = in(goodbinidx)
    for i in 1:nt
        position = spr.position[i]
        if isempty(position)
            continue
        end
        xpos = [pos[1] for pos in position]
        ypos = [pos[2] for pos in position]
        points = [(xp,yp, zero(T)) for (xp,yp) in zip(xpos,ypos)]
        bidx = mapto(mm, points)
        for k in bidx
            for k1 in k
                if f(k1) 
                    spatial_count[k1] += 1.0
                end
            end
        end
    end
    spoc_weight = zeros(T2, size(spoc.weight,1))
    spoc_weight[goodbinidx] .= dropdims(sum(spoc.weight[goodbinidx, :],dims=2),dims=2)
    SpatialMapNew(mm, spatial_count, spoc_weight), h
end

function SpatialMapNew(;redo=false, do_save=true, kwargs...)
    h = process_kwargs(SpatialMapNew;kwargs...)
    fname = DPHT.filename(SpatialMapNew)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
        spm = load_jld2(SpatialMap,fname)
    end
    if !redo && isfile(fname)
        spm = load_jld2(SpatialMapNew, fname)
    else
        spoc = cd(DPHT.process_level("session")) do
            SpatialOccupancyNew()
        end
        spr = SpatialRepresentation(;kwargs...)
        spm,h = SpatialMapNew(spr, spoc;kwargs...)
        if do_save
            save_jld2(spm, fname)
        end
    end
    spm
end

function get_rate_map(spm::AbstractMap;invalidate_unvisited=true)
    Z = spm.weight./spm.occupancy
    unvisited = spm.occupancy.==0
    if !invalidate_unvisited
        Z[unvisited] .= zero(eltype(Z))
    end
    Z
end

function get_rate_map(spm::SmoothedMap;filter_unvisited=true)
    Z = spm.weight./spm.occupancy
    if filter_unvisited
        Z[spm.unvisited] .= NaN
    end
    Z
end


function explore(sm::SmoothedMap;filter_unoccupied=true,goodbinidx::Union{Nothing, Vector{Int64}}=nothing, kwargs...)
    mm = sm.mm
    tcolor = sm.weight./sm.occupancy 
    alpha = fill(1.0, length(tcolor))
    alpha[tcolor.==0] .= 0.0
    if goodbinidx !== nothing
        badbinidx = setdiff(1:length(tcolor), goodbinidx)
        tcolor[badbinidx] .= 0.0
        alpha[badbinidx] .= 0.0
    elseif filter_unoccupied
        tcolor[sm.unvisited] .= 0.0 
        alpha[sm.unvisited] .= 0.0
    end
    explore(mm;color=tcolor, alpha=alpha, showsegments=true, kwargs...)
end

function explore(spm::AbstractSpatialMap;separate_occupancy=false,filter_zeros=false)
    mms = Shadow("xy")(spm.mm)
    offset = 0 
    figheight = 300 
    if separate_occupancy
        figwidth = 3*figheight 
    else
        figwidth = figheight
    end
    fig = Figure(size=(figwidth,figheight))
    if separate_occupancy
        ax1 = Axis(fig[1,1])
        tcolor1 = spm.weight
        alpha1 = fill(1.0, size(tcolor1)...) 
        fidx1 = (!isfinite).(tcolor1)
        alpha1[fidx1] .= 0.0
        tcolor1[fidx1] .= 0.0
        if filter_zeros
            alpha1[tcolor1.==0.0] .= 0.0
        end

        viz!(ax1, mms;color=tcolor1, alpha=alpha1,showsegments=true)
        Colorbar(fig[1,2], colorrange=extrema(tcolor1[alpha1.>0]),label="Spike count")

        ax2 = Axis(fig[1,3])
        tcolor2 = spm.occupancy
        alpha2 = fill(1.0, size(tcolor2)...) 
        fidx2 = (!isfinite).(tcolor2)
        alpha2[fidx1] .= 0.0
        tcolor2[fidx1] .= 0.0
        if filter_zeros
            alpha2[tcolor2.==0.0] .= 0.0
        end
        viz!(ax2, mms;color=tcolor2, alpha=alpha2,showsegments=true)
        Colorbar(fig[1,4], colorrange=extrema(tcolor2[alpha2.>0]),label="Duration [s]")
        offset = 4
    end
    tcolor = spm.weight./spm.occupancy
    alpha = fill(1.0, size(tcolor)...) 
    fidx = (!isfinite).(tcolor)
    alpha[fidx] .= 0.0
    tcolor[fidx] .= 0.0
    if filter_zeros
        alpha[tcolor.==0.0] .= 0.0
    end
    ax = Axis(fig[1,offset+1])
    viz!(ax, mms;color=tcolor, alpha=alpha,showsegments=true)
    Colorbar(fig[1,offset+2], colorrange=extrema(tcolor[alpha.>0]),label="Firing rate [Hz]")

    display(fig)
    fig
end

function SpatialMap(xbins, ybins;redo=false, do_save=false,kwargs...)
    sp = Spiketrain()
    rp = cd(DPHT.process_level(RippleData;kwargs...)) do
        RippleData()
    end
    udata = cd(DPHT.process_level(UnityData;kwargs...)) do
        UnityData()
    end
    spoc = SpatialOccupancy(udata, xbins, ybins)
    spr = SpatialRepresentation(sp, rp, udata,spoc;kwargs...);
    SpatialMap(spr, spoc)
end

function filter_occupancy(spm::SpatialMap{T}) where T <: Real
    findall(spm.occupancy .== 0)
end

function filter_occupancy(spm::SpatialMapNew{T}) where T <: Real
    findall(spm.occupancy .== 0)
end

function filter_occupancy(spm::SmoothedSpatialMap{T}) where T <: Real
    spm.unvisited
end

function adaptive_smoothing(spm::SpatialMap{T}, α=T(10000.0)^2;filter_unoccupied=true) where T <: Real
    unoccupied = findall(spm.occupancy.==0)
    Z,X,Y = adaptive_smoothing(spm.weight, spm.occupancy, α)
    SmoothedSpatialMap(spm.xbins, spm.ybins, X, Y, unoccupied,α)
end

function compute_sic(spm::AbstractSpatialMap)
    compute_sic(spm.weight, spm.occupancy)
end

function compute_sic(weight, occupancy)
    x = weight./occupancy
    idx = isfinite.(x)
    p = occupancy[idx]./sum(occupancy[idx])
    r = sum(p.*x[idx])
    xr = x[idx]./r
    ll = log2.(xr)
    lidx = isfinite.(ll)
    sic = sum((p.*xr.*ll)[lidx])
end

"""
Compute Skagg's SIC for per bin firing rate `λ` and per-bin `occupancy`
"""
function compute_skaggs_sic(λ::Vector{T}, occupancy::Vector{T}) where T <: Real
    ps = occupancy./sum(occupancy)
    l1 = ps.*λ
    ll = sum(filter(isfinite, l1))
    l1./= ll
    l2 = log2.(λ./ll)
    sum(filter(isfinite, l1.*l2))
end

function compute_skaggs_sic(sm::SmoothedMap)
    weight = sm.weight
    λ = weight./sm.occupancy
    occupancy = similar(sm.occupancy)
    occupancy .= sm.occupancy
    occupancy[sm.unvisited] .= 0
    λ[sm.unvisited] .= NaN
    compute_skaggs_sic(λ, occupancy)
end

function compute_skaggs_sic(sm::SpatialMapNew)
    weight = sm.weight
    λ = weight./sm.occupancy
    compute_skaggs_sic(λ, sm.occupancy)
end

function compute_entropy(sp::Union{SpatialMap, SpatialOccupancy})
    pp = sp.weight./sum(sp.weight)
    -sum(filter(isfinite, pp.*log2.(pp)))
end

function Makie.convert_arguments(::Type{<:AbstractPlot}, spm::AbstractSpatialMap,args::NamedTuple=(;))
    default_arguments = Dict(:normalize=>true, :filter_gaze => true)
    for k in keys(args)
        v = args[k] 
        default_arguments[k] = v
    end
    X = copy(spm.weight)
    if default_arguments[:normalize]
        X = spm.weight./spm.occupancy
        label = "Firing rate [Hz]"
    else
        X = spm.weight
        label = "Spike count"
    end
    if default_arguments[:filter_gaze]
        oidx = filter_occupancy(spm)
        X[oidx] .= eltype(X)(NaN)
    end
    h = S.Heatmap(spm.xbins, spm.ybins, rotr90(X), colormap=:turbo)
    ax1 = S.Axis(plots=[h])
    ll = S.Colorbar(h,label=label)
    S.GridLayout([ax1 ll])
end

function Makie.convert_arguments(T::Type{<:AbstractPlot}, spm::SpatialMap,args::Vector{<:NamedTuple})
    a = [convert_arguments(T, spm, _arg) for _arg in args]
    S.GridLayout(a)
end

"""
Combine the spatial representions for each of the cells represented by `cellidirs` and use
those combined responses to regress to position in the maze
"""
function regress_space(celldirs::Vector{String};kwargs...)
    spr = map(celldirs) do celldir
       cd(celldir) do
            Hippocampus.SpatialRepresentation(;min_speed=0.3,trial_start=2)
        end
    end
    regress_space(spr;kwargs...)
end

function regress_space(spr::Vector{SpatialRepresentation{T1,T2}};n_spatial_clusters=256, kwargs...) where T1 <: Real where T2 <: Real
    X, Y = get_population_representation(spr)
    lq,pca, km_results, X2 = regression_space(X,Y;n_spatial_clusters=n_spatial_clusters, kwargs...)
    lq, pca, km_results, X2,X,Y 
end

function fit_spatial_glm(nspikes::Vector{T}, pos::Matrix{T}) where T <: Real
    m = glm(pos, nspikes)
end

function fit_spatial_glm(spr, km_results,pos)
    d,np = size(km_results.centers)
    Σ = zeros(d,d,np)
    for i in 1:np
        Σ[:,:,i] = cov(pos[:,km_results.assignments.==i];dims=2)
    end
    Pm_xy = [MultivariateNormal(km_results.centers[:,i],Σ[:,:,i]) for i in 1:np]
    P_xy = km_results.counts./sum(km_results.counts)
    nspikes = zeros(np)
    for trialpos in spr.position
       for _pos in trialpos
           idx = argmin(dropdims(sum(abs2, _pos .- km_results.centers,dims=1),dims=1))
           nspikes[idx] += 1
        end
    end
    m = glm(cat(permutedims(Float64.(km_results.centers)),fill(1.0, np),dims=2), round.(Int64, nspikes), Poisson())
    # construct Poisson distribution
    P_nspikes_xy = 1.0 .- cdf.(Poisson.(m.rr.mu), nspikes)
    # fit the number of spikes regardless of position
    # TODO: This is probably not quite correct
    P_nspikes = fit(Poisson, Int64.(nspikes))
    function P_xy_nspikes(x::Real,y::Real,nsp::Integer)
        λ = exp(m.pp.beta0[1]*x + m.pp.beta0[2]*y + m.pp.beta0[3])
        pxy = sum(pdf.(Pm_xy, [[x,y]]).*P_xy)
        pp = (1.0 .- cdf(Poisson(λ),nsp))*pxy
        pp./(1.0 .- cdf(P_nspikes,nsp))
    end
end

function merge_responses(X::Matrix{T}, km_results::Vector{Clustering.KmeansResult{Matrix{T},T,Int64}}) where T <: Real
    assignments = [km.assignments for km in km_results]
    centers = [km.centers for km in km_results]
    dm = size.(centers,1)
    d,n = size(X)
    nq = maximum.(assignments)
    nn = prod(nq)
    Xq = zeros(T, d, nq...)
    ny = zeros(Int64, 1, nq...)
    idx = [CartesianIndex(i,j) for (i,j) in zip(assignments...)]
    for (k,ii) in enumerate(idx)
        Xq[:,ii] .+= X[:,k]
        ny[1,ii] += 1
    end
    Y = zeros(T, sum(dm), nq...)
    iidx = CartesianIndices(tuple([1:_n for _n in nq]...))
    for ii in iidx
        offset = 0
        for j in 1:length(dm)
            Y[offset+1:offset+dm[j],ii] .= centers[j][:,ii.I[j]]
            offset += dm[j]
        end
    end
    # filter out combinations that did not happen
    fidx = findall(ny[1,:,:] .> 0)
    Xq[:,fidx]./ny[:,fidx], Y[:,fidx]
end

function merge_responses(X::Matrix{T},assignments::Vector{Vector{Int64}}) where T <: Real
    d,n = size(X)
    nq = maximum.(assignments)
    nn = prod(nq)
    Xq = zeros(T, d, nq...)
    ny = zeros(Int64, 1, nq...)
    idx = [CartesianIndex(i,j) for (i,j) in zip(assignments...)]
    for (k,ii) in enumerate(idx)
        Xq[:,ii] .+= X[:,k]
        ny[1,ii] += 1
    end
    Xq./ny
end

function merge_responses(X::Matrix{T}, assignment::AbstractVector{Int64}, weight::AbstractVector{T2}=ones(T, maximum(assignment))) where T <: Real where T2 <: Real
    nc = maximum(assignment)
    X2 = zeros(T, size(X,1), nc)
    for (i,k) in enumerate(assignment)
       X2[:,k] .+= X[:,i]
    end
    X2 ./= reshape(weight,1,size(X2,2))
    X2
end

function merge_by_time(Y::Matrix{<:Real}, twin::Vector{<:Real},seed=1;tmax=0.1)
    d = dropdims(sum(abs2, Y[:,seed:seed] .- Y,dims=1),dims=1)
    didx = sortperm(d)
    t = twin[seed]
    y = Y[:,seed]
    k = 1
    for i in 2:length(didx)
        t2 = t + twin[didx[i]]
        if t2 >= tmax
            break
        end
        t = t2
        y += Y[:,didx[i]]
        k += 1
    end
    y./k, t, didx[1:1+k]
end

function merge_by_time(Y1::Matrix{<:Real}, Y2::Matrix{<:Real}, twin::Vector{<:Real},seed=1;tmax=0.1)
    d1 = dropdims(sum(abs2, Y1[:,seed:seed] .- Y1,dims=1),dims=1)
    d2 = dropdims(sum(abs2, Y2[:,seed:seed] .- Y2,dims=1),dims=1)
    # normalize first
    d1 ./= maximum(d1)
    d2 ./= maximum(d2)
    d = d1+d2
    didx = sortperm(d)
    t = twin[seed]
    y1 = Y1[:,seed]
    y2 = Y2[:,seed]
    k = 1
    for i in 2:length(didx)
        t2 = t + twin[didx[i]]
        if t2 >= tmax
            break
        end
        t = t2
        y1 += Y1[:,didx[i]]
        y2 += Y2[:,didx[i]]
        k += 1
    end
    y1./k, y2./k, t, didx[1:1+k]
end

function regress_space(X::Matrix{T}, Y::Matrix{T};n_spatial_clusters=256,kwargs...) where T <: Real
    km_results = kmeans(Y, n_spatial_clusters)
    # sum up responses in each of the spatial bins returned by the kmean algorithm
    X2 = zeros(T, size(X,1), n_spatial_clusters)
    for (i,k) in enumerate(km_results.assignments)
       X2[:,k] .+= X[:,i]
    end
    # normalize to counts in for each point
    X2 ./= reshape(km_results.counts,1,size(X2,2))
    # perform PCA to decorrelate the inputs
    pca = fit(PCA, X2)
    Z = predict(pca, X2)
    lq = LinearRegressionUtils.llsq_stats(permutedims(Z), permutedims(km_results.centers))
    lq, pca, km_results, X2,X,Y 
end

function plot_regression_reults(lq, km_results, position, X)
    # compute histogram of positions
    xbins = range(-12.5f0, stop=12.5f0, length=40)
    ybins = xbins
    h = fit(Histogram, (eachrow(position)...,), (xbins, ybins))
    
    # project X onto the regression space
    q,r = qr(lq.β[1:end-1,:])
    w = q[:,1:2]
    Y = w'*X
    # also do the actual regresson
    Yp = lq.β[1:end-1,:]'*X .+ lq.β[end,:]
    with_theme(plot_theme) do
        fig = Figure(size=(600,700))
        ax1 = Axis(fig[1,1])
        Label(fig[1,1,TopLeft()], "A")
        hh = heatmap!(ax1, xbins, ybins, h.weights)
        Colorbar(fig[1,2], hh, label="Count") 
        scatter!(ax1, Point2f.(eachcol(km_results.centers)), color=:white,
                markersize=5px)

        
        lg1 = GridLayout(fig[2,1:2])
        ax1_1 = Axis(lg1[1,1])
        Label(lg1[1,1,TopLeft()], "B")
        scatter!(ax1_1, km_results.centers[1,:], Y[1,:],color=km_results.centers[1,:],
                            colormap=:solar)
        ax1_1.xlabel = "True x-pos"
        ax1_1.ylabel = "Estimated x-pos"
        ax1_2 = Axis(lg1[1,2])
        Label(lg1[1,2,TopLeft()], "C")
        scatter!(ax1_2, km_results.centers[2,:], Y[2,:], color=km_results.centers[2,:],
                        colormap=:matter)
        ax1_2.xlabel = "True y-pos"
        ax1_2.ylabel = "Estimated y-pos"

        ax1_3 = Axis(lg1[1,3])
        Label(lg1[1,3,TopLeft()], "D")
        barplot!(ax1_3, [1], [lq.r²])
        ax1_3.ylabel = "r²"
        ax1_3.xticksvisible = false
        ax1_3.xticklabelsvisible = false
        colsize!(lg1, 3, Relative(0.15))

        lg2 = GridLayout(fig[3,1:2])
        ax2 = Axis(lg2[1,1])
        Label(lg2[1,1,TopLeft()],"E")
        sc = scatter!(ax2, Point2f.(eachcol(Y)), color=km_results.centers[1,:],colormap=:solar)
        #Colorbar(lg2[2,1], sc, label="x-pos", vertical=false, flipaxis=false)
        ax3 = Axis(lg2[1,2])
        Label(lg2[1,2,TopLeft()], "F")
        scy = scatter!(ax3, Point2f.(eachcol(Y)), color=km_results.centers[2,:], colormap=:matter)
        #Colorbar(lg2[2,2], scy, label="y-pos", vertical=false, flipaxis=false)
        #rowsize!(lg2, 2, Relative(0.25))
        fig
    end
end