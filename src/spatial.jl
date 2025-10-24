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
    event::Vector{Vector{T2}}
end

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

function SpatialRepresentation(spikes::Spiketrain, rp::RippleData, udata::UnityData;min_speed=0.0,trial_start=1, gidx::Union{Vector{Vector{Bool}}, Nothing}=nothing)
    nt = numtrials(udata)
    position = Vector{Vector{Point2f}}(undef, nt)
    timestamp = Vector{Vector{Float64}}(undef, nt)
    events = Vector{Vector{Float64}}(undef, nt)
    sp = spikes.timestamps/1000.0 #convert to seconds
    for i in 1:nt
        tp,posx,posy,_ = get_trial(udata,i;trial_start=trial_start)
        tp .-= tp[1]
        if gidx !== nothing
            use_bin = gidx[i]
        else
            use_bin = fill(true,length(tp))
        end
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
        position[i] = Point2f[]
        for j in 1:nspikes
            k = searchsortedlast(tp,sp_trial[j])
            if 0 < k < length(posx)
                if use_bin[k]
                    push!(position[i], Point2f(posx[k],posy[k]))
                    push!(events[i], sp_trial[j])
                    push!(timestamp[i], tp[k])
                end
            end
        end
    end
    SpatialRepresentation(position,timestamp, events)
end

function SpatialRepresentation(;kwargs...)
    sptrain = Spiketrain()
    rdata = cd(DPHT.process_level(level(RippleData))) do
        RippleData()
    end
    udata = cd(DPHT.process_level(level(UnityData))) do
        UnityData()
    end
    SpatialRepresentation(sptrain, rdata, udata;kwargs...)
end

each(X::AbstractVector{<:Any}) = X
each(X::AbstractMatrix{<:Any}) = eachcol(X)

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
    if typeof(pq) == Point{2,T1}
        d = 2
    elseif typeof(pq) == Point{3,T1}
        d = 3
    else
        d = length(pq)
    end
    position = zeros(T1, d, nq)
    time_window = zeros(T1, nq)
    X = zeros(T1, nr, nq)
    for (i,(q,v)) in enumerate(Q)
        position[:,i] = q
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
    ax = S.Axis(plots=[S.Scatter(points)])
    S.GridLayout(ax)
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

abstract type AbstractSpatialMap end

"""
Contains information about the total time spent in each spatial bin.
"""
struct SpatialOccupancy{T<:Real} <: AbstractSpatialMap
    xbins::AbstractVector{T}
    ybins::AbstractVector{T}
    weight::Matrix{T}
end

function SpatialOccupancy(udata::UnityData, xbins::AbstractVector{T}, ybins::AbstractVector{T};trial_start=1,min_speed=0.0, min_duration=0.0, min_num_observations=0,gidx::Union{Vector{Vector{Bool}},Nothing}=nothing) where T <: Real
    nt = numtrials(udata)
    weight = zeros(T, length(xbins)-1, length(ybins)-1)
    w = fill!(similar(weight), zero(T))
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
                    w[xidx,yidx] += Δt
                end
            end
        end
        # check number of observations per bin
        weight .+= w
    end
    SpatialOccupancy(xbins, ybins, weight)
end

function SpatialOccupancy(xbins,ybins=xbins;kwargs...)
    udata = UnityData()
    SpatialOccupancy(udata, xbins, ybins;kwargs...)
end

struct SpatialMap{T<:Real} <: AbstractSpatialMap
    xbins::AbstractVector{T}
    ybins::AbstractVector{T}
    weight::Matrix{T}
    occupancy::Matrix{T}
end

struct SmoothedSpatialMap{T<:Real} <: AbstractSpatialMap
    xbins::AbstractVector{T}
    ybins::AbstractVector{T}
    weight::Matrix{T}
    occupancy::Matrix{T}
    unvisited::Vector{CartesianIndex{2}}
    α::T
end

DPHT.level(::Type{<:AbstractSpatialMap}) = "cell"

function SpatialMap(spr::SpatialRepresentation{T,T2}, spoc::SpatialOccupancy{T};kwargs...) where T <: Real where T2 <: Real
    xbins = spoc.xbins
    ybins = spoc.ybins
    spatial_count = zeros(T, length(xbins)-1, length(ybins)-1)
    nt = numtrials(spr)
    for i in 1:nt
        position = spr.position[i]
        xpos = [pos[1] for pos in position]
        ypos = [pos[2] for pos in position]
        h = fit(Histogram, (xpos,ypos), (xbins, ybins))
        spatial_count .+= h.weights
    end
    SpatialMap(xbins,ybins, spatial_count, spoc.weight)
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
    spr = SpatialRepresentation(sp, rp, udata;kwargs...);
    SpatialMap(spr, xbins, ybins, spoc)
end

function filter_occupancy(spm::SpatialMap{T}) where T <: Real
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

    x = spm.weight./spm.occupancy
    idx = isfinite.(x)
    p = spm.occupancy[idx]./sum(spm.occupancy[idx])
    r = sum(p.*x[idx])
    xr = x[idx]./r
    ll = log2.(xr)
    lidx = isfinite.(ll)
    sic = sum((p.*xr.*ll)[lidx])
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

function merge_responses(X::Matrix{T}, assignment::AbstractVector{Int64}, weight::Union{AbstractVector{T2},Nothing}=nothing) where T <: Real where T2 <: Real
    nc = maximum(assignment)
    if weight === nothing
        weight == fill(one(T), nc)
    end
    X2 = zeros(T, size(X,1), nc)
    for (i,k) in enumerate(assignment)
       X2[:,k] .+= X[:,i]
    end
    X2 ./= reshape(weight,1,size(X2,2))
    X2
end

function regress_space(X::Matrix{T}, Y::Matrix{T};n_spatial_clusters=256,kwargs...) where T <: Real
    km_results = kmeans(Y, n_spatial_clusters)
    # sum up responses in each of the spatial bins returned by the kmean algorithm
    X2 = zeros(eltype(X), size(X,1), n_spatial_clusters)
    for (i,k) in enumerate(km_results.assignments)
       X2[:,k] .+= X[:,i]
    end
    
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