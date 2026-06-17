using LinearAlgebra

struct DirectionFiltered
    anglebins::AbstractVector{<:Real}
    weight::Vector{Dict{CartesianIndex{4}, Float64}}
    occupancy::Vector{Dict{CartesianIndex{4}, Float64}}
    index::Vector{Vector{CartesianIndex{5}}}
end

struct CardinalPlaceFieldDirectionality{T<:Real}
    directions::Vector{Symbol}
    weight::Array{T,3}
    occupancy::Array{T,3}
    ee::Vector{T}
    ees::Matrix{T}
    args::Dict{Symbol,Any}
end

function issignificant(card::CardinalPlaceFieldDirectionality{T};pv_threshold=0.05, tail=:right) where T <: Real
    res = fill(false, length(card.ee))
    for i in 1:length(res)
        vidx = isfinite.(card.ees[i,:])
        l,u = percentile(card.ees[i,vidx], [100*pv_threshold, 100*(1-pv_threshold)])
        if tail == :both
            res[i] = card.ee[i] < l || card.ee[i] > u
        elseif tail == :right
            res[i] = card.ee[i] > u
        elseif tail == :left
            res[i] = card.ee[i] < l
        else
            error("`tail` should `:left`, `:right`, or `:both`, noth $(tail)")
        end
    end
    res
end

function CardinalPlaceFieldDirectionality(gidx::DirectionFiltered;nshuffles=1000, kwargs...)
     m_floor = Shadow("xy")(Hippocampus.floor_topology3(;nrefinements=3));
    Ls = get_normalize_laplacian(m_floor)
    X,Y = get_joint_cardinal_prob(gidx)
    Xs = laplace_smoothing(X, Ls, 0.1;niter=50); 
    Ys = laplace_smoothing(Y, Ls, 0.1;niter=50); 

    ees = zeros(size(X,3),nshuffles)
    ee = zeros(size(X,3))
    for j in 1:size(X,3)
        ee[j]= get_conditional_information(Xs[:,:,j],Ys[:,:,j])
    end
    for i in 1:nshuffles
        X,Y = get_joint_cardinal_prob(gidx;do_shuffle=true)
        Xs = laplace_smoothing(X, Ls, 0.1;niter=50); 
        Ys = laplace_smoothing(Y, Ls, 0.1;niter=50); 
        for j in 1:size(X,3)
            ees[j,i] = get_conditional_information(Xs[:,:,j],Ys[:,:,j])
        end
    end
    CardinalPlaceFieldDirectionality{Float64}([:North, :South, :East, :West], X,Y, ee, ees, Dict{Symbol,Any}())
end

function CardinalPlaceFieldDirectionality(;redo=fname->false, do_save=true, nshuffles=1000, kwargs...)
    fname = "cardinal_place_field_directionality.jld2"
    h = process_kwargs(CardinalPlaceFieldDirectionality;nshuffles=nshuffles, kwargs...)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    if isfile(fname) && !redo(fname)
        card = load_jld2(CardinalPlaceFieldDirectionality, fname)
    else
        gidx = DirectionFiltered(;kwargs...)
        card = CardinalPlaceFieldDirectionality(gidx;nshuffles=nshuffles)
        if do_save
            save_jld2(card, fname)
        end
    end
    card
end

function issignificant_old(gidx::DirectionFiltered;smooth=true, α=0.1, niter=1,pv_threshold=0.01)
    λ = get_direction_tuning(gidx;smooth=smooth, niter=niter, α=α) 
    z = λ.*exp.(gidx.anglebins*im)
    pv = zeros(size(z,2))
    for i in axes(z,2)
        pv[i] = pvalue(RayleighTest(filter(x->abs(x)>0, z[:,i])))
    end
    pv .< pv_threshold
end

function get_pvalue(gidx::DirectionFiltered;smooth=true, α=0.1, niter=1,nruns=1000, kwargs...)
    μr0,ϕ0 = get_directional_tuning_strength(gidx;smooth=smooth, α=α,niter=niter,kwargs...)
    μr = zeros(length(μr0),nruns)
    for i in 1:nruns
        μr[:,i],_ = get_directional_tuning_strength(gidx;do_shuffle=true, smooth=smooth, α=α, niter=niter, kwargs...)
    end
    pv = zeros(length(μr0))
    for i in 1:length(pv)
        if sum(isfinite.(μr[i,:])) > 20
            k = Makie.KernelDensity.kde(μr[i,:],boundary=(0,1))
            xx = sort(μr[i,:])
            cc = sum(pdf(k,xx[xx.<=μr0[i]]))/sum(pdf(k,xx))
            pv[i] = 1-cc
        end
    end
    pv
end

function issignificant(gidx::DirectionFiltered;smooth=true, α=0.1, niter=1,nruns=10000,kwargs...)
    μr0,ϕ0 = get_directional_tuning_strength(gidx;smooth=smooth, α=α,niter=niter,kwargs...)
    μr = zeros(length(μr0),nruns)
    for i in 1:nruns
        μr[:,i],_ = get_directional_tuning_strength(gidx;do_shuffle=true, smooth=smooth, α=α, niter=niter, kwargs...)
    end
    res = fill(false, length(μr0))
    # TODO: Maybe do some kind of interpolation here so that we can evaluate
    #       this with arbitrary p-values
    for i in 1:length(res)
        pp = percentile(μr[i,:], 100*(1-α))
        res[i] = μr0[i] > pp
    end
    res
end

function issignificant_param(gidx::DirectionFiltered;pv_threshold=0.01,kwargs...)
    pv = get_pvalue(gidx;kwargs...)
    pv .< pv_threshold
end

function get_directional_tuning_strength(gidx;idx=1:length(gidx.anglebins), kwargs...)
    λ = get_direction_tuning(gidx;idx=idx, kwargs...)
    get_directional_tuning_strength(λ, gidx.anglebins)
end

function get_directional_tuning_strength(λ::Matrix{<:Real}, anglebins::AbstractVector{<:Real})
    z = λ.*exp.(anglebins*im)
    μr = zeros(size(λ,2))
    ϕ = zeros(size(λ,2))
    for i in axes(λ,2)
        _z = filter(x->abs(x)>0, z[:,i])
        zm = sum(_z)/sum(abs.(_z))
        ϕ[i] = angle(zm)
        μr[i] = abs(zm)
    end
    μr,ϕ
end

abstract type AbstractDirection end
struct East <: AbstractDirection end
Base.string(::Type{East}) = "East"
struct West <: AbstractDirection end
Base.string(::Type{West}) = "West"
struct North <: AbstractDirection end
Base.string(::Type{North}) = "North"
struct South <: AbstractDirection end
Base.string(::Type{South}) = "South"

function get_direction(gidx::DirectionFiltered, ::Type{East};Δθ=π/3)
    union(findall(gidx.anglebins .> π - Δθ/2), findall(gidx.anglebins .< -π+Δθ/2))
end

function get_arrow(mm::SimpleMesh, ::Type{East})
    (-2.5, 0.0,0.0), (5.0, 0.0, 0.0)
end

function get_direction(gidx, ::Type{West};Δθ=π/3)
    union(findall(0 .<= gidx.anglebins .< Δθ/2), findall(-Δθ/2 .<= gidx.anglebins .< 0))
end

function get_arrow(mm::SimpleMesh, ::Type{West})
    (2.5, 0.0,0.0), (-5.0, 0.0, 0.0)
end

function get_direction(gidx, ::Type{North};Δθ=π/3)
    findall( π/2 - Δθ/2 .<= gidx.anglebins .< π/2 + Δθ/2)
end

function get_arrow(mm::SimpleMesh, ::Type{North})
    (0.0, -2.5, 0.0), (0.0, 5.0, 0.0)
end

function get_direction(gidx, ::Type{South};Δθ=π/3)
    findall( -π/2 - Δθ/2 .< gidx.anglebins .<= -π/2 + Δθ/2)
end

function get_arrow(mm::SimpleMesh, ::Type{South})
    (0.0, 2.5, 0.0), (0.0, -5.0, 0.0)
end



function get_view_rate_map(dd::DirectionFiltered, nbins)
    X = zeros(nbins)
    Y = zeros(nbins)
    for k in keys(dd.occupancy)
        k1 = getindex(k,1)
        Y[k1] += dd.occupancy[k] 
        X[k1] += get(dd.weight, k, 0.0)
    end
    X, Y
end

function get_view_rate_map(dd::DirectionFiltered, ii::Integer, nbins, directions::Vector{Vector{Int64}})
    nq = length(directions)
    X = zeros(nbins, nq)
    Y = zeros(nbins, nq)
    for k in keys(dd.occupancy[ii])
        k1 = getindex(k,1)
        k4 = getindex(k,4)
        for (jj,kr) in enumerate(directions)
            if k4 in kr
                Y[k1,jj] += dd.occupancy[ii][k] 
                X[k1,jj] += get(dd.weight[ii],k,0.0)
                break
            end
        end
    end
    X, Y
end

function get_view_rate_map(occupancy::Dict{CartesianIndex{4},Float64}, weight::Dict{CartesianIndex{4},Float64}, view_place_dir_idx::Vector{NTuple{3,Int64}};smooth=false, smoothing_method=:laplace, α=0.1, niter=100)
   # TODO: Make this more general 
    mm = Hippocampus.get_maze_mesh(;nrefinements=2);
    X = zeros(nelements(mm))
    Y = zeros(nelements(mm))
    for (k,v) in occupancy
        (vidx,pidx,hidx,lidx) = Tuple(k)
        if (vidx,pidx,lidx) in view_place_dir_idx
            Y[vidx] += v
            X[vidx] += get(weight, k, 0.0)
        end
    end
    if smooth
        Ls = get_normalize_laplacian(mm)
        X = laplace_smoothing(X, Ls,α;niter=niter)
        Ys = laplace_smoothing(Y, Ls,α;niter=niter)
    else
        Ys = Y
    end
    λ = X./Ys
    λ[Y.==0] .= NaN
    λ[:]
end

"""
    get_direction(pos::Matrix{<:Real}, mm::SimpleMesh, binidx::Vector{<:Integer})

Csompute the directional vector through the field characteried by `binidx`
"""
function get_direction(pos::Matrix{<:Real}, mm::SimpleMesh, binidx::AbstractVector{<:Integer})
    kn = KNearestSearch(mm,1)
    entered = 0 
    exited = 0
    for (ii,_pos) in enumerate(eachcol(pos))
        _idx,dd = searchdists(Meshes.Point(_pos...), kn)
        idx = first(_idx)
        if entered == 0
            if idx in binidx
                entered = ii 
            end
        end
        if entered > 0 && exited == 0
            if (!in(binidx))(idx)
                exited = ii-1 # we want the last point inside the field
                break # we found the exist; stop processing
            end
        end
    end
    if entered > 0
        # if the field was not exited at all, just use the last point
        if exited == 0
            exited = size(pos,2)-1
        end
    end
    return entered,exited 
end

function get_direction(pos::Vector{<:Matrix{<:Real}}, args...)
    n = length(pos)
    idx0 = fill(0,n) 
    idx1 = fill(0,n) 
    for (ii,_pos) in enumerate(pos)
        idx0[ii], idx1[ii] = get_direction(_pos[1:2,:], args...)
    end
    idx0,idx1
end

function get_direction(udata::UnityData, args...)
    nt = numtrials(udata) 
    p0 = zeros(2,nt)
    p1 = zeros(2,nt)
    for i in 1:nt
        tu,posx,posy,_ = get_trial(udata, i;trial_start=2)
        p0[:,i], p1[:,i] = get_direction(permutedims([posx posy]), args...)
    end
    p0,p1
end

function DirectionFiltered(qdata::UnityRaytraceData, vpvrp::ViewAndPlaceRepresentationNew, jocc::JointOccupancy, mm::SimpleMesh, rf::SpatialResponseFields;trial_start=1,only_full_traversal=false, kwargs...)
    nt = numtrials(qdata)
    θbins = range(-π, stop=π, length=24)
    nn = zeros(Int64, nt)
    nm = zeros(Int64, nt)
    clusters = merge_fields(mm, rf.binidx)
    nclusters = Hippocampus.get_num_fields(rf)
    cidx = findall(dropdims(mean(nclusters,dims=2),dims=2) .< 0.001)
    gidx = Vector{Vector{CartesianIndex{5}}}(undef, length(cidx))
    deltaT = Vector{Vector{Float64}}(undef, length(cidx))
    nspikes = Vector{Vector{Float64}}(undef, length(cidx))
    weight = Vector{Dict{CartesianIndex{4}, Float64}}(undef, length(cidx))
    occupancy = Vector{Dict{CartesianIndex{4},Float64}}(undef, length(cidx))
    for i in 1:length(cidx)
        gidx[i] = CartesianIndex{5}[]
        weight[i] = Dict{CartesianIndex{4}, Float64}()
        occupancy[i] = Dict{CartesianIndex{4}, Float64}()
        deltaT[i] = Float64[]
        nspikes[i] = Float64[]
    end
    for i in 1:nt
        for (ll,idx) in enumerate(clusters[cidx])
            # find the first point at which the field is entered and when it is exited
            tg,gaze,pos, fixmask,fo = get_trial(qdata,i;trial_start=trial_start);
            idx0,idx1 = get_direction(pos[1:2,:], mm,rf.binidx[idx])
            if only_full_traversal && idx1==length(tg)-1
                # skip this trial if we did not see a full traversal, i.e. if idx1 is at the penultimate step
                continue
            end
            # debug: check the number of times the field is visited
            qmidx = jocc.index[i][vpvrp.placeviewidx[i]]
            nn[i] = length(filter(k->in(idx)(getindex(k,2)), qmidx))
            # TODO: Also get the gaze for these positions
        
            if idx1 >= idx0 > 0
                v = qdata.position[i][1:2,idx1] - qdata.position[i][1:2,idx0]
                θ = atan(v[2],v[1])
                l = searchsortedfirst(θbins, θ)
                # get the mesh bin indices
                # we actually do not want this; we want all the bins, not just the ones with spikes
                #this is the subset with spikes

                # occupancy
                # @show length(jocc.index[i]), idx0, idx1
                qpidx = jocc.index[i][idx0:idx1]
                tt = qdata.timestamps[i]
                for (jj,qp) in zip(idx0:idx1, qpidx)
                    if qp == CartesianIndex(0,0,0)
                        # skip invalid bins
                        continue
                    end
                    
                    Δt = tt[jj+1]-tt[jj]
                    # since we have filtered out invalid bins above, we do not need to worry about gaps here
                    kk = CartesianIndex(qp[1], qp[2], qp[3], l)
                    occupancy[ll][kk] = get(occupancy[ll], kk, 0.0) + Δt
                    #if time point jj had a spike
                    vv = findfirst(vpvrp.placeviewidx[i].==jj)
                    if vv !== nothing
                        weight[ll][kk] = get(weight[ll], kk, 0.0) + 1.0
                        push!(nspikes[ll], 1.0)
                    else
                        push!(nspikes[ll],0.0)
                    end
                    push!(gidx[ll], CartesianIndex(qp[1], qp[2], qp[3], l, i))
                    push!(deltaT[ll], Δt)
                end
            end
        end
    end
    obj = DirectionFiltered(θbins, weight, occupancy, gidx )
    return obj
end

function process_kwargs(::Type{DirectionFiltered},h::UInt32=zero(UInt32);only_full_traversal=false, kwargs...)
    h = process_kwargs(UnityRaytraceData;kwargs...)
    h = process_kwargs(ViewAndPlaceRepresentationNew,h;kwargs...)
    h = process_kwargs(JointOccupancy,h;kwargs...)
    h = process_kwargs(SpatialResponseFields,h;kwargs...)
    if only_full_traversal
        h = crc32c(string(:only_full_traversal=>true),h)
    end
    h
end

function process_kwargs(::Type{CardinalPlaceFieldDirectionality},h::UInt32=zero(UInt32);nshuffles=10000, kwargs...)
    h = process_kwargs(DirectionFiltered,h;kwargs...)
    h = crc32c(string(:nshuffles=>nshuffles),h)
    h
end

function DirectionFiltered(;redo=fname->false, do_save=true,kwargs...)
    h = process_kwargs(DirectionFiltered;kwargs...)
    fname = "direction_filtered_placefields.jld2"
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    if !redo(fname) && isfile(fname)
        obj = load_jld2(DirectionFiltered,fname)
        #hack: I have no idea why  this is needed
        if !isa(obj, DirectionFiltered)
            obj = first(obj)
        end
    else
        nrefinements = get(kwargs, :nrefinements, (p=3,g=2))
        m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
        sessiondir = DPHT.get_level_path("session")
        qdata, jocc = cd(sessiondir) do
            qdata = UnityRaytraceData(;kwargs...)
            jocc = JointOccupancy(;kwargs...)
            # TODO: Also deal with other filtering here
            qdata, jocc
        end
        vpvrp = ViewAndPlaceRepresentationNew(;kwargs...)
        rf_spatial = get_response_fields(SpatialResponseFields, 10_000;kwargs...)
        obj = DirectionFiltered(qdata, vpvrp, jocc, m_floor, rf_spatial;kwargs...)
        if do_save
            save_jld2(obj, fname)
        end
    end
    obj
end

function get_directionality(qdata::UnityRaytraceData, vpvrp::ViewAndPlaceRepresentationNew, jocc::JointOccupancy, rf::T) where T <: AbstractResponseFields
    mm = get_mesh(T, rf.args[:nrefinements])
    clusters = merge_fields(rf)
    nt = numtrials(qdata)
    λ = fill(NaN, nt, length(clusters))
    θ = fill(NaN, nt, length(clusters))
    gaze = Matrix{Matrix{Float64}}(undef, nt, length(clusters))
    for (ii,cluster) in enumerate(clusters)
        λ[:,ii], θ[:,ii],gg = get_directionality(qdata, vpvrp, mm, rf.binidx[cluster]) 
        for j in 1:nt
            if isassigned(gg, j, ii)
                gaze[j,ii] = gg[j,ii]
            end
        end
    end
    λ, θ, gaze
end

function get_direction_tuning(gidx::DirectionFiltered;idx=1:length(gidx.anglebins), do_shuffle=false, smooth=false, α=0.1, niter=1, spikes_only=false, occupancy_only=false,kwargs...)
    nc = length(gidx.index)
    X = zeros(length(gidx.anglebins),nc)
    Y = zeros(length(gidx.anglebins),nc)
    for j in 1:nc
        for (k,v) in gidx.occupancy[j]
            if do_shuffle
                l = rand(1:length(gidx.anglebins))
            else
                l = getindex(k,4)
            end
            Y[l,j] += v
            if k in keys(gidx.weight[j])
                X[l,j] += gidx.weight[j][k]
            end
        end
    end
    if smooth
        Ls = get_normalize_laplacian(size(X,1))
        X = permutedims(laplace_smoothing(permutedims(X[idx,:]), Ls, α;niter=niter))
        Y = permutedims(laplace_smoothing(permutedims(Y[idx,:]), Ls, α;niter=niter))
    else
        X = X[idx,:]
        Y = Y[idx,:]
    end
    if spikes_only
        return X
    elseif occupancy_only
        return Y
    end
    X ./ Y
end

function get_cardinal_direction_tuning(gidx::DirectionFiltered;do_shuffle=false)
    direction_label = Symbol.([North, South, East, West])
    directions = [get_direction(gidx, d) for d in [North, South, East, West]]
    # group spike counts by cardinal directions
    # TODO: Do not use cardinal directions here (necessarily); rather use the main axis of the field
    # But for now just do this
    spike_count = Dict{Symbol,Vector{Float64}}()
    occupancy = Dict{Symbol,Vector{Float64}}()
    λ = Dict{Symbol,Vector{Float64}}()
    i = 1
    for (k,v) in gidx.occupancy[i]
        qidx = getindex(k,4)
        lidx = findfirst([in(d)(qidx) for d in directions])
        if lidx === nothing
            continue
        end
        dl = direction_label[lidx]
        if !(dl in keys(occupancy))
            spike_count[dl] = Float64[]
            occupancy[dl] = Float64[]
            λ[dl] = Float64[]
        end 
        push!(occupancy[dl], v)
        if k in keys(gidx.weight[i])
            push!(spike_count[dl],gidx.weight[i][k])
            push!(λ[dl], gidx.weight[i][k]/v)
        end
    end
    spike_count, occupancy, λ
end

function get_major_axis_direction_tuning_alt(gidx::DirectionFiltered,rf::SpatialResponseFields;do_shuffle=false)
    v = get_major_axis(rf)
    spike_count = Vector{Dict{Symbol,Vector{Float64}}}(undef, size(v,2))
    occupancy = Vector{Dict{Symbol,Vector{Float64}}}(undef, size(v,2))
    λ = Vector{Dict{Symbol,Vector{Float64}}}(undef, size(v,2))
    for i in 1:length(gidx.occupancy)
        spike_count[i] = Dict{Symbol,Vector{Float64}}()
        occupancy[i] = Dict{Symbol,Vector{Float64}}()
        λ[i] = Dict{Symbol,Vector{Float64}}()
        θ = atan(v[2],v[1])
        # grab a 60 degree cone around θ
        bidx1 = findall(cos.(gidx.anglebins .- θ) .> cos(π/6))
        f1 = in(bidx1)
        # grab a 60 degree cone around the opposite direction
        bidx2 = findall(cos.(gidx.anglebins .- (θ-π)) .> cos(π/6))
        f2 = in(bidx2)
        for (k,v) in gidx.occupancy[i]
            qidx = getindex(k,4)
            if f1(qidx)
                dl = :forward
            elseif f2(qidx)
                dl = :backward
            else
                continue
            end
            if !(dl in keys(occupancy[i]))
                spike_count[i][dl] = Float64[]
                occupancy[i][dl] = Float64[]
                λ[i][dl] = Float64[]
            end 
            push!(occupancy[i][dl], v)
            if k in keys(gidx.weight[i])
                push!(spike_count[i][dl],gidx.weight[i][k])
                push!(λ[i][dl], gidx.weight[i][k]/v)
            else
                push!(λ[i][dl], 0.0)
            end
        end
    end
    spike_count, occupancy, λ
end

struct MajorAxisDirectionTuning{T<:Real}
    v::Matrix{T} # direciton for each place field
    spike_count_forward::Vector{Vector{T}}
    spike_count_reverse::Vector{Vector{T}}
    occupancy_forward::Vector{Vector{T}}
    occupancy_reverse::Vector{Vector{T}}
    trialidx_forward::Vector{Vector{Int64}}
    trialidx_reverse::Vector{Vector{Int64}}
end

function process_kwargs(::Type{MajorAxisDirectionTuning},h::UInt32=zero(UInt32);kwargs...)
    h = process_kwargs(SpatialResponseFields,h;kwargs...)
    h = process_kwargs(DirectionFiltered, h;kwargs...)
    h
end

function get_major_axis_direction_tuning(gidx::DirectionFiltered,rf::SpatialResponseFields,jocc::JointOccupancy, vpvrp::ViewAndPlaceRepresentationNew,qdata::UnityRaytraceData;do_shuffle=false,trial_start=2)
    v = get_major_axis(rf)
    spike_count = Float64[]
    spike_count_1 = Vector{Vector{Float64}}(undef, size(v,2))
    spike_count_2 = Vector{Vector{Float64}}(undef, size(v,2))
    occupancy_1 = Vector{Vector{Float64}}(undef, size(v,2))
    occupancy_2 = Vector{Vector{Float64}}(undef, size(v,2))
    trialidx_forward = Vector{Vector{Int64}}(undef, size(v,2))
    trialidx_reverse = Vector{Vector{Int64}}(undef, size(v,2))
    occupancy = Float64[]
    sindex = CartesianIndex{5}[]
    for i in 1:length(spike_count_1)
        spike_count_1[i] = Float64[]
        spike_count_2[i] = Float64[]
        occupancy_1[i] = Float64[]
        occupancy_2[i] = Float64[]
        trialidx_forward[i] = Int64[]
        trialidx_reverse[i] = Int64[]
        θ = atan(v[2,i],v[1,i])
        bidx1 = findall(cos.(gidx.anglebins .- θ) .> cos(π/6))
        f1 = in(bidx1)
        # grab a 60 degree cone around the opposite direction
        bidx2 = findall(cos.(gidx.anglebins .- (θ-π)) .> cos(π/6))
        f2 = in(bidx2)
        # get all the trials
        gindex = gidx.index[i]
        trialidx= sort(unique(getindex.(gindex, 5)))
        for tidx in trialidx

            timestamps,gaze,pos, fixmask,fo = get_trial(qdata,tidx;trial_start=trial_start);
            placeviewidx = vpvrp.placeviewidx[tidx]
            # find the corresponding indices in jocc
            idx1 = [CartesianIndex(k[1],k[2],k[3]) for k in filter(q->(q[5]==tidx)&(f1(q[4])), gindex)]
            vidx1 = findall(in(idx1), jocc.index[tidx])

            idx2 = [CartesianIndex(k[1],k[2],k[3]) for k in filter(q->(q[5]==tidx)&(f2(q[4])), gindex)]
            vidx2 = findall(in(idx2), jocc.index[tidx])
            # vidx1/2 now contains all points for this trial where the subject was within the spatial field and moving either along
            # the major axis (1) or opposite to the major axis (2)
            # for each of them, grap the occupancy and the spike count
            # TODO: Do we want a trial by trial measure here instead? 
            #       That is, the activity for each trial that the subject pass through
            #       the field either in the forward or reverse direction
            if !isempty(vidx1)
                push!(occupancy_1[i], sum(timestamps[vidx1.+1] - timestamps[vidx1]))
                push!(spike_count_1[i], sum(in(placeviewidx).(vidx1)))
                push!(trialidx_forward[i], tidx)
            end
            if !isempty(vidx2)
                push!(spike_count_2[i], sum(in(placeviewidx).(vidx2)))
                push!(occupancy_2[i], sum(timestamps[vidx2.+1] - timestamps[vidx2]))
                push!(trialidx_reverse[i], tidx)
            end

            for (kq, vidx) in enumerate([vidx1, vidx2])
                w = timestamps[vidx.+1] - timestamps[vidx]
                sp = in(placeviewidx).(vidx)
                append!(spike_count, Float64.(sp))
                append!(occupancy, w)
                for qq in jocc.index[tidx][vidx]
                    push!(sindex, CartesianIndex(qq[1], qq[2], qq[3], kq, tidx))
                end
            end
        end
    end
    # TODO: Aggregate over windows
    #spike_count, occupancy, sindex, 
    (spike_count_1, occupancy_1), (spike_count_2, occupancy_2), (trialidx_forward, trialidx_reverse)
end

function MajorAxisDirectionTuning(;redo=fname->false, do_save=true, kwargs...)
    fname = "major_axis_direction_tuning.jld2"
    h = process_kwargs(MajorAxisDirectionTuning;kwargs...)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    if !redo(fname) && isfile(fname)
        obj = load_jld2(MajorAxisDirectionTuning, fname)
    else
        sessiondir = DPHT.get_level_path("session")
        qdata, jocc = cd(sessiondir) do
            qdata = UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=redo)
            jocc = Hippocampus.JointOccupancy(;redo=redo, do_save=false, kwargs...)
            qdata, jocc
        end
        rf_spatial =  get_response_fields(SpatialResponseFields, get(kwargs, :nshuffles, 10_000);redo=redo,kwargs...)
        v = get_major_axis(rf_spatial)
        gidx = DirectionFiltered(;redo=redo,kwargs...)
        vpvrp = ViewAndPlaceRepresentationNew(;redo=redo,kwargs...)
        (spike_count_1, occupancy_1), (spike_count_2, occupancy_2),(trialidx_forward, trialidx_reverse) = get_major_axis_direction_tuning(gidx, rf_spatial, jocc, vpvrp,qdata;trial_start=get(kwargs, :trial_start,2))
        obj = MajorAxisDirectionTuning(v, spike_count_1, spike_count_2,occupancy_1, occupancy_2, trialidx_forward, trialidx_reverse)
        if do_save
            save_jld2(obj, fname)
        end
    end
    return obj
end

function issignificant(madt::MajorAxisDirectionTuning;pv_threshold=0.05)
    res = fill(false, length(madt.spike_count_forward))
    for i in 1:length(res)
        λ1 = madt.spike_count_forward[i]./madt.occupancy_forward[i]
        λ2 = madt.spike_count_reverse[i]./madt.occupancy_reverse[i]
        μ = mean(λ1) - mean(λ2)
        if all(length.((λ1, λ2)).>=5)
            q1,q2 = permutation_test(mean, λ1,λ2) 
            res[i] = μ > percentile(q1-q2, 100*(1-pv_threshold)) || μ < percentile(q1-q2, 100*pv_threshold)
        end
    end
    res
end

function aggregate(spike_count::Vector{T}, occupancy::Vector{T}, sindex::Vector{CartesianIndex{5}};window=0.05) where T <: Real
    tidx,vidx = (getindex(sindex[1],5), getindex(sindex[1], 4))
    midx = findall(q->(q[4]==vidx)&(q[5]==tidx), sindex)
    wsum = sum(occupancy[midx])
    mm = round(Int64,floor(wsum/window))
    _window = wsum/mm

    spike_count_new = T[]
    occupancy_new = T[]
    sindex_new = CartesianIndex{5}[]
    w = zero(T)
    ss = zero(T)
    # TODO Balance w so that it covers all points within a trial roughyl equally
    for (sp,occ,sidx) in zip(spike_count, occupancy, sindex)
        if (sidx[5] == tidx) && (sidx[4]==vidx) && (w <= _window)
            w += occ
            ss += sp
        else
            push!(spike_count_new, ss)
            push!(occupancy_new, w)
            push!(sindex_new, sidx)
            w = zero(T)
            ss = zero(T)
            tidx = sidx[5]
            vidx = sidx[4]
             midx = findall(q->(q[4]==vidx)&(q[5]==tidx), sindex)
            wsum = sum(occupancy[midx])
            mm = round(Int64,floor(wsum/window))
            _window = wsum/mm
        end
    end
    spike_count_new, occupancy_new, sindex_new
end

function get_cardinal_direction_tuning(gidx::DirectionFiltered,card::CardinalPlaceFieldDirectionality;do_shuffle=false)
    direction_label = Symbol.([North, South, East, West])
    directions = [get_direction(gidx, d) for d in [North, South, East, West]]
    # group spike counts by cardinal directions
    # TODO: Do not use cardinal directions here (necessarily); rather use the main axis of the field
    # But for now just do this
    spike_count = Dict{Symbol,Vector{Float64}}()
    occupancy = Dict{Symbol,Vector{Float64}}()
    λ = Dict{Symbol,Vector{Float64}}()
    i = 1
    for (k,v) in gidx.occupancy[i]
        pidx = getindex(k,2)
        qidx = getindex(k,4)
        lidx = findfirst([in(d)(qidx) for d in directions])
        if lidx === nothing
            continue
        end
        dl = direction_label[lidx]
        if !(dl in keys(occupancy))
            spike_count[dl] = Float64[]
            occupancy[dl] = Float64[]
            λ[dl] = Float64[]
        end 
        push!(occupancy[dl], card.occupancy[pidx,lidx,i])
        push!(spike_count[dl],card.weight[pidx,lidx,i])
        push!(λ[dl], card.weight[pidx,lidx,i]/card.occupancy[pidx,lidx,i])
    end
    spike_count, occupancy, λ
end

function get_cardinal_direction_tuning_old(gidx::DirectionFiltered, deltaT, nspikes;nt::Union{Int64,Nothing}=nothing,do_shuffle=false)
    qidx = gidx.index
    # we need to get nt from somehwere
    # get nt from the index for now
    nc = length(qidx)
    if nt === nothing
        nt = 0
        for i in 1:nc
            nt = max(nt, maximum(getindex.(qidx[i],5)))
        end
    end
    directions = [get_direction(gidx, d) for d in [North, South, East, West]]
    X = zeros(4,nt,nc)
    Y = zeros(4,nt,nc)
    for j in 1:nc
        if do_shuffle
            _nspikes = shuffle(nspikes[j])
        else
            _nspikes = nspikes[j]
        end
        for (k,w,x) in zip(qidx[j], deltaT[j], _nspikes)
            l = getindex(k,4)
            tidx = getindex(k,5)
            idx = findfirst(dd->in(dd)(l), directions)
            if idx !== nothing
                Y[idx,tidx,j] += w
                X[idx,tidx,j] += x
            end
        end
    end
    X,Y
end

function get_cardinal_direction_tuning_strength_old(gidx::DirectionFiltered, deltaT, nspikes;nt::Union{Int64,Nothing}=nothing,do_shuffle=false)
    X,Y = get_cardinal_direction_tuning(gidx, deltaT, nspikes;nt=nt)
    λ = X./(Y .+ sqrt(eps(T)))
    Hx = compute_skaggs_sic(λ[:], Y[:])
    Ys = dropdims(sum(Y,dims=2),dims=(2,3))
    Ps = Ys./sum(Ys)
    Pxs = λ.*deltaT
end

function get_cardinal_direction_tuning_strength(gidx;nshuffles=1000)
    m_floor = Shadow("xy")(Hippocampus.floor_topology3(;nrefinements=3));
    Ls = get_normalize_laplacian(m_floor)
    X,Y = get_joint_cardinal_prob(gidx)
    Xs = laplace_smoothing(X, Ls, 0.1;niter=50); 
    Ys = laplace_smoothing(Y, Ls, 0.1;niter=50); 

    ees = zeros(size(X,3),nshuffles)
    ee = zeros(size(X,3))
    for j in 1:size(X,3)
        ee[j]= get_conditional_information(Xs[:,:,j],Ys[:,:,j])
    end
    for i in 1:nshuffles
        X,Y = get_joint_cardinal_prob(gidx;do_shuffle=true)
        Xs = laplace_smoothing(X, Ls, 0.1;niter=50); 
        Ys = laplace_smoothing(Y, Ls, 0.1;niter=50); 
        for j in 1:size(X,3)
            ees[j,i] = get_conditional_information(Xs[:,:,j],Ys[:,:,j])
        end
    end
    ee, ees, X,Y
end

function get_joint_cardinal_prob(gidx;do_shuffle=false)
    directions = [get_direction(gidx, d) for d in [North, South, East, West]]
    nc = length(gidx.index)
    X = zeros(1344,4,nc)
    Y = zeros(1344,4,nc)
    for (j,(occupancy, weight)) in enumerate(zip(gidx.occupancy, gidx.weight))
        for (k,v) in occupancy
            if do_shuffle
                lidx = rand(1:length(gidx.anglebins))
            else
                lidx = getindex(k,4)
            end
            pidx = getindex(k,2)
            idx = findfirst(dd->in(dd)(lidx), directions)
            if idx !== nothing
                Y[pidx,idx,j] += v
                if k in keys(weight)
                    X[pidx,idx,j] += weight[k] 
                end
            end
        end
    end
    X,Y
end

"""
Compute the conditional information between the firing rate X and bins represented by the second
dimennsion of X, conditioned on the bins represented by the first dimension of X.
"""
function get_conditional_information(X::Matrix{T},Y::Matrix{T}) where T <: Real
    λ = X./(Y .+ eps(Float64)) # to avoid zeros
    Pdy = Y./sum(Y)
    Psdy = λ.*Pdy  
    Py = sum(Pdy,dims=2)  
    Psy = sum(Psdy,dims=2)
    λm = sum(filter(isfinite, Psdy))
    sum(filter(isfinite, Psdy.*log2.((Psdy.*Py)./(Psy.*Pdy))))/λm
end

"""
Return the bins in `mm` that wall within a 60 degree wedge centered on `pos`
"""
function get_view_bins(mm::SimpleMesh, pos, ϕ;Δϕ=π/3)
    midx = Int64[]
    # elevation in 120 degree wedge centered on π/2
    for ξ in range(π/2-π/3, stop=π/2+π/3, length=15)
        ϕ0 = ϕ-Δϕ/2
        ϕ1 = ϕ+Δϕ/2
        while ϕ0 < ϕ1 
            r0 = Ray(Tuple(pos), (sin(ξ)*cos(ϕ0),sin(ξ)*sin(ϕ0), cos(ξ)))
            # find all points where the ray intersects the maze
            iq0 = [Meshes.intersect(r0, m) for m in mm] 
            qidx0 = findall(iq0.!==nothing)
            idx0 = qidx0[argmin(norm.(centroid.(mm[qidx0]) .- pos))] 
            if !(idx0 in midx)
                push!(midx, idx0)
            end
            ϕ0 += Δϕ/15 # 10 steps
        end
    end
    midx
end

function get_view_bins(anglebins::AbstractVector{<:Real}, mm::SimpleMesh, m_floor::SimpleMesh, place_dir_idx::Vector{Tuple{Int64, Int64}};z=0.5, kwargs...)
    left_gaze_place_dir_idx = NTuple{3, Int64}[]
    right_gaze_place_dir_idx = NTuple{3, Int64}[]
    for (pidx, lidx) in place_dir_idx
        # player position
        pos = centroid.(m_floor[pidx]) + Meshes.Vec(0.0, 0.0, z)
        ϕ = anglebins[lidx]

        midx_left = get_view_bins(mm, pos, ϕ+π/2)
        for m in unique(midx_left)
            cc = (m,pidx,lidx)
            if !(cc in left_gaze_place_dir_idx)
                push!(left_gaze_place_dir_idx, (m, pidx, lidx))
            end
        end
        midx_right = get_view_bins(mm, pos, ϕ-π/2)
        for m in unique(midx_right)
            cc = (m,pidx,lidx)
            if !(cc in right_gaze_place_dir_idx)
                push!(right_gaze_place_dir_idx, (m, pidx, lidx))
            end
        end
    end 
    left_gaze_place_dir_idx, right_gaze_place_dir_idx
end

function get_egocentric_view_tuning(gidx, mm::SimpleMesh, m_floor::SimpleMesh;kwargs...)
    λright = Vector{Vector{Float64}}(undef, length(gidx.occupancy))
    λleft = Vector{Vector{Float64}}(undef, length(gidx.occupancy))
    for i in 1:length(gidx.occupancy)
        # get the place bins and traversal directions for this field
        vq = [(k[2], k[4]) for k in keys(gidx.occupancy[i])]
        left_gaze_place_dir_idx, right_gaze_place_dir_idx = get_view_bins(gidx.anglebins, mm, m_floor, unique(vq))
        # seperate the occupancy and weight into left and right
        λleft[i] = Float64[]
        λright[i] = Float64[]
        # create a view map first, by conditioning on place and direction
        # then sum up bins according to left_.. and right_.. above
        for (kk,vv) in gidx.occupancy[i]
            if (kk[1], kk[2], kk[4]) in left_gaze_place_dir_idx
                # TODO: Allow smoothing here
                #       This is the view map, 
                λ = get(gidx.weight[i], kk, 0.0)/vv
                push!(λleft[i], λ)
            elseif (kk[1], kk[2], kk[4]) in right_gaze_place_dir_idx
                λ = get(gidx.weight[i], kk, 0.0)/vv
                push!(λright[i], λ)
            end
        end
    end
    λleft, λright
end

"""
Estimate the preference for left vs right gaze relative to the direction of traverersal through the field
"""
function get_egocentric_gaze(gidx;pv_threshold=0.01,z=0.5)
    m_floor = (Hippocampus.floor_topology3(;nrefinements=3))
    mm = Hippocampus.get_maze_mesh(;nrefinements=2)
    # this only makes sense for a direction field, so test that first
    # actually, we can still get egocentric coding even without directionality
    # the cell could be coding left/right gaze for any direction of traversal
    μr,ϕ = get_directional_tuning_strength(gidx;do_shuffle=false, smooth=true,niter=1)
    midx = Vector{Vector{CartesianIndex{6}}}(undef, length(μr))
    for i in 1:length(μr)
        midx[i] = CartesianIndex{6}[]
        for (j,(θ ,direction)) in enumerate(zip([π, π/2, 0.0, -π/2], [West, North, East, South]))
                # traversal direction
            aidx = get_direction(gidx, direction)
            for k in gidx.index[i]
                if k[4] in aidx
                    # the player is located z distance from the floor
                    pos = centroid(m_floor[k[2]]) + Meshes.Vec(0.0, 0.0, z)
                    # in the correct direction
                    for (kk,Δϕ) in enumerate([π/2, -π/2])
                        # TODO: Do we need to use multiple rays here; 
                        ϕ0 = θ+Δϕ-π/6
                        ϕ1 = θ+Δϕ+π/6
                        # TODO Also go vertical
                        for ξ in range(-π/6, stop=π/6, length=10)
                            while ϕ0 < ϕ1 
                                r0 = Ray(Tuple(pos), (sin(ξ)*cos(ϕ0),sin(ξ)*sin(ϕ0), cos(ξ)))
                                # find all points where the ray intersects the maze
                                iq0 = [Meshes.intersect(r0, m) for m in mm] 
                                qidx0 = findall(iq0.!==nothing)
                                idx0 = qidx0[argmin(norm.(centroid.(mm[qidx0]) .- pos))] 
                                cc = CartesianIndex(idx0, k[2], k[3], k[4], kk, k[5] )
                                if !(cc in mix[i])
                                    push!(midx[i], cc)
                                end
                                ϕ0 += π/30 # 10 steps
                            end
                        end
                    end

                end
                
            end
        end
    end
    midx
end


## plots

function plot_field_directionality!(lg::GridLayout, λ::AbstractVector{<:Real}, θ::AbstractVector{<:Real};kwargs...)
    # compute circular mean
    qidx = isfinite.(λ)
    x = abs(sum(λ[qidx].*exp.(-im*θ[qidx]))/sum(λ[qidx]))
    θbins = range(-π, stop=π, length=24)
    h = fit(Histogram, θ, Weights(λ),θbins)
    ax = PolarAxis(lg[1,1],)
    hidedecorations!(ax)
    lines!(ax, θbins, [h.weights;h.weights[1]];linewidth=2.0, kwargs...)
    ax.thetaticklabelsvisible = false
    ax.rticksvisible = true
    ax.thetagridvisible=true
    ax.rgridvisible=true
    Label(lg[1,2], L"\mu_r=%$(round(x, sigdigits=2))", tellheight=false,rotation=-π/2)
end

function plot_field_directionality(args...)
    with_theme(plot_theme) do
        fig = Figure(size=(700,400))
        lg = GridLayout(fig[1,1])
        plot_field_directionality!(lg,args...)
        fig
    end
end

function plot_field_directionality!(lg, gidx::DirectionFiltered, rf::T) where T <: AbstractResponseFields
    λ = get_direction_tuning(gidx)
    plot_field_directionality!(lg, λ, gidx.anglebins, rf)
end

function plot_field_directionality!(lg, λ::Matrix{<:Real}, θ::AbstractVector{<:Real}, rf::T) where T <: AbstractResponseFields
    lg1 = GridLayout(lg[1,1])
    Label(lg1[1,1,TopLeft()], "A")
    plot_response_fields!(lg1, rf)
    lg2 = GridLayout(lg[1,2])
    Label(lg2[1,1,TopLeft()], "B")
    lgi = [GridLayout(lg2[i,1]) for i in 1:size(θ,2)]
    colors = Makie.wong_colors()
    for (ii,_lg) in enumerate(lgi)
        plot_field_directionality!(_lg, λ[:,ii], θ;color=colors[ii])
    end
    colsize!(lg, 1, Relative(0.7))
end

function plot_field_directionality!(lg, gidx::DirectionFiltered, rf::SpatialResponseFields;kwargs...)
    lg1 = GridLayout(lg[1,1])
    Label(lg1[1,1,TopLeft()], "A")
    plot_response_fields!(lg1, rf)
    lg2 = GridLayout(lg[1,2])
    plot_directional_tuning!(lg2, gidx;kwargs...)
    Label(lg2[0,1,TopLeft()], "B")
    Label(lg2[0,1, Top()], "Firing rate", tellwidth=false, fontsize=14)
    rowsize!(lg2, 0, 5)

    lg3 = GridLayout(lg[1,3])
    plot_directional_tuning!(lg3, gidx;occupancy_only=true, kwargs...)
    Label(lg3[0,1,TopLeft()], "C")
    Label(lg3[0,1, Top()], "Occupancy", tellwidth=false, fontsize=14)
    rowsize!(lg3, 0, 5)
    colsize!(lg, 1, Relative(0.6))
end

function plot_directional_view_field(dd::DirectionFiltered, args...;kwargs...)
    with_theme(plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_directional_view_field!(lg, dd, args...;kwargs...)
        link_cameras_lscene(fig)
        fig
    end
end

function plot_directional_view_field!(lg, dd::DirectionFiltered, place_field_idx::AbstractVector{<:Integer},idx::Integer;floor_offset=-30, directions=[North, South])
    # south to north vs north to south
    # TODO: Should be tailored to each place field
    mm = Hippocampus.get_maze_mesh(;nrefinements=2)
    binranges = [get_direction(dd, d) for d in directions]
    X,Y = get_view_rate_map(dd, idx, 1312,binranges)
    Ls = get_normalize_laplacian(mm)
    Xs = laplace_smoothing(permutedims(X), Ls, 0.1;niter=100);
    Ys = laplace_smoothing(permutedims(Y), Ls, 0.1;niter=100);
    λ = Xs./Ys
    cr = extrema(filter(isfinite, λ))
    m_floor = Translate(0.0, 0.0, floor_offset)(Hippocampus.floor_topology3(;nrefinements=3))
    cm = Meshes.Point(mean(to.(centroid.(m_floor[place_field_idx])))...)
    cmp = Point3f(ustrip(cm.coords.x), ustrip(cm.coords.y), ustrip(cm.coords.z))
    Z = zeros(nelements(m_floor))
    Z[place_field_idx] .= 1.0
   # Label(lg[1,1], "$(directions[1]) → $(directions[2])", tellwidth=false)
    lscene1 = LScene(lg[2,1], show_axis=false)
    plotmesh!(lscene1, mm;alpha=0.0, showsegments=true, segmentcolor=:darkgray, ceiling_offset=10, floor_offset=-15)
    plotmesh!(lscene1, mm;color=Xs[1,:]./Ys[1,:], ceiling_offset=10, floor_offset=-15, colorrange=cr)
    viz!(lscene1, m_floor;color=Z)
    # indicate directionality
    a,b = get_arrow(mm, directions[2])
    arrows3d!(lscene1, cmp-Point3f(a), Point3f(b))

   # Label(lg[1,2], "$(directions[2]) → $(directions[1])",tellwidth=false)
    lscene2 = LScene(lg[2,2], show_axis=false)
    plotmesh!(lscene2, mm;alpha=0.0, showsegments=true, segmentcolor=:darkgray, ceiling_offset=10, floor_offset=-15)
    plotmesh!(lscene2, mm;color=Xs[2,:]./Ys[2,:], ceiling_offset=10, floor_offset=-15, colorrange=cr)
    viz!(lscene2, m_floor;color=Z)
    a,b = get_arrow(mm, directions[1])
    arrows3d!(lscene2, cmp-Point3f(a), Point3f(b))

    # combined
    Xsc = laplace_smoothing(dropdims(sum(X,dims=2),dims=2), Ls, 0.1;niter=100);
    Ysc = laplace_smoothing(dropdims(sum(Y,dims=2),dims=2), Ls, 0.1;niter=100);
    λc = Xsc./Ysc
    Label(lg[1,3], "Combined", tellwidth=false)
    lscene3 = LScene(lg[2,3], show_axis=false)
    plotmesh!(lscene3, mm;alpha=0.0, showsegments=true, segmentcolor=:darkgray, ceiling_offset=10, floor_offset=-15)
    plotmesh!(lscene3, mm;color=λc, ceiling_offset=10, floor_offset=-15, colorrange=cr)
    viz!(lscene3, m_floor;color=Z)

    # occupancy only 
    Label(lg[1,4], "Occupancy", tellwidth=false)
    lscene4 = LScene(lg[2,4], show_axis=false)
    plotmesh!(lscene4, mm;alpha=0.0, showsegments=true, segmentcolor=:darkgray, ceiling_offset=10, floor_offset=-15)
    plotmesh!(lscene4, mm;color=Ysc, ceiling_offset=10, floor_offset=-15)
    viz!(lscene4, m_floor;color=Z)
end

function plot_directional_tuning!(lg, gidx::DirectionFiltered;kwargs...)
    # mild smoothing
    λ = get_direction_tuning(gidx;kwargs...)
    μr,ϕ = get_directional_tuning_strength(gidx;kwargs...)
    @show μr, ϕ
    pv = get_pvalue(gidx;kwargs...)
    for i in axes(λ,2)
        ax = PolarAxis(lg[i,1])
        hidedecorations!(ax)
        ax.thetagridvisible = true
        ax.rgridvisible = true
        ax.spinevisible = false
        ax.rticklabelsize = 10
        ax.thetaticklabelsize=10
        ax.rticksize = 0
        ax.thetaticksize = 0
        ax.thetaticklabelpad = 0
        lines!(ax, gidx.anglebins, λ[:,i], color=Cycled(i))
        ym = maximum(filter(isfinite, λ[:,i]))
        ym *= μr[i]
        linesegments!(ax,[ϕ[i],ϕ[i]], [0.0, ym], color=:red, linewidth=1)
        # indicate siginifance
        if pv[i] < 0.001
            lq = "**"
        elseif pv[i] < 0.01
            lq = "*"
        else
            lq = "ns"
        end
        ax.title = lq
        ax.titlesize = 12
        ax.titlegap = 0
        if i > 1
            rowgap!(lg, i-1, 0)
        end
        #Label(lg[i,2], lq, rotation=-π/2, tellheight=false,fontsize=14)
        #Label(lg[i,2], L"$μ_r = %$(round(μr[i], sigdigits=2))$ \\ $p = %$(round(pv[i], sigdigits=2))$", rotation=-π/2, tellheight=false, 
         #               fontsize=14)
    end
end

function plot_directional_tuning(gidx::DirectionFiltered;kwargs...)
    with_theme(plot_theme) do
        fig = Figure(size=(350,600))
        lg = GridLayout(fig[1,1])
        plot_directional_tuning!(lg, gidx;kwargs...)
        fig
    end
end

function plot_egocentric_gaze_tuning!(lg, λleft::Vector{<:Vector{<:Real}}, λright::Vector{<:Vector{<:Real}})
    for i in 1:length(λleft)
        ax = Axis(lg[i,1])
        fidx_left = findall(isfinite, λleft[i])
        fidx_right = findall(isfinite, λright[i])
        xx = [fill(1.0, length(fidx_left));fill(2.0, length(fidx_right))]
        yy = [λleft[i][fidx_left];λright[i][fidx_right]]
        boxplot!(ax, xx, yy, color=Cycled(i),show_outliers=false)
        ax.xticks = (1:2, ["Left","Right"])
        xlims!(ax, 0.5, 2.5)
        ax.ylabel = "Firing rate [Hz]"
    end
end

function plot_egocentric_gaze_tuning(λleft::Vector{<:Vector{<:Real}}, λright::Vector{<:Vector{<:Real}})
    with_theme(plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_egocentric_gaze_tuning!(lg, λleft, λright)
        fig
    end
end

function plot_field_directionality(rf_spatial, gidx, λleft, λright)
    with_theme(plot_theme) do
        fig = Figure(size=(1024, 600))
        lg1 = GridLayout(fig[1,1])  
        plot_field_directionality!(lg1, gidx, rf_spatial)
        lg2 = GridLayout(fig[1,2])
        plot_egocentric_gaze_tuning!(lg2, λleft, λright)
        Label(lg2[1,1,TopLeft()], "D")
        colsize!(fig.layout, 2, Relative(0.2))
        fig
    end
end

function plot_left_vs_right_rate(args...)
    with_theme(plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_left_vs_right_rate!(lg, args...)
        link_cameras_lscene(fig)
        fig
    end
end

function plot_left_vs_right_rate!(lg, gidx::DirectionFiltered, left_gaze_place_dir_idx::Vector{NTuple{3,Int64}}, right_gaze_place_dir_idx::Vector{NTuple{3,Int64}}, mm::SimpleMesh, m_floor::SimpleMesh)
    occupancy = gidx.occupancy[1]
    weight = gidx.weight[1]
    lscene_left = LScene(lg[1,1], show_axis=false)
    lscene_right = LScene(lg[1,2], show_axis=false)

     λ_left = Hippocampus.get_view_rate_map(occupancy, weight, left_gaze_place_dir_idx)
     λ_right = Hippocampus.get_view_rate_map(occupancy, weight, right_gaze_place_dir_idx)

    cr = extrema(filter(isfinite, [λ_left;λ_right]))

    plotmesh!(lscene_left, mm, color=:lightgray, showsegments=true, hide_ceiling=true, hide_floor=true)
    plotmesh!(lscene_left, mm, color=λ_left , colorrange=cr, showsegments=false, hide_ceiling=true, hide_floor=false)

    plotmesh!(lscene_right, mm, color=:lightgray, showsegments=true, hide_ceiling=true, hide_floor=true)
    plotmesh!(lscene_right, mm, color=λ_right , colorrange=cr, showsegments=false, hide_ceiling=true, hide_floor=false)
    Colorbar(lg[1,3], colorrange=cr, label="Firing rate [Hz]")
    Label(lg[1,1,Top()], "Left gaze", tellwidth=false)
    Label(lg[1,2,Top()], "right gaze", tellwidth=false)
    # TODO: Indicate place field and direction
    place_idx = unique(getindex.([left_gaze_place_dir_idx;right_gaze_place_dir_idx],2))
    angle_idx = mode(getindex.([left_gaze_place_dir_idx;right_gaze_place_dir_idx],3))
    ϕ = gidx.anglebins[angle_idx]
    #bb = find_boundary(m_floor[place_idx])
    m_floor2 = Translate(0.0, 0.0, -0.1)(m_floor)
    bb = m_floor2[place_idx]
    cm = mean(Point3f.(Tuple.(centroid.(bb))))
    viz!(lscene_left, m_floor2, color=:darkgray, showsegments=true)
    viz!(lscene_left, bb, color=:blue)
    arrows3d!(lscene_left, cm,2.0*Point3f(cos(ϕ), sin(ϕ), 0.0), color=:orange)
    viz!(lscene_right, m_floor2, color=:darkgray, showsegments=true)
    viz!(lscene_right, bb, color=:blue)
    arrows3d!(lscene_right, cm,2.0*Point3f(cos(ϕ), sin(ϕ), 0.0), color=:orange)
end

function plot_field_traversals!(lg, gidx::DirectionFiltered, udata::UnityData;floorcolor=:lightgray, kwargs...)
    # get the trajectories 
    trialidx = unique(getindex.(gidx.index[1], 5))
    trajectories = Point2f[]
    trajectory_color = Float64[]
    for tidx in trialidx
        tu,posx, posy,_ = get_trial(udata,tidx;trial_start=2);
        color_idx = range(0.0, stop=1.0, length=length(posx))
        for p in zip(posx, posy) 
            push!(trajectories, Point2f(p...))
        end
        push!(trajectories, Point2f(NaN))
        append!(trajectory_color, collect(color_idx))
        push!(trajectory_color, NaN)
    end
    ax = Axis(lg[1,1], aspect=1)
    m_floor = Shadow("xy")(floor_topology3(;nrefinements=3))
    viz!(ax, m_floor;color=floorcolor)
    hidedecorations!(ax)
    ax.bottomspinevisible = false
    ax.leftspinevisible = false
    colormap = get(kwargs, :colormap, :viridis)
    ll = lines!(ax, trajectories, color=trajectory_color, colormap=colormap)
    # color bar
    Colorbar(lg[2,1],ll, vertical=false, flipaxis=false,ticks=([0.05,0.95],["start","end"]),label="Subject position")
    ax
end

function plot_field_traversals(gidx, qdata;_plot_theme=plot_theme, kwargs...)
    with_theme(_plot_theme) do
        fig = Figure(size=(600,600))
        lg = GridLayout(fig[1,1])
        plot_field_traversals!(lg, gidx, qdata;kwargs...)
        fig
    end
end
