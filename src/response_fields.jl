using Meshes
using BetaKDE
abstract type AbstractResponseFields end

struct SpatialResponseFields <: AbstractResponseFields
    λ::Vector{Float64}
    binidx::Vector{Int64}
    gamma_params::Matrix{Float64} # gaamma parameter fit for the null distribution; two parameters per bin
    args::Dict{Symbol,Any}
end

abstract type AbstractResponseFieldsSimple <: AbstractResponseFields end

struct SpatialResponseFieldsSimple <: AbstractResponseFieldsSimple
    λ::Vector{Float64}
    fieldindices::Vector{Vector{Int64}}
    baseline_percentile_threshold::Float64
    peak_percentile_threshold::Float64
    peak_threshold::Float64 
    args::Dict{Symbol,Any}
end

SpatialResponseFieldsAll = Union{SpatialResponseFields, SpatialResponseFieldsSimple}

DPHT.filename(::Type{SpatialResponseFieldsSimple}) = "spatial_response_fields_simple.jld2"

struct GazeResponseFieldsSimple <: AbstractResponseFieldsSimple
    λ::Vector{Float64}
    fieldindices::Vector{Vector{Int64}}
    baseline_percentile_threshold::Float64
    peak_percentile_threshold::Float64
    peak_threshold::Float64 
    args::Dict{Symbol,Any}
end

DPHT.filename(::Type{GazeResponseFieldsSimple}) = "gaze_response_fields_simple.jld2"

function process_kwargs(::Type{<:AbstractResponseFieldsSimple},h::UInt32=zero(UInt32);baseline_percentile_threshold=10, peak_percentile_threshold=95, peak_threshold=0.5, kwargs...)
    h = process_kwargs(JointMap,h;kwargs...)
    h = CRC32c.crc32c(string(:baseline_percentile_threshold=>baseline_percentile_threshold),h)
    h = CRC32c.crc32c(string(:peak_percentile_threshold=>peak_percentile_threshold), h)
    h = CRC32c.crc32c(string(:peak_threshold=>peak_percentile_threshold),h)
    h
end

getfields(x::AbstractResponseFieldsSimple;kwargs...) = x.fieldindices

get_mesh(::Type{SpatialResponseFieldsSimple},nrefinements::NamedTuple) = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
maptype(::Type{SpatialResponseFieldsSimple}) = SpatialMapNew

get_mesh(::Type{GazeResponseFieldsSimple},nrefinements::NamedTuple) = get_maze_mesh(;nrefinements=nrefinements.g)
maptype(::Type{GazeResponseFieldsSimple}) = ViewMapNew


function get_colors(colormap::Symbol)
    if in([:rain , :navia])(colormap)
        ccolors = [:red, :gold, :orangered, :orange, :salmon, :coral3, :goldenrod1, :firebrick, :tan1, :sienna]
    elseif in([:jet])(colormap)
        ccolors = [:darkorchid4, :sandybrown, :antiquewhite2, :wheat, :lightsalmon]
    else
        ccolors = Makie.wong_colors()
    end 
    return ccolors
end

DPHT.filename(::Type{SpatialResponseFields}) = "spatial_response_fields.jld2"
get_mesh(::Type{SpatialResponseFields},nrefinements::NamedTuple) = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
get_mesh(::Type{SpatialResponseFields},nrefinements::Integer) = Shadow("xy")(floor_topology3(;nrefinements=nrefinements))
maptype(::Type{SpatialResponseFields}) = SpatialMapNew
get_sic_type(::Type{SpatialResponseFields}) = SpatialInformationContent

struct GazeResponseFields <: AbstractResponseFields
    λ::Vector{Float64}
    binidx::Vector{Int64}
    gamma_params::Matrix{Float64} # gaamma parameter fit for the null distribution; two parameters per bin
    args::Dict{Symbol,Any}
end

GazeResponseFieldsAll = Union{GazeResponseFields, GazeResponseFieldsSimple}

DPHT.filename(::Type{GazeResponseFields}) = "gaze_response_fields.jld2"
get_mesh(::Type{GazeResponseFields},nrefinements::NamedTuple) = get_maze_mesh(;nrefinements=nrefinements.g)
maptype(::Type{GazeResponseFields}) = ViewMapNew
get_sic_type(::Type{GazeResponseFields}) = GazeInformationContent

function process_kwargs(::Type{<:AbstractResponseFields},h::UInt32=zero(UInt32);nshuffles=10_000, nrefinements=(p=3,g=2),trial_start=2,smooth=false, smoothing_method=:gaussian, σ=3, α=1000.0^2,pv_threshold=0.05, use_trials=:all, kwargs...)
    h = process_kwargs(JointMap,h;kwargs...)
    if nshuffles != 10_000
        h = CRC32c.crc32c(string(:nshuffles=>nshuffles),h)
    end
    if nrefinements != (p=3,g=2)
        h = CRC32c.crc32c(string(:nrefinements=>nrefinements),h)
    end
    if trial_start != 2
        h = CRC32c.crc32c(string(:trial_start=>trial_start),h)
    end
    if pv_threshold != 0.05
        h = CRC32c.crc32c(string(:pv_threshold=>pv_threshold),h)
    end
    if smooth
        h = CRC32c.crc32c(string(:smoothing_method=>smoothing_method),h)
        if smoothing_method==:gaussian
            h = CRC32c.crc32c(string(:σ=>σ),h)
        elseif smoothing_method==:laplace
            α = get(kwargs, :α, 0.1)
            niter = get(kwargs, :niter, 100)
            h = CRC32c.crc32c(string(:α=>α),h)
            h = CRC32c.crc32c(string(:niter=>niter),h)
        elseif smoothing_method == :adaptive
            h = CRC32c.crc32c(string(:α=>α),h)
        else
            error("Unkown smoothing method $(smoothing_method)")
        end
    end
    if use_trials != :all
         h = CRC32c.crc32c(string(:use_trials=>use_trials),h)
    end
    h
end

function getfields(rf::Union{SpatialResponseFields,GazeResponseFields};cluster_threshold=0.001)
    clusters = get_num_fields(rf,cluster_threshold)
    [rf.binidx[c] for c in clusters]
end

function get_binindex(rf::AbstractResponseFields)
    rf.binidx
end

function get_binindex(rf::AbstractResponseFields, pv_threshold)
    λ = rf.λ
    exceeds = fill(false, size(λ,1))
    gamma_params = rf.gamma_params
    for (i,_λ) in enumerate(λ)
        _params = gamma_params[:,i]
        if all(isfinite.(_params))
            G = Gamma(_params...)
            threshold = percentile(G,100*(1-pv_threshold))
            exceeds[i] = _λ > threshold
        end
    end
    binidx = findall(exceeds) 
    binidx
end

function get_response_fields(::Type{T}, nshuffles::Integer;nrefinements=(p=3,g=2), trial_start=2, redo=fname->false, do_save=true, smooth=false, prog_offset=0, load_only=false, pv_threshold=0.05, use_fitted=false, kwargs...) where T <: Union{GazeResponseFields, SpatialResponseFields}
    h = process_kwargs(T;nshuffles=nshuffles, nrefinements=nrefinements,trial_start=trial_start,smooth=smooth,pv_threshold=pv_threshold, kwargs...)
    args = Dict(:nshuffles=>nshuffles, :nrefinements=>nrefinements, :trial_start=>trial_start, :smooth=>smooth,:pv_threshold=>pv_threshold,:use_trials=>get(kwargs, :use_trials, :all))
    @assert typeof(args) == fieldtype(T, :args)
    args[:dir] = pwd()
    if smooth
        args[:smoothing_method] = get(kwargs, :smoothing_method, :gaussian)
        smoothing_method = get(kwargs, :smoothing_method, :gaussian)
        if smoothing_method == :gaussian
            args[:σ] = get(kwargs, :σ, 3.0)
        elseif smoothing_method == :laplace
            args[:α] = get(kwargs, :α, 0.1)
            args[:niter] = get(kwargs, :niter, 100)
        elseif smoothing_method == :adaptive
            args[:α]= get(kwargs, :α, 1000.0^2)
        end
    end
    fname = DPHT.filename(T)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    do_compute = true 
    obj = nothing
    if !redo(fname) && isfile(fname)
        obj = load_jld2(T, fname)
        if isa(obj, JLD2.ReconstructedMutable)
            do_compute = true
        else
            do_compute = false
        end
    elseif use_fitted
        # load an object where the only difference is the pv_threshold and check whether
        # this object has a gamma fit for the surrogate distribution
        _fname = DPHT.filename(T)
        files = glob(replace(_fname, ".jld2"=>"*.jld2"))
        if !isempty(files)
            for f in files
                _obj = load_jld2(T,f)
            end

        end
        obj = nothing
        do_compute = false
    elseif load_only
        obj = nothing
        do_compute = false
    end
    if do_compute
        # TODO: Here we can actually check if we an object already computed and surrogates fitted
        sp = Spiketrain()

        rp = cd(DPHT.process_level(level(RippleData))) do
            RippleData()
        end
        sp_r = RandomlyShiftedSpiketrains(;nshifts=nshuffles, use_trials = args[:use_trials], trial_start=args[:trial_start])

        unity_gaze_data = cd(DPHT.process_level("session")) do
            UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=fname->false)
        end
        jocc = cd(DPHT.process_level(level(JointOccupancy))) do
            JointOccupancy(;redo=fname->false, nrefinements=nrefinements,trial_start=trial_start,kwargs...)
        end
        jocc_filtered = JointFilteredOccupancy(jocc, unity_gaze_data;kwargs...)
        sic = zeros(nshuffles)
        mm = get_mesh(T,nrefinements)
        vpvrpb = ViewAndPlaceRepresentationNew(sp,rp,unity_gaze_data;kwargs...)
        use_trials = get(kwargs, :use_trials, :all)
        jmb = JointMap(vpvrpb, jocc, jocc_filtered;use_trials=use_trials)
        spmb = maptype(T)(jmb,mm)
        if smooth
            if args[:smoothing_method] == :gaussian
                dmatrix = distancematrix(mm)
                smg = SmoothedMap(spmb;method=:gaussian, σ=args[:σ])
            elseif args[:smoothing_method] == :laplace
                Ls = get_normalize_laplacian(mm)
                @show "Laplace smooth"
                smg = SmoothedMap(spmb;method=:laplace, Ls=Ls, α=args[:α], niter=args[:niter])
            elseif args[:smoothing_method] == :adaptive
                smg = SmoothedMap(spmb;method=:adaptive, α=args[:α])
            else
                error("Unknown smoothing method $(args[:smoothing_method])")
            end
            λ = smg.weight./smg.occupancy
            λ[smg.unvisited] .= NaN
        else
            # get the rate per bin
            λ = spmb.weight./spmb.occupancy
        end
        λ_shuffled = zeros(size(λ,1), nshuffles)
        @showprogress "Computing shuffled firing rates..." offset=prog_offset for (i,sptrain) in enumerate(eachcol(sp_r.timestamps))
            vpvrp = ViewAndPlaceRepresentationNew(sptrain/1000.0,rp,unity_gaze_data;kwargs...)
            jm = JointMap(vpvrp, jocc, jocc_filtered;use_trials=use_trials)
            spm = maptype(T)(jm,mm)
            if smooth
                if args[:smoothing_method] == :gaussian
                    smg = SmoothedMap(spm;dmatrix=dmatrix, method=:gaussian, σ=args[:σ])
                elseif args[:smoothing_method] == :laplace
                    smg = SmoothedMap(spm;Ls=Ls, method=:laplace, α=args[:α], niter=args[:niter])
                elseif args[:smoothing_method] == :adaptive
                    smg = SmoothedMap(spmb;method=:adaptive, α=args[:α])
                end 
                λ_shuffled[:,i] .= smg.weight./smg.occupancy
                λ_shuffled[smg.unvisited,i] .= NaN
            else
                λ_shuffled[:,i] .= spm.weight./spm.occupancy
            end
        end
        exceeds = fill(false, size(λ,1))
        gamma_params = fill(NaN, 2,length(λ))
        for (i,_λ) in enumerate(λ)
            fidx = findall(isfinite, λ_shuffled[i,:])
            if ~isempty(fidx)
                G = fit(Gamma, λ_shuffled[i,fidx])
                gamma_params[:,i] .= params(G)
                threshold = percentile(λ_shuffled[i,fidx],100*(1-pv_threshold))
                exceeds[i] = _λ > threshold
            end
        end
        binidx = findall(exceeds)
        obj = T(λ,binidx, gamma_params, args)
        if do_save
            save_jld2(obj, fname)
        end
    end
    obj
end


function overlap(qq1::Quadrangle, qq2::Quadrangle)
    s1 = segments(boundary(qq1))
    s2 = segments(boundary(qq2))
    res = false
    for _s1 in s1
        for _s2 in s2
            iq = intersection(_s1,_s2)
            if type(iq) != NotIntersecting
                res = true
                break
            end
        end
        if res
            break
        end
    end
    res
end

function consolidate(qq::Vector{<:Quadrangle})
    n = length(qq)
    can_merge = fill(false, n, n)
    #can_merge[diagind(can_merge)] .= true
    for i in 1:n-1
        for j in i+1:n
            #can_merge[j,i] = overlap(qq[i], qq[j])
            iq = intersection(qq[i], qq[j])
            if type(iq) != NotIntersecting
                can_merge[j,i] = true
            end
        end
    end
    can_merge
end


function consolidate_sym(qq::Vector{<:Quadrangle})
    n = length(qq)
    can_merge = fill(false, n, n)
    #can_merge[diagind(can_merge)] .= true
    for i in 1:n
        for j in 1:n 
            if Meshes.intersects(qq[i], qq[j])
                can_merge[j,i] = true
            end
        end
    end
    can_merge
end

function consolidate(can_merge::Matrix{Bool})
    to_process = [1:size(can_merge,2);]
    clusters = Vector{Int64}[]
    while !isempty(to_process)
        vidx = consolidate(can_merge, first(to_process))
        push!(clusters, unique(vidx))
        to_process = setdiff(to_process, vidx)
    end
    clusters
end

function consolidate(can_merge::Matrix{Bool},idx::Int64,visited::Set{Int64}=Set{Int64}())
    idx in visited && return Int64[]
    push!(visited, idx)
    vidx = findall(can_merge[:,idx])
    for _idx in vidx
        qidx = consolidate(can_merge, _idx, visited)
        append!(vidx, qidx)
    end
    push!(vidx, idx)
    return vidx 
end

function set_refinements(::Type{SpatialResponseFields},nr)
    (p=nr,g=2)
end

function set_refinements(::Type{GazeResponseFields},nr)
    (p=3,g=nr)
end

sic_type(::Type{GazeResponseFields}) = GazeInformationContent
sic_type(::Type{SpatialResponseFields}) = SpatialInformationContent

"""
    find_number_of_response_fields(celldir::String;kwargs...)

Estimate the number of distinct response fields by looking across spatial scales.
The size of a response field is defined by the scale at which that response field disappears.
"""
function find_number_of_response_fields(::Type{T}, celldir::String;nrefinements=[0,1,2,3], nshuffles=10_000, kwargs...) where T <: AbstractResponseFields
    #look for intersecting boundary
    qq = Any[]
    cd(celldir) do
        for nr in sort(nrefinements,rev=true)
            _nrefinements = set_refinements(T, nr)
            mm = get_mesh(T,_nrefinements)
            rf = get_response_fields(T,nshuffles;nrefinements=_nrefinements,kwargs...) 
            sic = compute_skaggs_sic(sic_type(T), nshuffles;nrefinements=_nrefinements,kwargs...)
            #only include this level if it is significant overall
            if issignificant(sic)
                append!(qq, mm[rf.binidx])
            end
        end
    end
    qq = convert(Vector{typeof(qq[1])}, qq)
    qqa = merge_fields(qq) 
end

function merge_fields(qqf::Vector{<:Quadrangle})
    can_merge = consolidate_sym(qqf)
    #clusters = Hippocampus.consolidate(can_merge)
    g = SimpleGraph(can_merge)
    clusters = connected_components(g)
    aa = Any[]
    for cluster in clusters
        if length(cluster) == 1
            push!(aa, boundary(qqf[cluster[1]]))
        else
            _aa = boundary(qqf[cluster[1]])
            for k in cluster[2:end]
                _aa = merge(_aa,boundary(qqf[k]))
            end
            push!(aa, _aa)
        end
    end
    aa
end

function merge_fields(rf::T) where T <: AbstractResponseFields
    mm = get_mesh(T, rf.args[:nrefinements])
    merge_fields(mm, rf.binidx)
end

function merge_fields(rf::T) where T <: SpatialResponseFields
    mm = Shadow("xy")(Hippocampus.floor_topology3(;nrefinements=rf.args[:nrefinements].p))
    merge_fields(mm, rf.binidx)
end

function get_major_axis(rf::SpatialResponseFields)
    mm = Shadow("xy")(Hippocampus.floor_topology3(;nrefinements=rf.args[:nrefinements].p))
    clusters = merge_fields(rf)
    nclusters = get_num_fields(rf)
    cidx = findall(dropdims(mean(nclusters,dims=2),dims=2) .< 0.001)
    v = Matrix{Float64}(undef, 2, length(cidx))
    μ = Matrix{Float64}(undef, 2, length(cidx))
    ms = zeros(length(cidx)) # field size
    for (i,c) in enumerate(cidx)
        mq = centroid.(mm[rf.binidx[clusters[c]]])
        points = Tuple.(mq)
        XX = cat([[x,y] for (x,y) in points]...,dims=2)
        ss = svd(XX .- mean(XX, dims=2))
        # or just use the direction of maximum deviation?
        # major axis
        v[:,i] = ss.U[:,1]
        μ[:,i] = mean(Point2f.(Tuple.(mq)))
        ms[i] = ustrip(sum(measure.(mm[rf.binidx[clusters[c]]])))
    end
    v,μ,ms
end

"""
    get_num_fields(rf::T) where T <: AbstractResponseFields

Get the null distribution for the number of random fields with size at least 
as big as those in the actual `rf`.
"""
function get_num_fields(rf::T) where T <: AbstractResponseFields
    mm = get_mesh(T, rf.args[:nrefinements])
    nm = nelements(mm)
    clusters = merge_fields(rf)
    nc = length.(clusters)
    n_sig = sum(nc)
    # cluster analysis
    nclusters = fill(0, length(clusters), 10_000)
    for i in 1:10_000
        bidx = sort(shuffle(1:nm)[1:n_sig])
        fclusters = merge_fields(mm, bidx)
        for j in 1:length(clusters)
            nclusters[j,i] = sum(length.(fclusters) .>= nc[j])
        end
    end
    nclusters
end

function get_num_fields(rf::T,α::Real) where T <: AbstractResponseFields
    nclusters = get_num_fields(rf)
    clusters = merge_fields(rf)
    clusters[dropdims(mean(nclusters,dims=2),dims=2) .< α]
end

function merge_fields(mm::SimpleMesh, idx::Vector{<:Integer})
    # merge fields that are within diagonal distance
    # euclidean distance between centroids
    cm = centroid.(mm[idx])
    D = ustrip.(norm.(cm .- permutedims(cm)))
    ss = ustrip.(sqrt.(measure.(mm[idx])))
    # compare to the element size
    Sq = sqrt.(2*(0.5*ss .+ 0.5*permutedims(ss)).^2)
    # those elements that are within diagonal distance are connected
    Aq = D .<= Sq
    G = SimpleGraph(Aq)
    connected_components(G)
end

function merge_fields(xbins::AbstractVector{<:Real}, ybins::AbstractVector{<:Real}, idx::Vector{<:Integer})
    xy = [xbins ybins]
    D = pairwise(Euclidean(), permutedims(xy))
    ss = mean(sqrt.(diff(xbins).^2 + diff(ybins).^2))
    Aq = D .<= ss 
    G = SimpleGraph(Aq)
    connected_components(G)
end

function merge_fields(A::AbstractMatrix{<:Number}, idx::Vector{<:Integer})
    Aq = A[idx,idx]
    G = SimpleGraph(Aq)
    connected_components(G)
end

function run_cluster_analysis()
end

function find_boundaries(rf::T) where T <: AbstractResponseFields
    mm = get_mesh(T, rf.args[:nrefinements])
     clusters = merge_fields(rf)
     nclusters = get_num_fields(rf)
     cidx = findall(dropdims(mean(nclusters,dims=2),dims=2) .< 0.001)
     boundaries = map(cidx) do _cidx
        find_boundary(mm, rf.binidx[clusters[_cidx]])
     end
     boundaries
end

"""
    find_fields(spm::T; peak_threshold=0.5, baseline_percentile_threshold=10, peak_percentile_threshold=95) where T <: AbstractMap

Identify contiguous peaks where the rate exceeds `peak_percentile_threshold` of the overall rate, and extend these
peaks until the activity reaches `peak_threshold*(λ - b0) + b0` where b0 is `baseline_percentile_threshold` of the overall
firing rate.
"""
function find_fields(spm::T; peak_threshold=0.5, baseline_percentile_threshold=10, peak_percentile_threshold=95,kwargs...) where T <: AbstractMap
    mm = spm.mm
    A = adjacencymatrix(mm)
    # TODO: Make this a bit more data dependent
    # Hm, maybe compute SIC after removing peaks and see when the SIC is no longer significant?
    # what does significance mean? Maybe not feasible to shuffle after every peak?
    # Maybe just use the stats we computed from the full analysis?
    # find the larget peak
    λ = get_rate_map(spm)
    fidx = findall(isfinite, λ)
    b0 = percentile(λ[fidx], baseline_percentile_threshold)
    b1 = percentile(λ[fidx], peak_percentile_threshold)
    fields = Vector{Int64}[]
    while true
        mx,midx = findmax(λ[fidx])
        if mx < b1
            break
        end
        use_threshold = peak_threshold*(λ[fidx[midx]]-b0)
        if use_threshold <= 0
            break
        end
        use_threshold += b0
        g = grow_region(midx, λ[fidx], A[fidx,fidx], use_threshold)  
        push!(fields, fidx[g])
        fidx = setdiff(fidx, fields[end])
        if isempty(fidx)
            break
        end
    end
    fields
end

function get_response_fields(::Type{T},args...;redo=fname->false, do_save=true, kwargs...) where T <: AbstractResponseFieldsSimple
    fname = DPHT.filename(T)
    h = process_kwargs(T;kwargs...)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    do_compute = false
    if !isfile(fname) || redo(fname)
        do_compute=true
    end
    if do_compute
        nrefinements = get(kwargs, :nrefinements, (p=3,g=2))
        mm = get_mesh(T, nrefinements)
        jm = JointMap(;kwargs...)
        spm = maptype(T)(jm,mm)
        method = get(kwargs, :smoothing_method, :laplace)
        spml = SmoothedMap(spm;method=method, kwargs...) 
        ff = find_fields(spml;kwargs...)
        dd = Dict{Symbol,Any}(kwargs)
        baseline_percentile_threshold = pop!(dd, :baseline_percentile_threshold, 10)
        peak_percentile_threshold = pop!(dd, :peak_percentile_threshold, 95)
        peak_threshold = pop!(dd, :peak_threshold, 0.5)
        dd[:dir] = pwd()
        obj = T(get_rate_map(spml), ff, baseline_percentile_threshold, peak_percentile_threshold, peak_threshold,dd)
        if do_save
            save_jld2(obj, fname)
        end
    else
        # load it
        obj = load_jld2(T, fname)
    end
    obj
end

function find_view_field_intersections!(qq::Dict{Symbol,Any};kwargs...)
    nshuffles = get(kwargs, :nshuffles, 10_000)
    rf_gaze = get_response_fields(GazeResponseFields, nshuffles;kwargs...)
    find_view_field_intersections!(qq, rf_gaze;kwargs...)
end

function find_view_field_intersections!(qq::Dict{Symbol,Any}, rf_gaze::GazeResponseFieldsAll;kwargs...)
    # TODO: This is hardcoded and should be made to dependend on the particular session We
    #       are looking at
    poster_pillar_walls = [1, 3, 6,9, 11, 16]
    nrefinments = get(kwargs, :nrefinements, (p=3,g=2))
    mm = get_maze_mesh(;nrefinements=nrefinments.g)
    grouped_bins = group_bins(mm)
    vidx = reduce(vcat, getfields(rf_gaze;cluster_threshold=0.001))
    qq_pillars = fill(false, 4)
    fv = in(vidx)
    if !(:pillars in keys(qq))
        qq[:pillars] = fill(0, 4)
    end
    for i in 1:4
        idx0 = (i-1)*4+1
        idx1 = i*4
        qq[:pillars][i] += any(fv.(reduce(vcat, grouped_bins.pillar_idx[idx0:idx1])))
    end
    for (jj,ii) in enumerate(poster_pillar_walls)
        kk = Symbol("poster$jj")
        if !(kk in keys(qq))
            qq[kk] = 0
        end
        qq[kk] += any(fv.(grouped_bins.pillar_idx[ii]))
    end
    for (k1,k2) in zip([:west_wall_idx, :east_wall_idx, :north_wall_idx, :south_wall_idx, :ceiling_idx, :floor_idx],[:west_wall, :east_wall, :north_wall, :south_wall, :ceiling, :floor])
        if !(k2 in keys(qq))
            qq[k2] = 0
        end
        qq[k2] += any(fv.(get(grouped_bins, k1, 0)))
    end
    qq
end

function find_view_field_intersections(celldirs::Vector{String};kwargs...)
    qq = Dict{Symbol,Any}()
    for celldir in celldirs
        cd(celldir) do
            find_view_field_intersections!(qq;kwargs...)
        end
    end
    qq
end

## plots

function plot_n_fields(::Type{T}, celldirs::Vector{String};figsize=(900,500), kwargs...) where T <: AbstractResponseFields
    with_theme(plot_theme) do
        fig = Figure(size=figsize)
        lg = GridLayout(fig[1,1])
        plot_n_fields!(lg, T, celldirs;kwargs...)
        fig
    end
end

function plot_n_fields!(lg, ::Type{T}, celldirs::Vector{String};labels=["A","B","C","D"], redo=false, kwargs...) where T <: AbstractResponseFields
    h = process_kwargs(T;kwargs...)
    h = CRC32c.crc32c(string(celldirs),h)
    cluster_threshold = get(kwargs, :cluster_threshold, 0.001)
    if cluster_threshold != 0.01
        h = CRC32c.crc32c(string(:cluster_threshold => cluster_threshold),h)
    end
    hs = string(h, base=16)
    fname = joinpath(@__DIR__, "..","data","field_stats_$(hs).jld2")
    @show fname
    if !redo && isfile(fname)
        Z,nfields,field_size, field_boundaries, peak_firing_rate, args = JLD2.load(fname, "Z", "nfields","field_size", "field_boundaries","peak_firing_rate", "args")
    else
        rfs = process_dirs(celldirs) do
            get_response_fields(T, 10_000;load_only=true, kwargs...)
        end
        kk = filter(k->k[2]!==nothing, rfs)
        mm = get_mesh(T,rfs[first(keys(kk))].args[:nrefinements])
        nfields = Dict()
        field_size = Dict()
        # histogram of field sizes
        Z = zeros(nelements(mm))
        bb = Any[]
        peak_firing_rate = Dict()
        for (k,v) in rfs
            if v !== nothing
                peak_firing_rate[k] = maximum(filter(isfinite, v.λ))
                clusters = merge_fields(v)
                if length(clusters) > 0
                    nclusters = get_num_fields(v)
                    threshold = dropdims(sum(nclusters.>0,dims=2),dims=2)
                    cluster_idx = findall(threshold .<= 0.01)
                    nfields[k] = length(cluster_idx)
                    field_size[k] = fill(0.0, length(cluster_idx))
                    for (jj,ii) in enumerate(cluster_idx)
                        cluster = clusters[ii]
                        push!(bb, find_boundary(mm, v.binidx[cluster]))
                        field_size[k][jj] = ustrip(sum(measure.(mm[cluster])))
                        Z[v.binidx[cluster]] .+= 1.0
                    end
                else
                    nfields[k] = 0
                end
            end
        end
        Z[Z.==0.0] .= NaN
        args = rfs[first(keys(kk))].args 
        JLD2.save(fname, Dict("Z"=>Z, "nfields"=>nfields, "field_size"=>field_size, "args"=>args, "field_boundaries"=>bb,"peak_firing_rate"=>peak_firing_rate))
    end
    # normalize to density
    Z ./= sum(filter(isfinite, Z))
    mm = get_mesh(T,args[:nrefinements])
    cc = countmap(values(nfields))
    non_place_cells = collect(keys(filter(k->k[2]==0, nfields)))
    place_cells = collect(keys(filter(k->k[2]>0, nfields)))

    kk = sort(collect(keys(cc)))
    field_sizes = Float64[]
    for (k,v) in field_size
        append!(field_sizes, v)
    end
    colormap = get(kwargs, :colormap, :rain)
    with_theme(plot_theme) do
        lg1 = GridLayout(lg[1,1])
        ax = Axis(lg1[1,1])
        Label(lg1[1,1,TopLeft()], labels[1], padding=(0, 0, 10, 0))
        bcolors = fill(parse(Colorant, :gray), length(kk))
        bcolors[findall(k->k==0, kk)] .= parse(Colorant, :black)
        barplot!(ax, kk, [cc[k] for k in kk],color=bcolors)
        ax.xlabel = "No of fields"
        ax.ylabel = "Count"
        ax2 = Axis(lg1[2,1],xticks=WilkinsonTicks(3))
        Label(lg1[2,1, TopLeft()], labels[2], padding=(0,0,10,0))
        hist!(ax2, field_sizes,color=:gray)
        ax2.xlabel = "Field size [unit^2]"
        lg3 = GridLayout(lg[1,2])
        axf = Axis(lg3[1,1])
        Label(lg[1,2, TopLeft()], labels[3])
        yy_pp = [Float64(peak_firing_rate[k]) for k in place_cells]
        kde_pp = betakde(yy_pp;lower=0.0, upper=maximum(yy_pp))
        yy_np = [Float64(peak_firing_rate[k]) for k in non_place_cells]
        kde_np = betakde(yy_np;lower=0.0, upper=maximum(yy_np))
        xx_pp = range(minimum(yy_pp), stop=1.1*maximum(yy_pp), length=30)
        Δx = mean(diff(xx_pp))
        xx_np = range(minimum(yy_np), stop=1.1*maximum(yy_np),step=Δx) 

        @show KruskalWallisTest(yy_pp, yy_np)
        @show median(yy_pp), median(yy_np)
        hh_pp = normalize(fit(Histogram, yy_pp, xx_pp),mode=:pdf)
        hh_np = normalize(fit(Histogram, yy_np, xx_pp);mode=:pdf)
        
        #TODO: Instead of scatter, use bars, with the gray bars thicker than the black
        #      Alter draw order to avoid bars blocking
        hist_bins = [hh_pp.edges[1][1:end-1];hh_np.edges[1][1:end-1]]
        hist_weight = [hh_pp.weights;hh_np.weights]
        hist_color = [fill(:gray, length(hh_pp.weights));fill(:black, length(hh_np.weights))]
        hist_dodge = [fill(1, length(hh_pp.weights));fill(2, length(hh_np.weights))]
        barplot!(axf, hist_bins, hist_weight, dodge=hist_dodge, color=hist_color, direction=:x)
        #hist!(axf, yy_pp, bins=xx_pp, normalization=:pdf, color=:gray, direction=:x)
        #scatter!(axf, hh_pp.weights, hh_pp.edges[1][1:end-1], color=:gray)
        #barplot!(axf, hh_pp.edges[1][1:end-1], hh_pp.weights, color=:gray, direction=:x)
        #lines!(axf, kde_pp.density, kde_pp.x, color=:gray)
        #hist!(axf, yy_np, bins=xx_np, normalization=:pdf, color=:black, direction=:x)
        #barplot!(axf,hh_np.edges[1][1:end-1], hh_np.weights, color=:black, direction=:x, width=0.7*Δx)
        #lines!(axf, kde_np.density, kde_np.x, color=:black)
        axf.ylabel = "Peak firing rate [Hz]"
        axf.xticklabelsvisible = false
        axf.xticksvisible = true 
        axf.bottomspinevisible = true 
        axf.xlabel = "Density"
        lg2 = GridLayout(lg[1,3])
        floor_offset = get(kwargs, :floor_offset, -20)
        if embeddim(mm) == 2
            ax3 = Axis(lg2[1,1],aspect=1)
            hidedecorations!(ax3)
            ax3.topspinevisible = true
            ax3.rightspinevisible = true
            viz!(ax3, mm;color=:lightgray)
            viz!(ax3, mm;color=Z, showsegments=false, colormap=colormap)
        else
            ax3 = LScene(lg2[1,1], show_axis=false)
            hide_ceiling = get(kwargs, :hide_ceiling, true)
            plotmesh!(ax3, mm;color=Z, showsegments=true, indicate_north=false, segmentcolor=:lightgray, floor_offset=floor_offset, ceiling_offset=10, colormap=colormap, hide_ceiling=hide_ceiling)
        end
        # indicate where the pillars are
        plot_pillars!(ax3;floor_offset=floor_offset)
        Label(lg2[1,1,TopLeft()], labels[4])
        Colorbar(lg2[1,2],colorrange=extrema(filter(isfinite, Z)), ticksvisible=true, ticklabelsvisible=false, label="Spatial density", colormap=colormap)
        #rowsize!(lg, 1, Relative(0.4))
        colsize!(lg, 1, Relative(0.2))
        colsize!(lg, 2, Relative(0.2))
        [ax,ax2, ax3]
    end
end

function plot_response_fields(rf::GazeResponseFields,args...;_plot_theme=plot_theme, kwargs...)
    with_theme(_plot_theme) do
        fig = Figure(size=(579,639))
        lg = GridLayout(fig[1,1])
        lscene = plot_response_fields!(lg, rf, args...;kwargs...)
        fig
    end
end

function plot_response_fields(rf::SpatialResponseFields,args...;_plot_theme=plot_theme, kwargs...)
    with_theme(_plot_theme) do
        fig = Figure(size=(600,500))
        lg = GridLayout(fig[1,1])
        lscene = plot_response_fields!(lg, rf, args...;kwargs...)
        fig
    end
end

function plot_response_fields!(lg::GridLayout, rf::T, λ=rf.λ;filter_spurious=true, show_colorbar=true, show_points=true, label="Firing rate [Hz]",colorbar_below=false,  show_boundaries=false, segmentsize=2.0,showsegments=true, kwargs...) where T <: GazeResponseFieldsAll
    mm = get_mesh(T, get(rf.args,:nrefinements, (p=3,g=2)))
    m_floor, m_ceiling, m_middle = get_floor_and_ceiling(mm)
    lscene = LScene(lg[1,1],show_axis=false)
    floor_offset = get(kwargs, :floor_offset,-10)
    ceiling_offset = get(kwargs, :ceiling_offset, 10)
    colormap = get(kwargs, :colormap, :binary)
    mazecolor = get(kwargs, :mazecolor, :lightgray)
    colorrange = get(kwargs, :colorrange, extrema(filter(isfinite, λ)))
    pointsize = get(kwargs, :pointsize, 10)
    if mazecolor !== nothing
        plotmesh!(lscene, mm;color=mazecolor, showsegments=false, ceiling_offset=ceiling_offset, floor_offset=floor_offset, colormap=colormap,colorrange=colorrange, kwargs...)
    end
    plotmesh!(lscene, mm;color=λ, showsegments=showsegments, ceiling_offset=ceiling_offset, floor_offset=floor_offset, colormap=colormap,colorrange=colorrange, kwargs...)
    plot_pillars!(lscene;floor_offset=floor_offset)
    clusters = getfields(rf;cluster_threshold=0.001)

    ccolors = get_colors(colormap)
    if show_points || show_boundaries
        for (cc,cluster) in zip(ccolors[1:length(clusters)],clusters)
            pidx = cluster
            cpoints = centroid.(mm[pidx])
            floor_points = filter(Meshes.intersects(m_floor), cpoints)
            ceil_points = filter(Meshes.intersects(m_ceiling), cpoints)
            mid_points = setdiff(cpoints, union(floor_points, ceil_points))
            if !isempty(floor_points)
                if show_points
                    viz!(lscene, Translate(0.0, 0.0, floor_offset)(floor_points),color=cc, pointsize=pointsize)
                else
                    midx = findall(Meshes.intersects(m_floor), cpoints)
                    bb = find_boundary(mm[pidx[midx]])
                    viz!(lscene, Translate(0.0, 0.0, floor_offset)(bb),color=cc, pointsize=pointsize, segmentsize=segmentsize)
                end
            end
            if !isempty(ceil_points) && !hide_ceiling
                if show_points
                    viz!(lscene, Translate(0.0, 0.0, ceiling_offset)(ceil_points),color=cc, pointsize=pointsize)
                else
                    midx = findall(Meshes.intersects(m_ceil), cpoints)
                    bb = find_boundary(mm[pidx[midx]])
                    viz!(lscene, Translate(0.0, 0.0, ceiling_offset)(bb),color=cc, pointsize=pointsize, segmentsize=segmentsize)
                end
            end
            if !isempty(mid_points)
                if show_points
                viz!(lscene, mid_points,color=cc, pointsize=pointsize)
                else
                    midx = findall(Meshes.intersects(m_middle), cpoints)
                    bb = find_boundary(mm[pidx[midx]])
                    viz!(lscene, bb, color=cc, segmentsize=segmentsize)
                end
            end
        end
    end
    if show_colorbar
        ticks = WilkinsonTicks(3)
        if colorbar_below
            Colorbar(lg[2,1], colorrange=colorrange, colormap=get(kwargs, :colormap, :binary), label=label, vertical=false, flipaxis=false, ticksvisible=true,ticks=ticks)
        else
            Colorbar(lg[1,2], colorrange=colorrange, colormap=get(kwargs, :colormap, :binary), label=label)
        end
    end
    lscene
end

function plot_response_fields!(lscene::LScene, rf::GazeResponseFields,idx::Union{Nothing,Integer}=nothing, λ=rf.λ ;ceiling_offset=10, floor_offset=-20)
    mm = get_mesh(GazeResponseFields, rf.args[:nrefinements])
    m_floor, m_ceiling, m_middle = get_floor_and_ceiling(mm)
    clusters = merge_fields(rf)

    ccolors = Makie.wong_colors()
    for (jj,(cc,cluster)) in enumerate(zip(ccolors,clusters))
        if idx === nothing || jj == idx
            pidx = rf.binidx[cluster]
            cpoints = centroid.(mm[pidx])
            floor_points = filter(Meshes.intersects(m_floor), cpoints)
            ceil_points = filter(Meshes.intersects(m_ceiling), cpoints)
            mid_points = setdiff(cpoints, union(floor_points, ceil_points))
            if !isempty(floor_points)
                viz!(lscene, Translate(0.0, 0.0, floor_offset)(floor_points),color=cc)
            end
            if !isempty(ceil_points)
                viz!(lscene, Translate(0.0, 0.0, ceiling_offset)(ceil_points),color=cc)
            end
            if !isempty(mid_points)
                viz!(lscene, mid_points,color=cc)
            end
        end
    end
end

function plot_response_fields!(lg::GridLayout, rf::T, λ::AbstractVector{<:Real}=rf.λ;filter_spurious=true, show_points=true, show_boundaries=false, colorbar_below=false, show_colorbar=true, label="Firing rate [Hz]", kwargs...) where T <: SpatialResponseFieldsAll
    mm = Shadow("xy")(floor_topology3(;nrefinements=get(rf.args,:nrefinements, (p=3,g=2)).p))
    ax = Axis(lg[1,1],aspect=1)
    hidedecorations!(ax)
    ax.bottomspinevisible = false
    ax.leftspinevisible = false
    colormap = get(kwargs, :colormap, :binary)
    mazecolor = get(kwargs, :mazecolor, :lightgray)
    colorrange = get(kwargs, :colorrange, extrema(filter(isfinite, λ)))
    viz!(ax, mm;color=mazecolor)
    viz!(ax, mm;color=λ,colormap=colormap, colorrange=colorrange, showsegments=get(kwargs, :showsegments, false),segmentcolor=get(kwargs, :segmentcolor, :black))
    plot_pillars!(ax)
    clusters = getfields(rf;cluster_threshold=0.001)
    ccolors = get_colors(colormap)
    if show_points
        for (cc,cluster) in zip(ccolors[1:length(clusters)],clusters)
            pidx = cluster
            cpoints = centroid.(mm[pidx])
            viz!(ax, cpoints, color=cc)
        end
    elseif show_boundaries
        bb = find_boundaries(rf)
        for (cc,_bb) in zip(ccolors, bb)
            viz!(ax, _bb, color=cc)
        end
    end
    if show_colorbar
        ticks = WilkinsonTicks(3)
        if colorbar_below
            Colorbar(lg[2,1], colorrange=colorrange, colormap=colormap, label=label, vertical=false, flipaxis=false, ticksvisible=true,ticks=ticks)
        else
            Colorbar(lg[1,2], colorrange=colorrange, colormap=colormap, label=label,ticks=ticks)
        end
    end
    ax
end

function plot_response_fields!(lscene::LScene, rf::SpatialResponseFields,idx::Union{Nothing,Integer}=nothing,;offset=0.0)
    mm = Translate(0.0, 0.0, offset)(floor_topology3(;nrefinements=rf.args[:nrefinements].p))
    clusters = merge_fields(rf)

    ccolors = Makie.wong_colors()
    for (jj,(cc,cluster)) in enumerate(zip(ccolors,clusters))
        if idx === nothing || (jj == idx)
            pidx = rf.binidx[cluster]
            cpoints = centroid.(mm[pidx])
            viz!(lscene, cpoints, color=cc)
        end
    end
end

function plot_response_fields!(lg, rf_place::SpatialResponseFields, rf_gaze::GazeResponseFields;kwargs...)
    lgs = GridLayout(lg[1,1])
    plot_response_fields!(lgs, rf_place;kwargs...)
    lgg = GridLayout(lg[1,2])
    plot_response_fields!(lgg, rf_gaze;kwargs...)
end