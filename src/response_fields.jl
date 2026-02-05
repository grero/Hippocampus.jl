using Meshes
abstract type AbstractResponseFields end


struct SpatialResponseFields <: AbstractResponseFields
    binidx::Vector{Int64}
    args::Dict{Symbol,Any}
end

DPHT.filename(::Type{SpatialResponseFields}) = "spatial_response_fields.jld2"
get_mesh(::Type{SpatialResponseFields},nrefinements::NamedTuple) = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
maptype(::Type{SpatialResponseFields}) = SpatialMapNew

struct GazeResponseFields <: AbstractResponseFields
    binidx::Vector{Int64}
    args::Dict{Symbol,Any}
end

DPHT.filename(::Type{GazeResponseFields}) = "gaze_response_fields.jld2"
get_mesh(::Type{GazeResponseFields},nrefinements::NamedTuple) = get_maze_mesh(;nrefinements=nrefinements.g)
maptype(::Type{GazeResponseFields}) = ViewMapNew

function process_kwargs(::Type{<:AbstractResponseFields};nshuffles=10_000, nrefinements=(p=3,g=2),trial_start=2,smooth=false, smoothing_method=:gaussian, σ=3, α=1000.0^2,kwargs...)
    h = zero(UInt32)
    if nshuffles != 10_000
        h = CRC32c.crc32c(string(:nshuffles=>nshuffles),h)
    end
    if nrefinements != (p=3,g=2)
        h = CRC32c.crc32c(string(:nrefinements=>nrefinements),h)
    end
    if trial_start != 2
        h = CRC32c.crc32c(string(:trial_start=>trial_start),h)
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
        elseif smoothing_method == :adpative
            h = CRC32c.crc32c(string(:α=>α),h)
        else
            error("Unkown smoothing method $(smoothing_method)")
        end
    end
    h
end

function get_response_fields(::Type{T}, nshuffles::Integer;nrefinements=(p=3,g=2), trial_start=2, redo=fname->false, do_save=true, smooth=false, prog_offset=0, load_only=false, kwargs...) where T <: AbstractResponseFields
    h = process_kwargs(T;nshuffles=nshuffles, nrefinements=nrefinements,trial_start=trial_start,smooth=smooth,kwargs...)
    args = Dict(:nshuffles=>nshuffles, :nrefinements=>nrefinements, :trial_start=>trial_start, :smooth=>smooth)
    @assert typeof(args) == fieldtype(T, :args)
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
    if !redo(fname) && isfile(fname)
        obj = load_jld2(T, fname)
    elseif load_only
        return nothing
    else
        sp = Spiketrain()
        sp_r = RandomlyShiftedSpiketrains(sp;nshifts=nshuffles, kwargs...)

        rp = cd(DPHT.process_level(level(RippleData))) do
            RippleData()
        end
        unity_gaze_data = cd(DPHT.process_level(level(UnityRaytraceData))) do
            UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=false)
        end
        jocc = cd(DPHT.process_level(level(JointOccupancy))) do
            JointOccupancy(;redo=false, nrefinements=nrefinements,trial_start=trial_start)
        end
        jocc_filtered = JointFilteredOccupancy(jocc, unity_gaze_data;kwargs...)
        sic = zeros(nshuffles)
        mm = get_mesh(T,nrefinements)
        vpvrpb = ViewAndPlaceRepresentationNew(sp,rp,unity_gaze_data;kwargs...)
        jmb = JointMap(vpvrpb, jocc, jocc_filtered)
        spmb = maptype(T)(jmb,mm)
        if smooth
            if args[:smoothing_method] == :gaussian
                dmatrix = distancematrix(mm)
                smg = SmoothedMap(spmb;method=:gaussian, σ=args[:σ])
            elseif args[:smoothing_method] == :laplace
                Ls = get_normalize_laplacian(mm)
                @show "Laplace smooth"
                smg = SmoothedMap(spmb;method=:laplace, Ls=Ls, α=args[:α], niter=args[:niter])
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
            jm = JointMap(vpvrp, jocc, jocc_filtered)
            spm = maptype(T)(jm,mm)
            if smooth
                if args[:smoothing_method] == :gaussian
                    smg = SmoothedMap(spm;dmatrix=dmatrix, method=:gaussian, σ=args[:σ])
                else args[:smoothing_method] == :laplace
                    smg = SmoothedMap(spm;Ls=Ls, method=:laplace, α=args[:α], niter=args[:niter])
                end 
                λ_shuffled[:,i] .= smg.weight./smg.occupancy
                λ_shuffled[smg.unvisited,i] .= NaN
            else
                λ_shuffled[:,i] .= spm.weight./spm.occupancy
            end
        end
        exceeds = fill(false, size(λ,1))
        for (i,_λ) in enumerate(λ)
            fidx = findall(isfinite, λ_shuffled[i,:])
            if ~isempty(fidx)
                threshold = percentile(λ_shuffled[i,fidx],95)
                exceeds[i] = _λ > threshold
            end
        end
        binidx = findall(exceeds)
        obj = T(binidx, args)
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

function run_cluster_analysis()
end