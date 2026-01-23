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
        elseif smoothing_method == :adpative
            h = CRC32c.crc32c(string(:α=>α),h)
        end
    end
    h
end

function get_response_fields(::Type{T}, nshuffles::Integer;nrefinements=(p=3,g=2), trial_start=2, redo=fname->false, do_save=true, smooth=false, prog_offset=0, kwargs...) where T <: AbstractResponseFields
    h = process_kwargs(T;nshuffles=nshuffles, nrefinements=nrefinements,trial_start=trial_start,smooth=smooth,kwargs...)
    args = Dict(:nshuffles=>nshuffles, :nrefinements=>nrefinements, :trial_start=>trial_start, :smooth=>smooth)
    @assert typeof(args) == fieldtype(T, :args)
    if smooth
        args[:smoothing_method] = get(kwargs, :smoothing_method, :gaussian)
        smoothing_method = get(kwargs, :smoothing_method, :gaussian)
        if smoothing_method == :gaussian
            args[:σ] = get(kwargs, :σ, 3.0)
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
            dmatrix = distancematrix(mm)
            smg = SmoothedMap(spmb;kwargs...)
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
                smg = SmoothedMap(spm;dmatrix=dmatrix, kwargs...)
                λ_shuffled[:,i] .= spm.weight./spm.occupancy
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
            #can_merge[j,i] = overlap(qq[i], qq[j])
            iq = intersection(qq[i], qq[j])
            if type(iq) != NotIntersecting
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


"""
    find_number_of_response_fields(celldir::String;kwargs...)

Estimate the number of distinct response fields by looking across spatial scales.
The size of a response field is defined by the scale at which that response field disappears.
"""
function find_number_of_response_fields(celldir::String;nrefinements=[0,1,2,3], nshuffles=10_000, kwargs...)
    #look for intersecting boundary
    bb = Any[]
    qq = Any[]
    cd(celldir) do
        for nr in sort(nrefinements,rev=true)
            _nrefinements = (p=nr, g=2)
            mm = get_mesh(SpatialResponseFields,_nrefinements)
            rf = get_response_fields(SpatialResponseFields,nshuffles;nrefinements=_nrefinements,kwargs...) 
            sic = compute_skaggs_sic(SpatialInformationContent, nshuffles;nrefinements=_nrefinements,kwargs...)
            #only include this level if it is significant overall
            if issignificant(sic)
                A = adjacencymatrix(mm)
                idx = rf.binidx
                bbs = boundary.(mm[idx])
                push!(qq, mm[idx])
                bbq = Any[]
                do_include = fill(true, length(bbs)) 
                for i in 1:length(idx)-1
                    for j in i+1:length(idx)
                        if A[i,j] != 0
                            #TODO: Make this unto a single ring
                            #bq = merge(bbs[i],bbs[j])
                            #do_include[[i,j]] .= false 
                            #push!(bbq, bq)
                        end
                    end
                end
                append!(bbq, bbs[do_include])
                # get the positions
                push!(bb, bbq)
            end
        end
    end

    # start from the larget scale and absorb regions that intersect in the lower scales
    # TODO: Also merge fields that are adjacent
    # FIXME: For some reason this doesn't seem to work
    avail = [fill(true, length(_qq)) for _qq in qq]
    for k in 2:length(qq)
        for _qq in qq[k]
            c = centroid(_qq)
            m = sqrt(measure(_qq))/2
            for j in (k-1):-1:1
                lidx = findall(avail[j]) 
                #vidx = Meshes.intersects.(bb[j][lidx], _bb)
                for (l,_qq2) in enumerate(qq[j][lidx])
                    #vidx = any(func.(Meshes.vertices(_bb2)))
                    # contained within
                    v = centroid(_qq2) - c
                    vidx = all(abs.(v) .< m)
                    if vidx
                        avail[j][lidx[l]] = false
                    end
                    # if adjacent?
                end
                #avail[j][lidx[vidx]] .= false
            end
        end
    end
    qqf = [_qq[aa] for (_qq,aa) in zip(qq,avail)]
    # check for adjacent regions and merge them
    cat(qqf...,;dims=1)
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