using Meshes
abstract type AbstractResponseFields end


struct SpatialResponseFields <: AbstractResponseFields
    binidx::Vector{Int64}
    gamma_params::Matrix{Float64} # gaamma parameter fit for the null distribution; two parameters per bin
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

function process_kwargs(::Type{<:AbstractResponseFields},h::UInt32=zero(UInt32);nshuffles=10_000, nrefinements=(p=3,g=2),trial_start=2,smooth=false, smoothing_method=:gaussian, σ=3, α=1000.0^2,pv_threshold=0.05, kwargs...)
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
        elseif smoothing_method == :adpative
            h = CRC32c.crc32c(string(:α=>α),h)
        else
            error("Unkown smoothing method $(smoothing_method)")
        end
    end
    h
end

function get_response_fields(::Type{T}, nshuffles::Integer;nrefinements=(p=3,g=2), trial_start=2, redo=fname->false, do_save=true, smooth=false, prog_offset=0, load_only=false, pv_threshold=0.05, kwargs...) where T <: AbstractResponseFields
    h = process_kwargs(T;nshuffles=nshuffles, nrefinements=nrefinements,trial_start=trial_start,smooth=smooth,pv_threshold=pv_threshold, kwargs...)
    args = Dict(:nshuffles=>nshuffles, :nrefinements=>nrefinements, :trial_start=>trial_start, :smooth=>smooth,:pv_threshold=>pv_threshold)
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
        if isa(obj, JLD2.ReconstructedMutable)
            obj = T(obj.binidx, fill(NaN, 2, nshuffles), obj.args) 
        end
    elseif load_only
        return nothing
    else
        # TODO: Here we can actually check if we an object already computed and surrogates fitted
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
        obj = T(binidx, gamma_params, args)
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
    min_nc = minimum(nc)
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

function plot_n_fields(::Type{T}, celldirs::Vector{String};kwargs...) where T <: AbstractResponseFields
    with_theme(plot_theme) do
        fig = Figure(size=(600,600))
        lg = GridLayout(fig[1,1])
        plot_n_fields!(lg, T, celldirs;kwargs...)
        fig
    end
end

function plot_n_fields!(lg, ::Type{T}, celldirs::Vector{String};kwargs...) where T <: AbstractResponseFields
    rfs = process_dirs(celldirs) do
        get_response_fields(T, 10_000;load_only=true, kwargs...)
    end
    kk = filter(k->k[2]!==nothing, rfs)
    mm = get_mesh(T,rfs[first(keys(kk))].args[:nrefinements])
    nfields = Dict()
    field_size = Dict()
    # histogram of field sizes
    Z = zeros(nelements(mm))
    for (k,v) in rfs
        if v !== nothing
            clusters = merge_fields(v)
            nfields[k] = length(clusters)
            field_size[k] = fill(0.0, length(clusters))
            for (ii,cluster) in enumerate(clusters)
                field_size[k][ii] = ustrip(sum(measure.(mm[cluster])))
                Z[v.binidx[cluster]] .+= 1.0
            end
        end
    end
    Z[Z.==0.0] .= NaN
    cc = countmap(values(nfields))
    kk = sort(collect(keys(cc)))
    field_sizes = Float64[]
    for (k,v) in field_size
        append!(field_sizes, v)
    end
    with_theme(plot_theme) do
        ax = Axis(lg[1,1])
        barplot!(ax, kk, [cc[k] for k in kk])
        ax.xlabel = "No of fields"
        ax.ylabel = "Count"
        ax2 = Axis(lg[1,2])
        hist!(ax2, field_sizes)
        ax2.xlabel = "Field size [unit^2]"
        lg2 = GridLayout(lg[2,1:2])
        if embeddim(mm) == 2
            ax3 = Axis(lg2[1,1],aspect=1)
            hidedecorations!(ax3)
            ax3.topspinevisible = true
            ax3.rightspinevisible = true
            viz!(ax3, mm;color=:lightgray)
            viz!(ax3, mm;color=Z, showsegments=true, segmentcolor=:lightgray)
        else
            ax3 = LScene(lg2[1,1], show_axis=false)
            plotmesh!(ax3, mm;color=Z, showsegments=true, segmentcolor=:lightgray, floor_offset=-20, ceiling_offset=10)
        end
        Colorbar(lg2[1,2],colorrange=extrema(filter(isfinite, Z)), ticksvisible=true, label="Count")
        rowsize!(lg, 1, Relative(0.4))
        [ax,ax2, ax3]
    end
end