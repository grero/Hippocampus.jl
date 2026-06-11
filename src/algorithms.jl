abstract type AbstractInformationContent end

struct SpatialInformationContent <: AbstractInformationContent
    sic::Vector{Float64}
    sic0::Float64
    args::Dict{Symbol,Any}
end

DPHT.filename(::Type{SpatialInformationContent}) = "spatial_information_content.jld2"

function issignificant(sic::AbstractInformationContent;pv_threshold=0.05,kwargs...)
    sic_th = percentile(sic.sic, 100*(1-pv_threshold))
    sic.sic0 > sic_th
 end

function issignificant(::Type{T};pv_threshold=0.05,kwargs...) where T <: AbstractInformationContent
    sic = compute_skaggs_sic(T, get(kwargs, :nshuffles,10_000);kwargs...)
    if sic !== nothing
        return issignificant(sic;pv_threshold=pv_threshold)
    end
    return nothing

end

function get_sic(::Type{T}, celldirs::Vector{String};skip_error=false, kwargs...) where T <: AbstractInformationContent
    sic0 = zeros(length(celldirs))
    for (i,celldir) in enumerate(celldirs)
        try
            sic = cd(celldir) do 
                compute_skaggs_sic(T, get(kwargs, :nshuffles,10_000);kwargs...)
            end
            sic0[i] = sic.sic0
        catch ee
            if skip_error
                @show celldir
            else
                rethrow(ee)
            end
        end
    end
    sic0
end

function issignificant(::Type{T},celldirs::Vector{String};skip_error=false, kwargs...) where T <: AbstractInformationContent
    res = Vector{Union{Bool,Nothing}}(undef, length(celldirs))
    for (i,c) in enumerate(celldirs)
        try
            cd(c) do
                res[i] = issignificant(T;kwargs...)
            end
        catch ee
            if skip_error
                @show c
            else
                rethrow(ee)
            end
        end
    end
    res
end

struct JointInformationContent <: AbstractInformationContent
    sic::Vector{Float64}
    sic0::Float64
    args::Dict{Symbol,Any}
end

DPHT.filename(::Type{JointInformationContent}) = "joint_information_content.jld2"

function compute_skaggs_sic(jm::JointMap, nrefinements::@NamedTuple{p::Int64, g::Int64};do_smooth=false, smoothing_method=:laplace, α=0.01, niter=1000)
    # get the joint counts
    m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
    mm = get_maze_mesh(;nrefinements=nrefinements.g)
    ng = nelements(mm)
    np = nelements(m_floor)
    X = zeros(ng,np)
    Y = zeros(ng,np)
    for (idx,w,o) in zip(jm.index, jm.weight, jm.occupancy)
        gidx = getindex(idx,1)
        pidx = getindex(idx,2)
        X[gidx,pidx] += w
        Y[gidx,pidx] += o
    end
   
    if do_smooth
        if smoothing_method == :laplace
            A_floor = adjacencymatrix(m_floor)
            Dp_floor = Diagonal(vec(1.0./sqrt.(sum(A_floor,dims=2))))
            Ls_floor = I - Dp_floor*A_floor*Dp_floor
            A_mm = adjacencymatrix(mm)
            Dp_mm = Diagonal(vec(1.0./sqrt.(sum(A_mm,dims=2))))
            Ls_mm = I - Dp_mm*A_mm*Dp_mm
            X = laplace_smoothing(permutedims(X), Ls_floor, α;niter=niter)
            X = laplace_smoothing(permutedims(X), Ls_mm, α;niter=niter)

            Y = laplace_smoothing(permutedims(Y), Ls_floor, α;niter=niter)
            Y = laplace_smoothing(permutedims(Y), Ls_mm, α;niter=niter)
        end
    end
    compute_skaggs_sic(X[:], Y[:])
end

"""
    process_dirs(func::Function,celldirs::Vector{String}, args...;skip_error=false, kwargs...)

Run the function `func`, accepting positional arguments `args` and keyword arguments `kwargs` in each `celldirs`, collecting
the return values in a dictionary with `celldirs` as the keys.
"""
function process_dirs(func::Function,celldirs::Vector{String}, args...;skip_error=false, kwargs...)
    res = Dict{String,Any}()
    @showprogress "Processing...." for c in celldirs
        try
            cd(c) do
                res[c] = func(args...; kwargs...)
            end
        catch ee
            if skip_error
                @show c
            else
                rethrow(ee)
            end
        end
    end
    res
end

function process_kwargs(::Type{<:AbstractInformationContent};nshuffles=10_000, nrefinements=(p=3,g=2),trial_start=2,smooth=false, smoothing_method=:gaussian, σ=3, α=1000.0^2,rseed::UInt32=UInt32(1234), use_trials=:all, kwargs...)
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
        elseif smoothing_method == :laplace
            h = CRC32c.crc32c(string(:α=>α),h)
            niter = get(kwargs, :niter, 1000)
            h = CRC32c.crc32c(string(:α=>α),h)
            h = CRC32c.crc32c(string(:niter=>niter),h)
        end
    end
    if rseed != 1234
        h = CRC32c.crc32c(string(:rseed=>rseed),h)
    end
    if use_trials != :all
        h = CRC32c.crc32c(string(:use_trials>use_trials),h)
    end
    h
end

struct GazeInformationContent <: AbstractInformationContent
    sic::Vector{Float64}
    sic0::Float64
    args::Dict{Symbol,Any}
end

DPHT.filename(::Type{GazeInformationContent}) = "gaze_information_content.jld2"

maptype(::Type{SpatialInformationContent}) = SpatialMapNew
get_mesh(::Type{SpatialInformationContent},nrefinements::NamedTuple) = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
maptype(::Type{GazeInformationContent}) = ViewMapNew
get_mesh(::Type{GazeInformationContent},nrefinements::NamedTuple) = get_maze_mesh(;nrefinements=nrefinements.g)
maptype(::Type{JointInformationContent}) = JointMap

"""
Compute Skagg's SIC for the cell in the current working directory using `nshuffles`
random circular shuffles of the underlying spike train data
"""
function compute_skaggs_sic(::Type{T}, nshuffles::Integer;nrefinements=(p=3,g=2), load_only=fname->false, redo=fname->false, do_save=true, do_skip=fname->false, smooth=false, prog_offset=0, rseed=UInt32(1234), use_trials=:all, kwargs...) where T <: AbstractInformationContent
    h = process_kwargs(T;nshuffles=nshuffles, nrefinements=nrefinements,smooth=smooth,rseed=rseed, use_trials=use_trials, kwargs...)
    args = Dict(:nshuffles=>nshuffles, :nrefinements=>nrefinements, :trial_start=>get(kwargs, :trial_start, 1), :smooth=>smooth, :rseed=>rseed,:use_trials=>use_trials)
    @assert typeof(args) == fieldtype(T, :args)
    if smooth
        args[:smoothing_method] = get(kwargs, :smoothing_method, :gaussian)
        smoothing_method = get(kwargs, :smoothing_method, :gaussian)
        if smoothing_method == :gaussian
            args[:σ] = get(kwargs, :σ, 3.0)
        elseif smoothing_method == :laplace
            args[:α] = get(kwargs, :α, 0.01)
            args[:niter] = get(kwargs, :niter, 1000)
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
        sicobj = load_jld2(T, fname)
    elseif (load_only(fname) && !isfile(fname)) || do_skip(fname)
        return nothing
    else
        rp = cd(DPHT.process_level(level(RippleData))) do
            RippleData()
        end
        sp = Spiketrain()
        # TODO: shuffle only some spikes
        trialstart = args[:trial_start]
        nt = div(size(rp.timestamps,1), 2)
        if use_trials == :firstHalf
            tmin = 1000.0*rp.timestamps[1,trialstart]
            tmax = 1000.0*rp.timestamps[nt,3]
        elseif use_trials == :secondHalf
            tmin = 1000.0*rp.timestamps[nt+1,trialstart]
            tmax = 1000.0*rp.timestamps[end,3]
        else
            # TODO: This should probably come from rpdata as well
            tmin = 1000.0*rp.timestamps[1,trialstart]
            tmax = 1000.0*rp.timestamps[end,3]
        end
        sp_r = RandomlyShiftedSpiketrains(sp;nshifts=nshuffles, rseed=rseed, tmin=tmin, tmax=tmax, kwargs...)

        jm = JointMap(;use_trials=use_trials, kwargs...)
        unity_gaze_data = cd(DPHT.process_level("session")) do
            UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=redo)
        end
        jocc = cd(DPHT.process_level(level(JointOccupancy))) do
            JointOccupancy(;redo=redo, kwargs...)
        end
        jocc_filtered = JointFilteredOccupancy(jocc, unity_gaze_data;redo=redo,kwargs...)
        vpvrpb = ViewAndPlaceRepresentationNew(sp,rp,unity_gaze_data;redo=redo,kwargs...)
        jmb = JointMap(vpvrpb, jocc, jocc_filtered;use_trials=use_trials)
        jmb_fname = DPHT.filename(JointMap;use_trials=use_trials, kwargs...)
        save_jld2(jmb,jmb_fname)
        sic = zeros(nshuffles)
        if T <: JointInformationContent
            if smooth
                if smoothing_method == :laplace
                     jml = Hippocampus.JointSmoothedMap(jm;method=:laplace, α=args[:α], niter=args[:niter])
                     sic0 = compute_skaggs_sic(jml)
                else
                    error("Only laplace smoothing is implemented for joint maps")
                end
            else
                sic0 = compute_skaggs_sic(jm)
            end
        else
            mm = get_mesh(T,args[:nrefinements])
            spmb = maptype(T)(jmb,mm)
            if smooth
                if smoothing_method == :gaussian
                    dmatrix = distancematrix(mm)
                    smg = SmoothedMap(spmb;method=:gaussian, dmatrix=dmatrix, kwargs...)
                else
                    smoothing_method = :laplace
                    Ls = get_normalize_laplacian(mm)
                    smg = SmoothedMap(spmb;method=:laplace, Ls=Ls, kwargs...)
                end
                sic0 = compute_skaggs_sic(smg)
            else
                sic0 = compute_skaggs_sic(spmb)
            end
        end
        @showprogress "Computing sic...." offset=prog_offset for (i,sptrain) in enumerate(eachcol(sp_r.timestamps))
            vpvrp = ViewAndPlaceRepresentationNew(sptrain/1000.0,rp,unity_gaze_data;kwargs...)
            jm = JointMap(vpvrp, jocc, jocc_filtered;use_trials=use_trials)
            if T <: JointInformationContent
                if smooth
                    if smoothing_method == :laplace
                        jml = Hippocampus.JointSmoothedMap(jm;method=:laplace, α=args[:α], niter=args[:niter])
                        sic[i]= compute_skaggs_sic(jml)
                    end
                else
                    sic[i] = compute_skaggs_sic(jm)
                end
            else
                spm = maptype(T)(jm,mm)
                if smooth
                    if smoothing_method == :gaussian
                        smg = SmoothedMap(spm;method=:gaussian, dmatrix=dmatrix, kwargs...)
                    else
                        smg = SmoothedMap(spm;method=:laplace, Ls=Ls, kwargs...)
                    end
                    sic[i] = compute_skaggs_sic(smg)
                else
                    sic[i] = compute_skaggs_sic(spm)
                end
            end
        end
        sicobj = T(sic,sic0, args)
        if do_save
            save_jld2(sicobj, fname)
        end
    end
    sicobj
end

function shuffle_sic(;nshuffles=1000,α=0.1, niter=100)
    m_floor = Shadow("xy")(floor_topology3(;nrefinements=3))
    Lsf = Hippocampus.get_normalize_laplacian(m_floor)
    mm = get_maze_mesh(;nrefinements=2)
    Lsm = Hippocampus.get_normalize_laplacian(mm)
    jm = JointMap(;nrefinements=(p=3,g=2),redo=fname->false)
    weight,occupancy = get_maps(jm)
    jml = JointSmoothedMap(jm;method=:laplace, α=0.1, niter=100);
    Z = get_rate_map(jml)
    fidx = findall(isfinite, Z)
    sic_true = Hippocampus.compute_skaggs_sic(Z[fidx], jml.occupancy[fidx])
    sic_shuffle = zeros(nshuffles)
    X_shuffle = zeros(size(jml.weight)..., nshuffles)
    @show size(X_shuffle)
    for i in 1:length(sic_shuffle)
        # FIXME: This issue here is that we are also mixing space and view. Ideally, we should shuffle the spike trains
        X_shuffle[:,:,i] .= shuffle(weight)
    end
    # shuffle space first
    Xs = reshape(permutedims(X_shuffle, [3,1,2]), nshuffles*size(X_shuffle,1), size(X_shuffle,2))
    Xp = laplace_smoothing(Xs, Lsf, α;niter=niter)
    @show size(Xp)
    Xg = reshape(permutedims(reshape(Xp, nshuffles, size(X_shuffle,1), size(X_shuffle,2)),[1,3,2]), nshuffles*size(X_shuffle,2), size(X_shuffle,1))
    Xgs = laplace_smoothing(Xg, Lsm,α;niter=niter)
    @show size(Xgs)
    Xgs = permutedims(reshape(Xgs, nshuffles, size(X_shuffle,2), size(X_shuffle,1)), [3,2,1])
    @assert size(Xgs) == size(X_shuffle)
    Xgs ./= reshape(occupancy, size(occupancy)...,1)
    # smooth occupancy
    #Yg = laplace_smoothing(jml.occupancy, Lsf,α;niter=niter)
    #Yg = permutedims(laplace_smoothing(permutedims(Yg, [2,1]), Lsm,α;niter=niter))
    for i in 1:nshuffles
        Xgs[jml.unvisited,i] .= NaN
        midx = findall(isfinite, Xgs[:,:,i])
        sic_shuffle[i] = compute_skaggs_sic(Xgs[midx,i], jml.occupancy[midx])
    end
    sic_true, sic_shuffle, Xgs
end
