abstract type AbstractInformationContent end

struct SpatialInformationContent <: AbstractInformationContent
    sic::Vector{Float64}
    sic0::Float64
    args::Dict{Symbol,Any}
end

DPHT.filename(::Type{SpatialInformationContent}) = "spatial_information_content.jld2"

function issignificant(sic::AbstractInformationContent;α=0.05,kwargs...)
    sic.sic0 > percentile(sic.sic, 100*(1-α))
end

function issignificant(::Type{T};α=0.05,kwargs...) where T <: AbstractInformationContent
    sic = compute_skaggs_sic(T, get(kwargs, :nshuffles,10_000);kwargs...)
    issignificant(sic;α=α)
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
    res = fill(false, length(celldirs))
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
                res[c] = func(args...;prog_offset=1, kwargs...)
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

function process_kwargs(::Type{<:AbstractInformationContent};nshuffles=10_000, nrefinements=(p=3,g=2),trial_start=2,smooth=false, smoothing_method=:gaussian, σ=3, α=1000.0^2,kwargs...)
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

"""
Compute Skagg's SIC for the cell in the current working directory using `nshuffles`
random circular shuffles of the underlying spike train data
"""
function compute_skaggs_sic(::Type{T}, nshuffles::Integer;nrefinements=(p=3,g=2), trial_start=2, redo=fname->false, do_save=true, smooth=false, prog_offset=0, kwargs...) where T <: AbstractInformationContent
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
        sicobj = load_jld2(T, fname)
    else
        sp = Spiketrain()
        sp_r = RandomlyShiftedSpiketrains(sp;nshifts=nshuffles, kwargs...)

        jm = JointMap(;kwargs...)
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
        jmb_fname = DPHT.filename(JointMap;kwargs...)
        save_jld2(jmb,jmb_fname)
        spmb = maptype(T)(jmb,mm)
        if smooth
            dmatrix = distancematrix(mm)
            smg = SmoothedMap(spmb;kwargs...)
            sic0 = compute_skaggs_sic(smg)
        else
            sic0 = compute_skaggs_sic(spmb)
        end
        @showprogress "Computing sic...." offset=prog_offset for (i,sptrain) in enumerate(eachcol(sp_r.timestamps))
            vpvrp = ViewAndPlaceRepresentationNew(sptrain/1000.0,rp,unity_gaze_data;kwargs...)
            jm = JointMap(vpvrp, jocc, jocc_filtered)
            spm = maptype(T)(jm,mm)
            if smooth
                smg = SmoothedMap(spm;dmatrix=dmatrix, kwargs...)
                sic[i] = compute_skaggs_sic(smg)
            else
                sic[i] = compute_skaggs_sic(spm)
            end
        end
        sicobj = T(sic,sic0, args)
        if do_save
            save_jld2(sicobj, fname)
        end
    end
    sicobj
end
