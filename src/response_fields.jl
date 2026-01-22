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
