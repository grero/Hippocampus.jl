using Eyelink
using Makie

function zerounless(::Type{Eyelink.Event};kwargs...)
    args = Any[]
    for k in fieldnames(Eyelink.Event)
        if k in keys(kwargs)
            push!(args, kwargs[k])
        else
            tt = fieldtype(Eyelink.Event, k)
            if tt <: Real
                push!(args, zero(tt))
            else
                if tt <: Symbol
                    push!(args, :none)
                elseif tt <: String
                    push!(args, "")
                else
                    push!(args, nothing)
                end
            end
        end
    end
    Eyelink.Event(args...)
end

struct EyelinkData
    triggers::Matrix{Int64}
    timestamps::Matrix{UInt64}
    analogtime::Vector{UInt64}
    gazex::Matrix{Float32}
    gazey::Matrix{Float32}
    fixation_start::Vector{UInt64}
    fixation_end::Vector{UInt64}
    saccade_start_time::Vector{UInt64}
    saccade_end_time::Vector{UInt64}
    saccade_start_pos::Matrix{Float32}
    saccade_end_pos::Matrix{Float32}
    header::Dict
end

numtrials(x::EyelinkData) = size(x.triggers,1)

function Base.show(io::IO, x::EyelinkData)
    n = length(x.analogtime)
    m = length(x.saccade_start_time)
    nt = size(x.triggers,1)
    print(io, "Eyelink data with $n samples and $m saccade events over $(nt) trials")
end

DPHT.filename(::Type{EyelinkData}) = "eyelink.mat"
DPHT.level(::Type{EyelinkData}) = "day"

function DPHT.load(::Type{EyelinkData}, fname::String)
    qdata = MAT.matread(fname)
    if "metadata" in keys(qdata)
        meta = qdata["metadata"]
    end
    EyelinkData(qdata)
end

function DPHT.save(edata::EyelinkData)
    fname = DPHT.filename(EyelinkData)
    qdata = convert(Dict{String,Any}, edata)
    MAT.matwrite(fname, qdata)
end

function EyelinkData(qdata::Dict)
    args = Any[]
    for k in fieldnames(EyelinkData)
        push!(args, qdata[string(k)])
    end
    EyelinkData(args...)
end

function EyelinkData(;do_save=true,redo=false)
    outfile = DPHT.filename(EyelinkData)
    if !redo && isfile(outfile)
        return DPHT.load(EyelinkData, outfile)
    end

    # create the object
    edffiles = glob("*.edf")
    if isempty(edffiles)
        error("No EDF files found")
    end
    edata = EyelinkData(first(edffiles))
    if do_save
        DPHT.save(edata)
    end
    edata
end

"""
    get_markers(messages::Vector{Eyelink.Event})

Extract trial markers from Eyelink message events.
"""
function get_markers(messages::Vector{Eyelink.Event})
    triggers = Int64[]
    timestamps = UInt64[]
    for msg in messages
        if startswith(msg.message, "Start Trial") ||
        startswith(msg.message, "End Trial") ||
        startswith(msg.message, "Cue Offset") ||
        startswith(msg.message, "Timeout")
                trigger = parse(Int64, split(msg.message)[end])
                push!(triggers, trigger)
                push!(timestamps, msg.sttime)
        end
    end
    trial_markers, trial_timestamps = reshape_triggers(triggers, timestamps)
end

function EyelinkData(fname::String;do_save=true, redo=false, kvs...)
    eyelinkdata = Eyelink.load(fname)
    header = Dict()
    #get gaze coordinates
    for ee in eyelinkdata.events
        if ee.eventtype == :messageevent
            if startswith(ee.message, "GAZE_COORDS")
                qs = split(ee.message)
                gaze_coords = parse.(Float64, qs[2:end])
                header["gaze_coords"] = gaze_coords
                break
            end
        end
    end
    # figure out which eye we are recording from
    min_gaze_value = minimum(eyelinkdata.samples.gx,dims=2)
    eye_recorded = fill(true,2)
    for k in [1,2]
        if all(min_gaze_value[k] == typemin(Int16))
            eye_recorded[k] = false
        end
    end
    header["eye_recorded"] = findall(eye_recorded)
    screen_height = header["gaze_coords"][4]
    # process fixation events
    fixations = filter(ee->ee.eventtype==:endfix, eyelinkdata.events)
    nfix = length(fixations)
    fixation_start = Vector{UInt64}(undef, nfix)
    fixation_end = Vector{UInt64}(undef, nfix)
    for (i,fixation) in enumerate(fixations)
        fixation_start[i] = fixation.sttime
        fixation_end[i] = fixation.entime
    end
    #process saccade events
    saccades = filter(ee->ee.eventtype==:endsacc, eyelinkdata.events)
    nsacc = length(saccades)
    saccade_start_time = Vector{UInt64}(undef, nsacc)
    saccade_end_time = Vector{UInt64}(undef, nsacc)
    saccade_start_pos = Matrix{Float32}(undef, 2, nsacc)
    saccade_end_pos = Matrix{Float32}(undef, 2, nsacc)

    for (i,saccade) in enumerate(saccades)
        saccade_start_time[i] = saccade.sttime
        saccade_end_time[i] = saccade.entime
        saccade_start_pos[:,i] = [saccade.gstx, screen_height-saccade.gsty]
        saccade_end_pos[:,i] = [saccade.genx, screen_height-saccade.geny]
    end

    # get the messages
    messages = filter(ee->ee.eventtype==:messageevent, eyelinkdata.events)

    trial_markers, trial_timestamps = get_markers(messages)
    qdata = Dict{String,Any}()
    merge!(qdata, Dict("triggers"=>trial_markers, "timestamps"=>trial_timestamps, "analogtime"=>eyelinkdata.samples.time,
            "gazex"=>eyelinkdata.samples.gx,"gazey"=>screen_height .- eyelinkdata.samples.gy, "fixation_start"=>fixation_start,
            "fixation_end"=>fixation_end, "saccade_start_time"=>saccade_start_time,
            "saccade_end_time"=>saccade_end_time, "saccade_start_pos"=>saccade_start_pos, "saccade_end_pos"=>saccade_end_pos))
    qdata["header"] = header
    EyelinkData(qdata)
end

function Base.convert(::Type{Dict{String, Any}}, edata::EyelinkData)
    qdata = Dict{String,Any}()
    qdata["triggers"] = edata.triggers
    qdata["timestamps"] = edata.timestamps
    qdata["analogtime"] = edata.analogtime
    qdata["gazex"] = edata.gazex
    qdata["gazey"] = edata.gazey
    qdata["fixation_start"] = edata.fixation_start
    qdata["fixation_end"] = edata.fixation_end
    qdata["saccade_start_time"] = edata.saccade_start_time
    qdata["saccade_end_time"] = edata.saccade_end_time
    qdata["saccade_start_pos"] = edata.saccade_start_pos
    qdata["saccade_end_pos"] = edata.saccade_end_pos
    qdata["header"] = edata.header
    meta = Dict{String,Any}()
    tag!(meta;storepatch=true)
    qdata["metadata"] = meta
    qdata
end

function get_trial(edata::EyelinkData, i;trial_start=1)
    idx0 = searchsortedfirst(edata.analogtime, edata.timestamps[i,trial_start])
    idx1 = searchsortedfirst(edata.analogtime, edata.timestamps[i,3])
    trial_time = edata.analogtime[idx0:idx1]
    er = edata.header["eye_recorded"]

    fixation_mask = fill(false, idx1-idx0+1)
    for i in 1:length(edata.fixation_start)
        fix_start = edata.fixation_start[i]
        fix_end = edata.fixation_end[i]
        fidx0 = searchsortedfirst(trial_time, fix_start)
        fidx1 = searchsortedfirst(trial_time, fix_end)
        fidx1 = min(fidx1, length(trial_time))
        fixation_mask[fidx0:fidx1] .= true
        if fidx1 >= length(trial_time)
            break
        end
    end
    gx = edata.gazex[er,idx0:idx1]
    gx[gx.==1.0f8] .= NaN32
    gy = edata.gazey[er,idx0:idx1]
    gy[gy.==1.0f8] .= NaN32

    trial_time, gx, gy,fixation_mask
end

function Makie.convert_arguments(::Type{<:AbstractPlot}, x::EyelinkData)
    gx = x.gazex[1,:]
    gy = x.gazey[1,:]
    idx = findall((abs.(gx) .>= 32768).|(abs.(gy) .>= 32768))
    gx[idx] .= NaN
    gy[idx] .= NaN
    [PlotSpec(Lines, x.analogtime, gx),
    PlotSpec(Lines, x.analogtime, gy)]
end

function Makie.convert_arguments(::Type{<:AbstractPlot}, x::EyelinkData, trial::Trial)
    t,gx,gy = get_trial(x, trial.i)
    idx = findall((abs.(gx[1,:]) .>= 32768).|(abs.(gy[1,:]) .>= 32768))
    gx[1,idx] .= NaN
    gy[1,idx] .= NaN
    PlotSpec(Lines, gx[1,:], gy[1,:])
end

function visualize!(ax::Makie.Block, edata::EyelinkData;trial::Observable{Trial}=Observable(Trial(1)),current_time::Observable{Float64}=Observable(0.0),kwargs...)

    edata_trial = lift(trial) do _trial
        te, gx, gy,fixmask = get_trial(edata, _trial.i)
    end
    # set the limits
    x0,y0,x1,y1 = edata.header["gaze_coords"]
    xlims!(ax, x0,x1)
    ylims!(ax, y0,y1)
    ax.xgridvisible = false
    ax.ygridvisible = false
    gaze_pos = Observable([Point2f(NaN)])
    current_j = 1
    onany(edata_trial, current_time) do _edt, ct
        te = _edt[1]
        tef = (te .- te[1])/1000.0
        j = searchsortedfirst(tef, ct)
        if 0 < j <= length(tef)
            j0 = min(j,current_j)
            j1 = max(j, current_j)
            gaze = permutedims([_edt[2][:] _edt[3][:]])
            gaze_pos[] = Point2f.(eachcol(gaze[:,j0:j1]))
            current_j = j
        end
    end
    if isa(ax, Axis)
        scatter!(ax, gaze_pos;color=:red)
    else
        # plot into the near plane
    end
end
