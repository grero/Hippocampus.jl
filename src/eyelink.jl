using Eyelink
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

function Base.show(io::IO, x::EyelinkData)
    n = length(x.analogtime)
    m = length(x.saccade_start_time)
    nt = size(x.triggers,1)
    print(io, "Eyelink data with $n samples and $m saccade events over $(nt) trials")
end

DPHT.filename(::Type{EyelinkData}) = "eyelink.mat"

function EyelinkData(fname::String;do_save=true, redo=false, kvs...)
    outfile = DPHT.filename(EyelinkData)
    if isfile(outfile) && !redo
        qdata = MAT.matread(outfile)
        meta = qdata["metadata"]
        # TODO: Check code version here
    else
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
        meta = Dict{String,Any}()
        tag!(meta;storepatch=true)
        qdata = Dict{String,Any}()
        merge!(qdata, Dict("triggers"=>trial_markers, "timestamps"=>trial_timestamps, "analogtime"=>eyelinkdata.samples.time,
                "gazex"=>eyelinkdata.samples.gx,"gazey"=>screen_height .- eyelinkdata.samples.gy, "fixation_start"=>fixation_start,
                "fixation_end"=>fixation_end, "saccade_start_time"=>saccade_start_time,
                "saccade_end_time"=>saccade_end_time, "saccade_start_pos"=>saccade_start_pos, "saccade_end_pos"=>saccade_end_pos))
        qdata["metadata"] = meta
        qdata["header"] = header 
        if do_save
            MAT.matwrite(outfile, qdata)
        end
    end
    args = Any[]
    for k in fieldnames(EyelinkData)
        push!(args, qdata[string(k)])
    end
    EyelinkData(args...)
end

function get_trial(edata::EyelinkData, i)
    idx0 = searchsortedfirst(edata.analogtime, edata.timestamps[i,1])
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
