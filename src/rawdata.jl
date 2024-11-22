using Makie

struct NeuralData
    analogtime::AbstractVector{UInt64}
    data::Matrix{Float64}
    sampling_rate::Float64
end

function NeuralData(fname::String)
    pth,ext = splitext(fname)
    if ext == ".mda"
        data = MDAFiles.load(fname)
        return NeuralData(1:length(data), data, 30_000.0)
    end
    error("Unknown data format $(ext)")
end


function plot_data(data::NeuralData;kvs...)
    fig = Figure()
    ax = Axis(fig[1,1])
    plot_data!(ax, data;kvs...)
    fig
end

function plot_data!(ax, data::NeuralData;window=1.0)
    #TODO: Check if this is registered first
    deregister_interaction!(ax, :scrollzoom)

    tmin,tmax = extrema(data.analogtime)
    tmin /= data.sampling_rate
    tmax /= data.sampling_rate
    current_time = Observable(data.analogtime[1]/data.sampling_rate)

    on(events(ax).scroll) do (dx,dy)
        #@show tmin,current_time[]+dx, tmax-window
        if tmin < current_time[]+dx < tmax - window
            current_time[] = current_time[] + dx
        end
    end
    
    _data = lift(current_time) do ct
        _ct = round(UInt64, ct*data.sampling_rate)
        idx0 = searchsortedfirst(data.analogtime, _ct)
        idx1 = searchsortedlast(data.analogtime,(ct+window)*data.sampling_rate) 
        [Point2f(tmin + idx/data.sampling_rate,data.data[idx]) for idx in idx0:idx1]
    end

    on(_data) do _d
        autolimits!(ax)
        tightlimits!(ax)
    end

    lines!(ax, _data)
    current_time
end

struct RippleData
    triggers::Matrix{Int64}
    timestamps::Matrix{Float64}
    header::Dict
end

numtrials(rpdata::RippleData) = size(rpdata.triggers,1)
DPHT.level(::Type{RippleData}) = "session"

DPHT.filename(::Type{RippleData}) = "rplparallel.mat"

function RippleData(;do_save=true, redo=false, kvs...)
    outfile = DPHT.filename(RippleData)
    if isfile(outfile) && !redo
        qdata = MAT.matread(outfile)
        # hack
        if "df" in keys(qdata)
            trial_markers = round.(Int64, qdata["df"]["data"]["markers"])
            trial_timestamps = qdata["df"]["data"]["timeStamps"]
            header = Dict()
        else
            trial_markers, trial_timestamps, header = (qdata["triggers"], qdata["timestamps"], get(qdata, "header", Dict()))
            meta = qdata["metadata"]
        end
        rp = RippleData(trial_markers, trial_timestamps, header)
    else
        # extract the markers from the raw file
        nev_files = glob("*.nev")
        if isempty(nev_files)
            error("No nev files found")
        end
        fname = first(nev_files)
        markers,timestamps = RippleTools.extract_markers(fname)
        rp = RippleData(markers, timestamps)
        if do_save
            meta = Dict{String,Any}()
            tag!(meta;storepatch=true)
            qdata = Dict("triggers"=>rp.triggers, "timestamps"=>rp.timestamps, "header"=>rp.header, "metadata"=>meta)
            MAT.matwrite(outfile, qdata)
        end
    end
    rp
end

function RippleData(markers::Vector{T}, timestamps::Vector{Float64}) where T <: Real
    idx = markers.>0
    trial_markers, trial_timestamps = reshape_triggers(Int64.(markers[idx]), timestamps[idx])
    RippleData(trial_markers,trial_timestamps,Dict())
end


