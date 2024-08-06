using MakieCore

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


function MakieCore.convert_arguments(::Type{<:AbstractPlot},data::NeuralData)
    
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