using Makie
"""
A spatial representation of events
"""
struct SpatialRepresentation
    position::Vector{Vector{Point2f}}
    events::Vector{Vector{Float64}}
end

function SpatialRepresentation(spikes::Spiketrain, rp::RippleData, udata::UnityData)
    nt = numtrials(udata)
    position = Vector{Vector{Point2f}}(undef, nt)
    events = Vector{Vector{Float64}}(undef, nt)
    sp = spikes.timestamps/1000.0 #convert to seconds
    for i in 1:nt
        tp,posx,posy,_ = get_trial(udata,i)
        tp .-= tp[1]
        timestamps = rp.timestamps[i,:]
        idx0 = searchsortedfirst(sp, timestamps[1])
        idx1 = searchsortedlast(sp, timestamps[3])
        # align to trial start
        sp_trial = sp[idx0:idx1] .- timestamps[1]
        nspikes = idx1-idx0+1
        events[i] = sp_trial
        position[i] = Vector{Point2f}(undef, nspikes)
        for j in 1:nspikes
            k = searchsortedfirst(tp,sp_trial[j])
            if 0 < k <= length(posx)
                position[i][j] = Point2f(posx[k],posy[k])
            end
        end
    end
    SpatialRepresentation(position,events)
end

function SpatialRepresentation()
    sptrain = Spiketrain()
    rdata = cd(DPHT.process_level(level(RippleData))) do
        RippleData()
    end
    udata = cd(DPHT.process_level(level(UnityData))) do
        UnityDate()
    end
    SpatialRepresentation(sptrain, rdata, udata)
end

numtrials(spr::SpatialRepresentation) = length(spr.position)

function Makie.convert_arguments(::Type{<:AbstractPlot}, spr::SpatialRepresentation)
    nt = numtrials(spr)
    points = Point2f[]
    for pos in spr.position
        append!(points, Point2f.(pos))
    end
    ax = S.Axis(plots=[S.Scatter(points)])
    S.GridLayout(ax)
end

function visualize!(lscene, spr::SpatialRepresentation;trial::Observable{Trial}=Observable(Trial(1)), current_time::Observable{Float64}=Observable(0.0))
    nt = numtrials(spr)
    trial_events = lift(trial) do _trial
        if 0 < _trial.i <= nt
            return spr.position[_trial.i]
        else
            return [Point2f(NaN)]
        end
    end
    scatter!(lscene, trial_events)
end