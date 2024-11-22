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

"""
Contains information about the total time spent in each spatial bin.
"""
struct SpatialOccupancy{T<:Real}
    xbins::AbstractVector{T}
    ybins::AbstractVector{T}
    weight::Matrix{T}
end

function SpatialOccupancy(udata::UnityData, xbins::AbstractVector{T}, ybins::AbstractVector{T};trial_start=1) where T <: Real
    nt = numtrials(udata)
    weight = fill(0.0, length(xbins)-1, length(ybins)-1)
    for i in 1:nt
        tu, posx, posy, _ = get_trial(udata, i;trial_start=trial_start)
        for j in 2:length(tu)
            Δt = tu[j]-tu[j-1]
            xidx = searchsortedlast(xbins, posx[j-1])
            yidx = searchsortedlast(ybins, posy[j-1])
            if 0 < xidx <= size(weight,1) && 0 < yidx <= size(weight,2)
                weight[xidx,yidx] += Δt
            end
        end
    end
    SpatialOccupancy(xbins, ybins, weight)
end

function SpatialOccupancy(xbins,ybins=xbins;kwargs...)
    udata = UnitData()
    SpatialOccupancy(udata, xbins, ybins;kwargs...)
end

struct SpatialMap{T<:Real}
    xbins::AbstractVector{T}
    ybins::AbstractVector{T}
    weight::Matrix{T}
    occupancy::Matrix{T}
end

DPHT.level(::Type{SpatialMap}) = "cell"

function SpatialMap(spr::SpatialRepresentation, xbins::AbstractVector{T}, ybins::AbstractVector{T},spoc::SpatialOccupancy;kwargs...) where T <: Real
    spatial_count = fill(0.0, length(xbins)-1, length(ybins)-1)
    nt = numtrials(spr)
    for i in 1:nt
        position = spr.position[i]
        xpos = [pos[1] for pos in position]
        ypos = [pos[2] for pos in position]
        h = fit(Histogram, (xpos,ypos), (xbins, ybins))
        spatial_count .+= h.weights
    end
    SpatialMap(xbins,ybins, spatial_count, spoc.weight)
end

function SpatialMap(xbins, ybins;redo=false, do_save=false,kwargs...)
    sp = Spiketrain()
    rp = cd(DPHT.process_level(RippleData;kwargs...)) do
        RippleData()
    end
    udata = cd(DPHT.process_level(UnityData;kwargs...)) do
        UnityData()
    end
    spoc = SpatialOccupancy(udata, xbins, ybins)
    spr = SpatialRepresentation(sp, rp, udata;kwargs...);
    SpatialMap(spr, xbins, ybins, spoc)
end

function compute_sic(spm::SpatialMap)

    x = spm.weight./spm.occupancy
    idx = isfinite.(x)
    p = spm.occupancy[idx]./sum(spm.occupancy[idx])
    r = sum(p.*x[idx])
    xr = x[idx]./r
    ll = log2.(xr)
    lidx = isfinite.(ll)
    sic = sum((p.*xr.*ll)[lidx])
end

function compute_entropy(sp::Union{SpatialMap, SpatialOccupancy})
    pp = sp.weight./sum(sp.weight)
    -sum(filter(isfinite, pp.*log2.(pp)))
end

function Makie.convert_arguments(::Type{<:AbstractPlot}, spm::SpatialMap;normalize=true)
    if normalize
        X = spm.weight./spm.occupancy
        label = "Firing rate"
    else
        X = spm.weight
        label = "Spike count"
    end
    h = S.Heatmap(spm.xbins, spm.ybins, X)
    ax1 = S.Axis(plots=[h])
    ll = S.Colorbar(h,label=label)
    S.GridLayout([ax1 ll])
end

