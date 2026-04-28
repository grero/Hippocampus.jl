 using Random
 function random_shift(sptimes::AbstractVector{T}, args...;kwargs...) where T <: Real
    spnew = similar(sptimes)
    random_shift!(spnew, sptimes, args...;kwargs...)
    spnew
 end

 function random_shift!(new_sptimes::AbstractVector{T}, sptimes::AbstractVector{T},min_shift::Real,max_shift::Real;tmin=minimum(sptimes), tmax=maximum(sptimes)) where T <: Real
    Δt = (max_shift-min_shift)*rand() + min_shift
    shift_spiketimes!(new_sptimes, sptimes, Δt;tmin=tmin, tmax=tmax)
 end

 function shift_spiketimes(sptimes::AbstractVector{T},Δt::T;tmin=minimum(sptimes), tmax=maximum(sptimes)) where T <: Real
    new_sptimes = similar(sptimes)
    shift_spiketimes!(new_sptimes, sptimes, Δt;tmin=tmin, tmax=tmax)
 end

 function shift_spiketimes!(new_sptimes::AbstractVector{T}, sptimes::AbstractVector{T},Δt::T;tmin=minimum(sptimes), tmax=maximum(sptimes)) where T <: Real
    new_sptimes .= sptimes .+ Δt
    fidx = new_sptimes .< tmin
    new_sptimes[fidx] .= tmax .- (tmin .- new_sptimes[fidx])
    bidx = new_sptimes .> tmax
    new_sptimes[bidx] .= tmin .+ (new_sptimes[bidx] .- tmax)
    sort!(new_sptimes)
    new_sptimes
 end

struct Spiketrain
    timestamps::Vector{Float64}
    reference::String
    components::Int64
end

function Spiketrain(fname::String)
    q = MAT.matread(fname)
    reference = get(q, "reference", "None")
    components = get(q, "components", 0)
    Spiketrain(q["timestamps"][:], reference, components)
end

function Spiketrain()
    if isfile("spiketrain.mat")
        return Spiketrain("spiketrain.mat")
    elseif isfile("spiketrain.csv")
        timestamps = parse.(Float64, readlines(open("spiketrain.csv")))
        return Spiketrain(timestamps, "", 0)
    else
        return Spiketrain(Float64[], "",0)
    end
end

function Base.show(io::IO, x::Spiketrain)
    nspikes = length(x.timestamps)
    print(io, "Spiketrain with $(nspikes) spikes")
end

struct TrialAlignedSpiketrain
    spiketimes::Vector{Vector{Float64}}
    trigger_timestamps::Matrix{Float64}
    alignto::Int64
end

TrialAlignedSpiketrain(spiketimes, trigger_timestamps) = TrialAlignedSpiketrain(spiketimes, trigger_timestamps,1)

numtrials(sp::TrialAlignedSpiketrain) = length(sp.spiketimes)

function TrialAlignedSpiketrain(sp::Spiketrain, rp::RippleData;kwargs...)
    sptimes =  sp.timestamps/1000.0 
    TrialAlignedSpiketrain(sptimes, rp;kwargs...)
end

function TrialAlignedSpiketrain(sptimes::AbstractVector{T}, rp::RippleData;trial_start=1,alignto=1,Δt0=-0.5,Δt1=0.5) where T <: Real
    nt = numtrials(rp)
    spiketimes = Vector{Vector{Float64}}(undef, nt)
    for i in 1:nt
        timestamps = rp.timestamps[i,:]
        idx0 = searchsortedfirst(sptimes, timestamps[trial_start]+Δt0)
        idx1 = searchsortedlast(sptimes, timestamps[3]+Δt1)
        sp_trial = sptimes[idx0:idx1] .- timestamps[alignto]
        spiketimes[i] = sp_trial
    end
    TrialAlignedSpiketrain(spiketimes, rp.timestamps,alignto)
end

function TrialAlignedSpiketrain(;kwargs...)
    rp = cd(DPHT.process_level(RippleData)) do 
        RippleData()
    end
    sp = Spiketrain()        
    TrialAlignedSpiketrain(sp, rp;kwargs...)
end

function compute_psth(sp::TrialAlignedSpiketrain, binsize::Float64,w=1;tmax=Inf, kwargs...)
    nt = numtrials(sp)
    # find the bins
    nspikes = length.(sp.spiketimes)
    fidx = nspikes .> 0
    tmin = minimum(minimum.(sp.spiketimes[fidx]))
    tmax = min(tmax, maximum(maximum.(sp.spiketimes[fidx])))
    bins = range(tmin, stop=tmax, step=binsize)
    compute_psth(sp, bins,w)
end

function compute_psth(sp::TrialAlignedSpiketrain, bins::AbstractVector{Float64},w=1;kwargs...)
    nt = numtrials(sp)
    hh = StatsBase.fit(Histogram, sp.spiketimes[1], bins;kwargs...)
    weight = zeros(length(hh.weights), nt)
    weight[:,1] .= hh.weights
    for i in 2:nt
        _hh = StatsBase.fit(Histogram, sp.spiketimes[i],bins;kwargs...)
        weight[:,i] .= _hh.weights
    end
    if w > 1
        for j in 1:size(weight,1)-w
            weight[j,:] .= dropdims(sum(weight[j:j+w-1,:],dims=1),dims=1)
        end
        weight = weight[1:size(weight,1)-w,:]
        bins = bins[1:end-w]
    end
    weight, bins
end

function compute_psth(sp::Vector{TrialAlignedSpiketrain}, binsize::Float64,w=1;tmax=Inf, kwargs...)
    # find the bins
    _tmin = Inf
    _tmax = -Inf
    nt = 0
    for _sp in sp
        nt = max(nt, numtrials(_sp))
        nspikes = length.(_sp.spiketimes)
        fidx = nspikes .> 0
        _tmin = min(_tmin, minimum(minimum.(_sp.spiketimes[fidx])))
        _tmax = max(_tmax, maximum(maximum.(_sp.spiketimes[fidx])))
    end
    tmax = min(tmax, _tmax)
    tmin = _tmin
    bins = range(tmin, stop=tmax, step=binsize)
    weight = zeros(length(bins)-w-1, nt, length(sp))
    for (ii,_sp) in enumerate(sp)
        weight[:,:,ii],_ = compute_psth(_sp,bins,w)
    end
    weight, bins
end

struct RandomlyShiftedSpiketrains
    Δt::Vector{Float64} # shifts
    timestamps::Matrix{Float64}
end

function RandomlyShiftedSpiketrains(sp::Spiketrain;nshifts::Integer=10_000, shift_min=0.1, shift_max=0.9, kwargs...)
    n = length(sp.timestamps)
    tmin,tmax = extrema(sp.timestamps)
    dur = tmax-tmin

    sptrains = zeros(n, nshifts) 
    Δt = zeros(nshifts)
    Δs = shift_max - shift_min
    for i in 1:nshifts
        Δt[i] = Δs*dur*rand() + shift_min*dur
        shift_spiketimes!(view(sptrains, :, i), sp.timestamps, Δt[i];tmax=tmax)
    end
    RandomlyShiftedSpiketrains(Δt, sptrains)
end