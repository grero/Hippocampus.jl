using Random
using Distributions

"""
Simulate a simple place field neuron using behavioural data in `udata`
"""
function model_place_field(udata::UnityData, rpdata::RippleData;λmin=0.1, λmax=3.0,dt=0.01,σ=1.0,μ=[3.5,1.2],rng=Random.default_rng())
    Σ = [σ 0.0;0.0 σ]
    G = MvNormal(μ, Σ)
    G0 = pdf(G, μ)
    nt = numtrials(udata)

    spikes = Float64[]
    q = -log(rand(rng))
    r = 0.0
    for i in 1:nt
        t,mposx,mposy = get_trial(udata, i;trial_start=2)
        # get the time from ripple
        trp = rpdata.timestamps[i,2]
        # reference to ripple time
        t .= t .- t[1] .+ trp[1]
        Δt = t[2]-t[1]
        for (t0,posx,posy) in zip(t,mposx, mposy)
            λ = pdf(G, [posx,posy]) # firing rate based on place field
            # scale firing rate
            λ = λmax*λ/G0 + λmin
            _t = t0
            while _t < t0+Δt
                r += λ*dt
                if r >= q
                    push!(spikes, _t)
                    r = 0.0
                    q = -log(rand(rng))
                end
                _t += dt
            end
        end
    end
    spikes
end
