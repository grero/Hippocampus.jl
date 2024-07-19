using GLMakie
using ContinuousWavelets


function get_time_slice(x::AbstractVector{T}, idx) where T
    x[idx]
end

function get_observable(x::AbstractVector{T},idx) where T
    Observable(get_time_slice(x,idx))
end

function get_observable(X::Matrix{T},idx) where T
    Observable(get_time_slice(X,idx))
end

function get_time_slice(X::Matrix{T},idx) where T
    X[idx,:]
end

function plot_heatmap(x::Vector{T1}, X::Matrix{T2},M::Integer;t=1:size(X,1),freqs=1:size(X,2),freq_bands::Union{Nothing, Vector{Tuple{Float64,Float64}}}=nothing) where T1 where T2
    fig = Figure()
    axes = [Axis(fig[i,1]) for i in 1:3]
    rowsize!(fig.layout, 1, Relative(0.1))
    linkxaxes!(axes[2:end]...)
    n = size(X,1)
    idx = 1
    indices = idx:idx+M-1
    obs = [get_observable(_x, indices) for _x in [x,X]]
    tt = get_observable(t, indices)
    tend = Observable([t[idx+M-1]])
    on(events(fig.scene).keyboardbutton) do event
        if event.action == Keyboard.press
            do_update = false
            if event.key == Keyboard.left
            if idx > M
                    idx -= M
                    do_update = true
            end
            elseif event.key == Keyboard.right
                if idx < n-2*M
                    idx += M
                    do_update = true
                end
            end
            if do_update
                #img[] = X[idx:idx+M-1,:]
                for (_x,_ob) in zip([x,X], obs)
                    _ob[] = get_time_slice(_x, idx:idx+M-1) 
                end
                tt[] = t[idx:idx+M-1]
                tend[] = [t[idx+M-1]]
                for ax in axes[2:end]
                    autolimits!(ax)
                end
            end
        end
    end
    for ax in axes[1:end-1]
        ax.xticklabelsvisible = false
    end
    barplot!(axes[1], [1.0], tend, direction=:x)
    xlims!(axes[1], t[1], t[end])
    axes[1].yticklabelsvisible = false
    axes[1].yticksvisible = false
    lines!(axes[2], tt, obs[1])
    heatmap!(axes[3], tt, freqs, obs[2])
    if freq_bands !== nothing
        # highlight these bands
        for fb in freq_bands
            hlines!(axes[3], fb[1],fb[2],color=:white)
        end
    end
    ax = axes[3]
    ax.xlabel = "Time"
    ax.ylabel = "Frequency"
    fig
end

function get_spectrum(x::Vector{T},fs::T2;β=2) where T <: Real where T2 <: Real
    n = length(x)
    c = wavelet(Morlet(π), averagingType=NoAve(), β=β)
    daughers, ω = computeWavelets(n,c)
    freqs = getMeanFreq(daughers,fs)
    res = cwt(x, c, daughters)
    res, freqs
end

function get_reward_aligned_spectrum(;window=(-500,500), β=1.5, maxfreq=100.0)
    data = MAT.matread("vmlfp.mat")
    fs = data["vp"]["data"]["analogInfo"]["SampleRate"]
    dd = data["vp"]["data"]["analogData"][:]
    # findall rewarded trials
    ridx = findall(30 .< data["vp"]["data"]["markers"][:,3] .< 40)
    # reward onsets in units of ms
    reward_onsets = round.(Int64, ceil.(data["vp"]["data"]["timeStamps"][ridx,3]*1000))
   
    #this is probably massively inefficient

    # compute wavelet-bases power-spectrum
    n = window[2]-window[1]
    c = wavelet(Morlet(π), averagingType=NoAve(), β=β)
    daughters, ω = computeWavelets(n,c)
    freqs = getMeanFreq(daughters,fs) 

    # Average reward aligned power spectrum
    P = fill(0.0, 1000, length(freqs))
    Ps = fill!(similar(P), 0.0)
    P2 = fill(0.0, 1000, length(freqs))
    for onset in reward_onsets
        x = dd[onset+window[1]:onset+window[2]-1]
        res = cwt(x, c, daughters) 
        res_s = cwt(circshift(x, rand(-n+1:n-1)),c,daughters)
        ll = log.(abs.(res))
        ll_s = log.(abs.(res_s))
        P .+= ll
        Ps .+= ll_s
        P2 .+= ll.*ll
    end
    P ./= length(reward_onsets)
    Ps ./= length(reward_onsets)
    P2 ./= length(reward_onsets)
    P2 .-= P.*P
    P2 .= sqrt.(P2)
    P, Ps, range(window[1], stop=window[2], length=n), freqs
end

function get_trial_data(unity_data::Matrix{T}, lfp_data::Vector{T2}, unity_triggers::Matrix{T}, lfp_triggers::Matrix{T2};fs=1000.0, β=1.5) where T <: Real where T2 <: Real
    ntrials = size(lfp_triggers,1)
    nu = size(unity_data,1)
    nl = length(lfp_data)
    pos = Vector{Matrix{Float64}}(undef, ntrials)
    lfp = Vector{Vector{Float64}}(undef, ntrials)
    spec = Vector{Matrix{ComplexF64}}(undef, ntrials)
    t = Vector{Vector{Float64}}(undef, ntrials)
    freqs = Vector{Vector{Float64}}(undef, ntrials)
    for tidx in 1:ntrials
        uidx0 = max(round(Int64, unity_triggers[tidx, 1]), -1, 1)
        uidx1 = min(round(Int64, unity_triggers[tidx, 3])+1, nu)
        pos[tidx] = unity_data[uidx0:uidx1,3:4] 

        idx0 = max(round(Int64, lfp_triggers[tidx, 1]*1000)-1000,1)
        idx1 = min(round(Int64, lfp_triggers[tidx, 3]*1000)+1000,nl)
        lfp[tidx] = lfp_data[idx0:idx1]
        t[tidx] = range(idx0/1000.0, stop=idx1/1000, length=idx1-idx0+1)
        # compute power spectrum using wavelets
        spec[tidx],freqs[tidx] = get_spectrum(lfp[tidx],fs;β=β)
    end
    lfp,spec,t,pos,freqs
end

function plot_trial(unity_data::Matrix{T}, lfp_data::Vector{T2}, unity_triggers::Matrix{T}, lfp_triggers::Matrix{T2}) where T <: Real where T2 <: Real
    ntrials = size(lfp_triggers,1)
    nu = size(unity_data,1)
    nl = length(lfp_data)
    xmin,xmax = extrema(unity_data[:,3])
    ymin,ymax = extrema(unity_data[:,4])
    tidx = Observable(1)
    pos = Observable(Point2f[])
    lfp = Observable(Point2f[])
    trial_start_lfp = Observable(0)
    trial_start_unity = Observable(Point2f[])
    on(tidx) do _tidx
        uidx0 = max(round(Int64, unity_triggers[tidx[], 1]), -1, 1)
        uidx1 = min(round(Int64, unity_triggers[tidx[], 3])+1, nu)
        pos[] = [Point2f(unity_data[i,3:4]...) for i in uidx0:uidx1]
        trial_start_unity[] = [Point2f(unity_data[uidx0+1,3], unity_data[uidx0+1,4])]

        idx0 = max(round(Int64, lfp_triggers[tidx[], 1]*1000)-1000,1)
        idx1 = min(round(Int64, lfp_triggers[tidx[], 3]*1000)+1000,nl)
        lfp[] = [Point2f.(i, lfp_data[i]) for i in idx0:idx1]
        trial_start_lfp[] = idx0+1000
    end

    fig = Figure()

    on(events(fig.scene).keyboardbutton) do event
        if event.action == Keyboard.press
            do_upate = false
            if event.key == Keyboard.left
                if tidx[] > 0 
                    tidx[] = tidx[] - 1 
                    do_update = true
                end
            elseif event.key == Keyboard.right
                if tidx[] < ntrials 
                    tidx[] = tidx[] + 1 
                    do_update = true
                end
            end
            if do_update
                autolimits!(axes[2])
            end
        end
    end

    axes = [Axis(fig[i,1]) for i in 1:2]
    lines!(axes[1], pos)
    scatter!(axes[1], trial_start_unity, color=:red)
    xlims!(axes[1], xmin, xmax)
    ylims!(axes[1], ymin, ymax)

    lines!(axes[2], lfp)
    vlines!(axes[2],trial_start_lfp, color=:red)
    tidx[] = 1
    fig
end

function plot_trial(pos::Vector{Matrix{T}}, lfp::Vector{Vector{T}}, spec::Vector{Matrix{ComplexF64}}, freqs::Vector{Vector{T}},t::Vector{Vector{T}};
                                                        arena_size::Union{Nothing, NTuple{4,Float64}}=nothing) where T <: Real
    ntrials = length(pos)
    tidx = Observable(1)
    p_pos = Observable([Point2f(pos[tidx[]][i,:]...) for i in 1:size(pos[tidx[]],1)])
    p_lfp = Observable([Point2f(_t,_lfp) for (_t,_lfp) in zip(t[tidx[]], lfp[tidx[]])])
    p_spec = Observable(abs.(spec[tidx[]]))
    p_freqs = Observable(freqs[tidx[]])
    p_t = Observable(t[tidx[]])
    trial_start_unity = Observable(Point2f[])

    fig = Figure()
    axes = [Axis(fig[i,1]) for i in 1:3]
    linkxaxes!(axes[2], axes[3])
    if arena_size !== nothing
        xlims!(axes[1], arena_size[1:2]...)
        ylims!(axes[1], arena_size[3:4]...)
    end

    on(tidx) do _tidx
        # hack
        _t = t[_tidx]
        if length(_t) > 16383
            qidx = round.(Int64,range(1, stop=length(_t), length=16383))
        else
            qidx = 1:length(_t)
        end
        p_pos[] =[Point2f(pos[_tidx][i,:]...) for i in 1:size(pos[_tidx],1)] 
        trial_start_unity[] = p_pos[][1:1]
        p_lfp[] = [Point2f(_t,_lfp) for (_t,_lfp) in zip(t[tidx[]][qidx], lfp[tidx[]][qidx])]
        p_t[] = t[_tidx][qidx]
        p_freqs[] = freqs[_tidx]
        p_spec[] = abs.(spec[_tidx][qidx,:])
    end

    on(events(fig.scene).keyboardbutton) do event
        if event.action == Keyboard.press
            do_update = false
            if event.key == Keyboard.left
                if tidx[] > 1 
                    tidx[] = tidx[] - 1 
                    do_update = true
                end
            elseif event.key == Keyboard.right
                if tidx[] < ntrials 
                    tidx[] = tidx[] + 1 
                    do_update = true
                end
            end
            if do_update
                autolimits!(axes[2])
                autolimits!(axes[3])
            end
        end
    end

    lines!(axes[1], p_pos)
    scatter!(axes[1], trial_start_unity, color=:red)
    heatmap!(axes[2], p_t, p_freqs, p_spec)
    axes[2].xticklabelsvisible = false
    lines!(axes[3], p_lfp)

    tidx[] = 1
    fig
end