module SpatialAnalyses
using MakieCore
using GLMakie
using Meshes
import Makie.SpecApi as S
using Colors
using ContinuousWavelets
using DataProcessingHierarchyTools
const DPHT = DataProcessingHierarchyTools
using ProgressMeter
using MAT
using DrWatson
using DelimitedFiles
import StatsBase
using RippleTools
using Eyelink
using Distributions

include("unity.jl")

import Base.show

export numtrials

abstract type AbstractSpatialPreferenceMap end

struct SpatialPreferenceMap <: AbstractSpatialPreferenceMap
    xbins::AbstractVector{Float64}
    ybins::AbstractVector{Float64}
    occupancy::Matrix{Float64}
    preference::Matrix{Float64}
end

struct LFPBandSpatialPreferenceMap <: AbstractSpatialPreferenceMap 
    xbins::AbstractVector{Float64}
    ybins::AbstractVector{Float64}
    occupancy::Matrix{Float64}
    preference::Matrix{Float64}
    central_freq::Float64
end

MakieCore.convert_arguments(::Type{<:AbstractPlot}, x::SpatialPreferenceMap) = PlotSpec(Heatmap, x.xbins, x.ybins, x.preference./x.occupancy)

"""
Spectrogram aligned to trial events, computed using Morlet wavelets with the specified β
"""
struct TrialAlignedSpectrogram
    t::Vector{Vector{Float64}}
    spec::Vector{Matrix{ComplexF64}}
    freqs::Vector{Vector{Float64}}
    sampling_rate::Float64
    β::Float64
end

get_sampling_rate(x::TrialAlignedSpectrogram) = x.sampling_rate

function Base.show(io::IO, spec::TrialAlignedSpectrogram)
    print(io, "Morlet-based spectrogram")
end

Base.length(x::TrialAlignedSpectrogram) = length(x.spec) 

#function Base.show(io,  ::MIME"text/plain", spec::TrialAlignedSpectrogram)
 #   print(io, "Morlet-based spectrogram\n")
#end

VectorOrMatrix{T} = Union{Matrix{T}, Vector{T}}
is_int64(x) = round(Int64, x) == x

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

struct Trial
    i::UInt64
end


struct RippleData
    triggers::Matrix{Int64}
    timestamps::Matrix{Float64}
    header::Dict
end

DPHT.filename(::Type{RippleData}) = "rplparallel.mat"

function RippleData(fname::String;do_save=true, redo=false, kvs...)
    outfile = DPHT.filename(RippleData)
    if isfile(outfile) && !redo
        qdata = MAT.matread(outfile)
        trial_markers, trial_timestamps, header = (qdata["triggers"], qdata["timestamps"], get(qdata, "header", Dict()))
        rp = RippleData(trial_markers, trial_timestamps, header)
        meta = qdata["metadata"]
    else
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





struct Spiketrain
    timestamps::Vector{Float64}
    reference::String
    components::Int64
end

function Spiketrain(fname::String)
    q = MAT.matread(fname)
    Spiketrain(q["timestamps"][:], q["reference"], q["components"])
end

function Spiketrain()
    if isfile("spiketrain.mat")
        return Spiketrain("spiketrain.mat")
    else
        return Spiketrain(Float64[], "",0)
    end
end

function Base.show(io::IO, x::Spiketrain)
    nspikes = length(x.timestamps)
    print(io, "Spiketrain with $(nspikes) spikes")
end

function MakieCore.convert_arguments(::Type{<:AbstractPlot}, x::EyelinkData) 
    gx = x.gazex[1,:]
    gy = x.gazey[1,:]
    idx = findall((abs.(gx) .>= 32768).|(abs.(gy) .>= 32768))
    gx[idx] .= NaN
    gy[idx] .= NaN
    [PlotSpec(Lines, x.analogtime, gx),
    PlotSpec(Lines, x.analogtime, gy)]
end

function MakieCore.convert_arguments(::Type{<:AbstractPlot}, x::EyelinkData, trial::Trial)
    t,gx,gy = get_trial(x, trial.i)
    idx = findall((abs.(gx[1,:]) .>= 32768).|(abs.(gy[1,:]) .>= 32768))
    gx[1,idx] .= NaN
    gy[1,idx] .= NaN
    PlotSpec(Lines, gx[1,:], gy[1,:])
end

function MakieCore.convert_arguments(::Type{<:AbstractPlot}, x::UnityData) 
    PlotSpec(Lines, x.position[:,1], x.position[:,2])
end

function MakieCore.convert_arguments(::Type{<:AbstractPlot}, x::UnityData, trial::Trial) 
    t,posx,posy = get_trial(x, trial.i)
    PlotSpec(Lines, posx, posy)
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
    daughters, ω = computeWavelets(n,c)
    freqs = getMeanFreq(daughters,fs)
    res = cwt(x, c, daughters)
    res, freqs
end

function get_spectrum(x::Vector{Vector{Float64}},fs;offset=0, kvs...)
    n = length(x)
    res = Vector{Matrix{ComplexF64}}(undef, n)
    freqs = Vector{Vector{Float64}}(undef, n)
    prog = Progress(n, desc="Computing spectrogram...",offset=offset)
    for (i,_x) in enumerate(x)
        res[i], freqs[i] = get_spectrum(_x, fs;kvs...)
        next!(prog)
    end
    res, freqs
end

"""
    get_frequency_band(X::Vector{Matrix{ComplexF64}},freqs::Vector{Vector{Float64}}, freq::Float64)

Extract the real part of `X` asociated with the specified frequency
"""
function get_frequency_band(X::Vector{Matrix{ComplexF64}},freqs::Vector{Vector{Float64}}, freq::Float64)
    nn = length(X)
    y = Vector{Vector{Float64}}(undef,nn)
    for (i,(x,ff)) in enumerate(zip(X,freqs))
        fidx = searchsortedfirst(ff, freq)
        y[i] = real.(x[:,fidx])
    end
    y
end

function TrialAlignedSpectrogram(x::Vector{Vector{Float64}}, t::Vector{Vector{Float64}},fs;β=1.5)
    res,freqs = get_spectrum(x, fs;β=β)
    TrialAlignedSpectrogram(t, res, freqs, fs,β)
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

"""
    get_trial_data()

Get trial aligned LFP, Spectrogram and trajectories for the current channel
"""
function get_trial_data(;do_save=true, redo=false,kvs...)
    fname = "trial_aligned_lfp.mat"
    args = Dict(kvs)
    if isfile(fname) && !redo
        qdata = MAT.matread(fname)
        aligned_lfp = convert(Vector{Vector{Float64}} ,qdata["aligned_lfp"])
        aligned_lfp_time = convert(Vector{Vector{Float64}}, qdata["aligned_lfp_time"])
        aligned_res = convert(Vector{Matrix{ComplexF64}}, qdata["aligned_spec"])
        aligned_freqs = convert(Vector{Vector{Float64}}, qdata["aligned_freqs"])
        aligned_eyepos = convert(Vector{Matrix{Float32}}, qdata["aligned_eye_pos"])
        aligned_eye_t = convert(Vector{Vector{Float64}}, qdata["aligned_eye_time"])
        aligned_maze_pos = convert(Vector{Matrix{Float64}}, qdata["aligned_maze_pos"])
        aligned_maze_time = convert(Vector{Vector{Float64}}, qdata["aligned_maze_time"])
    else
        # go down to the session directory to get the unity maze triggers
        udata,edata,rdata = cd(DataProcessingHierarchyTools.process_level("session")) do
            MAT.matread("unitymaze.mat"), MAT.matread("eyelink.mat"), MAT.matread("rplparallel.mat")
        end
        aligned_maze_pos_and_head_direction, aligned_maze_time = get_trial_data(udata["um"]["data"]["unityData"][:,3:5], udata["um"]["data"]["unityTime"][:], round.(Int64,udata["um"]["data"]["unityTriggers"]))
        aligned_maze_pos = [a[:,1:2] for a in aligned_maze_pos_and_head_direction]
        aligned_maze_head_direction = [a[:,3] for a in aligned_maze_pos_and_head_direction]

        #eye link data are in units of miliseconds
        aligned_eyepos, aligned_eye_t = get_trial_data(edata["el"]["data"]["eye_pos"], edata["el"]["data"]["timestamps"][:], edata["el"]["data"]["trial_timestamps"];
        pre_buffer=1000.0, post_buffer=1000.0)

        vdata = MAT.matread("rpllfp.mat")

        #lfp data are in units of seconds
        aligned_lfp, aligned_lfp_time = get_trial_data(vdata["df"]["data"]["analogData"][:], vdata["df"]["data"]["analogTime"][:], rdata["rp"]["data"]["timeStamps"])
        fs = vdata["df"]["data"]["analogInfo"]["SampleRate"]
        aligned_res, aligned_freqs = get_spectrum(aligned_lfp, fs;β=1.5)

        if do_save
            DrWatson.tag!(args, storepatch=true)
            MAT.matwrite(fname, Dict("aligned_lfp"=>aligned_lfp, "aligned_spec"=>aligned_res, "aligned_lfp_time"=>aligned_lfp_time,
                                    "aligned_maze_pos"=>aligned_maze_pos, "aligned_maze_time" => aligned_maze_time,
                                    "aligned_maze_head_direction"=>aligned_maze_head_direction[:,3], "aligned_freqs"=>aligned_freqs,
                                    "aligned_eye_pos"=>aligned_eyepos, "aligned_eye_time"=>aligned_eye_t, "args"=>args))
        end
    end
    aligned_lfp, aligned_lfp_time, aligned_res, aligned_freqs, aligned_eyepos, aligned_eye_t, aligned_maze_pos, aligned_maze_head_direction, aligned_maze_time
end

"""
    get_trial_data(data::Vector{T}, t::AbstractVector{T2}, triggers::Matrix{T2})

Align `data` to trials using triggers. Both `t` and `triggers` should be in the same units
"""
function get_trial_data(data::T4, t::AbstractVector{T2}, triggers::Matrix{T3};pre_buffer=1.0, post_buffer=1.0) where T4 <: VectorOrMatrix{T} where T <: Real where T2 <: Real where T3 <: Real 
    ntrials = size(triggers,1)
    aligned_data = Vector{T4}(undef, ntrials)
    aligned_t = Vector{Vector{Float64}}(undef, ntrials)
    for tidx in 1:ntrials
        idx0 = searchsortedfirst(t, triggers[tidx,1]-pre_buffer)
        idx1 = searchsortedlast(t, triggers[tidx,end]+post_buffer)
        aligned_data[tidx] = get_time_slice(data,idx0:idx1)
        aligned_t[tidx] = Float64.(t[idx0:idx1]) .- Float64(triggers[tidx,1])
    end
    aligned_data, aligned_t
end

function get_trial_data(data::VectorOrMatrix{T}, t::AbstractVector{T2}, triggers::Matrix{Int64};pre_buffer=1.0, post_buffer=1.0) where T <: Real where T2 <: Real
    ntrials = size(triggers,1)
    aligned_data = Vector{typeof(data)}(undef, ntrials)
    aligned_t = Vector{Vector{Float64}}(undef, ntrials)
    for tidx in 1:ntrials
        idx0 = triggers[tidx,1]
        idx0 = searchsortedfirst(t, t[idx0]-pre_buffer)
        idx1 = triggers[tidx,end]
        idx1 = searchsortedfirst(t, t[idx1]+post_buffer)
        aligned_data[tidx] = get_time_slice(data,idx0:idx1)
        aligned_t[tidx] = Float64.(t[idx0:idx1]) .- Float64(t[triggers[tidx,1]])
    end
    aligned_data, aligned_t    
end

function get_trial_data(unity_data::Matrix{T}, lfp_data::Vector{T2}, eyelink_data::Matrix{Float32}, unity_triggers::Matrix{T}, lfp_triggers::Matrix{T2}, eyelink_triggers::Matrix{T3};fs=1000.0, β=1.5) where T <: Real where T2 <: Real where T3 <: Real
    ntrials = size(lfp_triggers,1)
    nu = size(unity_data,1)
    nl = length(lfp_data)
    pos = Vector{Matrix{Float64}}(undef, ntrials)
    eye_pos = Vector{Matrix{Float64}}(undef, ntrials)
    lfp = Vector{Vector{Float64}}(undef, ntrials)
    spec = Vector{Matrix{ComplexF64}}(undef, ntrials)
    t = Vector{Vector{Float64}}(undef, ntrials)
    freqs = Vector{Vector{Float64}}(undef, ntrials)
    prog = Progress(ntrials, desc="Computing per trial spectrogram...")
    for tidx in 1:ntrials
        # unity
        uidx0 = max(round(Int64, unity_triggers[tidx, 1]), -1, 1)
        uidx1 = min(round(Int64, unity_triggers[tidx, 3])+1, nu)
        pos[tidx] = unity_data[uidx0:uidx1,3:4] 

        # eyelink
        uidx0 = max(round(Int64, unity_triggers[tidx, 1]), -1, 1)
        uidx1 = min(round(Int64, unity_triggers[tidx, 3])+1, nu)

        idx0 = max(round(Int64, lfp_triggers[tidx, 1]*1000)-1000,1)
        idx1 = min(round(Int64, lfp_triggers[tidx, 3]*1000)+1000,nl)
        lfp[tidx] = lfp_data[idx0:idx1]
        t[tidx] = range(idx0/1000.0, stop=idx1/1000, length=idx1-idx0+1)
        # compute power spectrum using wavelets
        spec[tidx],freqs[tidx] = get_spectrum(lfp[tidx],fs;β=β)
        next!(prog)
    end
    lfp,spec,t,pos,freqs
end


function plot_trial(udata::UnityData, edata::EyelinkData)
    ntrials = numtrials(udata)
    fig = Figure(size=(1000,750))
    ax = Axis(fig[1,1])

    trial = Observable(1)
    current_time = plot_trial!(ax, udata, edata, trial)
    ax2 = Axis(fig[1,2])
    plot_trial!(ax2, edata, trial, current_time)

    # trial selection via keys
    on(events(fig.scene).keyboardbutton) do event
        if event.action == Keyboard.press
            do_upate = false
            if event.key == Keyboard.left
                if trial[] > 0 
                    trial[] = trial[] - 1 
                    do_update = true
                end
            elseif event.key == Keyboard.right
                if trial[] < ntrials 
                    trial[] = trial[] + 1 
                    do_update = true
                end
            end
        end
    end
    fig
end

function plot_trial(udata::UnityData, edata::EyelinkData, trial)
    fig = Figure(size=(1000,750))
    ax = Axis(fig[1,1])
    current_time = plot_trial!(ax, udata, edata, trial)
    ax2 = Axis(fig[1,2])
    plot_trial!(ax2, edata, trial, current_time)
    fig
end

function plot_trial!(ax, edata::EyelinkData, trial::Int64, current_time)
    plot_trial!(ax, edata,Observable(trial),current_time)
end
"""
Plot eye link positions
"""
function plot_trial!(ax, edata::EyelinkData, trial::Observable{Int64}, current_time)
    gx0, gy0, gxm, gym = edata.header["gaze_coords"]
    ax.aspect = gxm/gym
    xlims!(ax, gx0, gxm)
    ylims!(ax, gy0, gym)
    _edata = lift(trial) do _trial
        te, gazex, gazey = get_trial(edata, trial[])
        te = (te .- te[1])/1000.0 # s
        gazey = gym .- gazey
        te,gazex, gazey 
    end

    tidx_e = lift(current_time) do ct
        te = _edata[][1]
        idx = searchsortedfirst(te, ct)
        idx
    end
    current_eyepos = lift(tidx_e) do _tidxe
        gazex = _edata[][2]
        gazey = _edata[][3]
        gx = gazex[_tidxe]
        gy = gazey[_tidxe]
        [Point2f(gx,gy)]
    end
    scatter!(ax, current_eyepos, color=:black)
end

function get_focal_plane_width(fov, focal_length)
    w = 2*focal_length*tan(fov/2)
    w
end

"""
Convert `x` from pixel coordinates to sensor coordinates
"""
function scale_to_camera(x, sensor_width, screen_width)
    x = sensor_width*(x -0.5*screen_width)/screen_width
end

"""
Trace a ray from position `x,y` through the camera with focal length `focal_length` until it
impacts something in the arena
"""
function raytrace(x, y, pos,direction, focal_length;camera_height=2.5)
    # find the angle of the point
    θ = atan(x, focal_length)
    ϕ = atan(y, focal_length)
    θ += direction
    # now cast along the line θ

    dl = 0.01
    xp,yp,zp = (pos[1], pos[2],camera_height)
    sθ = sin(θ)
    cθ = cos(θ)
    sϕ = sin(ϕ)
    while true
        dx,dy,dz = (dl*cθ, dl*sθ,dl*sϕ)
        xp += dx
        yp += dy 
        zp += dz 
        if impacts([xp,yp,zp])
            xp -= dx
            yp -= dy
            zp -= dz
            break
        end
    end
    xp,yp,zp
end

function plot_trial!(ax, udata::UnityData, edata::EyelinkData, trial::Int64)
    plot_trial!(ax, udata, edata, Observable(trial))
end

function plot_trial!(ax, udata::UnityData, edata::EyelinkData, trial::Observable{Int64})
    # TODO: Get this from data
    xmin,ymin, xmax, ymax = (-13.0, -13.0, 13.0, 13.0)
    room_height = 5.0
    camera_height = 2.5
    # 
    screen_width = edata.header["gaze_coords"][3] - edata.header["gaze_coords"][1]
    screen_height = edata.header["gaze_coords"][4] - edata.header["gaze_coords"][2]
    gx0 = edata.header["gaze_coords"][1]
    gy0 = edata.header["gaze_coords"][2]
    t, posx, posy,direction = get_trial(udata, trial[])
    _udata = lift(trial) do _trial
        t, posx, posy,direction = get_trial(udata, _trial)
        t .-= t[1] # start at zero
        t, posx, posy, direction
    end
    _edata = lift(trial) do _trial
        te, gazex, gazey = get_trial(edata, _trial)
        # TODO: Fix at source, eyelink flips the y-coordinate
        gazey = screen_height .- gazey
        te = (te .- te[1])./1000.0 # convert to seconds start at zero
        te, gazex, gazey
    end

    tlims = lift(_udata) do _ud
        extrema(_ud[1])
    end

    xlims!(ax, xmin, xmax)
    ylims!(ax, ymin, ymax)
    deregister_interaction!(ax, :scrollzoom)

    current_time = lift(_udata) do _ud
        _ud[1][1]
    end

    tidx = lift(current_time) do ct
        t = _udata[][1]
        idx = searchsortedfirst(t, ct)
        idx
    end

    tidx_e = lift(current_time) do ct
        te = _edata[][1]
        idx = searchsortedfirst(te, ct)
        idx
    end

    current_pos = lift(tidx) do _tidx
        posx,posy = _udata[][2:3]
        [Point2f(posx[_tidx], posy[_tidx])]
    end
    current_pos2 = lift(tidx) do _tidx
        posx,posy = _udata[][2:3]
        repeat([Point2f(posx[_tidx], posy[_tidx])],2)
    end

    current_dir = lift(tidx) do _tidx
        direction = _udata[][4]
        θ = π*direction[_tidx]/180 # convert to radians
        θ
    end

    current_arrow = lift(current_dir) do θ
        [Point2f(cos(θ), sin(θ))]
    end

    current_eyepos = lift(tidx_e) do _tidxe
        gazex = _edata[][2]
        gazey = _edata[][3]
        gx = gazex[_tidxe]
        if gx < 0.0 || gx > screen_width
            gx = NaN
        end
        gy = gazey[_tidxe]
        if gy < 0.0 || gy > screen_height
            gy = NaN
        end
        [Point2f(gx,gy)]
    end

    _fov = π*50.0/180
    l = 30.0
    fov = lift(current_dir, current_pos) do θ, _pos
        pos = _pos[1]
        θ0 = θ-0.5*_fov
        θ1 = θ+0.5*_fov

        # brute force raycasting
        # extent lines until they hit something
        l0 = 0.0
        dx = 0.01
        while true
            l0 += dx 
            x0 = pos[1] + l0*cos(θ0)
            y0 = pos[2] + l0*sin(θ0)
            if impacts([x0,y0, 0.0])
                # the ray hit something, so we stop here
                break
            end
        end

        l1 = 0.0
        while true
            l1 += dx 
            x1 = pos[1] + l1*cos(θ1)
            y1 = pos[2] + l1*sin(θ1)
            if impacts([x1,y1, 0.0])
                break
            end
        end

        [Point2f(l0*cos(θ0), l0*sin(θ0)), Point2f(l1*cos(θ1), l1*sin(θ1))]
    end

    projected_eyepos = lift(current_dir, current_pos, current_eyepos) do _dir, _pos, _epos
        # convert to sensor coordinates
        gx = _epos[1][1]
        gy = _epos[1][2]
        if isfinite(gx) && isfinite(gy)
            # Unity camera sensor size: (36,24)
            x = scale_to_camera(gx-gx0, 36.0, screen_width)
            y = scale_to_camera(gy-gy0, 24.0, screen_height)

            #TODO: We need to know the height of the camera (and the distance to the ceiling)
            # to correctly project onto the floor or the ceiling
            #ray trace
            xp,yp = raytrace(-x, -y, _pos[1], _dir, 50.0)
        else
            xp,yp = (NaN,NaN)
        end
        [Point2(xp,yp)]
    end

    #handle scroll event
    on(events(ax).scroll) do (dx,dy)
        tmin,tmax = tlims[]
        if tmin < current_time[]+dx < tmax 
            current_time[] = current_time[] + dx
        end
    end
    ax.aspect = 1.0
    plot_arena!(ax)
    scatter!(ax, current_pos)
    arrows!(ax, current_pos, current_arrow)
    arrows!(ax, current_pos2, fov)
    scatter!(ax, projected_eyepos, color=:black)
    current_time[] = t[1]
    current_time
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

function plot_trial(pos::Vector{Matrix{T}}, eye_pos::Vector{Matrix{T2}}, lfp::Vector{Vector{T}}, spec::Vector{Matrix{ComplexF64}}, freqs::Vector{Vector{T}},t::Vector{Vector{T}};
                                                        arena_size::Union{Nothing, NTuple{4,Float64}}=nothing,max_freq=Inf) where T <: Real where T2 <: Real
    ntrials = length(pos)
    if max_freq < Inf
        _spec = Vector{Matrix{ComplexF64}}(undef, length(spec))
        _freqs = Vector{Vector{T}}(undef, length(spec))
        for i in 1:ntrials
            fidx = searchsortedlast(freqs[i], max_freq)
            _spec[i] = spec[i][:,1:fidx]
            _freqs[i] = freqs[i][1:fidx]
        end
    else
        _spec = spec
        _freqs = freqs
    end

    tidx = Observable(1)
    fidx = Observable(0)
    p_pos = Observable([Point2f(pos[tidx[]][i,:]...) for i in 1:size(pos[tidx[]],1)])
    p_eyepos = Observable([Point2f(eye_pos[tidx[]][i,:]...) for i in 1:size(eye_pos[tidx[]],1)])

    p_lfp = Observable([Point2f(_t,_lfp) for (_t,_lfp) in zip(t[tidx[]], lfp[tidx[]])])

    p_spec = Observable(abs.(_spec[tidx[]]))
    p_freqs = Observable(_freqs[tidx[]])
    p_phase = Observable(fill(Point2f(NaN), length(t[tidx[]])))

    p_t = Observable(t[tidx[]])
    trial_start_unity = Observable(Point2f[])
    current_time = Observable(p_t[][1])
    current_pos = Observable([Point2f(p_pos[][1])])
    current_eyepos = Observable([Point2f(p_eyepos[][1])])
    current_freq = Observable(NaN)


    fig = Figure()
    lg = GridLayout(fig[1,1])
    ax1 = Axis(lg[1,1])
    ax2 = Axis(lg[1,2])
    axes = [Axis(fig[i+1,1]) for i in 1:3]
    rowsize!(fig.layout, 3, Relative(0.05))
    hidedecorations!(axes[2])
    axes[2].bottomspinevisible = false
    axes[2].ylabelvisible = true
    axes[2].ylabel = "Phase"
    for ax in axes[1:end]
        ax.xgridvisible = false
        ax.ygridvisible = false
        ax.topspinevisible = false
        ax.rightspinevisible = false
    end
    # override scroll zoom
        deregister_interaction!(ax, :scrollzoom)
    end
    linkxaxes!(axes...)
    if arena_size !== nothing
        xlims!(ax1, arena_size[1:2]...)
        ylims!(ax1, arena_size[3:4]...)
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
        p_eyepos[] =[Point2f(eye_pos[_tidx][i,:]...) for i in 1:size(eye_pos[_tidx],1)] 
        current_pos[] = [p_pos[][1]]
        current_eyepos[] = [p_eyepos[][1]]
        trial_start_unity[] = p_pos[][1:1]
        if fidx[] == 0
            # show entire lfp
            p_phase[] = fill(Point2f(NaN), length(t[tidx[]]))
            p_lfp[] = [Point2f(_t,_lfp) for (_t,_lfp) in zip(t[tidx[]][qidx], lfp[_tidx][qidx])]
        else
            p_phase[] = [Point2f(_t,a) for (_t,a) in zip(t[tidx[]][qidx], angle.(spec[_tidx][qidx,fidx[]]))]
            p_lfp[] = [Point2f(_t,_lfp) for (_t,_lfp) in zip(t[tidx[]][qidx], real.(spec[_tidx][qidx,fidx[]]))]
        end
        p_t[] = t[_tidx][qidx]
        current_time[] = p_t[][1]

        p_freqs[] = _freqs[_tidx]
        p_spec[] = abs.(_spec[_tidx][qidx,:])
    end

    on(fidx) do _fidx
        _tidx = tidx[]
        _t = t[_tidx]
        if length(_t) > 16383
            qidx = round.(Int64,range(1, stop=length(_t), length=16383))
        else
            qidx = 1:length(_t)
        end
        if _fidx > 0 
            #TODO: Fix the code dupication between this and the tidx update
            p_phase[] = [Point2f(_t,a) for (_t,a) in zip(t[tidx[]][qidx], angle.(spec[_tidx][qidx,fidx[]]))]
            p_lfp[] = [Point2f(_t,_lfp) for (_t,_lfp) in zip(t[_tidx][qidx], real.(spec[_tidx][qidx,_fidx]))]
            current_freq[] = freqs[_tidx][_fidx]
        else
            p_phase[] = fill(Point2f(NaN), length(t[tidx[]]))
            p_lfp[] = [Point2f(_t,_lfp) for (_t,_lfp) in zip(t[_tidx][qidx], real.(lfp[_tidx]))]
            current_freq[] = NaN 
        end
        autolimits!(axes[2])
        autolimits!(axes[3])
    end

    function handle_scroll(dx)
        tmax = p_t[][end]
        tmin = p_t[][1]
        Δt = tmax - tmin
        n_pos = length(p_pos[])
        n_eyepos = length(p_eyepos[])
        if tmin < current_time[]+dx < tmax
            current_time[] = current_time[] + dx
            pidx = min(max(1, round(Int64,n_pos*(current_time[]-tmin)/Δt)),n_pos) 
            current_pos[] = [p_pos[][pidx]]
            pidx = min(max(1, round(Int64,n_eyepos*(current_time[]-tmin)/Δt)),n_eyepos) 
            current_eyepos[] = [p_eyepos[][pidx]]
        end
    end

    #handle scroll event
    on(events(fig.scene).scroll) do (dx,dy)
        handle_scroll(dx)
    end

    on(events(fig.scene).keyboardbutton) do event
        if event.action == Keyboard.press
            do_scroll = false
            if Keyboard.left_shift in events(fig.scene).keyboardstate
                do_scroll = true
            end
            do_update = false
            if event.key == Keyboard.left
                if do_scroll
                    handle_scroll(-0.1)
                else
                    if tidx[] > 1 
                        tidx[] = tidx[] - 1 
                        do_update = true
                    end
                end
            elseif event.key == Keyboard.right
                if do_scroll
                    handle_scroll(0.1)
                else
                    if tidx[] < ntrials 
                        tidx[] = tidx[] + 1 
                        do_update = true
                    end
                end
            end
            if do_update
                autolimits!(axes[1])
                autolimits!(axes[2])
                autolimits!(axes[3])
            end
        end
    end

    lines!(ax1, p_pos)
    lines!(ax2, p_eyepos;color=:lightgray)
    scatter!(ax1, trial_start_unity, color=:red)
    scatter!(ax1, current_pos, color=:green)
    scatter!(ax2, current_eyepos, color=:green)
    h = heatmap!(axes[1], p_t, p_freqs, p_spec)
    vlines!(axes[1], current_time, color=:green)
    hlines!(axes[1], current_freq, color=:white)
    axes[1].ylabel = "Frequency [Hz]"
    axes[1].xticklabelsvisible = false
    lines!(axes[2], p_phase)
    lines!(axes[3], p_lfp)
    vlines!(axes[2], current_time, color=:green)
    vlines!(axes[3], current_time, color=:green)
    axes[3].xlabel = "Time [s]"
    axes[3].ylabel = "LFP"

    mouseevents = addmouseevents!(fig.scene, h)

    onmouseleftdoubleclick(mouseevents, priority=4) do event
        plt, i = pick(fig.scene)
        fidx[] = 0
        Consume(true)
    end

     # handle click
     onmouseleftclick(mouseevents, priority=3) do event
        plt, i = pick(fig.scene)
        _pos = mouseposition(axes[1])
        # find the clicked frequency
        fidx[] = searchsortedfirst(p_freqs[], _pos[2])
        Consume(true)
    end

    tidx[] = 1
    fig
end

function pan_through(X::Vector{Matrix{T}};label::Union{Vector{String},Nothing}=nothing) where T
    fig = Figure()
    # TODO: Ensure that all elements of `X` have the same dimensions
    axes = [Axis(fig[i,j]) for i in 1:size(X[1],1), j in 1:size(X[1],2)]

    ntrials = length(X)
    tidx = Observable(1)
    xs = [Observable(X[tidx[]][i])  for i in 1:length(X[1])]
    header = Label(fig[0, 1:size(axes,2)], "Index $(tidx[])")
    on(tidx) do _tidx
        for i in 1:length(xs)
            xs[i][] = X[_tidx][i]
        end
        #axes[1].title[] = "Index $(_tidx)"
        header.text[] = "Index $(_tidx)"
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
                for ax in axes
                    autolimits!(ax)
                end
            end
        end
    end
    if label !== nothing
        _label = label
    else
        _label = ["" for i in 1:length(xs)]
    end
    for (ax,x,l) in zip(axes, xs,_label)
        ax.ylabel = l
        plot!(ax, x)
        ax.xgridvisible = false
        ax.ygridvisible = false
        ax.xticklabelsvisible = false
        ax.xticksvisible = false
        ax.yticklabelsvisible = false
        ax.yticksvisible = false
    end
    tidx[] = 1
    fig
end

"""
Create a spatial map weighted by LFP power. This is essentially a map showing which parts of the arena elicit that highest
LFP response.
"""
function get_spatial_map(maze_pos::Vector{Matrix{Float64}}, maze_time::Vector{Vector{Float64}}, lfp::Vector{Vector{Float64}}, lfp_time::Vector{Vector{Float64}})
    xmin,xmax,ymin,ymax = (Inf,-Inf,Inf,-Inf)
    for mpos in maze_pos
        for pos in eachrow(mpos)
            xmin = min(xmin, pos[1])
            xmax = max(xmax, pos[1])
            ymin = min(ymin, pos[2])
            ymax = max(ymax, pos[2])
        end
    end

    xbins = range(xmin, stop=xmax, length=20)
    ybins = range(ymin, stop=ymax, length=20)

    # get occupancy
    nn = 0
    img = fill(0.0, length(xbins), length(ybins))
    for mpos in maze_pos
        for pos in eachrow(mpos)
            xidx = searchsortedfirst(xbins, pos[1])
            yidx = searchsortedfirst(ybins, pos[2])
            img[xidx,yidx] += 1.0
            nn += 1.0
        end
    end
    img ./= nn

    # get lfp power weighted spatial map
    pp = 0.0
    img2 = fill(0.0, length(xbins), length(ybins))
    for i in 1:length(lfp_time)
        u_t = maze_time[i]
        u_pos = maze_pos[i]
        for (_t,d) in zip(lfp_time[i], lfp[i])
            # match time between unity and lfp
            tidx = searchsortedfirst(u_t, _t)
            if 0 < tidx <= length(u_t)
                pos = u_pos[tidx,:] 
                xidx = searchsortedfirst(xbins, pos[1])
                yidx = searchsortedfirst(ybins, pos[2])
                img2[xidx, yidx] += abs(d)
                pp += abs(d)
            end
        end
    end
    img2 ./= pp
    SpatialPreferenceMap(xbins, ybins, img, img2)
end

function get_spatial_map(udata::UnityData; grid_limits=(-12.0, -12.0, 12.0, 12.0), nbins=(20,20))
    nt = size(udata.timestamps,1)

    xmin,ymin,xmax,ymax = grid_limits
    xbins = range(xmin, stop=xmax, length=nbins[1])
    ybins = range(ymin, stop=ymax, length=nbins[2])
    nn = 0
    img = fill(0.0, length(xbins), length(ybins))
    for i in 1:nt
        _,mposx,mposy = get_trial(udata, i;trial_start=2) # start of navigation
        for (posx,posy) in zip(mposx, mposy)
            xidx = searchsortedfirst(xbins, posx)
            yidx = searchsortedfirst(ybins, posy)
            img[xidx,yidx] += 1.0
            nn += 1.0
        end
    end
    img ./= nn
    xbins,ybins, img
end

function get_spatial_map(udata::UnityData, spiketrain::Spiketrain, rpdata::RippleData; grid_limits=(-12.0, -12.0, 12.0, 12.0), nbins=(20,20), 
                                                                                       normalize_by_occupancy=true)
    xbins,ybins,img = get_spatial_map(udata;grid_limits=grid_limits, nbins=nbins)
    # use the ripple markers to get trial structure
    timestamps_rp = rpdata.timestamps
    nt = size(timestamps_rp,1)
    nn = 0
    img2 = fill!(similar(img), 0.0)
    for i in 1:nt
        idx0 = searchsortedfirst(spiketrain.timestamps, 1000*timestamps_rp[i,2]) # start of navigation
        idx1 = searchsortedlast(spiketrain.timestamps, 1000*timestamps_rp[i,3])
        _spiketrain = spiketrain.timestamps[idx0:idx1]/1000.0 # convert to seconds
        # get the corresponding position
        t,mposx, mposy = get_trial(udata,i;trial_start=2)
        # re-reference to ripple time
        trp = timestamps_rp[i,2]
        t = t .- t[1] .+ trp[1]
        for sp in _spiketrain
            idx = searchsortedfirst(t, sp)
            if idx < length(t)
                _posx = mposx[idx]
                _posy = mposy[idx]
                xidx = searchsortedfirst(xbins, _posx)
                yidx = searchsortedfirst(ybins, _posy)
                img2[xidx, yidx] += 1.0
                nn += 1.0
            end
        end
    end
    img2 ./= nn
    if normalize_by_occupancy
        img2 ./= img
    end
    img2, xbins, ybins
end

"""
Simulate a simple place field neuron using behavioural data in `udata`
"""
function model_place_field(udata::UnityData, rpdata::RippleData;λmin=0.1, λmax=3.0,dt=0.01,σ=1.0,μ=[3.5,1.2])
    Σ = [σ 0.0;0.0 σ]
    G = MvNormal(μ, Σ)
    G0 = pdf(G, μ)
    nt = numtrials(udata)

    spikes = Float64[]
    q = -log(rand())
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
                    q = -log(rand())
                end
                _t += dt
            end
        end
    end
    spikes
end

function get_spatial_map(;freq::Union{Nothing, Float64}=nothing, redo=false)
    aligned_lfp, aligned_lfp_time, aligned_spec,aligned_freqs,~,~,aligned_maze_pos, aligned_maze_time = SpatialAnalyses.get_trial_data(;redo=redo)
    if freq === nothing
        spatial_map = get_spatial_map(aligned_maze_pos, aligned_maze_time, aligned_lfp, aligned_lfp_time)
    else
        X = get_frequency_band(aligned_spec, aligned_freqs, freq)
        spatial_map = get_spatial_map(aligned_maze_pos, aligned_maze_time, X, aligned_lfp_time)
    end
    spatial_map
end

end