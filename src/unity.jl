using Makie
using Colors
using Meshes
using LinearAlgebra
using DelimitedFiles

xBound = [-12.5, 12.5, 12.5, -12.5, -12.5]
zBound = [12.5, 12.5, -12.5, -12.5, 12.5]
x1Bound = [-7.5, -2.5, -2.5, -7.5, -7.5]  # yellow pillar
z1Bound = [7.5, 7.5, 2.5, 2.5, 7.5]
x2Bound = [2.5, 7.5, 7.5, 2.5, 2.5]  # red pillar
z2Bound = [7.5, 7.5, 2.5, 2.5, 7.5]
x3Bound = [-7.5, -2.5, -2.5, -7.5, -7.5]  # blue pillar
z3Bound = [-2.5, -2.5, -7.5, -7.5, -2.5]
x4Bound = [2.5, 7.5, 7.5, 2.5, 2.5]  # green pillar
z4Bound = [-2.5, -2.5, -7.5, -7.5, -2.5]

# TODO: Use actual values here
camera_height = 2.5
ceiling_height = 5.0

struct UnityData
    time::Vector{Float64}
    position::Matrix{Float64}
    head_direction::Vector{Float64}
    triggers::Matrix{Int64}
    timestamps::Matrix{Float64}
    header::Dict
end

function UnityData(fname::String)
    data, header, column_names = read_unity_file(fname)
    _time = cumsum(data[:,2])

    # organize triggers into a trial structure
    # triggers 1x indicate trial start, where x indicates the poster number
    trial_start_idx = findall(10.0 .< data[:,1] .< 20.0)

    # an additional marker indicates when cue for the monkey to start navigating
    trial_start_nav = findall(20 .<= data[:,1] .< 30)
    # trial end is either 3x, or 4x, where 3 indicates success and 4 indicates time out
    trial_end_idx = findall(30.0 .< data[:,1] .< 50.0)
    length(trial_start_idx) == length(trial_start_nav) == length(trial_end_idx) || error("Inconsitent number of triggers")
    nt = length(trial_start_idx)
    triggers = [data[trial_start_idx,1] data[trial_start_nav,1] data[trial_end_idx,1]]
    timestamps = [_time[trial_start_idx] _time[trial_start_nav] _time[trial_end_idx]]
    # the unity coordinate system swaps x and y
    UnityData(_time, data[:,[4,3]], data[:,5], triggers, timestamps, header)
end

numtrials(x::UnityData) = size(x.triggers,1)

function get_trial(data::UnityData, i;trial_start=1)
    idx0 = searchsortedfirst(data.time, data.timestamps[i,trial_start])
    idx1 = searchsortedfirst(data.time, data.timestamps[i,3])
    data.time[idx0:idx1], data.position[idx0:idx1, 1], data.position[idx0:idx1, 2], data.head_direction[idx0:idx1]
end

"""
Read the file `fname`, assuming each column is separated by a single space, and
the first 14 rows contain header information
"""
function read_unity_file(fname::String;header_rows=14)
    # TODO: Return the config as well
    column_names = ["marker","Δt","xpos","ypos","direction"]  
    data = readdlm(fname, ' ', skipstart=header_rows) 
    # check if first column should be skipped
    if size(data,2) == 6 # should only be 5 columns
        data = data[:,2:end]
    end
    header = Dict()
    open(fname) do fid
        for i in 1:header_rows
            _line = readline(fid)
            k,v = split(_line, ':')
            header[k] = v
        end
    end
    data = convert(Matrix{Float64}, data)
    data, header, column_names
end

function plot_arena()
    fig = Figure()
    ax = Axis(fig[1,1])
    plot_arena!(ax)
    fig,ax
end

function plot_arena!(ax)
    poly!(ax, Point2f.(zip(zBound, xBound)),color=:grey)
    poly!(ax, Point2f.(zip(z1Bound, x1Bound)), color=:yellow)
    poly!(ax, Point2f.(zip(z2Bound, x2Bound)), color=:red)
    poly!(ax, Point2f.(zip(z3Bound, x3Bound)), color=:blue)
    poly!(ax, Point2f.(zip(z4Bound, x4Bound)), color=:green)
end

function soft_range(start::T, stop::T,step::T) where T <: Real
    nn = round(Int64,(start-stop)/step)
    Δ = (start-stop)/nn
    range(start,stop=stop,step=Δ)
end

"""
Check whether the point `pos` intersects with any of the pillars, or the walls
of the arena.
"""
function impacts(pos)
    x,y,z = pos
    b1 = -7.5 <= x <= -2.5 && 2.5 <= y <= 7.5
    b2 = -7.5 <= x <= -2.5 && -7.5 <= y <= -2.5
    b3 = 2.5 <= x <= 7.5 && -7.5 <= y <= -2.5
    b4 = 2.5 <= x <= 7.5 && 2.5 <= y <= 7.5
    b5 = x < -12.5 || x > 12.5
    b6 = y < -12.5 || y > 12.5
    b7 = z > ceiling_height || z < 0.0
    b1 || b2 || b3 || b4 || b5 || b6 || b7
end

"""
Return the meshes representing the maze
"""
function create_maze(;kvs...)
    # unity uses 40x40 bins on the floor
    floor_bins = range(-12.5, stop=12.5, length=40)
    Δb = step(floor_bins)
    @show Δb
    bins = Matrix{NTuple{3,AbstractVector{Float64}}}(undef, 4,4)
    wall_bins = Vector{NTuple{3,AbstractVector{Float64}}}(undef, 4)
    zbins = range(0.0, stop=5.0, length=10)
    normals_bins = Matrix{Vector{Float64}}(undef, 4,4)
    normals_walls = Vector{Vector{Float64}}(undef, 4)
    Δ = get(kvs, :Δz, 0.1)
    # pillar 1
    xbins = soft_range(-7.5, -2.5, Δb)
    y0 = 2.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[1,1] = (xbins, ybins, zbins)
    normals_bins[1,1] = [0.0, -1.0, 0.0]

    y0 = 7.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[2,1] = (xbins, ybins, zbins)
    normals_bins[2,1] = [0.0, 1.0, 0.0]

    ybins = soft_range(2.5, 7.5,Δb) 
    x0 = -7.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[3,1] = (xbins, ybins, zbins)
    normals_bins[3,1] = [-1.0, 0.0, 0.0]
    x0 = -2.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[4,1] = (xbins, ybins, zbins)
    normals_bins[4,1] = [1.0, 0.0, 0.0]

    # pillar 2
    xbins = soft_range(2.5, 7.5, Δb)
    y0 = 2.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[1,2] = (xbins,ybins, zbins)
    normals_bins[1,2] = [0.0, -1.0, 0.0]
    y0 = 7.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[2,2] = (xbins, ybins, zbins)
    normals_bins[2,2] = [0.0, 1.0, 0.0]

    ybins = soft_range(2.5, 7.5,Δb)
    x0 = 2.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[3,2] = (xbins, ybins, zbins)
    normals_bins[3,2] = [-1.0, 0.0, 0.0]
    x0 = 7.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[4,2] = (xbins, ybins, zbins)
    normals_bins[4,2] = [1.0, 0.0, 0.0]

    # pillar 3
    xbins = soft_range(2.5, 7.5, Δb)
    y0 = -2.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[1,3] = (xbins, ybins, zbins)
    normals_bins[1,3] = [0.0, 1.0, 0.0]
    y0 = -7.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[2,3] = (xbins, ybins, zbins)
    normals_bins[2,3] = [0.0 , -1.0, 0.0] 
    ybins = soft_range(-7.5, -2.5,Δb)
    x0 = 2.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[3,3] = (xbins, ybins, zbins)
    normals_bins[3,3] = [-1.0, 0.0, 0.0]
    x0 = 7.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[4,3] = (xbins, ybins, zbins)
    normals_bins[4,3] = [1.0, 0.0, 0.0]

    # pillar 4
    xbins = soft_range(-7.5, -2.5, Δb)
    y0 = -2.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[1,4] = (xbins, ybins, zbins)
    normals_bins[1,4] = [0.0, 1.0, 0.0]
    y0 = -7.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[2,4] = (xbins, ybins, zbins)
    normals_bins[2,4] = [0.0, -1.0, 0.0]
    ybins = soft_range(-7.5, -2.5,Δb)
    x0 = -2.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[3,4] = (xbins, ybins, zbins)
    normals_bins[3,4] = [1.0, 0.0, 0.0]
    x0 = -7.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[4,4] = (xbins, ybins, zbins)
    normals_bins[4,4] = [-1.0, 0.0, 0.0]

    # walls
    xbins = soft_range(-12.5, 12.5, Δb)
    y0 = -12.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    wall_bins[1] = (xbins,ybins,zbins)
    normals_walls[1] = [0.0, 1.0,0.0]
    y0 = 12.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    wall_bins[2] = (xbins,ybins,zbins)
    normals_walls[2] = [0.0, -1.0, 0.0]

    ybins = soft_range(-12.5, 12.5, Δb)
    x0 = -12.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    wall_bins[3] = (xbins, ybins, zbins)
    normals_walls[3] = [1.0, 0.0, 0.0]
    x0 = 12.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    wall_bins[4] = (xbins, ybins, zbins)
    normals_walls[4] = [-1.0, 0.0, 0.0]
   (pillars=bins, walls=wall_bins), (pillars=normals_bins, walls=normals_walls)
end

"""
Compute a histogram of `pos` projected onto the plane at `z0`
"""
function compute_histogram(pos::Matrix{Float64}, xbins,ybins,z0=0.0;Δz=0.1)
    counts = fill(0.0, length(xbins), length(ybins))
    nn = 0.0
    for (x,y,z) in eachcol(pos)
        idx_x = searchsortedlast(xbins, x)
        idx_y = searchsortedlast(ybins, y)
        if idx_x > 0 && idx_y > 0
            if z0 - Δz <= z <= z0 + Δz
                counts[idx_x, idx_y] += 1.0
                nn += 1.0 
            end
        end
    end
    counts
end


function compute_histogram(pos::Matrix{Float64},bins)
    counts = [fill(0.0, length.(bin)) for bin in bins]
    compute_histogram!(counts, pos,bins)
end

function compute_histogram(pos::Vector{Matrix{Float64}},bins)
    counts = [fill(0.0, length.(bin)) for bin in bins]
    for _pos in pos
        compute_histogram!(counts, _pos, bins)
    end
    counts
end

function compute_histogram!(counts, pos::Matrix{Float64},bins)
    qpos = ([pos[i,:] for i in 1:size(pos,1)]...,)
    for (bin,count) in zip(bins,counts)
        # hackish; add one bin to the end
        Δs = [step(b) for b in bin]
        h = fit(Histogram, qpos, ([[b;b[end]+Δ] for (b,Δ) in zip(bin,Δs)]...,))
        count .+= h.weights
    end
    counts
end

