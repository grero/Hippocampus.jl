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

pillar_color = Dict(:pillar_1=>:yellow, :pillar_2=>:red, :pillar_3=>:blue, :pillar_4=>:green)

poster_pos = [[-5, -7.55], [-7.55, 5], [7.55, -5], [5, 7.55], [5, 2.45], [-5, -2.45]]
# for some reason x and y appear to be flipped
poster_pos = reverse.(poster_pos)
poster_img = joinpath.("artefacts",  ["camel 1.png","cat 1.png","crocodile.png","donkey 1.png","pig 1.png","rabbit 1.png"])

# TODO: Use actual values here
camera_height = 2.5
ceiling_height = 5.0

struct UnityData
    time::Vector{Float64}
    position::Matrix{Float64}
    head_direction::Vector{Float64}
    triggers::Matrix{Int64}
    timestamps::Matrix{Float64}
    poster_pos::Vector{Vector{Float64}}
    header::Dict
end

DPHT.level(::Type{UnityData}) = "session"

function UnityData()
    # attempt to find data file
    # assume we are at the sesison level
    if isdir("RawData_T1-100")
        ff = glob("session_*.txt", "RawData_T1-100")
        if !isempty(ff)
            return UnityData(first(ff))
        end
    end
    error("No unity data found")
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

    # parse the poster locations
    poster_position = header["PosterLocations"]
    loc = split(poster_position)
    rr = r"P([\d+.])\(([-\d\.]+),[-\d.]+,([-\d.]+)\)"
    poster_loc = Vector{Vector{Float64}}(undef, length(loc))
    for (i,_loc) in enumerate(loc)
        m = match(rr, _loc)
        if m !== nothing
            # first is the poster number
            j = parse(Int64, m.captures[1])
            poster_loc[j] = parse.(Float64, reverse(m.captures[2:end]))
        end
    end
    UnityData(_time, data[:,[4,3]], data[:,5], triggers, timestamps, poster_loc, header)
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

    bins = Dict{Symbol,Vector{NTuple{3,AbstractVector{Float64}}}}()
    normals = Dict{Symbol,Vector{Vector{Float64}}}()

    zbins = range(0.0, stop=5.0, length=10)
    Δ = get(kvs, :Δz, 0.1)
    # pillar 1
    bins[:pillar_1] = Vector{AbstractVector{Float64}}(undef, 4)
    normals[:pillar_1] = Vector{Vector{Float64}}(undef, 4)
    xbins = soft_range(-7.5, -2.5, Δb)
    y0 = 2.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[:pillar_1][1] = (xbins,ybins,zbins)
    normals[:pillar_1][1] = [0.0, -1.0,0.0]

    y0 = 7.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[:pillar_1][2] = (xbins,ybins,zbins)
    normals[:pillar_1][2] = [0.0, 1.0, 0.0]

    ybins = soft_range(2.5, 7.5,Δb) 
    x0 = -7.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[:pillar_1][3] = (xbins,ybins,zbins)
    normals[:pillar_1][3] = [-1.0, 0.0, 0.0]
    x0 = -2.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[:pillar_1][4] = (xbins,ybins,zbins)
    normals[:pillar_1][4] = [1.0, 0.0, 0.0]

    # pillar 2
    bins[:pillar_2] = Vector{AbstractVector{Float64}}(undef, 4)
    normals[:pillar_2] = Vector{Vector{Float64}}(undef, 4)
    xbins = soft_range(2.5, 7.5, Δb)
    y0 = 2.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[:pillar_2][1] = (xbins,ybins,zbins)
    normals[:pillar_2][1] = [0.0, -1.0, 0.0]
    y0 = 7.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[:pillar_2][2] = (xbins,ybins,zbins)
    normals[:pillar_2][2] = [0.0, 1.0, 0.0]

    ybins = soft_range(2.5, 7.5,Δb)
    x0 = 2.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[:pillar_2][3] = (xbins,ybins,zbins)
    normals[:pillar_2][3] = [-1.0, 0.0, 0.0]
    x0 = 7.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[:pillar_2][4] = (xbins, ybins,zbins)
    normals[:pillar_2][4] = [1.0, 0.0, 0.0]

    # pillar 3
    bins[:pillar_3] = Vector{AbstractVector{Float64}}(undef, 4)
    normals[:pillar_3] = Vector{Vector{Float64}}(undef, 4)
    xbins = soft_range(2.5, 7.5, Δb)
    y0 = -2.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[:pillar_3][1] = (xbins,ybins,zbins)
    normals[:pillar_3][1] = [0.0, 1.0, 0.0]
    y0 = -7.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[:pillar_3][2] = (xbins, ybins,zbins)
    normals[:pillar_3][2] = [0.0, -1.0, 0.0]
    ybins = soft_range(-7.5, -2.5,Δb)
    x0 = 2.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[:pillar_3][3] = (xbins, ybins,zbins)
    normals[:pillar_3][3] = [-1.0, 0.0, 0.0]
    x0 = 7.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[:pillar_3][4] = (xbins, ybins, zbins)
    normals[:pillar_3][4] = [1.0, 0.0, 0.0]

    # pillar 4
    bins[:pillar_4] = Vector{AbstractVector{Float64}}(undef,4)
    normals[:pillar_4] = Vector{Vector{Float64}}(undef, 4)
    xbins = soft_range(-7.5, -2.5, Δb)
    y0 = -2.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[:pillar_4][1] = (xbins, ybins, zbins)
    normals[:pillar_4][1] = [0.0, 1.0, 0.0]
    y0 = -7.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[:pillar_4][2] = (xbins, ybins, zbins)
    normals[:pillar_4][2] = [0.0, -1.0, 0.0]
    ybins = soft_range(-7.5, -2.5,Δb)
    x0 = -2.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[:pillar_4][3] = (xbins, ybins, zbins)
    normals[:pillar_4][3] = [1.0, 0.0, 0.0]
    x0 = -7.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[:pillar_4][4] = (xbins, ybins, zbins)
    normals[:pillar_4][4] = [-1.0, 0.0, 0.0]

    # walls
    bins[:walls] = Vector{AbstractVector{Float64}}(undef, 4)
    normals[:walls] = Vector{Vector{Float64}}(undef, 4)
    xbins = soft_range(-12.5, 12.5, Δb)
    y0 = -12.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[:walls][1] = (xbins,ybins,zbins)
    normals[:walls][1] = [0.0, 1.0,0.0]
    y0 = 12.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[:walls][2] = (xbins,ybins,zbins)
    normals[:walls][2] = [0.0, -1.0, 0.0]

    ybins = soft_range(-12.5, 12.5, Δb)
    x0 = -12.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[:walls][3] = (xbins, ybins, zbins)
    normals[:walls][3] = [1.0, 0.0, 0.0]
    x0 = 12.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[:walls][4] = (xbins, ybins, zbins)
    normals[:walls][4] = [-1.0, 0.0, 0.0]

    # floor
    xbins = floor_bins
    ybins = floor_bins
    z0 = 0.0
    zbins = range(z0-Δ, stop=z0+Δ, length=2)
    bins[:floor] = [(xbins, ybins, zbins)]
    normals[:floor] = [[0.0, 0.0, 1.0]]

    # ceiling
    xbins = range(-12.5, stop=12.5, step=Δb)
    ybins = range(-12.5, stop=12.5, step=Δb)
    z0 = 5.0
    zbins = range(z0-Δ, stop=z0+Δ, length=2)
    bins[:ceiling] = [(xbins, ybins, zbins)]
    normals[:ceiling]=[[0.0, 0.0 -1.0]]

    bins, normals
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

