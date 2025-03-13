using Makie
using Colors
using Meshes
using LinearAlgebra
using DelimitedFiles
using FileIO
using StatsBase
using ImageFiltering

# TODO: Unclear if these are the latest values. Perhaps update?
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

pillar_positions = Dict()
pillar_positions[:blue] = Dict(:lower_left => (-7.5, -7.5),
                                :upper_right => (-2.5, -2.5))
pillar_positions[:yellow] = Dict(:lower_left => (-7.5, 2.5),
                                 :upper_right => (-2.5, 7.5))
pillar_positions[:green] = Dict(:lower_left => (2.5, -7.5),
                                :upper_right => (7.5, -2.5))
pillar_positions[:red] = Dict(:lower_left => (2.5, 2.5),
                              :upper_right => (7.5, 7.5))
poster_pos = Dict{Symbol,NTuple{2,Float64}}()
poster_pos[:camel] = (-7.567, 5.0)
poster_pos[:cat] = (-5.0, -7.6)
poster_pos[:pig] = (-5.0, -2.44)
poster_pos[:donkey] = (5.168, 7.561)
poster_pos[:croc] = (5.0, 2.433)
poster_pos[:rabbit] = (7.561, -5.0)


#poster_pos = [[-5, -7.55], [-7.55, 5], [7.55, -5], [5, 7.55], [-5, 2.45], [5, -2.45]]
# for some reason x and y appear to be flipped
#poster_pos = reverse.(poster_pos)
poster_img = Dict(zip([:camel,:cat,:croc, :donkey,:pig,:rabbit], joinpath.(@__DIR__, "..","artefacts",  ["camel 1.png","cat 1.png","crocodile.png","donkey 1.png","pig 1.png","rabbit 1.png"])))

# TODO: Use actual values here
camera_height = 2.5
ceiling_height = 4.93
pillar_height = 3.11

struct UnityData
    time::Vector{Float64}
    position::Matrix{Float64}
    head_direction::Vector{Float64}
    triggers::Matrix{Int64}
    timestamps::Matrix{Float64}
    poster_pos::Dict{Symbol, NTuple{3,Float64}}
    header::Dict
end

DPHT.level(::Type{UnityData}) = "session"

function UnityData()
    # attempt to find data file
    # assume we are at the sesison level
    _datadir = glob("RawData[_-]T*")
    if !isempty(_datadir)
        datadir = first(_datadir)
        if isdir(datadir)
            ff = glob("session_*.txt", datadir)
            if !isempty(ff)
                return UnityData(first(ff))
            end
        end
    end
    error("No unity data found")
end

function UnityData(fname::String)
    data, header, column_names = read_unity_file(fname)
    UnityData(data, header)
end

function UnityData(data, header)
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

    UnityData(_time, data[:,[3,4]], data[:,5], triggers, timestamps, header["PosterLocations"], header)
end


"""
    parse_poster_location(loc::String)

Parse poster location from a string like "P1(-5,1.5,-7.55)"
"""
function parse_poster_location(loc::AbstractString)
    rr = r"P([\d+.])\(([-\d\.]+),([-\d.]+),([-\d\.]+)\)"
    m = match(rr, loc)
    if m !== nothing
        # first is the poster number
        j = parse(Int64, m.captures[1])
        poster_loc = tuple(parse.(Float64, m.captures[2:end])...)
        return j,poster_loc
    end
    return 0, (NaN,NaN)
end

function Base.parse(::Type{NTuple{3,Float64}},s::AbstractString)
    rr = r"\(([-\d\.]+),[ ]*([-\d\.]+),[ ]*([-\d\.]+)\)"
    m = match(rr,s)
    if m !== nothing
        tt = tuple(parse.(Float64, m.captures)...)
        return tt
    end
    throw(ArgumentError("Unable for parse NTuple{3,Float64} for $s"))
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
function read_unity_file(fname::String)
    # TODO: Return the config as well
    column_names = ["marker","Δt","xpos","ypos","direction"]
    # check if first column should be skipped
    header = Dict()
    header_rows = 0
    lines = open(fname) do fid
        _lines = String[]
        while true
            _line = readline(fid)
            if isempty(_line)
                continue
            end
            if isdigit(first(strip(_line)))
                break
            end
            push!(_lines, _line)
            header_rows += 1
        end
        _lines
    end
    _header = header
    poster_idx = 0
    for (i,_line) in enumerate(lines) 
        if occursin(':', _line)
            k,v = split(_line, ':')
            if k == "PosterLocations"
                # check the next line
                if startswith(lines[i+1],"name: Poster")
                    # different format: start 
                    header["PosterLocations"] = Dict()
                    _header = header["PosterLocations"]
                else
                    # older format
                    header["PosterLocations"] = Dict()
                    loc_strings = split(v)
                    _keys = sort(collect(keys(poster_img)),by=k->string(k))
                    header["PosterLocations"] = Dict(_keys[ll[1]]=>ll[2] for ll in parse_poster_location.(loc_strings))
                end
            elseif k == "name" || k == "posterPosition"
                if k == "name" 
                    # grab the next line which shoud be the position
                    k2,v2 = split(lines[i+1], ':')
                    _header[poster_idx] = parse(NTuple{3,Float64},v2)
                elseif k == "posterPosition"
                    # do nothing since we've already processed this line above
                end
            else
                if k == "MazeName"
                    poster_idx = get_poster_index(v)
                end
                _header = header
                _header[k] = v
            end
        end
    end
    # post processing
    if "MazeName" in keys(header)
        # figure out the poster index
        poster_idx = get_poster_index(header["MazeName"])
    end
    data = readdlm(fname, ' ', skipstart=header_rows) 
    if size(data,2) == 6 # should only be 5 columns
        data = data[:,2:end]
    end
    data = convert(Matrix{Float64}, data)
    data, header, column_names
end

function get_poster_index(mazename::AbstractString)
    pp = split(mazename)
    poster_idx = findfirst(k->occursin(lowercase(pp[end]),k),poster_img) 
end

function Makie.convert_arguments(::Type{<:AbstractPlot}, x::UnityData) 
    PlotSpec(Lines, x.position[:,1], x.position[:,2])
end

function Makie.convert_arguments(::Type{<:AbstractPlot}, x::UnityData, trial::Trial) 
    t,posx,posy = get_trial(x, trial.i)
    PlotSpec(Lines, posx, posy)
end

function angle2arrow(a::Float64)
    θ = π*a/180
    sin(θ), cos(θ)
end

function visualize!(lscene, udata::UnityData;trial::Observable{Trial}=Observable(Trial(1)), current_time::Observable{Float64}=Observable(0.0), show_maze=false,show_north=false, kwargs...)
    _ntrials = numtrials(udata)
    udata_trial = lift(trial) do _trial
        if 0 < _trial.i <= _ntrials
            tp,posx,posy,dir = get_trial(udata, _trial.i)
            position = [Point3f(px,py, 0.1) for (px,py) in zip(posx,posy)]
            return tp .- tp[1], position, dir
        else
            return Float64[], fill(0.0, 3, 0), Float64
        end
    end

    position = lift(udata_trial) do _udt
        _udt[2]
    end

    head_direction = lift(udata_trial) do _udt
        _udt[3]
    end

    current_pos = Observable(position[][1:1])
    current_head_direction = Observable([Point3f(angle2arrow(head_direction[][1])...,0.0)])

    north_direction = Observable([Point3f(0.0, 1.0,0.0)])
    # trigger curent_pos both on trial and time change
    onany(position, current_time) do _position, _ct
        # figure out which position to use
        tp = udata_trial[][1]
        j = searchsortedfirst(tp, _ct)
        if 0 < j <= length(tp)
            current_pos[] = _position[j:j]
            aa = angle2arrow(head_direction[][j])
            current_head_direction[] = [Point3f(aa...,0.0)]
        end
    end
    if show_maze
        mm = MazeModel()
        visualize!(lscene,mm;kwargs...)
    end
    lines!(lscene, position, color=:black)
    scatter!(lscene, current_pos, color=:blue)
    arrows!(lscene, current_pos, current_head_direction, color=:blue)
    if show_north
        arrows!(lscene, current_pos, north_direction, color=:red)
    end
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

"""
    soft_range(start::T, stop::T,step::T) where T <: Real

Range where the step is adjusted to match the start and stop
"""
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
    b1 = -7.5 <= x <= -2.5 && 2.5 <= y <= 7.5 && z < pillar_height
    b2 = -7.5 <= x <= -2.5 && -7.5 <= y <= -2.5 && z < pillar_height
    b3 = 2.5 <= x <= 7.5 && -7.5 <= y <= -2.5 && z < pillar_height
    b4 = 2.5 <= x <= 7.5 && 2.5 <= y <= 7.5 && z < pillar_height
    b5 = x < -12.5 || x > 12.5
    b6 = y < -12.5 || y > 12.5
    b7 = z > ceiling_height || z < 0.0
    b1 || b2 || b3 || b4 || b5 || b6 || b7
end

"""
A mesh with an attached normal
"""
struct OrientedMesh{T<:AbstractVector{<: Real}}
    bins::NTuple{3,T}
    normal::Vec3f
end

OrientedMesh(bins, normal::AbstractVector{T}) where T <: Real = OrientedMesh(bins, Vec3f(normal))

struct MazeModel{T<:AbstractVector{<: Real}}
    walls::Vector{OrientedMesh{T}}
    pillars::Vector{Vector{OrientedMesh{T}}}
    floor::OrientedMesh{T}
    ceiling::OrientedMesh{T}
end

function get_maze_colors(mm::MazeModel)
    colors = Dict{Symbol, Any}()
    colors[:walls] = fill(RGB(0.3f0, 0.21470589f0, 0.21470589f0), 4)
    colors[:floor] = to_color(:gray)
    colors[:ceiling] = to_color(:gray)
    _colors = [RGB(255/255,195/255,0.0), # yellow
               RGB(255/255,111/255,120/255),# red
               RGB(0.0, 118/255,1.0), # blue
               RGB(5/255,223/255,0.0) # green
            ]
    colors[:pillars] = [[to_color(c) for _ in  1:length(mm.pillars[i])] for (i,c) in enumerate(_colors)]
    colors
end

"""
Convert the histogram counts on each surface of `mm` to a color
"""
function get_maze_colors(mm::MazeModel, counts::Dict{Symbol, Vector{Array{T,3}}};kernel=nothing, invert_ceiling=false,kwargs...) where T <: Real
    # TODO: This can probably be done more efficiently and robustly
    base_colors = get_maze_colors(mm)
    colors = Dict{Symbol,Union{Vector{Vector{T}}, Vector{Vector{Vector{T}}}}}()
    normals = get_normals(mm)
    for k in keys(normals)
        colors[k] = Vector{Vector{T}}(undef, length(normals[k]))
        if k == :ceiling
            invert = invert_ceiling
        else
            invert = false
        end
        for (i,n) in enumerate(normals[k])
            # we want to color only the inside
            c = counts[k][i]
            _color = fill!(similar(c), 0.0)
            for d in 1:length(n)
                if n[d] == 0
                    continue
                end
                if (n[d] < 0) && (!invert)
                    idx = ntuple(dim->dim==d ? 1 : axes(c,dim), 3)
                else
                    idx = ntuple(dim->dim==d ? size(c,d) : axes(c,dim), 3)
                end
                cq = dropdims(sum(_c->isfinite(_c) ? _c : zero(_c),c,dims=d),dims=d)
                if kernel !== nothing
                    cq = imfilter(cq, kernel)
                end
                _color[idx...] .= cq
            end
            colors[k][i] = _color[:]
        end
        if k == :pillars
            # reshape
            npillars = length(mm.pillars)
            n_pillar_walls = length.(mm.pillars)
            nn = length(colors[k])
            _colors = colors[k]
            colors[k] = Vector{Vector{Vector{T}}}(undef, npillars)
            offset = 0
            for (i,j) in enumerate(n_pillar_walls)
                colors[k][i] = _colors[offset+1:offset+j]
                offset += j
            end
        end
    end
    colors
end

function MazeModel(bins::Dict{Symbol,T}, normals) where T
    walls = [OrientedMesh(m,n) for (m,n) in zip(bins[:walls], normals[:walls])]
    pillars = [[OrientedMesh(m,n) for (m,n) in zip(bins[k],normals[k])] for k in [:pillar_1, :pillar_2, :pillar_3, :pillar_4]]
    _floor = OrientedMesh(first(bins[:floor]), first(normals[:floor]))
    _ceiling = OrientedMesh(first(bins[:ceiling]), first(normals[:ceiling]))
    MazeModel(walls, pillars, _floor, _ceiling)
end

function MazeModel(;kvs...)
    bins,normals = create_maze(;kvs...)
    MazeModel(bins,normals)
end

"""
Return a dictionary containg the 3D bins for each element of the maze `mm`.
"""
function get_bins(mm::MazeModel{T}) where T <: AbstractVector{T2} where T2 <: Real
    bins = Dict{Symbol, Vector{NTuple{3,T}}}()
    bins[:walls] = [w.bins for w in mm.walls]
    bins[:pillars] = cat([[p.bins for p in pp] for pp in mm.pillars]...,dims=1)
    bins[:ceiling] = [mm.ceiling.bins]
    bins[:floor] = [mm.floor.bins]
    bins
end

"""
Return a dictionary containg the normal vector for each element of the maze `mm`.
"""
function get_normals(mm::MazeModel{T}) where T <: AbstractVector{T2} where T2 <: Real
    normals = Dict{Symbol, Vector{Vector{Float64}}}()
    normals[:walls]  = [w.normal for w in mm.walls]
    normals[:pillars] = cat([[p.normal for p in pp] for pp in mm.pillars]...,dims=1)
    normals[:ceiling] = [mm.ceiling.normal]
    normals[:floor] = [mm.floor.normal]
    normals
end

function get_physical_size(mm::MazeModel)
   # floor size 
   (xmin,xmax),(ymin,ymax),_ = extrema.(mm.floor.bins)
   _,_,(zmin,zmax) = extrema.(mm.walls[1].bins)
   Δx = xmax-xmin
   Δy = ymax-ymin
   Δz = zmax-zmin
   Δx,Δy,Δz
end

function Base.show(io::IO, mm::MazeModel)
    npillars = length(mm.pillars)
    Δx,Δy,Δz = get_physical_size(mm)
    print(io, "$(Δx) by $(Δy) by $(Δz) maze with $npillars pillars")
end

function visualize!(lscene, mm::MazeModel;color::Dict{Symbol,<:Any}=get_maze_colors(mm), offsets::Union{Nothing, Dict{Symbol, Vector{Vector{Float64}}}}=nothing, show_ceiling=false, show_floor=true, show_walls=true, show_pillars=true, show_normals=false, flip_ceiling=false, colormap=:Blues, kwargs...)
    if show_floor
        #floor
        bin = mm.floor.bins
        m = CartesianGrid(first.(bin), last.(bin);dims=length.(bin))
        #hackish
        # only do this if color[:floor] is a vector of something
        if typeof(color[:floor]) <: AbstractVector
            _color = color[:floor][1][:]
        else
            _color = color[:floor]
        end
        viz!(lscene, m, color=_color,colormap=colormap)
    end

    if show_ceiling
        if typeof(color[:ceiling]) <: AbstractVector
            _color = color[:ceiling][1][:]
        else
            _color = color[:ceiling]
        end
        bin = mm.ceiling.bins
        m = CartesianGrid(first.(bin), last.(bin);dims=length.(bin))
        if flip_ceiling
            pp = Meshes.Vec(mean.(bin))
            m = Translate(-pp...)(m)
            m = Rotate(RotX(π))(m)
            m = Translate(pp...)(m)
        end
        if offsets !== nothing && :ceiling in keys(offsets)
            m = Translate(offsets[:ceiling][1]...)(m)
        end
        viz!(lscene, m, color=_color,colormap=colormap)
    end

    #pillars
    if show_pillars
        for (oms,ccs) in zip(mm.pillars,color[:pillars])
            for (om,cc) in zip(oms,ccs)
                bin = om.bins
                m = CartesianGrid(first.(bin), last.(bin);dims=length.(bin))
                viz!(lscene, m, color=cc,colormap=colormap)
                if show_normals
                    arrows!(lscene, [Point3f(mean.(om.bins))], [Point3f(om.normal)])
                end
            end
        end
    end
    if show_walls
        #walls
        for (ii,om) in enumerate(mm.walls)
            bin = om.bins
            m = CartesianGrid(first.(bin), last.(bin);dims=length.(bin))
            viz!(lscene, m, color=color[:walls][ii],colormap=colormap)
        end
    end
end

struct Posters{T<:RGB,T2<:Integer,T3<:Point3, T4<:Point2,T5<:Vec3}
    sprite::Vector{Sprite{T, T2, T3, T4, T5}}
end

function Posters(mm::MazeModel,udata::UnityData;kvs...)
    Posters(mm, udata.header["PosterLocations"];kvs...)
end

function Posters(mm::MazeModel, _poster_pos::Dict{Symbol, NTuple{3, Float64}};kvs...)
    _poster_pos = Dict(k=>(p[1],p[3]) for (k,p) in _poster_pos)
    Posters(mm, _poster_pos;kvs...)
end

function Posters(mm::MazeModel,_poster_pos=poster_pos;z=2.5)
    wall_pillar_idx = assign_posters(mm,_poster_pos)
    wall_idx = wall_pillar_idx.pillar_wall_idx
    pillar_idx = wall_pillar_idx.pillar_idx
    rot = LinearMap(RotX(3π/2))
    images = Dict(k=>load(v) for (k,v) in poster_img)
    # hack just to figure out the type
    sp = sprite(first(images)[2], Rect2(-1.25, -2.5/1.2/2, 2.5, 2.5/1.2))
    sprites = Vector{typeof(sp)}(undef, length(_poster_pos))
    for (ii,pk) in enumerate(keys(_poster_pos))
        pp = _poster_pos[pk]
        img = images[pk]
        sp = sprite(img, Rect2(-1.25, -2.5/1.2/2, 2.5, 2.5/1.2))
        sp2 = rot(sp)
        μ = mean(sp2.points) 
        # trans is relative
        trans = LinearMap(Translation(pp[1]-μ[1],pp[2]-μ[2], z))
        nn = mm.pillars[pillar_idx[pk]][wall_idx[pk]].normal
        θ = acos(sp2.normals[1]'*nn)
        rot2 = LinearMap(RotZ(θ))
        sp3 = trans(rot2(sp2))
        sprites[ii] = sp3
    end
    Posters(sprites)
end

function show_posters(args...;kwargs...)
    fig = Figure()
    lscene = LScene(fig[1,1])
    show_posters!(lscene, args...;kwargs...)
    fig
end

function visualize!(lscene, posters::Posters;kwargs...)
    for sp3 in posters.sprite
        plot!(lscene, sp3)
    end
end


"""
    create_mesh(lower_left::T, upper_right::T,Δb::Float64) where T <: NTuple{2,Float64}

Create a mesh representing an object (pillar,wall) with the specified corners, using the specified bin width
"""
function create_mesh(lower_left::T, upper_right::T,Δb::Float64, height::Float64=5.0, Δ::Float64=0.1;flip_normals=false) where T <: NTuple{2,Float64}
    bins = Vector{NTuple{3,Vector{Float64}}}(undef, 4)
    normals = Vector{Vector{Float64}}(undef, 4)
    
    zbins = range(0.0, stop=height, length=10)

    xbins = soft_range(lower_left[1], upper_right[1], Δb)
    y0 = lower_left[2] 
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[1] = (xbins,ybins,zbins)
    normals[1] = [0.0, -1.0,0.0]

    y0 = upper_right[2] 
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[2] = (xbins,ybins,zbins)
    normals[2] = [0.0, 1.0, 0.0]

    ybins = soft_range(lower_left[2], upper_right[2],Δb)
    x0 = lower_left[1] 
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[3] = (xbins,ybins,zbins)
    normals[3] = [-1.0, 0.0, 0.0]
    x0 = upper_right[1] 
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[4] = (xbins,ybins,zbins)
    normals[4] = [1.0, 0.0, 0.0]
    if flip_normals
        normals .*= -1.0
    end
    bins, normals
end

"""
Return the meshes representing the maze
"""
function create_maze(;xmin=-12.5, xmax=12.5, ymin=xmin,ymax=xmax, Δ=0.01, height=5.0, pillar_height=3.0,kvs...)
    # unity uses 40x40 bins on the floor
    floor_bins = range(xmin, stop=xmax, length=40)
    Δb = step(floor_bins)

    bins = Dict{Symbol,Vector{NTuple{3,Vector{Float64}}}}()
    normals = Dict{Symbol,Vector{Vector{Float64}}}()

    zbins = range(0.0, stop=height, length=10)

    # pillar 1; yellow pillar
    #lower_left = (-7.5, 2.5)
    #upper_right = (-2.5, 7.5)
    pos = pillar_positions[:yellow]
    bins[:pillar_1], normals[:pillar_1] = create_mesh(pos[:lower_left], pos[:upper_right],Δb,pillar_height, Δ)

    # pillar 2;red 
    lower_left = (2.5, 2.5)
    upper_right = (7.5, 7.5)
    pos = pillar_positions[:red]
    bins[:pillar_2], normals[:pillar_2] = create_mesh(pos[:lower_left], pos[:upper_right],Δb,pillar_height,Δ)

    # pillar 3; blue
    lower_left = (-7.5, -7.5)
    upper_right = (-2.5, -2.5)
    pos = pillar_positions[:blue]
    bins[:pillar_3], normals[:pillar_3] = create_mesh(pos[:lower_left], pos[:upper_right],Δb,pillar_height,Δ)

    # pillar 4
    lower_left = (2.5, -7.5)
    upper_right = (7.5, -2.5)
    pos = pillar_positions[:green]
    bins[:pillar_4], normals[:pillar_4] = create_mesh(pos[:lower_left], pos[:upper_right],Δb,pillar_height, Δ)


    zbins = range(0.0, stop=height, length=10)
    # walls
    bins[:walls] = Vector{Vector{Float64}}(undef, 4)
    normals[:walls] = Vector{Vector{Float64}}(undef, 4)
    xbins = soft_range(xmin, xmax, Δb)
    y0 = ymin
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[:walls][1] = (xbins,ybins,zbins)
    normals[:walls][1] = [0.0, 1.0,0.0]
    y0 = ymax 
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[:walls][2] = (xbins,ybins,zbins)
    normals[:walls][2] = [0.0, -1.0, 0.0]

    ybins = soft_range(ymin, ymax, Δb)
    x0 = xmin
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[:walls][3] = (xbins, ybins, zbins)
    normals[:walls][3] = [1.0, 0.0, 0.0]
    x0 = xmax 
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
    xbins = range(xmin, stop=xmax, step=Δb)
    ybins = range(ymin, stop=ymax, step=Δb)
    z0 = height 
    zbins = range(z0-Δ, stop=z0+Δ, length=2)
    bins[:ceiling] = [(xbins, ybins, zbins)]
    normals[:ceiling]=[[0.0, 0.0, -1.0]]

    bins, normals
end

function assign_posters_old(bins, normals)
    wall_idx = Vector{Tuple{Symbol,Int64}}(undef, length(poster_pos))
    for (i,pp) in enumerate(poster_pos)
        d = Inf
        for pillar in [:pillar_1, :pillar_2, :pillar_3, :pillar_4]
            for (j,pq) in enumerate(bins[pillar])
                _d = (pp[1] - mean(pq[1]))^2 + (pp[2] - mean(pq[2]))^2
                if _d < d
                    d = _d
                    wall_idx[i] = (pillar,j)
                end
            end
        end
    end
    wall_idx
end

function assign_posters(mm::MazeModel, _poster_pos::Dict{Symbol,NTuple{3,Float64}})
    # ignore the second coordinate
    assign_posters(mm, Dict(k=>(p[1],p[3]) for (k,p) in _poster_pos))
end

function assign_posters(mm::MazeModel, _poster_pos::Dict{Symbol,NTuple{2,Float64}}=poster_pos)
    pillar_idx = Dict{Symbol,Int64}()
    wall_idx = Dict{Symbol,Int64}()
    for (kp,pp) in _poster_pos
        d = Inf
        for (k,pillar) in enumerate(mm.pillars)
            for (j,_wall) in enumerate(pillar)
                pq = _wall.bins
                _d = (pp[1] - mean(pq[1]))^2 + (pp[2] - mean(pq[2]))^2
                if _d < d
                    d = _d
                    pillar_idx[kp] = k
                    wall_idx[kp] = j
                end
            end
        end
    end
    (pillar_idx=pillar_idx, pillar_wall_idx=wall_idx)
end

# TODO: The Unity raytracer uses 40x40 bins on the floor as the baseline
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

function compute_histogram(pos::Vector{Matrix{Float64}},bins,weight::Union{Nothing,Vector{Vector{Float64}}}=nothing)
    counts = [fill(0.0, length.(bin)) for bin in bins]
    if weight === nothing
        weight = [fill(1.0, size(_pos,2)) for _pos in pos]
    end
    for (ii,_pos) in enumerate(pos)
        compute_histogram!(counts, _pos, bins;weight=weight[ii])
    end
    counts
end

fstep(aa::Vector{T}) where T <: Real = mean(diff(aa))
fstep(aa::StepRangeLen) = step(aa)

function compute_histogram!(counts, pos::Matrix{Float64},bins;weight=fill(1.0, size(pos,2)))
    qpos = ([pos[i,:] for i in 1:size(pos,1)]...,)
    w = aweights(weight)
    for (bin,count) in zip(bins,counts)
        # hackish; add one bin to the end
        Δs = [fstep(b) for b in bin]
        h = fit(Histogram, qpos, w, ([[b;b[end]+Δ] for (b,Δ) in zip(bin,Δs)]...,))
        count .+= h.weights
    end
    counts
end

function compute_histogram(pos::Matrix{Float64}, bins::NTuple{N, T}) where T <: AbstractVector{T2} where T2 <: Real where N
    counts = fill(0.0, length.(bins))
    compute_histogram!(counts, pos, bins)
    counts
end

function compute_histogram!(counts::Array{Float64,3}, pos::Matrix{Float64}, bins::NTuple{N, T}) where T <: AbstractVector{T2} where T2 <: Real where N
    qpos = eachcol(pos)
    Δs = [step(b) for b in bins]
    h = fit(Histogram, qpos, ([[b;b[end]+Δ] for (b,Δ) in zip(bin,Δs)]...,))
    count .+= h.weights
end

