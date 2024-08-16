using Makie
using Colors
using Meshes

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

function show_maze(bins,normals,counts::Union{Vector{Array{T,3}},Nothing}=nothing;explore=false, position::Union{Nothing, Matrix{Float64}}=nothing, head_direction::Union{Nothing, Vector{Float64}}=nothing) where T <: Real
    fig = Figure()
    #ax = Axis3(fig[1,1],aspect=:data)
    lscene = LScene(fig[1,1], show_axis=false)
    for (i,(bin,n)) in enumerate(zip(bins,normals))
        m = CartesianGrid(first.(bin), last.(bin);dims=length.(bin))
        if counts !== nothing
            # we want to color only the inside
            c = counts[i]
            _color = fill!(similar(c), 0.0)
            for d in 1:length(n)
                if n[d] < 0
                    idx = ntuple(dim->dim==d ? 1 : axes(c,dim), 3)
                    _color[idx...] .= dropdims(sum(c,dims=d),dims=d)
                elseif n[d] > 0
                    idx = ntuple(dim->dim==d ? size(c,d) : axes(c,dim), 3)
                    _color[idx...] .= dropdims(sum(c,dims=d),dims=d)
                end
            end
            _color = _color[:]
        else
            _color = RGB(0.8, 0.8, 0.8) 
        end
        viz!(lscene, m, color=_color,colormap=:Blues)
    end

    lookat = Point3f(1.0, 0.0, 2.5)
    if position !== nothing
        # replay experiment with the supplied position
        #cc = Makie.Camera3D(lscene.scene, projectiontype = Makie.Perspective, rotation_center=:eyeposition, center=false)
        cc = cameracontrols(lscene.scene)
        ii = Observable(1)

        on(ii) do i
            pos = Point3f(position[1,i], position[2,i], 2.5)
            θ = π*head_direction[i]/180
            cc.lookat[] = Point3f(cos(θ), sin(θ), 0.0) + pos
            cc.eyeposition[] = pos
            update_cam!(lscene.scene, cc)
        end
        θ = π*head_direction[1]/180
        _pos  = Point3f(position[1,1], position[2,1], 2.5)
        _lookat = Point3f(cos(θ),sin(θ), 0.0) + _pos
        update_cam!(lscene.scene, _pos, _lookat)

        @async for j in 1:size(position,2)
            ii[] = j
            yield()
            sleep(0.01)
        end
    elseif explore
        # set up camera inside of the maze
        cc = Makie.Camera3D(lscene.scene, projectiontype = Makie.Perspective, rotation_center=:eyeposition, center=false)
        #cc.eyeposition[] = Point3f(0.0, 0.0, 2.5)
        eyepos = Point3f(0.0, 0.0, 2.5)
        v = lookat - eyepos 
        v = v./norm(v)
        #translate_cam!(lscene.scene, cc, Point3f(0.0, 0.0,2.5))
        update_cam!(lscene.scene, eyepos, lookat)
        on(events(lscene.scene).keyboardbutton, priority=20) do event
            if ispressed(lscene.scene, Keyboard.up)
                pos = cc.eyeposition[]
                dx = 0.1*v
                npos = pos + dx
                translate_cam!(lscene.scene, cc, Point3f(0.0, 0.0, -0.1))
                if impacts(cc.eyeposition[])
                    # move back
                    # TODO: This doesn't quite work, but maybe we don't care
                    translate_cam!(lscene.scene, cc, Point3f(0.0, 0.0, 0.1))
                    # last coordinate if foward movement (for some inexplicable reason))
                end
            end
            if ispressed(lscene.scene, Keyboard.right)
                rotate_cam!(lscene.scene, cc, Point3f(0.0, -0.1, 0.0))
                return Consume()
            end
            if ispressed(lscene.scene, Keyboard.left)
                rotate_cam!(lscene.scene,cc, Point3f(0.0, 0.1, 0.0))
                return Consume()
            end
        end
    end
    fig
end