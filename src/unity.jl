using Makie
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
    bins = Matrix{NTuple{3,AbstractVector{Float64}}}(undef, 4,4)
    wall_bins = Vector{NTuple{3,AbstractVector{Float64}}}(undef, 4)
    zbins = range(0.0, stop=5.0, length=10)
    Δ = get(kvs, :Δz, 0.1)
    # pillar 1
    xbins = range(-7.5, stop=-2.5, length=10)
    y0 = 2.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[1,1] = (xbins, ybins, zbins)
    y0 = 7.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[2,1] = (xbins, ybins, zbins)
    ybins = range(2.5, stop=7.5, length=10)
    x0 = -7.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[3,1] = (xbins, ybins, zbins)
    x0 = -2.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[4,1] = (xbins, ybins, zbins)

    # pillar 2
    xbins = range(2.5, stop=7.5, length=10)
    y0 = 2.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[1,2] = (xbins,ybins, zbins)
    y0 = 7.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[2,2] = (xbins, ybins, zbins)

    ybins = range(2.5, stop=7.5,length=10)
    x0 = 2.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[3,2] = (xbins, ybins, zbins)
    x0 = 7.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[4,2] = (xbins, ybins, zbins)

    # pillar 3
    xbins = range(2.5, stop=7.5, length=10)
    y0 = -2.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[1,3] = (xbins, ybins, zbins)
    y0 = -7.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[2,3] = (xbins, ybins, zbins)
    ybins = range(-7.5, stop=-2.5,length=10)
    x0 = 2.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[3,3] = (xbins, ybins, zbins)
    x0 = 7.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[4,3] = (xbins, ybins, zbins)

    # pillar 4
    xbins = range(-7.5, stop=-2.5, length=10)
    y0 = -2.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[1,4] = (xbins, ybins, zbins)
    y0 = -7.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    bins[2,4] = (xbins, ybins, zbins)
    ybins = range(-7.5, stop=-2.5,length=10)
    x0 = -2.5
    xbins = range(x0-Δ, stop=x0+Δ,length=3)
    bins[3,4] = (xbins, ybins, zbins)
    x0 = -7.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    bins[4,4] = (xbins, ybins, zbins)

    # walls
    xbins = range(-12.5, stop=12.5, length=10)
    y0 = -12.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    wall_bins[1] = (xbins,ybins,zbins)
    y0 = 12.5
    ybins = range(y0-Δ, stop=y0+Δ,length=2)
    wall_bins[2] = (xbins,ybins,zbins)

    ybins = range(-12.5, stop=12.5, length=10)
    x0 = -12.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    wall_bins[3] = (xbins, ybins, zbins)
    x0 = 12.5
    xbins = range(x0-Δ, stop=x0+Δ,length=2)
    wall_bins[4] = (xbins, ybins, zbins)
   (pillars=bins, walls=wall_bins)
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

# TODO: Currently the counts are sometimes shown on the wrong side of the wall
function show_maze(bins,counts)
    fig = Figure()
    ax = Axis3(fig[1,1])
    for (c,bin) in zip(counts,bins)
        m = CartesianGrid(first.(bin), last.(bin);dims=length.(bin))
        viz!(ax, m, color=c[:],colormap=:Blues)
    end
    fig
end