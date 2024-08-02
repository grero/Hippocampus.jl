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

function plot_arena()
    fig,ax = poly(Point2f.(zip(xBound, zBound)),color=:grey)
    poly!(ax, Point2f.(zip(x1Bound, z1Bound)), color=:yellow)
    poly!(ax, Point2f.(zip(x2Bound, z2Bound)), color=:red)
    poly!(ax, Point2f.(zip(x3Bound, z3Bound)), color=:blue)
    poly!(ax, Point2f.(zip(x4Bound, z4Bound)), color=:green)
    fig
end