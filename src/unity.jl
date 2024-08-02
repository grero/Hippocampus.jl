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
    fig,ax = poly(Point2f.(zip(zBound, xBound)),color=:grey)
    ax.aspect = 1.0 
    poly!(ax, Point2f.(zip(z1Bound, x1Bound)), color=:yellow)
    poly!(ax, Point2f.(zip(z2Bound, x2Bound)), color=:red)
    poly!(ax, Point2f.(zip(z3Bound, x3Bound)), color=:blue)
    poly!(ax, Point2f.(zip(z4Bound, x4Bound)), color=:green)
    fig,ax
end

function impacts(pos)
    x,y = pos
    b1 = -7.5 <= x <= -2.5 && 2.5 <= y <= 7.5
    b2 = -7.5 <= x <= -2.5 && -7.5 <= y <= -2.5
    b3 = 2.5 <= x <= 7.5 && -7.5 <= y <= -2.5
    b4 = 2.5 <= x <= 7.5 && 2.5 <= y <= 7.5
    b5 = x < -12.5 || x > 12.5
    b6 = y < -12.5 || y > 12.5
    b1 || b2 || b3 || b4 || b5 || b6
end