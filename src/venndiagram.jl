using Optim

function intersection_area(r1, r2, d)
    a1 = r1^2*acos((d^2+r1^2-r2^2)/(2*d*r1))
    a2 = r2^2*acos((d^2+r2^2-r1^2)/(2*d*r2))
    a3 = 0.5*sqrt.((-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2))
    a1+a2-a3
end

function find_d(r1, r2, A)
    func(d) = abs2(A-intersection_area(r1,r2,first(d)))

    q = optimize(func, [0.3*r1])
end

function plot_venn!(ax, r1,r2,d;Delta=0.0)
    x = (r1^2-r2^2 + d^2)/(2*d)
    @show x
    y1 = -sqrt(r1^2-x^2)
    y2 = sqrt(r1^2-x^2)
    theta1 = asin(y1/r1)
    theta2 = asin(abs(y2)/r1)
    @show x d
    if x > d
        phi1 = asin(abs(y2)/r2)
        phi2 = 2pi - asin(abs(y1)/r2)
    else
        phi1 = pi - asin(abs(y1)/r2)
        phi2 = pi + asin(abs(y2)/r2)
    end
    # left wedge
    bp1 = BezierPath([MoveTo(Point(x,y2)),
                      EllipticalArc(Point(0.0,0.0), r1,r1,0,theta2, 2pi+theta1),
                      EllipticalArc(Point(d, 0.0), r2,r2,0,phi2,phi1),
                      ClosePath()])
    pp1 = poly!(ax, bp1)
    rr1 = sqrt((r1+0.5*Delta)^2/r1^2)
    rr2 = sqrt((r2+0.5*Delta)^2/r2^2)
    translate!(pp1, -0.5*Delta, 0, 0)
    scale!(Accum, pp1, rr1,rr1,1.0)

    # central wedge
    bp = BezierPath([MoveTo(Point(x, y2)),
                    EllipticalArc(Point(0.0, 0.0), r1,r1,0,theta1,theta2),
                    EllipticalArc(Point(d, 0.0),r2,r2,0,phi2,phi1),
                    ClosePath()]
    )
    poly!(ax, bp)

    bp2 = BezierPath([MoveTo(Point(x,y2)),
                      EllipticalArc(Point(0.0,0.0), r1, r1, 0, theta2, theta1),
                      EllipticalArc(Point(d,0.0), r2, r2, 0, phi2, 2pi+phi1),
                      ClosePath()])
    pp2 = poly!(ax, bp2)
    translate!(pp2, 0.5*Delta, 0, 0)

end

function plot_venn(r1,r2,d;kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1],aspect=1)
    plot_venn!(ax, r1, r2, d;kwargs...)
    hidedecorations!(ax)
    display(fig)
    fig,ax
end