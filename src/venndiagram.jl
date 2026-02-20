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
    theta2 = asin(y2/r1)
    @show x d
    if x > d
        phi1 = asin(abs(y2)/r2)
        phi2 = 2pi - asin(abs(y1)/r2)
    else
        phi1 = pi - asin(abs(y1)/r2)
        phi2 = pi + asin(abs(y2)/r2)
    end
    #arc!(ax, Point2f(0.0, 0.0), r1, 0, 2pi;color=:black)
    # left wedge
    arc!(ax, Point2f(-Delta, 0.0), r1+Delta, theta2, 2pi+theta1, color=:green)
    if Delta > 0
       arc!(ax, Point2f(d+Delta, 0.0), r2, phi1,phi2, color=:green)
    end

    #arc!(ax, Point2f(d-Delta, 0.0), r2+Delta, phi1, phi2, color=:green)

    # central wedge
    arc!(ax, Point2f(0.0, 0.0), r1, theta1, theta2, color=:red)
    arc!(ax, Point2f(d, 0.0), r2, phi1, phi2, color=:red)

    #arc!(ax, Point2f(d,0.0), r2, 0, 2pi;color=:black)
    arc!(ax, Point2f(d+Delta, 0.0), r2+Delta, phi2, 2pi+phi1, color=:blue)
    if Delta > 0 
        arc!(ax, Point2f(Delta, 0.0), r1+Delta, theta1, theta2, color=:blue)
    end

end

function plot_venn(r1,r2,d;kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1])
    plot_venn!(ax, r1, r2, d;kwargs...)
    display(fig)
    fig,ax
end