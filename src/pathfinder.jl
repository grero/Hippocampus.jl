using GLMakie
using DifferentialEquations
using LinearAlgebra

function wall_potential(x,x0, x1;β=5.0)
    (1/(1+exp(β*(x-x0))))+(1/(1+exp(-β*(x-x1))))
end

function pillar_potential(x,y,x0, x1, y0, y1;β=5.0)
    (1/(1+exp(-β*(x-x0))))*(1/(1+exp(-β*(y-y0))))*(1/(1+exp(β*(x-x1))))*(1/(1+exp(β*(y-y1))))
end

function potential(x,y;β=5.0)
     z1 = pillar_potential(x,y, 2.5, 7.5, 2.5, 7.5;β=β)
     z2 = pillar_potential(x,y,-7.5, -2.5, 2.5, 7.5;β=β)
     z3 = pillar_potential(x,y,-7.5, -2.5, -7.5, -2.5;β=β)
     z4 = pillar_potential(x,y, 2.5, 7.5, -7.5, -2.5;β=β)
     wp_x = wall_potential(x,-12.5, 12.5;β=β)
     wp_y = wall_potential(y,-12.5, 12.5;β=β)
     z1+z2+z3+z4 + wp_x + wp_y
end

function pillar_potential_grad(x,y, x0, x1, y0, y1;β=5.0)
    A1 = (1/(1+exp(-β*(x-x0))))
    A2 = (1/(1+exp(β*(x-x1))))
    B1 = (1/(1+exp(-β*(y-y0))))
    B2 = (1/(1+exp(β*(y-y1))))

    [+β*A1*(1-A1)*A2*B1*B2 - β*A1*A2*(1-A2)*B1*B2,
    +β*A1*A2*B1*(1-B1)*B2-β*A1*A2*B1*B2*(1-B2)]
end

function wall_potential_grad(x, x0, x1;β=5.0)
    A1 = (1/(1+exp(β*(x-x0))))
    A2 = (1/(1+exp(-β*(x-x1))))
    -β*A1*(1-A1)+A2 + β*A1+A2*(1-A2)
end

function potential_grad(x,y;β=5.0)
    d1 = pillar_potential_grad(x,y, 2.5, 7.5, 2.5, 7.5;β=β)
    d2 = pillar_potential_grad(x,y, -7.5, -2.5, 2.5, 7.5;β=β)
    d3 = pillar_potential_grad(x,y, -7.5, -2.5, -7.5, -2.5;β=β)
    d4 = pillar_potential_grad(x,y, 2.5, 7.5, -7.5, -2.5;β=β)
    wp_grad_x = wall_potential_grad(x,-12.5, 12.5;β=β)
    wp_grad_y = wall_potential_grad(y,-12.5, 12.5;β=β)
    d1+d2+d3+d4 + [wp_grad_x, wp_grad_y]
end

function check_boundaries(xy,v,dt)
    xy_new = xy+v*dt
    (x_new,y_new) = xy_new
    if x_new < -12.5 || x_new > 12.5
        v = Vec2f(0) # reflect
    end
    if y_new < -12.5 || y_new > 12.5
        v = Vec2f(0)
    end
    # check pillars
    for x0 in [-7.5, 2.5]
        for y0 in [-7.5, 2.5]
            if (x0 < x_new < x0+5.0) && (y0 < y_new < y0 + 5.0)
                v = Vec2f(0)
            end
        end
    end
    v
end

function integrate!(start::Observable{Point2f}, goal,v::Observable{Vec2f},g::Observable{Vec2f},goal_vector::Observable{Vec2f},avo::Observable{Vec2f};dt=0.01,a1=0.1, a2=0.5, a3=1.0,k=0.01,k0=2.0, β=5.0,kv=0.01, kτ=0.1, kω=2*sqrt(kτ), vmax=1.0, interactive=true)
    xy = start[]
    pth = Point2f[]
    ω = 0.0 #angular speed
    up = [0.0, 0.0, 1.0]
    Ry = [0.0 -1.0; 1.0 0.0]
    while norm(xy-goal) > 0.5
        _v = v[]
        # compute current direction to goal
        vn = _v/(norm(_v) + 1e-6)
        # orthogonal to v
        vq = Ry*vn
        vn3 = Vec3f(vn...,0.0)
        vq3 = Vec3f(vq...,0.0)

        dg = goal - xy
        goal_vector[] = dg
        # normalize
        ndg = norm(dg)
        ng = dg/(ndg+1.e-6)
        dg3 = Vec3f(ng...,0.0)
        # project current velocity onto vector 
        dg_parallel = (vn'*ng)*vn
        dg_normal = ng - dg_parallel
        # get the orthogonal component
        #cosθ = clamp(dg_normal'*vn/(norm(dg_normal) + 1e-6), -1.0, 1.0)
        # stering acceleration
        # actually, I think it makes sense to separate steering from acceleation here
        # to steer, just apply a torque porportional to the angle
        #τ = -kτ*dg_normal
        τ = sign(vn3'*dg3)*kτ*cross(vn3,dg3)'*up
        # apply acceleration along the current speed,
        # but proportional to the spring force towards the goal projected onto v
        # should try and counter the ortogonal component, 
        # need to do some kind of PID here
        gg = -potential_grad(xy...;β=β)
        gg3 = Vec3f(gg...,0.0)
        gg_parallel = (gg'*vn)*vn
        # get the component orthogonal to the velocity
        gg_normal = gg - gg_parallel 
        g[] = k0*gg_normal
        gg_normal_nn = gg_normal/(norm(gg_normal)+1e-6)

        #τ += -k0*norm(gg_normal)* - kω*ω
        #τ += -cross(vn, gg)
        τ += k0*sign(vq3'*gg3)*cross(vq3, gg3)'*up - kω*ω
        #τ = clamp(τ, -1.0, 1.0)
        ω += τ*dt
        #ω = clamp(ω, -2.0, 2.0) 
        R = [cos(ω*dt) -sin(ω*dt);sin(ω*dt) cos(ω*dt)]
        #av = k*ndg*ng - 0.8*2*sqrt(k)*v
        av = k*ndg*(ng'*vn)*vn #+ #0.1*k0*(-gg'*_v)*vn
        avn = norm(av)
        av = (av./avn)*(min(1.0, avn))
        # TODO: Add a force that slows us down as we approach the boundary
        #       This should be a force that only works along the component
        #       anti-parallel to the negative gradient
        aq = gg'*_v
        if aq < 0
           av += aq*(gg/norm(gg)) 
        end
        avo[] = av
        #hack 
        if norm(_v) == 0
            av = 0.1*k*ng
        end
        _v = _v+av*dt
        vq = max(norm(_v), 1e-6)
        _v = (_v./vq).*min(vmax, vq)
        # rotate v
        _v = R*_v
        # check actual collision
        _v = check_boundaries(xy,_v,dt)
        xy = xy + _v*dt
        start[] = xy
        v[] = _v
        if interactive
            sleep(dt)
            yield()
        end
        push!(pth, xy)
    end
    pth
end

function setup_slowness_field(;xrange=range(-12.5, stop=12.5, length=40), yrange=xrange)
    ϕ = ones(length(xrange), length(yrange))
    # slowness is zero in the interior, 1 inside the pillars
     for x0 in [-8.0, 2.0]
        idx_x0 = searchsortedfirst(xrange, x0)
        idx_x1 = searchsortedfirst(xrange, x0+6.0)
        for y0 = [-8.0, 2.0]
            idx_y0 = searchsortedfirst(yrange, y0)
            idx_y1 = searchsortedfirst(xrange, y0+6.0)
            ϕ[idx_x0:idx_x1, idx_y0:idx_y1] .= 100 
        end
    end
    for x0 in [-7.5, 2.5]
        idx_x0 = searchsortedfirst(xrange, x0)
        idx_x1 = searchsortedfirst(xrange, x0+5.0)
        for y0 = [-7.5, 2.5]
            idx_y0 = searchsortedfirst(yrange, y0)
            idx_y1 = searchsortedfirst(xrange, y0+5.0)
            ϕ[idx_x0:idx_x1, idx_y0:idx_y1] .= Inf 
        end
    end
    ϕ
end

function geodesic_equation!(ddu, du,u, p,t)
    x,y = u
    vx,vy = du
    β = first(p)
    f = 1 + potential(x,y;kwargs...)
    (df_dx,df_dy) = k*potential_grad(x,y;β=β)
    γ_xxx = df_dx/(2*f)
    γ_xxy = df_dy/(2*f)
    γ_xyy = -df_dx/(2*f)
    γ_yyy = df_dy/(2*f)
    γ_xyy = df_dx/(2*f)
    γ_yxx = df_dy/(2*f)

    d2x_dt2 = γ_xxx*vx^2 + 2*γ_xxy*vx*vy + γ_xyy*vy^2
    d2y_dt2 = γ_yxx*vx^2 + 2*γ_xyy*vx*vy + γ_yyy*vy^2
    du .= [vx,vy]
    ddu .= [d2x_dt2,d2y_dt2]
    nothing
end

function solve_geodesic_equation(start, goal;kwargs...)
    v0 = goal - start
    v0 = v0/norm(v0)
    y0 = [start[1], start[2], v0[1],v0[2]]
    ODE.SecondOrderODEProblem(geodesic_equation!, v0, start, tspan, ω)
end

# Define the metric tensor g_μν at a point (x, y)
function metric_tensor(x, y;kwargs...)
    f = 1 .+ potential(x,y;kwargs...)
    return [f 0.0; 0.0 f]  # Diagonal metric tensor
end

function geodesic_energy(path;kwargs...)
    E = 0.0
    for i in 1:length(path)-1
        # Current and next points
        γ_i = path[i]
        γ_i1 = path[i+1]

        # Difference vector
        Δγ = γ_i1 - γ_i

        # Metric tensor at γ_i (midpoint rule could also be used)
        g = metric_tensor(γ_i[1], γ_i[2])

        # Compute the Riemannian distance: sqrt(Δγ^T * g * Δγ)
        distance = sqrt(Δγ' * g * Δγ)
        E += distance
    end
    return E
end

function relax_path!(path::Vector{Point2f};niter=1000, λ=0.01, ϵ=0.01, kwargs...)
    ee = zeros(niter)
    for j in 1:niter
        for i in 2:length(path)-1
            # tangent vector
            v = path[i+1] - path[i-1]
            g_smooth = 2 * (2 * path[i] - path[i-1] - path[i+1])
            v = v/norm(v)
            # find normal vector
            g = potential_grad(path[i]...,;kwargs...)
            g_total = g + λ*g_smooth
            # move point a small distance in the negative gradient direction
            g_normalized = g_total / (norm(g_total) + 1e-6)
            gv = -ϵ*(g_normalized - (g_normalized'*v)*v)
            path[i] = path[i] + Vec2f(gv...)
            # compute energy
            f = 1 .+ potential(path[i]...;kwargs...) 
            ds = path[i] - path[i-1]
            ee[j] += f*ds'*ds
        end
    end
    path,ee
end
