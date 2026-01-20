using Random
using Distributions
using Distances

"""
Simulate a simple place field neuron using behavioural data in `udata`
"""
function model_place_field(udata::UnityData, rpdata::RippleData;λmin=0.1, λmax=3.0,dt=0.01,σ=1.0,μ=[3.5,1.2],rng=Random.default_rng())
    Σ = [σ 0.0;0.0 σ]
    G = MvNormal(μ, Σ)
    G0 = pdf(G, μ)
    nt = numtrials(udata)

    spikes = Float64[]
    q = -log(rand(rng))
    r = 0.0
    for i in 1:nt
        t,mposx,mposy = get_trial(udata, i;trial_start=2)
        # get the time from ripple
        trp = rpdata.timestamps[i,2]
        # reference to ripple time
        t .= t .- t[1] .+ trp[1]
        Δt = t[2]-t[1]
        for (t0,posx,posy) in zip(t,mposx, mposy)
            λ = pdf(G, [posx,posy]) # firing rate based on place field
            # scale firing rate
            λ = λmax*λ/G0 + λmin
            _t = t0
            while _t < t0+Δt
                r += λ*dt
                if r >= q
                    push!(spikes, _t)
                    r = 0.0
                    q = -log(rand(rng))
                end
                _t += dt
            end
        end
    end
    spikes
end

"""
Model population response
"""
function population_place_model(udata::UnityData, rp::RippleData;kwargs...)
    xbins = range(-12.5f0, stop=12.5f0, length=40)
    ybins = range(-12.5f0, stop=12.5f0, length=40)
    μ1 = [-5.0, 0.0]
    μ2 = [0.0, 5.0]
    μ3 = [-5.0, 0.0]
    μ4 = [0.0, -5.0]
    μ5 = [0.0, 0.0]
    μ6 = [-9.0, 9.0]
    μ7 = [-9.0, 9.0]
    μ8 = [9.0, 9.0]
    μ9 = [9.0, -9.0]
    # create spike trains
    spq = [model_place_field(udata, rp;μ=μ,kwargs...) for μ in [μ1, μ2, μ3, μ4, μ5, μ6, μ7, μ8, μ9]]
    spoc = SpatialOccupancy(udata, xbins, ybins)
    #create spatial representations
    spr = [SpatialRepresentation(sp,rp, udata) for sp in spq]
    # create spatial maps
    spm = [SpatialMap(_spr, spoc) for _spr in spr]
    # create smoothed spatial maps
    spm_adsm = [adaptive_smoothing(_spm) for _spm in spm]
    spm_adsm, spm, spq
end

function model_view_field(gdata::GazeOnMaze, rpdata::RippleData;λmin=0.1, λmax=3.0,dt=0.01,σ=[1.0, 0.01, 1.0],μ=[-5.0, 2.5, 1.5],rng=Random.default_rng())
    Σ = Diagonal(σ)
    G = MvNormal(μ, Σ)
    G0 = pdf(G, μ)
    nt = numtrials(gdata)

    spikes = Float64[]

    for i in 1:nt
        q = -log(rand(rng))
        r = 0.0
        t,gaze,fixmask = get_trial(gdata, i;trial_start=1)
        # get the time from ripple
        trp = rpdata.timestamps[i,1]
        # reference to ripple time
        t .= t .- t[1] .+ trp[1]
        Δt = t[2]-t[1]
        for (t0,pos) in zip(t,eachcol(gaze))
            λ = pdf(G, pos) # firing rate based on view field
            # scale firing rate
            λ = λmax*λ/G0 + λmin
            _t = t0
            while _t < t0+Δt
                r += λ*dt
                if r >= q
                    push!(spikes, _t)
                    r = 0.0
                    q = -log(rand(rng))
                end
                _t += dt
            end
        end
    end
    spikes
end


function view_field2(μ::Vector{T}, Σ::AbstractMatrix{T}, θp::T) where T <: Real
    G = MvNormal(μ, Σ)
    # internal representation must be a combination of place, view and head direction?
    # or an allocentric representation where the viewed location in actual 3D space is remembered
    function vf(x::Vector{T}, θ::T)
        #shift by θ 
        Δ = θp - θ
    end
end

"""
    view_field(x,y,θ)

Return the 2D guassian for the view field at position (x,y) facing in the direction θ
"""
function view_field(μp::Vector{T}, Σp::AbstractMatrix{T}, θp::T, σ::T,r::T, μv::AbstractVector{T}, Σv::AbstractMatrix{T}, rvp::T, rva::T) where T <: Real
    # Circular gaussian prior on θ
    # independent; just the product between a spatial guassian for 
    # x and y and a circular gaussin from θ
    # more complicated dependence?
    Σt = zeros(T, 5,5)
    μt = [μp; zero(T);μv]
    Σt[1:2,1:2] .= Σp
    Σt[3,3] = σ
    Σt[1,3] = sqrt(Σp[1,1])*σ*r
    Σt[2,3] = sqrt(Σp[2,2])*σ*r
    Σt[3,1] = Σt[1,3]
    Σt[3,2] = Σt[2,3]
    Σt[4:5, 4:5] .= Σv
    Σt[1:2, 4:5] .= sqrt.(Σp).*sqrt.(Σv)*rvp
    Σt[4:5,1:2] .= Σt[1:2,4:5]
    Σt[4,3] = sqrt.(Σv[1,1])*rva*σ
    Σt[3,4] = Σt[4,3]
    pxyθ = MvNormal(μt, Σt)
    function vf(x::Vector{T}, θ::T,xv::Vector{T})
        # probably an inefficient way of doing modulo 2π
        Δ = θ-θp
        Δ = sign(Δ)*acos(cos(θ-θp))
        pdf(pxyθ, [x;Δ;xv])
    end
end

function camera_project(x,y,z,θh, fov, z_near)
    # normal vector along the direction of θh
    v = [cos(θh), sin(θh)]

end

function embed(X::Matrix{T}, Xp::Matrix{T},Δ::T) where T <: Real
    d = pairwise(Euclidean(), Xp, X, dims=2)
    mm, midx = findmin(d, dims=1)
    idx = findall(mm[:] .<= Δ)
    qidx = [ii.I[1] for ii in midx[idx]]
    X[:,qidx], qidx
end

function model_view_field(rdata::RippleData, gdata::UnityRaytraceData;μ=[0.5, 0.5], Σ=diagm([1.0, 1.0]), λ_min=0.1, λ_max=3.0, dt=0.01, rng=Random.default_rng())
    G = Normal(μ, Σ)
    nt = numtrials(gdata)

    spikes = Float64[]

    for i in 1:nt
        q = -log(rand(rng))
        r = 0.0
        tg,gaze,fixmas = get_trial(gdata, i;trial_start=1)
        tg = gdata.timestamps[i]
        gaze = gdata.gaze[i]
        pos = gdata.pos[i]
        fixmask = gdata.fixating[i]
        # get the time from ripple
        trp = rdata.timestamps[i,1]
        # reference to ripple time
        t .= t .- t[1] .+ trp[1]
        Δt = t[2]-t[1]
        for (t0,pos,gpos) in zip(t,eachcol(pos), eachcol(gaze))
            # need to project gpos onto the view plane
            # just use a normalized view
            # actually, we could just as well use the eyelink data here, after normalzing
            # that is infact what the view field looks like.
            # generate probability from place field
            # project the gaze
            z = pdf(G, pos)
            λ = pdf(G, pos) # firing rate based on view field
            # scale firing rate
            λ = λmax*λ/G0 + λmin
            _t = t0
            while _t < t0+Δt
                r += λ*dt
                if r >= q
                    push!(spikes, _t)
                    r = 0.0
                    q = -log(rand(rng))
                end
                _t += dt
            end
        end
    end
    spikes
end

function Makie.convert_arguments(::Type{<:AbstractPlot}, cam::DummyCam)
    fwidth = tan(cam.fov/2)*cam.z_near
    fheight = fwidth/cam.frustrum_ratio
    b = [cam.dir nullspace(permutedims(cam.dir))]
    zpos = cam.pos - cam.z_near*b[:,1]
    near_plane = [zpos-fwidth*b[:,2]-fheight*b[:,3],
                  zpos-fwidth*b[:,2]+fheight*b[:,3],
                  zpos+fwidth*b[:,2]+fheight*b[:,3],
                  zpos+fwidth*b[:,2]-fheight*b[:,3]]

    S.LineSegments([cam.pos=>cam.pos+cam.dir, near_plane[1]=>near_plane[2], near_plane[2]=>near_plane[3], near_plane[3]=>near_plane[4],near_plane[4]=>near_plane[1]])
end

function plot_camera!(ax, cam::DummyCam)
    fwidth = tan(cam.fov/2)*cam.z_near
    fheight = fwidth/cam.frustrum_ratio
    b = [cam.dir nullspace(permutedims(cam.dir))]
    zpos = cam.pos - cam.z_near*b[:,1]
    near_plane = [zpos-fwidth*b[:,2]-fheight*b[:,3],
                  zpos-fwidth*b[:,2]+fheight*b[:,3],
                  zpos+fwidth*b[:,2]+fheight*b[:,3],
                  zpos+fwidth*b[:,2]-fheight*b[:,3]]

    linesegments!(ax, [cam.pos=>cam.pos+cam.dir, near_plane[1]=>near_plane[2], near_plane[2]=>near_plane[3], near_plane[3]=>near_plane[4],near_plane[4]=>near_plane[1]])
end

function visualize_camera(cam::DummyCam, pos::Vector{Float64})
    fwidth = tan(cam.fov/2)*cam.z_near
    fheight = fwidth/cam.frustrum_ratio
    fig = Figure()
    ax1 = Axis3(fig[1,1])
    ax2 = Axis(fig[1,2], aspect=cam.frustrum_ratio)

    b = [cam.dir nullspace(permutedims(cam.dir))]
    zpos = cam.pos - cam.z_near*b[:,1]
    near_plane = [zpos-fwidth*b[:,2]-fheight*b[:,3],
                  zpos-fwidth*b[:,2]+fheight*b[:,3],
                  zpos+fwidth*b[:,2]+fheight*b[:,3],
                  zpos+fwidth*b[:,2]-fheight*b[:,3]]
    ppos = projecto(cam, pos)
    @show ppos
    scatter!(ax2, [ppos])
    limits!(ax2, -1.0, 1.0, -1.0, 1.0)
    apos = b[:,2:3]*([ppos...].*[fwidth, fwidth/cam.frustrum_ratio]) + cam.pos - cam.z_near*b[:,1]
    scatter!(ax1, [Point3f(pos), Point3f(apos)])
    linesegments!(ax1, [cam.pos=>cam.pos+cam.dir])
    #linesegments!(ax1, [cam.pos - cam.z_near*b[:,1]-fwidth*b[:,2]=>cam.pos-cam.z_near*b[:,1]+fwidth*b[:,2]])
    linesegments!(ax1, [near_plane[1]=>near_plane[2], near_plane[2]=>near_plane[3], near_plane[3]=>near_plane[4],near_plane[4]=>near_plane[1]])

    linesegments!(ax1, [Point3f(pos)=>Point3f(cam.pos)])
    linesegments!(ax1, [Point3f(pos)=>Point3f(apos)])
    fig
end

"""
    get_cube_normals(x::AbstractVector{T},w::T,h::T,d::T) where T <: Real

The the normals of a cube at position `x`
"""
function get_cube_side(x::AbstractVector{T},w::T,h::T,d::T) where T <: Real
    #figure out the side
    s = 0
    if (0 < x[1] < w) && ( 0 < x[3] < h)
        if x[2] == 0.0
            s = 4
        else
            s = 3 
        end
    elseif (0 < x[2] < d) && (0 < x[3]< h)
        if x[1] == 0.0
            s = 2
        else
            s = 1 
        end
    elseif (0 < x[3] < h) && (0 < x[1] < w)
        if x[2] == 0
            s = 6 
        else
            s = 5
        end
    end
    return s
end

function get_cube_side2(x::AbstractVector{T},w::T,h::T,d::T) where T <: Real
    #figure out the side
    s = 0
    if x[2] == 0.0
        if (0 < x[1] < w) && ( 0 < x[3] < h)
            s = 4
        else
            if x[1] <= 0
                # wrap around to the adjacent surface
                s = 2
            elseif x[1] >= w
                s = 1
            end
            if x[3] <= 0
                s = 6
            elseif x[3] >= h
                s = 5
            end
        end
    elseif x[2] == d
        if (0 < x[1] < w) && ( 0 < x[3] < h)
            s = 3
        else
            if x[1] <= 0
                s = 2
            elseif x[1] >= w
                s = 1
            end
            if x[3] <= 0
                s = 6
            elseif x[3] >= h
                s = 5
            end
        end
    elseif x[1] == 0.0
        if (0 < x[2] < d) && (0 < x[3]< h)
            s = 2 
        else
            if x[2] <= 0
                s = 4
            elseif x[2] >= d
                s = 3
            end
            if x[3] <= 0
                s = 6
            elseif x[3] >= h
                s = 5
            end
        end
    elseif x[1] == w
        if (0 < x[2] < d) && (0 < x[3]< h)
            s = 1
        else
            if x[2] <= 0
                s = 4
            elseif x[2] >= d
                s = 3
            end
            if x[3] <= 0
                s = 6
            elseif x[3] >= h
                s = 5
            end
        end
    elseif x[3] == 0
        if (0 < x[1] < w)  && (0 < x[2] < d)
            s = 6
        else
            if x[1] <= 0
                s = 2
            elseif x[1] >= w
                s = 1
            end
            if x[2] <= 0
                s = 4
            elseif x[2] >= d
                s = 3
            end
        end
    elseif x[3] == h
        if (0 < x[1] < w)  && (0 < x[2] < d)
            s = 5
        else
            if x[1] <= 0
                s = 2
            elseif x[1] >= w
                s = 1
            end
            if x[2] <= 0
                s = 4
            elseif x[2] >= d
                s = 3
            end
        end
    end
    s
end


function get_cube_normals(s::Int64,::Type{T}) where T <: Real
    n = zeros(T, 6, 3)
    n[1,:] = [1.0, 0.0, 0.0]
    n[2,:] = [-1.0, 0.0, 0.0]
    n[3,:] = [0.0, 1.0, 0.0]
    n[4,:] = [0.0, -1.0, 0.0]
    n[5,:] = [0.0, 0.0, 1.0]
    n[6,:] = [0.0, 0.0, -1.0]
    n[s,:]
end

function get_cube_normals(x::AbstractVector{T},w::T,h::T,d::T) where T <: Real

    s = get_cube_side2(x, w, h,d)
    s > 0 || error("Point $(x) is not on a cube surface")
    return get_cube_normals(s,T)
end

function wander_cube(w::T,h::T,d::T,x0::Vector{T}=[0.1, 0.0, 0.1];nn=100,σ=0.01) where T <: Real
    x = x0
    X = zeros(T, 3, nn)
    X[:,1] = x0
    s0 = get_cube_side2(x0, w, h,d)
    n = get_cube_normals(s0,T) 
    b = nullspace(permutedims(n))
    for i in 2:nn
        v = σ*b*randn(2)
        x = X[:,i-1] + v
        s1 = get_cube_side2(x, w, h, d)
        if s0 != s1
            s0 = s1
            # we changed sides
            n = get_cube_normals(s1,T) 
            b = nullspace(permutedims(n))
            # also reset the component parallel to the new normal
            xx = n'*X[:,i] 
            if xx < 0.0
                xx = 0.0
            elseif xx > 1.0
                xx = 1.0
            end
            
            v = σ*b*randn(2)
            X[:,i] = X[:,i-1] +v
            X[:,i]  = (b*b'*X[:,i]) + n*xx
        else
            X[:,i] = x
        end
    end
    X
end