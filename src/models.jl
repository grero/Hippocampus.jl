using Random
using Distributions
using Distances
using Dierckx

struct ViewField
    μ::Vector{Float64}
    σ::Float64
    λ0::Float64
    mm::SimpleMesh
    D::Vector{Float64} # distance from the mean to all other points
end

function ViewField(μ, σ, λ0, mm::SimpleMesh)
    kn = KNearestSearch(mm, 1)
    idx = search(Meshes.Point(μ...), kn)
    D = distancematrix(mm;between_centroids=false)
    ViewField(μ, σ,λ0, mm,D[:,first(idx)])
end

function ViewField(μ::Vector{Float64}, σ::Vector{Float64}, ρ::Matrix{Float64},mm::SimpleMesh)
    d = length(x)
    Σ = zeros(d,d)
    for i in 1:d
        for j in 1:d
            Σ[j,i] = vf.σ[i]*vf.σ[j]*ρ[i,j] 
        end
    end
    D = distancematrix(mm;between_centroids=false)
    # gaussian based on geodesic distance
    # find the location the mean on the mesh
    kn = KNearestSearcv(vf.mm, 1)
    idx = search(Meshes.Point(μ...), kn)
    ViewField(μ, Σ,mm,D[:,first(idx)])
end

function (vf::ViewField)(x::Vector{Float64})
    # need to project onto the mesh
    kn = KNearestSearch(vf.mm, 1)
    idx = search(Meshes.Point(x...), kn)
    vf(first(idx))
end

function (vf::ViewField)(pidx::Int64)
    d = vf.D[pidx]
    vf.λ0*exp(-d^2/(2*vf.σ^2))
end

sigmoid(x, x0, a) = 1.0/(1+exp(a*(x-x0)))
"""
Simulate a simple place field neuron using behavioural data in `udata`
"""
function model_place_field(udata::Union{UnityData,UnityRaytraceData}, rpdata::RippleData;λmin=0.1, λmax=3.0,dt=0.01,σ1=1.0,σ2=σ1, μ=[3.5,1.2],ρ=1.0, fd=0.0, temporal_factor=0.0, sigmoid_params=(Inf,1.0), view_field::Union{ViewField,Nothing}=nothing, rng=Random.default_rng())
    Σ = [σ1^2 ρ*σ1*σ2;ρ*σ1*σ2 σ2^2]
    G = MvNormal(μ, Σ)
    G0 = pdf(G, μ)
    nt = numtrials(udata)
    # vector along the major axis
    aa = (σ1^2+σ2^2)/2 + sqrt(((σ1^2-σ2^2)/2)^2 + ρ*σ1*σ2)
    θ = atan(aa-σ1^2, ρ*σ1*σ2)
    vm = [cos(θ), sin(θ)]
    spikes = Float64[]
    q = -log(rand(rng))
    r = 0.0
    fd = min(1.0, max(0.0, fd))
    tw = 0.5
    qt = 1.0
    qt_sigmoid = 1.0
    @showprogress for i in 1:nt
        if isa(udata, UnityData)
            t,mposx,mposy = get_trial(udata, i;trial_start=2)
        else
            t,gaze,pos,_,_ = get_trial(udata, i;trial_start=2)
            mposx = pos[1,:]
            mposy = pos[2,:]
        end
        # get the time from ripple
        trp = rpdata.timestamps[i,2]
        # reference to ripple time
        t .= t .- t[1] .+ trp
        # create a paramtric spline of the path
        Δ = permutedims([diff(mposx) diff(mposy)])
        # first find points where the gradient is non-zero
        idx = findall(norm.(eachcol(Δ)).>0)
        pos = permutedims([mposx mposy])
        spl = ParametricSpline(t[idx], pos[:,idx])

        Δt = t[2]-t[1]
        for (jk,(t0,posx,posy)) in enumerate(zip(t,mposx, mposy))
            λ = pdf(G, [posx,posy]) # firing rate based on place field
            # What do we do if there is no movement
            if fd > 0
                v = Dierckx.derivative(spl, t0)
                v ./= norm(v)
                f = 0.5*(v'*vm + 1.0) # from 0 to 1
                if isfinite(f)
                    f = (1.0 - fd) + fd*f
                else
                    f = 1.0
                end
            else
                f = 1.0
            end
            # scale firing rate
            if !isnothing(view_field)
                tv = view_field(gaze[:,jk]) 
            else
                tv = 1.0
            end
            λ = tv*qt_sigmoid*qt*f*λmax*λ/G0 + λmin
            _t = t0
            while _t < t0+Δt
                r += λ*dt
                if r >= q
                    push!(spikes, _t)
                    r = 0.0
                    q = -log(rand(rng))
                end
                _t += dt
                qt += temporal_factor*qt*dt
                qt_sigmoid = sigmoid(_t, sigmoid_params[1], sigmoid_params[2])
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


"""
Simulate a case in which apparent directional tuning comes simply from traversing a place field slightly 
differently when traversing in two different directions
"""
function simulate_fake_directionality(;numtrials=100, additive=1.0, multiplicative=1.5)
    m_floor = Shadow("xy")(floor_topology3(;nrefinements=3))

    place_field_center = (0.0, 0.0)
    place_field_size = 4.0
    λ_max = 5.0
    λ_min = 0.1

    # generate firing rate
    # TODO: Maybe not use cityblock here because it looks a bit weird
    D = distancematrix(m_floor)
    # find the centroid
    kn = KNearestSearch(m_floor,1)
    ii = first(search(Meshes.Point(place_field_center...), kn))
    λ = zeros(nelements(m_floor))
    for i in 1:nelements(m_floor)
        d = D[i,ii]
        λ[i] = λ_min + (λ_max-λ_min)*exp(-d^2/(2*place_field_size^2))
    end

    # now simulate trajectories, where we go through the field when going from south to north, off-center
    # when going the other way
    X = zeros(length(λ),2)
    X1 = zeros(size(X)...)
    X2 = zeros(size(X)...)
    Y = zeros(size(X)...)
    A = adjacencymatrix(m_floor)
    G = SimpleGraph(A)
    # northward
    starting_point = first(search(Meshes.Point(0.0, -11.0),kn))
    ending_point = first(search(Meshes.Point(0.0, 11), kn))
    dj = dijkstra_shortest_paths(G, starting_point)
    northward_trajectory =  get_path(dj, ending_point)

    #southward
    starting_point = first(search(Meshes.Point(1.5, 11.0),kn))
    ending_point = first(search(Meshes.Point(1.5, -10.0), kn))
    dj = dijkstra_shortest_paths(G, starting_point)
    southward_trajectory =  get_path(dj, ending_point)

    for i in 1:numtrials
        X[northward_trajectory,1] .+= rand.(Poisson.(λ[northward_trajectory]))
        X1[northward_trajectory,1] .+= rand.(Poisson.(λ[northward_trajectory] .+ additive)) # pure directional
        X2[northward_trajectory,1] .+= rand.(Poisson.(multiplicative*λ[northward_trajectory])) # interaction
        Y[northward_trajectory, 1] .+= 0.05
        X[southward_trajectory,2] .+= rand.(Poisson.(λ[southward_trajectory]))
        X1[southward_trajectory,2] .+= rand.(Poisson.(λ[southward_trajectory]))
        X2[southward_trajectory,2] .+= rand.(Poisson.(λ[southward_trajectory]))
        Y[southward_trajectory, 2] .+= 0.05
    end
    X, X1, X2,Y, λ, northward_trajectory, southward_trajectory
end