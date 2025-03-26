using Printf

function reshape_triggers(markers, timestamps)
    # the first marker is a session start; the remaining come in trios
    nn = length(markers)
    if markers[1] == 84
        nn -= 1
        _markers = markers[2:end]
        _timestamps = timestamps[2:end]
    else
        _markers = markers
        _timestamps = timestamps
    end
    rem(nn,3) == 0 || error("Inconsistent number of markers")
    nt = div(nn,3)
    trial_markers = permutedims(reshape(_markers, 3, nt))
    trial_timestamps = permutedims(reshape(_timestamps,3,nt))

    # sanity check; make sure that the last number digit is the same for each trial
    # and that the succession is 1,2,3 or 1,2,4.
    main_marker = floor.(trial_markers/10.0)
    p1 =sum(sum(main_marker .≈ [1.0 2.0 3.0],dims=2).==3)
    p2 =sum(sum(main_marker .≈ [1.0 2.0 4.0],dims=2).==3)

    p1+p2 == nt || error("Inconsistent main markers")

    # check the the minor markers are the same
    minor_marker = trial_markers - 10.0*main_marker
    p3 = sum(sum(minor_marker .== minor_marker[:,1:1],dims=2).==3)
    p3 == nt || error("Inconsistent cue markers")
    trial_markers, trial_timestamps
end

struct Trial
    i::UInt64
end

function parse_cellname(cellname::String)
    re = r"(\d{8,8})ch(\d{1,3})c(\d+)"
    m = match(re, cellname)
    channel = @sprintf("%03d", parse(Int64, m[2]))
    cell = @sprintf("%02d", parse(Int64, m[3]))
    pth = glob(joinpath(m[1], "session*","array*","channel$(channel)","cell$(cell)"))
    if !isempty(pth)
        return first(pth)
    else
        return nothing
    end
end

"""
Return a disk filter
"""
function disk(r)
    rn = round(r)
    n = 2*rn+1

    f = fill(0.0, n, n)
    center = div(n,2)+1
    for j in axes(f,2)
        for i in axes(f,1)
            d = sqrt((i - center)^2 + (j-center)^2)
            if d <= r
                f[i,j] = 1.0
            end
        end
    end
    f ./= sum(f)
end

function isinside(p::Point{N,T}, r::Rect{N,T}) where T where N
    origin(r) <= p <= origin(r) + widths(r)
end

"""
Compute the distance between the points `p0` and `p1` along the manifold `m`
"""
#function distance(p0::Point{N,T}, p1::Point{N,T}, m) where N where T <: Real
#end

struct ParametrizedManifold{N<:Any, N2<:Any, T3<:Real, T<:Vec{N,T3},T2<:NgonFace{N2,Int64}, T4 <: Point{N,T3}}
    normals::Vector{T}
    faces::Vector{T2}
    μ::Vector{T} 
    base::Vector{Matrix{T3}}
    points::Vector{T4}
end

function ParametrizedManifold(m::HyperRectangle{N,T}) where T <: Real where N
    nn = -decompose(GeometryBasics.Normal(Vec{N,T}),m).data
    points = decompose(Point{N,T},m)
    ff = faces(m)
    μ = Vector{Vec{N,T}}(undef, length(nn))
    bb = Vector{Matrix{T}}(undef, length(nn))
    for kk in 1:length(bb)
        bb[kk] = abs.(nullspace(permutedims(nn[kk])))
        μ[kk] = mean(points[ff[kk]])
    end
    ParametrizedManifold(nn, ff, μ, bb, points)
end

function distance(p0::T4, p1::T4, pm::T5;visited=fill(false, length(pm.faces)))  where T5<:ParametrizedManifold{N, M, T, T2, T3, T4} where T4 <: Point{N,T} where T2 <: Vec{N,T} where T3 <: NgonFace{M, <:Integer} where T <: Real  where N where M
    nn = pm.normals
    ff = pm.faces
    μ = pm.μ
    bb = pm.base
    sidx = assign_to_surface(p0, nn, ff, bb, μ, pm.points;visited=visited)
    visited[sidx] = true
    # project onto the surface
    p0p = Point{N,T}(bb[sidx]*bb[sidx]'*p0) + (μ[sidx]'*nn[sidx]).*nn[sidx]

    # now do the actual projection
    d = p1 - p0p

    # create a a boundingbox
    rf = Rect([pm.points[ff[sidx]].points...])
    dp = bb[sidx]*(bb[sidx]'*d)
    @debug "Some" p0 p0p d p0 + dp
    if p0p + dp ≈ p1
        #return norm(dp),[p1]
        return abs.(bb[sidx]'*dp),[p1]
    end
    # if the projected distance is zero and we reached here, 
    # that means that we need to travel to the end of the manifold
    @debug "norm" norm(dp) sidx

   # if norm(dp) == zero(T)
   #     # TODO: Choose the base vector most aligned with the distance vector
   #     dp = bb[sidx][:,2]
   # end
    # TODO: Find the edge point closest to p1
    p1p = Point{N,T}(bb[sidx]*bb[sidx]'*p1) + (μ[sidx]'*nn[sidx]).*nn[sidx]

    # TODO: This does not work for the cloor and ceiling. The pillars current sit on top of the
    # floor, but does not intersect it.
    dx = displacement_to_edge(rf, p1p)
    p0n,dp = move_to_edge(rf, p0p, dx)

    # project dp onto the base
    Δp = abs.(bb[sidx]'*dp)
    @debug p0n p0 p0p Δp p1p dx dp sidx

    @debug "Show after" p0 sidx #rf#dpn #visited #sidx[1] Δp rf bb[sidx[1]] fm dq nn[sidx[1]] dpn
    pth = [p0n]
    Δp1, pth1 = distance(p0n, p1, pm;visited=visited)
    # project the distance Δp1 onto its own 2D space again
    append!(pth, pth1)
    return Δp + Δp1,pth 

end

function assign_to_surface(p0::Point{N,T}, normals, _faces, bases, μ,points;visited=fill(false, length(normals))) where T <: Real where N
    nn = normals
    bb = bases
    ff = _faces
    i0 = 0 
    d0 = typemax(T)
    for (ii,_nn) in enumerate(nn)
        if visited[ii]
            continue
        end
        # project onto the plane
        # to make sure we are within the bounds of this surface
        pp = Point{N,T}(bb[ii]*bb[ii]'*p0) + (μ[ii]'*nn[ii]).*nn[ii]
        rf = Rect([points[ff[ii]].points...])
        if true
        #if in(rf)(pp)
            _d0 = norm(((p0-μ[ii])'*nn[ii]))
            if _d0 < d0
                i0 = ii
                d0 = _d0
            end
        end
    end
    i0
end

function find_closest_edgepoint(rf::HyperRectangle{N,T}, p1::Point{N,T}) where T <: Real where N
    _p1 = p1 - origin(rf)
    w = [widths(rf)...]
    x = zeros(T,N)
    iim = 0
    Δ = typemax(T)
    is_edgepoint = false
    for ii in eachindex(x)
        if w[ii] > 0 
            if _p1[ii] >= w[ii]
                x[ii] = w[ii] 
                is_edgepoint = true
                iim == ii
                Δ = 1.0
            elseif _p1[ii] <= 0.0
                x[ii] = zero(T)
                is_edgepoint = true
                iim = ii
                Δ = -1.0
            else
                x[ii] = _p1[ii]
                _Δ,iq = findmin(_p1[ii], w[ii]-_p1[ii])
                if _Δ < Δ
                    Δ = _Δ
                    iim = ii
                end
            end
        else
            x[ii] = _p1[ii]
        end
    end
    Δp = zeros(T, N)
    if is_edgepoint
        Δ = zero(T)
    else
        x[iim] -= Δ
        Δp[iim] = -1.0 
    end

    Point{N,T}(x) + origin(rf), Δp
end

function displacement_to_edge(rf::HyperRectangle{N,T}, p1::Point{N,T}) where T <: Real where N
    _p1 = p1 - origin(rf)
    w = [widths(rf)...]
    x = zeros(T,N)
    Δ = typemax(T)
    qq = 0
    iim = 0
    is_edgepoint = false
    for ii in eachindex(x)
        if w[ii] > 0 
            if _p1[ii] >= w[ii]
                x[ii] = w[ii] 
                is_edgepoint = true
                iim == ii
                x[ii]  = 1.0
            elseif _p1[ii] <= 0.0
                x[ii] = zero(T)
                is_edgepoint = true
                iim = ii
                x[ii] = -1.0
            else
                _Δ,iq = findmin([_p1[ii], w[ii]-_p1[ii]])
                if _Δ < Δ 
                    iim = ii
                    Δ = _Δ
                    if iq == 1
                        qq = -1
                    else
                        qq = 1
                    end
                end
            end
        end
    end
    if !is_edgepoint
        x[iim] = qq
    end
    Vec{N,T}(x)
end

function direction_to_edge(rf::HyperRectangle{N,T}, p1::Point{N,T}) where T <: Real where N

    _p1 = p1 - origin(rf)
    w = [widths(rf)...]
    x = zeros(T,N)
    Δ = typemax(T)
    qq = 0
    iim = 0
    is_edgepoint = false
    for ii in eachindex(x)
        if w[ii] > 0 
            if _p1[ii] >= w[ii]
                x[ii] = w[ii] 
                is_edgepoint = true
                iim == ii
                x[ii]  = 1.0
            elseif _p1[ii] <= 0.0
                x[ii] = zero(T)
                is_edgepoint = true
                iim = ii
                x[ii] = -1.0
            else
                _Δ,iq = findmin([_p1[ii], w[ii]-_p1[ii]])
                if _Δ < Δ 
                    iim = ii
                    Δ = _Δ
                    if iq == 1
                        qq = -1
                    else
                        qq = 1
                    end
                end
            end
        end
    end
    if !is_edgepoint
        x[iim] = qq
    end
    x
end

function move_to_edge(rf::HyperRectangle{N,T},p0::Point{N,T}, dp::Vec{N,T}) where T <: Real where N
    dpn = dp./norm(dp)
    w = widths(rf)
    Δ = fill(zero(T), length(w))
    _p0 = p0 - origin(rf)
    w = [widths(rf)...]
    Δ = typemax(T)
    # project the point onto the vector
    v = zeros(T,N)
    for ii in 1:N
        if _p0[ii] <= 0.0
            # already at the left edge
            if dpn[ii] > 0.0
                # need to move to the other edge
                v[ii] = w[ii]
            end # if not we need not move
        elseif _p0[ii] >= w[ii]
            if dpn[ii] < 0.0
                # need to move to the other edge
                v[ii] = -w[ii]
            end # if not no need to move
        else
            if dpn[ii] < 0.0
                v[ii] = -_p0[ii]
            elseif dpn[ii] > 0.0
                v[ii] = w[ii] - _p0[ii]
            end
        end
    end
    Point{N,T}(origin(rf) + _p0 + v), v
end

function move_to_edge_old(rf::HyperRectangle{N,T},p0::Point{N,T}, dp::Vec{N,T};incr=T(0.01)) where T <: Real where N
    dpn = dp./norm(dp)
    Δp = zero(T) 
    incr = T(0.01)
    # this doesn't include the interior
    inrf = in(rf)
    # figure out how far
    #p0 = p0 + incr*dpn 
    while true
        if !inrf(p0+incr*dpn)
            break
        end
        p0 = p0 + incr*dpn
        Δp += incr 
    end
    Δp 
end

function distance(p0::Point{N,T}, p1::Point{N,T}, m::T2;visited=fill(false, length(faces(m))))  where T2 <: Rect{N,T} where T <: Real where N
    # normals pointing to the space
    nn = -decompose(GeometryBasics.Normal(Vec3f),m).data
    points = decompose(Point{N,T},m)
    ff = faces(m)
    μ = Vector{Vec{N,T}}(undef, length(nn))

    # project points onto their nearest surface
    pp0 = Point{N,T}(zero(T))
    pp1 = Point{N,T}(zero(T))
    bb = Vector{Matrix{T}}(undef, length(nn))
    for kk in 1:length(bb)
        bb[kk] = abs.(nullspace(permutedims(nn[kk])))
        μ[kk] = mean(points[ff[kk]])
    end

    p = [p0,p1]
    pp = [pp0,pp1]
    sidx = fill(0, 2)
    d0 = typemax(T) 
    i0 = 0 
    for (ii,_nn) in enumerate(nn)
        if visited[ii]
            continue
        end
        #_d0 = norm((p0-μ[ii])'*_nn)
        _d0 = norm((p0-μ[ii]))
        if _d0 < d0
            i0 = ii
            d0 = _d0
        end
    end
    @debug "Some" μ[i0]
    sidx[1] = i0
    visited[i0] = true
    # project the point
    # FIXME: This does not work
    p0 = Point{N,T}(bb[sidx[1]]*bb[sidx[1]]'*p0) + (μ[sidx[1]]'*nn[sidx[1]]).*nn[sidx[1]]

    Q = diagm(fill(one(T),N))
    # now do the actual projection
    d = p1 - p0
    # project onto 
    # create scaling matrix
    fm = Mat{N,N}(Q - diagm(abs.(nn[sidx[1]])))
    # TODO: This does not appear to work
    dq = (μ[sidx[1]]'*nn[sidx[1]]).*nn[sidx[1]]
    rf = fm*(m - dq) + dq
    dp = bb[sidx[1]]*(bb[sidx[1]]'*d)
    @debug "Some" p0 d p0 + dp
    if p0 + dp ≈ p1
        return norm(dp)
    end
    # if the projected distance is zero and we reached here, 
    # that means that we need to travel to the end of the manifold
    if norm(dp) == zero(T)
        dp = bb[sidx[1]][:,1]
    end
    dpn = dp/norm(dp)
    Δp = zero(T) 
    incr = T(0.01)
    inrf = in(rf)
    while true
        if !inrf(p0+incr*dpn)
            break
        end
        p0 = p0 + incr*dpn
        Δp += incr 
    end
    # we need a way to indicate that we are the edge. Perhaps a flag to indicate which manifold we've alrady traversed?
    # we need another condition here; this will repeat as long as pp[1] is inside the current rectangle
    @debug "Show" p0 #rf#dpn #visited #sidx[1] Δp rf bb[sidx[1]] fm dq nn[sidx[1]] dpn
    return Δp + distance(p0, p1, m;visited=visited)
end