using Printf
using GeometryBasics
using GeometryBasics: Point, Rect, Vec, faces,Mat
using LinearAlgebra

function reshape_triggers(markers::AbstractVector{T1}, timestamps::AbstractVector{T2},session_start::Vector{UInt64}=UInt64[];perform_fix=false) where T1 <: Real where T2 <: Real
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
    p1 = 0
    p2 = 0
    nt = 1
    if rem(nn,3) == 0
        nt = div(nn,3)
        trial_markers = permutedims(reshape(_markers, 3, nt))
        trial_timestamps = permutedims(reshape(_timestamps,3,nt))

        # sanity check; make sure that the last number digit is the same for each trial
        # and that the succession is 1,2,3 or 1,2,4.
        main_marker = floor.(trial_markers/10.0)
        p1 =sum(sum(main_marker .≈ [1.0 2.0 3.0],dims=2).==3)
        p2 =sum(sum(main_marker .≈ [1.0 2.0 4.0],dims=2).==3)
    else
        if perform_fix
            main_marker = round.(Int64, floor.(_markers/10.0))
        else
            error("Inconsistent number of markers")
        end
    end

    if p1+p2 != nt 
        if perform_fix
            # try to recover by removing trials that break the pattern
            rmarkers = fix_markers(round.(Int64, permutedims(main_marker)[:]))
            rmarkers = reshape(rmarkers, 3, div(length(rmarkers),3))
            midx = ((!ismissing)).(rmarkers)

            trial_markers_r = Matrix{Union{Missing, T1}}(undef, size(rmarkers)...)
            trial_markers_r[midx] .= _markers 
            trial_timestamps_r = Matrix{Union{Missing, T2}}(undef, size(rmarkers)...)
            trial_timestamps_r[midx] .= _timestamps

            trial_markers = permutedims(trial_markers_r)
            trial_timestamps = permutedims(trial_timestamps_r)
            main_marker = floor.(trial_markers/10.0)
        else
            error("Inconsistent main markers")
        end
    end
    nt = size(trial_markers,1)
    # check the the minor markers are the same
    minor_marker = trial_markers - 10.0*main_marker
    p3 = 0
    for x in eachrow(minor_marker)
        q = filter(!ismissing, x)
        p3 += all(q.==q[1])
    end
    p3 == nt || error("Inconsistent cue markers")
    trial_markers, trial_timestamps
end

function fix_markers(markers)
    recovered_markers = Union{Missing, eltype(markers)}[]
    i = 1
    valid_marker = [2,(3,4),1,1]
    while i <= length(markers)
        # skip markers that do not follow the valid progression
        if i > 1 && !(markers[i] in valid_marker[markers[i-1]])
            push!(recovered_markers, missing)
        end
        push!(recovered_markers, markers[i])
        i += 1
    end 
    # we should end on a 3 or 4
    if !(recovered_markers[end] in [3,4])
        push!(recovered_markers, missing)
    end
    recovered_markers
end

"""
    field_outline(f::Matrix{T};t=2) where T <: Real

Extract contiguous patches of activity from `f` where the activity exceeds μ+t σ
where μ is the overall mean and σ is the overall standard deviation.
"""
function field_outline(f::Matrix{T};t=2) where T <: Real
    lidx = CartesianIndices((1:size(f,1), 1:size(f,2)))
    fidx = findall(isfinite, f)
    σ = std(f[fidx])
    μ = mean(f[fidx])
    fm,fi = findmax(f[fidx])
    avail = fill(true, length(fidx))
    patches = [[fidx[fi]]]
    avail[fi] = false
    pidx = 1
    skip = false
    D = [norm(Tuple(i) .- Tuple(j)) for i in lidx, j in lidx]
    while (sum(avail)>0) && (fm > μ+t*σ)
        skip = false
        # fill until we find it drops
        # grab the nearest point
        _fidx = fidx[avail]
        qq = sortperm(fm .- f[_fidx])
        did_change = false
        qidx = findall(avail)
        for _qq in qq
            _dd = minimum(D[_fidx[_qq], patches[pidx]])
            if _dd > 1
                continue
            end
            did_change = true
            if f[_fidx[_qq]]>= μ+t*σ
                push!(patches[pidx],_fidx[_qq])
            end
            avail[qidx[_qq]] = false
        end
        if !did_change 
            _fidx = fidx[avail]
            qidx = findall(avail)
            fm,fi = findmax(f[_fidx])
            if fm > μ+t*σ
                push!(patches, [_fidx[fi]])
                pidx += 1
            end
            avail[qidx[fi]] = false
        end
    end
    patches
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
    # label to indicate
    label::Vector{Tuple{Symbol,Int64,Int64}}
end

function smooth(counts::Dict{Symbol,Vector{T}}, D::Matrix{T2}, pidx::Vector{Tuple{Symbol, Int64, Int64, Int64}};σ=5.5) where T <: Array{T2,3}  where T2 <: Real
    Z = zeros(T2, size(D,1))
    for (k,vv) in counts
        # hack; we should get this from an associated MazeModel object
        if k == :pillars
            # 4 pillar with 4 walls
            idx = permutedims(CartesianIndices((1:4,1:4)))
        elseif k == :walls
            idx = CartesianIndices((1:1, 1:4))
        elseif k == :floor
            idx = CartesianIndices((1:1, 1:9))
        else #ceiling
            idx = CartesianIndices((1:1, 1:1))
        end
        for (qidx,v) in zip(idx,vv)
            # hackish; grap the side with the most counts
            _size = size(v)
            # dummy side has width 2
            didx = findfirst(_size.==2)
            udims = setdiff(1:3, didx)
            Xp = dropdims(sum(abs.(v),dims=(udims...,)),dims=(udims...,))
            sidx = argmax(Xp)
            vidx = ntuple(p->ifelse(p==didx, sidx, 1:_size[p]), 3)
            X = v[vidx...] 
            # find the points in the distance matrix
            idx = findall(x->(x[1]==k)&(x[2]==qidx[1])&(x[3]==qidx[2]), pidx)
            if length(idx) != length(X)
                @show qidx k
            end
            Z .+= exp.(-D[:,idx].^2/(2*σ^2))*X[:]
        end
    end
    Z
end

function visualize!(lscene, pm::ParametrizedManifold{N,N2,T3,T, T2,T4},color::AbstractVector{T5}, color_points::AbstractVector{T4}, pidx::AbstractVector{T7},sizes::Vector{T6};include_ceiling=true, kwargs...) where T4 <: Point{N,T3} where T <: Vec{N,T3} where T2 <: QuadFace{Int64} where T3 <: Real where N2 where T5 <: Real where T6 <: NTuple{3,Int64} where T7 <: Tuple{Symbol, Int64, Int64, Int64} where N
    # each face is implemented as a separate mesh
    n = length(pm.faces)
    # TODO Make sure the colorscale is the same
    cl = extrema(color)
    for ii in 1:n
        lidx = pm.label[ii]
        if lidx[1] == :ceiling && include_ceiling == false
            continue
        end
        offset = (ii-1)*N2
        points = pm.points[(ii-1)*N2+1:ii*N2]
        _ff = pm.faces[ii]
        fpoints = pm.points[_ff].points
        # rescale 
        ff = _ff .- offset
        nn = pm.normals[ii]
        # also a bit hackish

        # this is a bit clunky;
        # find the portion of the color vector pertaining to this surface
        # this is brittle as it depends on an exact match
        cidx = findall(x->(x[1]==lidx[1])&(x[2]==lidx[2])&&(x[3]==lidx[3]), pidx)
        @debug lidx length(cidx)
        if isempty(cidx)
            continue
        end
        _size = filter(x->x>1, sizes[ii])
        fq = reshape(color[cidx],_size)
        # I have no idea why this is necessary
        fq = rotl90(fq)
        pq = reshape(color_points[cidx], _size)

        # use the rect wrapper to find the origin
        rf = Rect([fpoints...])
        # identify the origin
        oidx = findfirst(_pq->_pq==origin(rf), pq)
        # identify the corners
        idx1 = findfirst(fp->fp==pq[oidx],points)
        idx2 = findfirst(fp->fp==pq[oidx.I[1],end],points)
        idx3 = findfirst(fp->fp==pq[end,oidx.I[2]],points)
        idx4 = findfirst(fp->fp==pq[end,end],points)
        uv = Vector{Vec2f}(undef,4)
        uv[idx1] = Vec2f(0,0)
        uv[idx2] = Vec2f(0,1)
        uv[idx3] = Vec2f(1,0)
        uv[idx4] = Vec2f(1,1)
        gb_mesh = GeometryBasics.Mesh(Point3f.(points), [ff];normal=Vec3f.([nn for _ in 1:N2]),uv=uv)
        # debug; plot the base
        #arrows!(lscene, Point3f.([pm.μ[ii], pm.μ[ii]]),Point3f.(eachcol(pm.base[ii])),color=:white)
        # indicate the uv points
        #scatter!(lscene, points[[idx1,idx2,idx3,idx4]],color=Makie.wong_colors()[1:4])
        mesh!(lscene, gb_mesh, color=fq,colorrange=cl)
    end
end

function visualize!(lscene, pm::ParametrizedManifold{N,N2,T3,T, T2,T4};kwargs...) where T4 <: Point{N,T3} where T <: Vec{N,T3} where T2 <: QuadFace{Int64} where T3 <: Real where N2 <: Val{4} where N
    # need to replicate normals to every point
    nn = Vector{Vec{N,T}}(undef, length(pm.points))
    for ii in 1:length(pm.normals)
        for k in 1:4
            nn[(ii-1)*4+k] = pm.normals[ii]
        end
    end
    gb_mesh = GeometryBasics.Mesh(pm.points, pm.faces;normal=nn)
    wireframe!(lscene, gb_mesh)
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
    sidx = assign_to_surface(p0, nn, μ;visited=visited)
    visited[sidx] = true
    bb2 = bb[sidx]*bb[sidx]'
    # project onto the surface
    p0p = Point{N,T}(bb2*p0) + (μ[sidx]'*nn[sidx]).*nn[sidx]

    # now do the actual projection
    d = p1 - p0p

    # create a a boundingbox
    rf = Rect([pm.points[ff[sidx]].points...])
    dp = bb2'*d
    @debug "Some" p0 p0p d p0 + dp
    if p0p + dp ≈ p1
        #return norm(dp),[p1]
        return abs.(bb[sidx]'*dp),[p1]
    end
    # if the projected distance is zero and we reached here, 
    # that means that we need to travel to the end of the manifold
    @debug "norm" norm(dp) sidx

    p1p = Point{N,T}(bb2'*p1) + (μ[sidx]'*nn[sidx]).*nn[sidx]

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

function assign_to_surface(p0::Point{N,T}, normals, μ;visited=fill(false, length(normals))) where T <: Real where N
    nn = normals
    i0 = 0 
    d0 = typemax(T)
    for (ii,_nn) in enumerate(nn)
        if visited[ii]
            continue
        end
        _d0 = norm(((p0-μ[ii])'*nn[ii]))
        if _d0 < d0
            i0 = ii
            d0 = _d0
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

"""
Return a discrete disc of radius `r` centered on `i`
"""
function disc(p0::CartesianIndex{2}, r::Int64,n1::Int64, n2::Int64)
    idx = Vector{CartesianIndex{2}}()
    # first make a square
    idx1 = max(1, p0.I[1]-r):min(p0.I[1]+r, n1)
    idx2 = max(1, p0.I[2]-r):min(p0.I[2]+r, n2)
    for i in idx1
        for j in idx2
            if (i-p0.I[1])^2 + (j-p0.I[2])^2 <= r^2
                push!(idx, CartesianIndex(i,j))
            end
        end
    end
    idx
end

function adaptive_smoothing(X::Matrix{T}, Y::Matrix{T}, α::T;stop_at_nan=true) where T <: Real
    n1,n2 = size(X)
    Z = fill!(similar(X), zero(T))
    for ii in CartesianIndices(size(X))
        nsp = X[ii]
        nocc = Y[ii]
        r = 1
        while nsp < α/(nocc^2*r^2) 
            idx = disc(ii, r, n1, n2)
            # stop expanding the kernel if we hit boundary
            if stop_at_nan && any(isnan.(Y[idx]))
                break
            end
            nsp = sum(X[idx])
            nocc = sum(Y[idx])
            r += 1
        end
        Z[ii] = nsp/nocc 
    end
    Z
end