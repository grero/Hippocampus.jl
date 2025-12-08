using Meshes
using Meshes: connect
using Graphs
using Graphs: DijkstraState
using Unitful

using GeometryBasics

function get_path(dj::DijkstraState, v::Integer)
    u = v
    path = [u]
    while dj.parents[u] != 0
        u = dj.parents[u]
        pushfirst!(path, u)
    end
    return path
end

function box_topology()
    _points = decompose(Point3f, Rect3(0.0, 0.0, 0.0, 1.0, 1.0, 1.0))
    points = [(p...,) for p in _points]
    cnx = [connect((1,2,4,3)), # x == 0
           connect((5,6,8,7)), # x == 1
           connect((1,2,6,5)), # y == 0
           connect((3,4,8,7)), # y == 1
           connect((1,3,7,5)), # z == 0
           connect((2,4,8,6))] # z == 1
    mesh = SimpleMesh(points, cnx;relations=true)
end

function box_topology2(n::Integer=2)
    box_mesh = boundary(Meshes.Box((0.0, 0.0, 0.0),(1.0, 1.0, 1.0)))
    for _ in 1:n
        box_mesh = refine(box_mesh,QuadRefinement())
    end
    # convert to using HalfEdgeTopology
    ref2 = SimpleMesh(box_mesh.vertices, convert(HalfEdgeTopology, box_mesh.topology))
end

function nearestneighbors(mm::SimpleMesh)
    A = adjacencymatrix(mm;rank=0)
    G = SimpleGraph(A)
    function nn(ii::Integer, n::Integer)
        dj = dijkstra_shortest_paths(G, ii, A;trackvertices=true)
        dj.closest_vertices[2:n+1]
    end
end

function ballsearch(mm::SimpleMesh;rank=paramdim(mm))
    A = adjacencymatrix(mm;rank=rank)
    G = SimpleGraph(A)
    function searchball(ii::Integer, r::Integer)
        dj = dijkstra_shortest_paths(G, ii, A;trackvertices=false)
        idx = findall(dj.dists .<= r)
    end
end

function get_circular_adjancency(n::Integer)
    A = diagm(fill(1.0,n))
    A[1,end] = 1.0
    A[end,1] = 1.0

    # each point is connected to two neighbors, in front and behind
    for i in 2:n-1
        A[i,i+1] = 1.0
        A[i+1,i] = 1.0
        A[i,i-1] = 1.0
        A[i-1,i] = 1.0
    end
    A
end

function distancematrix(mm::SimpleMesh;rank=paramdim(mm))
    nn = length(mm.vertices)
    A = adjacencymatrix(mm;rank=rank)
    if issymmetric(A)
        G = SimpleGraph(A)
    else
        G = SimpleDiGraph(A)
    end
    D = zeros(nn,nn)
    for ii in 1:nn
        dj = dijkstra_shortest_paths(G, ii;trackvertices=false)
        D[:,ii] = dj.dists 
    end
    D
end

function distancematrix(mm::Matrix{T}) where T <: Real
    lidx = CartesianIndices((1:size(mm,1), 1:size(mm,2)))
    D = [norm(Tuple(i) .- Tuple(j)) for i in lidx, j in lidx]
    D
end

function maze_topology(xmin=-12.5, xmax=12.5, ymin=xmin, ymax=xmax, zmin=0.0, zmax=5.0)
    # pillars
    x1Bound = [-7.5, -2.5, -2.5, -7.5, -7.5]  # yellow pillar
    z1Bound = [7.5, 7.5, 2.5, 2.5, 7.5]

    x2Bound = [2.5, 7.5, 7.5, 2.5, 2.5]  # red pillar
    z2Bound = [7.5, 7.5, 2.5, 2.5, 7.5]

    x3Bound = [-7.5, -2.5, -2.5, -7.5, -7.5]  # blue pillar
    z3Bound = [-2.5, -2.5, -7.5, -7.5, -2.5]

    x4Bound = [2.5, 7.5, 7.5, 2.5, 2.5]  # green pillar
    z4Bound = [-2.5, -2.5, -7.5, -7.5, -2.5]

    #first create the walls and the ceiling
    # these can be derived from a box
    points = decompose(Point3f, Rect3(xmin, ymin, zmin, xmax-xmin, ymax-ymin, zmax-zmin))
    #anchor points for the pillars
    push!(points, xmin, -7.5, 0.0)
    push!(points, xmin, -2.5, 0.0)
    push!(points, xmin, 2.5, 0.0)
    push!(points, xmin, 7.5, 0.0)

    push!(points, xmax, -7.5, 0.0)
    push!(points, xmax, -2.5, 0.0)
    push!(points, xmax, 2.5, 0.0)
    push!(points, xmax, 7.5, 0.0)

    push!(points, -7.5, ymin, 0.0)
    push!(points, -2.5, ymin, 0.0)
    push!(points, 2.5, ymin, 0.0)
    push!(points, 7.5, ymin, 0.0)

    push!(points, -7.5, ymin, 0.0)
    push!(points, -2.5, ymin, 0.0)
    push!(points, 2.5, ymin, 0.0)
    push!(points, 7.5, ymin, 0.0)

    pillar_points1 = decompose(Point3f, Rect3(-7.5, 2.5, 0.0, 5.0, 5.0, 1.5))
    pillar_points2 = decompose(Point3f, Rect3(2.5, 2.5, 0.0, 5.0, 5.0, 1.5))
    pillar_points3 = decompose(Point3f, Rect3(-7.5, -7.5, 0.0, 5.0, 5.0, 1.5))
    pillar_points4 = decompose(Point3f, Rect3(2.5, -7.5, 0.0, 5.0, 5.0, 1.5))
    append!(points, pillar_points1)
    append!(points, pillar_points2)
    append!(points, pillar_points3)
    append!(points, pillar_points4)

    walls = [connect((1,2,4,3)),
             connect(())]

    walls_and_ceiling = [connect((1,2,4,3)), # x == 0
                         connect((5,6,8,7)), # x == 1
                         connect((1,2,6,5)), # y == 0
                         connect((3,4,8,7)), # y == 1
                         connect((2,4,8,6))] # z == 1

    # floor is more complicated since we need to incorporate the pillars
    # connect each corner to each corner of the pillar
    offset = 8
    pillars = [connect(1,2*offset+1)]

end

function maze_topology2(xmin=-12.5, xmax=12.5, ymin=xmin, ymax=xmax, zmin=0.0, zmax=5.0)

    outer_walls = boundary(Meshes.Box((xmin,ymin,zmin),(ymax,ymax,zmax)))

    pillar_1 = boundary(Meshes.Box((-7.5, -7.5, 0.0),(-2.5, -2.5,1.5)))
    pillar_2 = boundary(Meshes.Box((-7.5, 2.5, 0.0),(-2.5, 7.5,1.5)))
    pillar_3 = boundary(Meshes.Box((2.5, -7.5, 0.0),(7.5, -2.5,1.5)))
    pillar_4 = boundary(Meshes.Box((2.5, 2.5, 0.0),(7.5, 7.5,1.5)))

    @show intersection(outer_walls, boundary(Meshes.Box((-7.5, -7.5, 0.0), (-2.5, -2.5, 0.0))))
    #merge geomtries one by one
    obj = Base.merge(outer_walls, pillar_1)
end

function connect_points(points, allpoints)
    idx = [findfirst(p->p==point, allpoints) for point in points]
    (idx...,)
end

function fix_connection(cnx::NTuple{4,<:Integer}, points::Vector{T};clockwise=true) where T <: NTuple{3,<:Real}
    cpoints = points[[cnx...]]
    # use consistent order
    # start from lower left
    idx = sortperm(cpoints)
    if clockwise
        return cnx[idx[[1,2,4,3]]]
    end
    cnx[idx[[1,3,4,2]]]
end

function has_edges(cc::NTuple{4, <:Integer}, cnx::Vector{NTuple{4,Int64}})
    eef = collect(zip(cc[1:end-1], cc[2:end]))
    do_reverse = false
    for _cnx in cnx
        _eef = collect(zip(_cnx[1:end-1], _cnx[2:end]))
        _eer = reverse.(_eef)
        if !isempty(intersect(eef, _eer))
            do_reverse = true
            break
        end
    end
    do_reverse
end

function maze_topology3(xmin=-12.5, xmax=12.5, ymin=xmin, ymax=xmax, zmin=0.0, zmax=5.0)
    points = [(xmin, ymin, 0.0), (-7.5, ymin,0.0),(-7.5, -7.5,0.0),(xmin, -7.5, 0.0)] 
    cnx = [fix_connection((1,2,3,4),points)]
    append!(points, [(-7.5, -2.5, 0.0),(xmin, -2.5, 0.0)])
    _cnx = (4,3,5,6)
    _do_reverse = has_edges(_cnx, cnx)
    push!(cnx, fix_connection(_cnx,points;clockwise=!_do_reverse))
    append!(points, [(xmin, 2.5, 0.0),(-7.5, 2.5, 0.0)])
    _cnx = (5,6,7,8)
    _do_reverse = has_edges(_cnx, cnx)
    push!(cnx, fix_connection(_cnx,points;clockwise=!_do_reverse))
    append!(points, [(-7.5, 7.5, 0.0), (xmin, 7.5, 0.0)])
    _cnx = (7,8,9,10)
    _do_reverse = has_edges(_cnx, cnx)
    push!(cnx, fix_connection(_cnx,points;clockwise=!_do_reverse))
    append!(points, [(xmin,ymax, 0.0),(-7.5, ymax, 0.0)])
    _cnx = (9,10,11,12)
    _do_reverse = has_edges(_cnx, cnx)
    push!(cnx, fix_connection(_cnx,points;clockwise=!_do_reverse))
    append!(points, [(-2.5, 7.5, 0.0),(-2.5, ymax, 0.0)])
    _cnx = (9,13,14,12)
    _do_reverse = has_edges(_cnx, cnx)
    push!(cnx, fix_connection(_cnx,points;clockwise=!_do_reverse))
    append!(points, [(2.5, 7.5, 0.0),(2.5, ymax, 0.0)])
    _cnx = (13,15,16,14)
    _do_reverse = has_edges(_cnx, cnx)
    push!(cnx, fix_connection(_cnx,points;clockwise=!_do_reverse))
    append!(points, [(7.5, 7.5, 0.0), (7.5, ymax, 0.0)])
    _cnx = (15,17,18,16)
    _do_reverse = has_edges(_cnx, cnx)
    push!(cnx, fix_connection(_cnx,points;clockwise=!_do_reverse))
    append!(points, [(xmax, 7.5, 0.0), (xmax, ymax, 0.0)])
    _cnx = (17,19,20,18)
    _do_reverse = has_edges(_cnx, cnx)
    push!(cnx, fix_connection(_cnx,points;clockwise=!_do_reverse))
    append!(points, [(7.5, 2.5, 0.0), (xmax, 2.5, 0.0)])
    _cnx = (21,22,19,17)
    _do_reverse = has_edges(_cnx, cnx)
    push!(cnx, fix_connection(_cnx,points;clockwise=!_do_reverse))
    append!(points, [(7.5, -2.5, 0.0), (xmax, -2.5, 0.0)])
    _cnx = (23,24,22,21)
    _do_reverse = has_edges(_cnx, cnx)
    push!(cnx, fix_connection(_cnx,points;clockwise=!_do_reverse))
    append!(points, [(7.5, -7.5, 0.0), (xmax, -7.5, 0.0)])
    _cnx = (25,26,24,23)
    _do_reverse = has_edges(_cnx, cnx)
    push!(cnx, fix_connection(_cnx,points;clockwise=!_do_reverse))
    append!(points, [(7.5, ymin, 0.0),(xmax, ymin, 0.0)])
    _cnx = (27,28,26,25)
    _do_reverse = has_edges(_cnx, cnx)
    push!(cnx, fix_connection(_cnx,points;clockwise=!_do_reverse))
    append!(points, [(2.5, ymin, 0.0), (2.5, -7.5, 0.0)])
    #push!(cnx, (29,27, 25,30))
    idx = connect_points([(2.5, ymin, 0.0), (2.5, -7.5, 0.0), (7.5,-7.5, 0.0),(7.5,ymin, 0.0)], points)
    _do_reverse = has_edges(idx, cnx)

    push!(cnx, fix_connection(idx,points;clockwise=!_do_reverse))
    append!(points, [(7.5, ymin, 0.0),(7.5, -7.5,0.0)])
    _cnx = (27,29,30,25)
    _do_reverse = has_edges(idx, cnx)
    push!(cnx, fix_connection(_cnx,points;clockwise=!_do_reverse))
    append!(points, [(-2.5, ymin, 0.0), (-2.5, -7.5, 0.0)])
    _cnx = (33,34,30,29)
    _do_reverse = has_edges(_cnx, cnx)
    _cnx = fix_connection(_cnx, points;clockwise=!_do_reverse)
    @show _cnx
    push!(cnx, fix_connection(_cnx,points;clockwise=!_do_reverse))
    push!(cnx, fix_connection((2,33,34,3),points))
    append!(points, [(-2.5, -2.5, 0.0), (2.5, -2.5, 0.0)])
    push!(cnx, fix_connection((34, 35, 36,30),points))
    append!(points, [(-2.5, 2.5, 0.0), (2.5, 2.5, 0.0)])
    push!(cnx, fix_connection((37, 13, 15,38),points)) 
    #append!(points, [()])
    #connect the inner points
    push!(cnx, fix_connection((35, 36,38,37),points))
    push!(cnx, fix_connection((5, 35, 37,8),points))
    push!(cnx, fix_connection((36, 23, 21, 38),points))
    floor_idx = 1:length(points)
    #add pillars
    # let's connect the pillar tops counter-clockwise
    # first pillar
    append!(points, [(-7.5, -7.5, 2.5), (-7.5, -2.5, 2.5)])
    append!(points, [(-2.5, -7.5, 2.5), (-2.5, -2.5, 2.5)])
    push!(cnx, fix_connection((3, 39, 40, 5),points))
    push!(cnx, fix_connection((34, 41, 42,35),points))
    push!(cnx, fix_connection((3, 39, 41, 34),points))
    push!(cnx, fix_connection((5, 40, 42,35),points))
    pillar_1_idx = (floor_idx[end]+1):length(points)

    #second pillar
    append!(points, [(2.5, -7.5, 2.5), (2.5, -2.5, 2.5)]) 
    append!(points, [(7.5, -7.5, 2.5), (7.5, -2.5, 2.5)]) 
    push!(cnx, fix_connection((30, 43, 44,36),points))
    push!(cnx, fix_connection((25, 45, 46, 23),points))
    push!(cnx, fix_connection((30,43, 45, 25),points))
    push!(cnx, fix_connection((36,44,46,23), points))
    pillar_2_idx = (pillar_1_idx[end]+1):length(points)

    #third pillar
    append!(points, [(2.5, 2.5, 2.5), (2.5, 7.5, 2.5)]) 
    append!(points, [(7.5, 2.5, 2.5), (7.5, 7.5, 2.5)]) 
    push!(cnx, fix_connection((38,47,48,15), points))
    push!(cnx, fix_connection((21,49,50,17),points))
    push!(cnx, fix_connection((38,47,49,21),points))
    push!(cnx, fix_connection((15, 48,50,17),points))
    pillar_3_idx = (pillar_2_idx[end]+1):length(points)


    #fourth pillar
    append!(points, [(-7.5, 2.5, 2.5),(-7.5, 7.5, 2.5)])
    append!(points, [(-2.5, 2.5, 2.5),(-2.5, 7.5,2.5)])
    push!(cnx, fix_connection((8,51,52,9),points))
    push!(cnx, fix_connection((37,53,54,13),points))
    push!(cnx, fix_connection((8,51,53,37),points))
    push!(cnx, fix_connection((9,52,54,13),points))
    pillar_4_idx = (pillar_3_idx[end]+1):length(points)


    #now the walls
    upper_wall_points = [(xmin,ymin, zmax),(xmin, -7.5, zmax),(xmin, -2.5, zmax), (xmin, 2.5, zmax),
                     (xmin, 7.5, zmax), (xmin, ymax, zmax), (-7.5, ymax,zmax), (-2.5, ymax,zmax),
                     (2.5, ymax,zmax),(7.5, ymax,zmax), (xmax, ymax,zmax), (xmax, 7.5,zmax),
                     (xmax, 2.5,zmax),(xmax, -2.5,zmax),(xmax, -7.5,zmax), (xmax, ymin,zmax),
                     (7.5, ymin,zmax),(2.5, ymin,zmax), (-2.5, ymin,zmax),(-7.5, ymin,zmax)]
    append!(points, upper_wall_points)
    #maybe be a bit smarter here
    for (p1,p2) in zip(upper_wall_points, circshift(upper_wall_points,1))
        idx0 = findfirst(p->p==(p2[1], p2[2],zmin), points)
        idx1 = findfirst(p->p==(p2[1], p2[2],zmax), points)
        idx2 = findfirst(p->p==(p1[1], p1[2],zmax), points)
        idx3 = findfirst(p->p==(p1[1], p1[2],zmin), points)
        push!(cnx, (idx0, idx1, idx2,idx3))
    end

    # ceiling
    above_pillar_points = Tuple{Float64,Float64,Float64}[]
    for x = [-7.5, -2.5, 2.5, 7.5]
        for y in [-7.5, -2.5, 2.5, 7.5]
            push!(above_pillar_points, (x,y,zmax))
        end
    end
    append!(points, above_pillar_points)
    # now connect them
    xpoints = [xmin, -7.5, -2.5, 2.5, 7.5, xmax]
    ypoints = [ymin, -7.5, -2.5, 2.5, 7.5, ymax]

    for (x1,x2) in zip(xpoints[1:end-1], xpoints[2:end]) 
        for (y1,y2) in zip(ypoints[1:end-1], ypoints[2:end])
            idx0 = findfirst(p->p==(x1,y1,zmax), points)
            idx1 = findfirst(p->p==(x1,y2,zmax), points)
            idx2 = findfirst(p->p==(x2,y2,zmax), points)
            idx3 = findfirst(p->p==(x2,y1,zmax), points)
            push!(cnx, (idx0,idx1,idx2,idx3))
        end
    end

    #hackish
    cnx = unique(cnx)
    points, cnx
end

function get_maze_mesh(args...;nrefinements=3, kwargs...)
    points,cnx = maze_topology3(args...;kwargs...)
    mm = SimpleMesh(points, connect.(cnx))
    for i in 1:nrefinements
        mm = refine(mm, QuadRefinement())
    end
    #mm2 = refine(refine(refine(mm, QuadRefinement()), QuadRefinement()),QuadRefinement())
    mm2 = SimpleMesh(mm.vertices, convert(HalfEdgeTopology, mm.topology))
    mm2
end

function get_floor_and_ceiling(mm::SimpleMesh)
    zmax = coords(maximum(mm.vertices)).z.val
    ppred(p1,p2) = ((coords(p1).z.val==0.0)&&(coords(p2).z.val==0.0))||((zmax > coords(p1).z.val > 0.0)&&(zmax > coords(p2).z.val>0.0))||((coords(p1).z.val==zmax)&&(coords(p2).z.val==zmax))
    parts = partition(mm, PointPredicatePartition(ppred))
    # the order is not consistent, but floor has the least number of elements, followed by the ceiling, and then the middle
    midx = sortperm(nelements.(parts))
    m_floor, m_ceiling, m_middle = parts[midx]
end

function count_on_manifold(mm::SimpleMesh, X::Matrix{T},w::AbstractVector{T}=ones(T,size(X,2))) where T <: Real
    # look for the nearest element
    kn = KNearestSearch(mm, 1)
    Z = zeros(nelements(mm))
    for (j,v) in enumerate(eachcol(X))
        idx = search(Meshes.Point(v...),kn)
        Z[idx] .+= w[j]
    end
    Z
end

function mapto(mm::SimpleMesh, x::AbstractVector{T}, y::AbstractVector{T},z::AbstractVector{T}) where T <: Real
    #get the size of a single element 
    idx = Int64[]
    kn = KNearestSearch(mm, 1)
    for _x in x
        for _y in y
            for _z in z
               _idx,dd = searchdists(Meshes.Point(_x,_y,_z), kn)
               _mm = mm[first(_idx)]
               Δ=mean(norm.(_mm.vertices .- centroid(_mm)))
               if all(dd .<= Δ)
                    push!(idx, first(_idx))
               else
                    push!(idx, 0)
               end
            end
        end
    end
    idx
end

function explore(mm::SimpleMesh;kwargs...)
    fig = Figure()
    lg = GridLayout(fig[1,1])
    explore!(lg, mm;kwargs...)
    fig
end

function explore(mm::SimpleMesh,tcolor::Vector{Vector{T}};kwargs...) where T <: Real
    fig = Figure()
    lgs = [GridLayout(fig[1,i]) for i in 1:length(tcolor)]
    for (lg,_tcolor) in zip(lgs, tcolor)
        alpha = get_alpha(_tcolor)
        explore!(lg, mm;color=_tcolor, alpha=alpha,kwargs...)
    end
    fig
end

"""
For use with MakiePan to get the layout
"""
function plotmesh(lg;kwargs...)
    lscene = LScene(lg[1,1];kwargs...)
    cb = Colorbar(lg[1,2];kwargs...)
    lscene,cb
end

function plotmesh(mm::SimpleMesh;kwargs...)
    with_theme(plot_theme) do
        fig = Figure()
        lscene = LScene(fig[1,1])
        plotmesh!(lscene, mm;kwargs...)
        display(fig)
        fig,lscene
    end
end

function plotmesh(lg, ii::Observable{Int64}, mm::SimpleMesh, color::Matrix{T};kwargs...) where T <: Real

end

function plotmesh!(lscene, mm::SimpleMesh;floor_offset=0.0, ceiling_offset=0.0,hide_ceiling=false, indicate_north=true, kwargs...)
    if (floor_offset != 0 || ceiling_offset != 0 || hide_ceiling)
        m_floor, m_ceiling, m_middle = get_floor_and_ceiling(mm)
         if floor_offset != 0
            m_floor2 = Translate(0.0, 0.0, floor_offset)(m_floor)
        else
            m_floor2 = m_floor
        end
        if ceiling_offset != 0
            m_ceiling2 = Translate(0.0, 0.0, ceiling_offset)(m_ceiling)
        else
            m_ceiling2 = m_ceiling
        end
        tqcolor = get(kwargs, :color,:lightgray) 
        if !isa(tqcolor,Observable)
            tcolor = Observable(tqcolor)
        else
            tcolor =tqcolor
        end
        kwargs = filter(k->k[1]!=:color, kwargs)
        use_color = Dict{Symbol,Observable{Any}}()
        for k in [:ceiling, :middle, :floor]
            use_color[k] = Observable(:lightgray)
        end
        cr = Observable(nothing)
        nanidx = nothing
        on(tcolor) do _tcolor
            if isa(_tcolor, AbstractArray{<:Any})
                # first filter out nans
                # need to separate into floor, middle, and ceiling
                nanidx = isnan.(_tcolor)
                qcolor = zero(_tcolor)
                qcolor .= _tcolor
                qcolor[nanidx] .= zero(eltype(qcolor))
                if eltype(tcolor) <: Real
                    cr[] = extrema(qcolor)
                else
                    cr[] = nothing
                end
                use_color[:middle][] = qcolor[m_middle.inds]
                use_color[:floor][] = qcolor[m_floor2.inds]
                use_color[:ceiling][] = qcolor[m_ceiling2.inds]
            else
                nanidx = nothing
                use_color[:middle][] = _tcolor 
                use_color[:floor][] = _tcolor 
                use_color[:ceiling][] = _tcolor 
                cr[] = nothing
            end
        end
        # kind of dumb;trigger a change
        tcolor[] = tcolor.val
        #talpha = get(kwargs, :alpha, fill(1.0, length(tcolor[])))
        talpha = get(kwargs, :alpha, nothing)
        if talpha === nothing
            if isa(tcolor[], AbstractVector{<:Real})
                talpha = fill(1.0, length(tcolor[]))
            end
        end
        kwargs = filter(k->k[1]!=:alpha, kwargs)
        use_alpha = Dict{Symbol,Any}()
        if isa(talpha, AbstractVector{<:Real})
            if nanidx !== nothing
                qalpha = fill!(similar(talpha), one(eltype(talpha)))
                qalpha .= talpha
                qalpha[nanidx] .= zero(eltype(qalpha)) 
            else
                qalpha = talpha
            end
            use_alpha[:middle] = qalpha[m_middle.inds]
            use_alpha[:floor] = qalpha[m_floor2.inds]
            use_alpha[:ceiling] = qalpha[m_ceiling2.inds]
        else
            use_alpha[:middle] = talpha
            use_alpha[:floor] = talpha
            use_alpha[:ceiling] = talpha
        end
        viz!(lscene, m_middle;color=use_color[:middle],alpha=use_alpha[:middle], colorrange=cr, kwargs...)
        viz!(lscene, m_floor2;color=use_color[:floor],alpha=use_alpha[:floor], colorrange=cr, kwargs...)
        if !hide_ceiling
            viz!(lscene, m_ceiling2;color=use_color[:ceiling],alpha=use_alpha[:ceiling], colorrange=cr, kwargs...)
        end
        if indicate_north
            arrows3d!(lscene, Point3f(0.0, 10.0, 7.0+ceiling_offset), Point3f(0.0, 5.0, 0.0), color=:black)
        end
    else
       viz!(lscene, mm;kwargs...) 
    end
end

function explore!(fig, mm::SimpleMesh;floor_offset=0.0, ceiling_offset=0.0, color=fill(0.0, nelements(mm)),alpha=fill(1.0, nelements(mm)),label::String="", show_axis=true, show_colorbar=true, kwargs...)
    # used for collision
    tcolor = zeros(length(color))
    tcolor .= color
    qidx = (!isfinite).(tcolor)
    alpha[qidx] .= 0.0
    tcolor[qidx] .= 0.0
    kn = KNearestSearch(mm, 1)
    lscene = LScene(fig[1,1],show_axis=show_axis)
    cm=get(Dict(kwargs), :colormap, :viridis)
    if show_colorbar
        cb = Colorbar(fig[1,2]; limits=extrema(filter(isfinite, tcolor)), colormap=cm,label=label)
    end
    plotmesh!(lscene, mm;color=color, alpha=alpha,floor_offset=floor_offset, ceiling_offset=ceiling_offset,kwargs...)
    #set up camera
    lookat = Point3f(1.0, 0.0, 0.7)
    cc = cameracontrols(lscene.scene)
    lookat0 = cc.lookat[]
    eyepos0 = cc.eyeposition[]
    upvector0 = cc.upvector[]
    near = cc.near[]
    far = cc.far[]

    eyepos = Point3f(0.0, 0.0, 0.7)
    v = lookat - eyepos 
    v = v./norm(v)
    #translate_cam!(lscene.scene, cc, Point3f(0.0, 0.0,2.5))
    update_cam!(lscene.scene, eyepos, lookat)
    fp = false
    on(events(lscene.scene).keyboardbutton, priority=20) do event
        if ispressed(lscene.scene, Keyboard.c)
                fp = ~fp
                if fp
                    cc.eyeposition[] = eyepos 
                    cc.lookat[] = eyepos + Makie.Vec(1.0, 0.0, 0.0)
                    cc.upvector[] = Makie.Vec(0.0, 0.0, 1.0)
                else
                    cc.eyeposition[] = eyepos0 
                    cc.lookat[] = lookat0
                    cc.upvector[] = upvector0
                    cc.near[] = near
                    cc.far[] = far
                end
                update_cam!(lscene.scene, cc)
        end
        if ispressed(lscene.scene, Keyboard.up) || ispressed(lscene.scene, Keyboard.down)
            if ispressed(lscene.scene, Keyboard.up)
                dx = Point3f(0.0, 0.0, -0.1)
            else
                dx = Point3f(0.0, 0.0, 0.1)
            end
            if fp
                translate_cam!(lscene.scene, cc, dx)
                #check for collision
                eyepos = cc.eyeposition[]
                qq = cc.eyeposition[]
                idx,dd = searchdists(Meshes.Point(qq...), kn)
                if dd[1] <= 0.1*Unitful.m
                    # move back
                    # TODO: This doesn't quite work, but maybe we don't care
                    translate_cam!(lscene.scene, cc, -dx)
                    # last coordinate if foward movement (for some inexplicable reason))
                end
            end
        end
        if ispressed(lscene.scene, Keyboard.right)
            rotate_cam!(lscene.scene, cc, Point3f(0.0, -0.1, 0.0))
            return Consume()
        end
        if ispressed(lscene.scene, Keyboard.left)
            rotate_cam!(lscene.scene,cc, Point3f(0.0, 0.1, 0.0))
            return Consume()
        end
    end
    lscene
end

"""
Scatter the points represented by `X` onto the mesh represented by `mm`
"""
function Makie.convert_arguments(::Type{<:Scatter}, X::Matrix{T}, mm::SimpleMesh) where T <: Real
    m_floor, m_ceiling, m_middle = get_floor_and_ceiling(mm)
    points = Tuple.(eachcol(X))
    # map to manifold
    kidx = mapto(mm, points, (zero(T), zero(T), zero(T));Δmax=100*Unitful.m)
    # find ceiling and floor points
    ceiling_idx = findall(in(m_ceiling.inds), first.(kidx))
    floor_idx = findall(in(m_floor.inds), first.(kidx))
    # offset ceiling and floor
    floor_offset = -20
    ceiling_offset = 10
    # copy X
    X2 = zeros(T, size(X)...)
    X2 .= X
    X2[3,ceiling_idx] .+= ceiling_offset
    X2[3,floor_idx] .+= floor_offset
    S.Scatter(Point3f.(eachcol(X2)))
end

function map_from_matlab(pillar_height=2.5f0)
    xbins = range(-12.5f0, stop=12.5f0, length=40)
    ybins = xbins
    zbins = range(0.0f0, stop=5.0f0, length=8)

    # each wall is 40 × 8
    #start at  bottom left corner, move clockwise
    wall_idx = reverse(permutedims(reshape(3203:3203+1280-1, 40*4,8)),dims=1)
    wall_points = NTuple{3, Float32}[]
    # move counter clockwise
    for zb in zbins 
        for xb in ybins
            push!(wall_points, (xb, ybins[1], zb))
        end
        for yb in ybins
            push!(wall_points, (xbins[end], yb, zb))
        end
        for xb in reverse(xbins) 
            push!(wall_points, (xb, ybins[end], zb))
        end
        for yb in reverse(ybins)
            push!(wall_points, (xbins[1], yb, zb))
        end
    end

    ceiling_idx = reverse(permutedims(reshape(1603:1603+1600-1, 40, 40)),dims=1)
    ceiling_points = NTuple{3, Float32}[]
    for xb in xbins
        for yb in ybins
            push!(ceiling_points, (xb,yb, last(zbins)))
        end
    end
    floor_idx = reverse(permutedims(reshape(3:3+1600-1, 40, 40)))
    floor_points = NTuple{3, Float32}[]
    for xb in xbins
        for yb in ybins
            push!(floor_points, (xb,yb, first(zbins)))
        end
    end

    # pillars
    pillar_points = Vector{Vector{NTuple{3, Float32}}}(undef, 4)
    p1_br_idx = reverse(permutedims(reshape(4483:4483+160-1, 8*4, 5)),dims=1) 
    p2_bl_idx = reverse(permutedims(reshape(4643:4643+160-1, 8*4, 5)),dims=1)
    p3_tr_idx = reverse(permutedims(reshape(4803:4803+160-1, 8*4, 5)),dims=1)
    p4_tl_idx = reverse(permutedims(reshape(4963:4963+160-1, 8*4, 5)),dims=1)

    lower_left = [(-7.5, -7.5),(-7.5, 2.5), (2.5, 2.5), (2.5, -7.5)]
    for (i,ll) in enumerate(lower_left)
        _points = NTuple{3,Float32}[]
        _ybins = range(ll[2], stop=ll[2]+5.0f0, length=8)
        _xbins = range(ll[1], stop=ll[1]+5.0f0, length=8)
        for zb in range(0.0f0, stop=pillar_height, length=5)
            for xb in _xbins
                push!(_points, (xb, _ybins[1],zb))
            end
            for yb in _ybins
                push!(_points, (_xbins[end], yb, zb))
            end
            for xb in reverse(_xbins) 
                push!(_points, (xb, _ybins[end], zb))
            end
            for yb in reverse(_ybins)
                push!(_points, (xbins[1], yb, zb))
            end
        end
        pillar_points[i] = _points
    end
     
    nn = length(wall_idx) + length(ceiling_idx) + length(floor_idx) + length(p1_br_idx) + length(p2_bl_idx) + length(p3_tr_idx) + length(p4_tl_idx)
    wall_points, ceiling_points, floor_points, pillar_points, nn
    allidx = [vec(floor_idx);vec(ceiling_idx);vec(wall_idx);vec(p1_br_idx);vec(p2_bl_idx);vec(p3_tr_idx);vec(p4_tl_idx)]
    [floor_points;ceiling_points;wall_points;pillar_points...], allidx
end