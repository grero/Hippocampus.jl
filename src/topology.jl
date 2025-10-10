using Meshes
using Meshes: connect
using Graphs
using Graphs: DijkstraState

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

function ballsearch(mm::SimpleMesh)
    A = adjacencymatrix(mm;rank=0)
    G = SimpleGraph(A)
    function searchball(ii::Integer, r::Integer)
        dj = dijkstra_shortest_paths(G, ii, A;trackvertices=false)
        idx = findall(dj.dists .<= r)
    end
end

function distancematrix(mm::SimpleMesh)
    nn = length(mm.vertices)
    A = adjacencymatrix(mm;rank=0)
    G = SimpleGraph(A)
    D = zeros(nn,nn)
    for ii in 1:nn
        dj = dijkstra_shortest_paths(G, ii;trackvertices=false)
        D[:,ii] = dj.dists 
    end
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
    append!(points, [(7.5, ymin, 0.0), (7.5, -7.5, 0.0)])
    _cnx = (31,29,30,32)
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
    points, cnx
end

function count_on_manifold(mm::SimpleMesh, X::Matrix{T}) where T <: Real
    # look for the nearest element
    kn = KNearestSearch(mm, 1)
    Z = zeros(nelements(mm))
    for v in eachcol(X)
        idx = search(Meshes.Point(v...),kn)
        Z[idx] .+= 1.0
    end
    Z
end