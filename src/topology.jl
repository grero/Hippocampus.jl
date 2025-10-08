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

function maze_topology3(xmin=-12.5, xmax=12.5, ymin=xmin, ymax=xmax, zmin=0.0, zmax=5.0)
    points = [(xmin, ymin, 0.0), (-7.5, ymin,0.0),(-7.5, -7.5,0.0),(xmin, -7.5, 0.0)] 
    cnx = [(1,2,3,4)]
    append!(points, [(-7.5, -2.5, 0.0),(xmin, -2.5, 0.0)])
    push!(cnx, (4,3,5,6))
    append!(points, [(xmin, 2.5, 0.0),(-7.5, 2.5, 0.0)])
    push!(cnx, (5,6,7,8))
    append!(points, [(-7.5, 7.5, 0.0), (xmin, 7.5, 0.0)])
    push!(cnx, (7,8,9,10))
    append!(points, [(xmin,ymax, 0.0),(-7.5, ymax, 0.0)])
    push!(cnx, (9,10,11,12))
    append!(points, [(-2.5, 7.5, 0.0),(-2.5, ymax, 0.0)])
    push!(cnx, (9,13, 14, 11))
    append!(points, [(2.5, 7.5, 0.0),(2.5, ymax, 0.0)])
    push!(cnx, (13,15,16,14))
    append!(points, [(7.5, 7.5, 0.0), (7.5, ymax, 0.0)])
    push!(cnx, (15,17,18,16))
    append!(points, [(xmax, 7.5, 0.0), (xmax, ymax, 0.0)])
    push!(cnx, (17,19,20,18))
    append!(points, [(7.5, 2.5, 0.0), (xmax, 2.5, 0.0)])
    push!(cnx, (21,22,19,17))
    append!(points, [(7.5, -2.5, 0.0), (xmax, -2.5, 0.0)])
    push!(cnx, (23,24,22,21))
    append!(points, [(7.5, -7.5, 0.0), (xmax, -7.5, 0.0)])
    push!(cnx, (25,26,24,23))
    append!(points, [(7.5, ymin, 0.0),(xmax, ymin, 0.0)])
    push!(cnx, (27,28,26,25))
    append!(points, [(2.5, ymin, 0.0), (2.5, -7.5, 0.0)])
    push!(cnx, (29,27, 25,30))
    append!(points, [(7.5, ymin, 0.0), (7.5, -7.5, 0.0)])
    push!(cnx, (31,29,30,32))
    append!(points, [(-2.5, ymin, 0.0), (-2.5, -7.5, 0.0)])
    push!(cnx, (33, 29, 32,34))
    push!(cnx, (2,33,34,3))
    append!(points, [(-2.5, -2.5, 0.0), (2.5, -2.5, 0.0)])
    push!(cnx, (34, 35, 36,30))
    append!(points, [(-2.5, 2.5, 0.0), (2.5, 2.5, 0.0)])
    push!(cnx, (37, 13, 15,38)) 
    #append!(points, [()])
    #connect the inner points
    push!(cnx, (35, 36,38,37))
    push!(cnx, (5, 35, 37,8))
    push!(cnx, (36, 23, 21, 38))

    #add pillars
    # first pillar
    append!(points, [(-7.5, -7.5, 2.5), (-7.5, -2.5, 2.5)])
    append!(points, [(-2.5, -7.5, 2.5), (-2.5, -2.5, 2.5)])
    push!(cnx, (3, 39, 40, 5))
    push!(cnx, (34, 41, 42,35))
    push!(cnx, (3, 39, 41, 34))
    push!(cnx, (5, 40, 42,35))

    #second pillar
    append!(points, [(2.5, -7.5, 2.5), (2.5, -2.5, 2.5)]) 
    append!(points, [(7.5, -7.5, 2.5), (7.5, -2.5, 2.5)]) 
    push!(cnx, (30, 43, 44,36))
    push!(cnx, (25, 45, 46, 23))
    push!(cnx, (30,43, 45, 25))
    push!(cnx, (36,44,46,23))

    #third pillar
    append!(points, [(2.5, 2.5, 2.5), (2.5, 7.5, 2.5)]) 
    append!(points, [(7.5, 2.5, 2.5), (7.5, 7.5, 2.5)]) 
    push!(cnx, (38,47,48,15))
    push!(cnx, (21,49,50,17))
    push!(cnx, (38,47,49,21))
    push!(cnx, (15, 48,50,17))

    #fourth pillar
    append!(points, [(-7.5, 2.5, 2.5),(-7.5, 7.5, 2.5)])
    append!(points, [(-2.5, 2.5, 2.5),(-2.5, 7.5,2.5)])
    push!(cnx, (8,51,52,9))
    push!(cnx, (37,53,54,13)) 
    push!(cnx, (8,51,53,37))
    push!(cnx, (9,52,54,13))

    #now the walls
    upper_wall_points = [(xmin,ymin, zmax),(xmin, -7.5, zmax),(xmin, -2.5, zmax), (xmin, 2.5, zmax),
                     (xmin, 7.5, zmax), (xmin, ymax, zmax), (-7.5, ymax,zmax), (-2.5, ymax,zmax),
                     (2.5, ymax,zmax),(7.5, ymax,zmax), (xmax, ymax,zmax), (xmax, 7.5,zmax),
                     (xmax, 2.5,zmax),(xmax, -2.5,zmax),(xmax, -7.5,zmax), (xmax, ymin,zmax),
                     (7.5, ymin,zmax),(2.5, ymin,zmax), (-2.5, ymin,zmax),(-7.5, ymin,zmax)]
    append!(points, upper_wall_points)
    #maybe be a bit smarter here
    for (p1,p2) in zip(upper_wall_points, circshift(upper_wall_points,1))
        idx0 = findfirst(p->p==(p1[1], p1[2],zmin), points)
        idx1 = findfirst(p->p==(p1[1], p1[2],zmax), points)
        idx2 = findfirst(p->p==(p2[1], p2[2],zmax), points)
        idx3 = findfirst(p->p==(p2[1], p2[2],zmin), points)
        push!(cnx, (idx0, idx1, idx2,idx3))
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