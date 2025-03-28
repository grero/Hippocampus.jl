"""
Return all unique paths between `n` fully connected vertices
"""
function get_all_unique_paths(n::Int64,s,t)
    path = Tuple{Int64,Int64}[]
    full = false
    added = false
    while !full
        added = false
        candidates = shuffle(setdiff(1:n,i))
        for j in candidates
            pp = (i,j)
            if !in(pp,path)
                push!(path, pp)
                added = true
                i = j
                break
            end
        end
        full = !added
    end

    # we don't want to end up at the same vertex as we started out from
    if path[end][2] == path[1][1]
        pop!(path)
    end
    path
end

"""
    find_path(n::T, s::T2, t::T2;path_length::Union{Int64,Nothing}=nothing,max_npaths=1000) where T <: Int64 where T2 <: Integer

The a maximum of `max_npaths` of length `path_length` paths with unique edges starting from `s` end ending up in `t`. The path can
visited all vertices, including `s` and `t` more than once.
"""
function find_path(n::T, s::T2, t::T2;path_length::Union{Int64,Nothing}=nothing,max_npaths=1) where T <: Int64 where T2 <: Integer
    result = Any[]
    if path_length === nothing
        # all possible pairs
        # the number of edges is path_length -1
        path_length = binomial(n,2)*2
        # we might not get exactly path_length vertices
        path_length = round(Int64, 0.9*path_length)
    end

    function dfs(current, path, visited_edges)
        if length(result) >= max_npaths
            return
        end
        if current == t && length(path) >= path_length
            push!(result, deepcopy(path))
            return
        end
        added = false
        for j in 1:n
            if j == current
                continue
            end
            edge = (current,j)
            if !in(edge, visited_edges)
                added = true
                push!(visited_edges, edge)
                push!(path, j)
                dfs(j, path, visited_edges)
                #free up this edge for the next path
                pop!(visited_edges, edge)
                # to get back to `current`
                pop!(path)
            end
        end
        if added == false
            # no more vertices to explore
            return
        end
    end
    visited_edges = Set{Tuple{T2,T2}}()
    dfs(s, [s], visited_edges)
    return result
end