"""
    outneighbors_and_weights(g::SimpleWeightedDiGraph, u::Integer)

Efficiently list the tuples `(v, w)` where `v` is an outneighbor of `u` and `w` is the weight of the corresponding edge.
"""
function outneighbors_and_weights(g::SimpleWeightedDiGraph, u::Integer)
    w = g.weights
    (; rowval, colptr, nzval) = w
    interval = colptr[u]:(colptr[u + 1] - 1)
    return zip(view(rowval, interval), view(nzval, interval))
end

weight_type(::SimpleWeightedDiGraph{T,W}) where {T,W} = W

struct PathNotFoundError <: Exception
    s::Int
    t::Int
end

function Base.showerror(io::IO, e::PathNotFoundError)
    return print(io, "PathNotFoundError: No path found from $(e.s) to $(e.t)")
end

function dijkstra!(
    heap::BinaryHeap,
    parents::AbstractVector,
    dists::AbstractVector,
    g::SimpleWeightedDiGraph,
    s::Integer,
)
    T, W = eltype(g), eltype(dists)
    empty!(heap.valtree)  # internal, will be released in DataStructures v0.19
    fill!(parents, zero(T))
    fill!(dists, typemax(W))
    # Add source
    push!(heap, s => zero(W))
    # Main loop
    while !isempty(heap)
        u, Δu = pop!(heap)
        if Δu <= dists[u]
            dists[u] = Δu
            for (v, w_uv) in outneighbors_and_weights(g, u)
                if Δu + w_uv < dists[v]
                    parents[v] = u
                    dists[v] = Δu + w_uv
                    push!(heap, v => Δu + w_uv)
                end
            end
        end
    end
end

function astar!(
    heap::BinaryHeap,
    parents::AbstractVector,
    dists::AbstractVector,
    path::AbstractVector,
    g::SimpleWeightedDiGraph,
    s::Integer,
    t::Integer,
    heuristic::AbstractVector,
)
    T, W = eltype(g), eltype(dists)
    empty!(heap.valtree)  # internal, will be released in DataStructures v0.19
    fill!(parents, zero(T))
    fill!(dists, typemax(W))
    fill!(path, zero(T))
    # Add source
    dists[s] = zero(W)
    push!(heap, s => heuristic[s])
    # Main loop
    while !isempty(heap)
        u, _ = pop!(heap)
        Δu = dists[u]
        if u == t
            k = 1
            path[k] = t
            while parents[u] != zero(T)
                u = parents[u]
                k += 1
                path[k] = u
            end
            return nothing
        else
            for (v, w_uv) in outneighbors_and_weights(g, u)
                heuristic[v] == typemax(W) && continue
                if Δu + w_uv < dists[v]
                    parents[v] = u
                    dists[v] = Δu + w_uv
                    push!(heap, v => Δu + w_uv + heuristic[v])
                end
            end
        end
    end
    throw(PathNotFoundError(s, t))
end

function init_dijkstra(g::SimpleWeightedDiGraph)
    T = eltype(g)
    W = weight_type(g)
    heap = BinaryHeap(Base.By(last), Pair{T,W}[])
    sizehint!(heap, nv(g))
    parents = Vector{T}(undef, nv(g))
    dists = Vector{W}(undef, nv(g))
    return (; heap, parents, dists)
end

function init_astar(g::SimpleWeightedDiGraph)
    T = eltype(g)
    W = weight_type(g)
    heap = BinaryHeap(Base.By(last), Pair{T,W}[])
    sizehint!(heap, nv(g))
    parents = Vector{T}(undef, nv(g))
    dists = Vector{W}(undef, nv(g))
    path = Vector{T}(undef, nv(g))
    return (; heap, parents, dists, path)
end

"""
    dijkstra(g::SimpleWeightedDiGraph, s::Integer)

Apply Dijkstra's algorithm to graph `g` starting from source vertex `s`.

Return a tuple of vectors `(parents, dists)`, where

  - `parents` maps each vertex to its predecessor in the shortest path tree
  - `dists` contains the distances from the source
"""
function dijkstra(g::SimpleWeightedDiGraph, s::Integer)
    (; heap, parents, dists) = init_dijkstra(g)
    dijkstra!(heap, parents, dists, g, s)
    return parents, dists
end

"""
    astar(g::SimpleWeightedDiGraph, s::Integer)

Apply the A* algorithm to graph `g` for a path from `s` to `t`, with a given `heuristic`.

Return a tuple `(path, dist)` where

  - `path` is a vector containing the vertices on a shortest path, ordered from `t` to `s`
  - `dist` is a scalar indicating the cost of that path
"""
function astar(g::SimpleWeightedDiGraph, s::Integer, t::Integer, heuristic::AbstractVector)
    (; heap, parents, dists, path) = init_astar(g)
    astar!(heap, parents, dists, path, g, s, t, heuristic)
    k = 1
    while path[k] != s
        k += 1
    end
    return view(path, 1:k), dists[t]
end
