function outneighbors_and_weights(g::SimpleWeightedDiGraph, u::Integer)
    w = g.weights
    (; rowval, colptr, nzval) = w
    interval = colptr[u]:(colptr[u + 1] - 1)
    return zip(view(rowval, interval), view(nzval, interval))
end

weight_type(::SimpleWeightedDiGraph{T,W}) where {T,W} = W

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
    empty!(heap.valtree)  # internal, will be released in DataStructures v0.19
    fill!(parents, zero(eltype(parents)))
    fill!(dists, typemax(eltype(dists)))
    fill!(path, zero(eltype(path)))
    T, W = eltype(g), eltype(dists)
    # Add source
    dists[s] = zero(W)
    push!(heap, s => heuristic[s])
    # Main loop
    while !isempty(heap)
        u, _ = pop!(heap)
        Δ_u = dists[u]
        if u == t
            k = 1
            path[k] = t
            while parents[u] != zero(T)
                u = parents[u]
                k += 1
                path[k] = u
            end
            break
        else
            for (v, Δ_uv) in outneighbors_and_weights(g, u)
                heuristic[v] == typemax(W) && continue
                if Δ_u + Δ_uv < dists[v]
                    parents[v] = u
                    dists[v] = Δ_u + Δ_uv
                    h_v = Δ_u + Δ_uv + heuristic[v]
                    push!(heap, v => h_v)
                end
            end
        end
    end
    if first(path) == 0
        error("A* could not find a path from $s to $t")
    end
end

function init_astar(g::SimpleWeightedDiGraph)
    T = eltype(g)
    W = weight_type(g)
    heap = BinaryHeap(Base.By(last), Pair{T,W}[])
    sizehint!(heap, ne(g))
    parents = Vector{T}(undef, nv(g))
    dists = Vector{W}(undef, nv(g))
    path = Vector{T}(undef, nv(g))
    return (; heap, parents, dists, path)
end

function astar(g::SimpleWeightedDiGraph, s::Integer, t::Integer, heuristic::AbstractVector)
    (; heap, parents, dists, path) = init_astar(g)
    astar!(heap, parents, dists, path, g, s, t, heuristic)
    return path
end
