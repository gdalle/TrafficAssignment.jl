@kwdef struct DijkstraStorage{W,T,H<:BinaryHeap}
    heap::H
    parents::Vector{T}
    edge_ids::Vector{T}
    dists::Vector{W}
end

function DijkstraStorage(g::SimpleWeightedDiGraph)
    T, W = vertex_type(g), weight_type(g)
    heap = BinaryHeap(Base.By(last), Pair{T,W}[])
    sizehint!(heap, nv(g))
    parents = Vector{T}(undef, nv(g))
    edge_ids = Vector{T}(undef, nv(g))
    dists = Vector{W}(undef, nv(g))
    return DijkstraStorage(; heap, parents, edge_ids, dists)
end

function reset!(storage::DijkstraStorage{W,T}) where {W,T}
    (; heap, parents, edge_ids, dists) = storage
    empty!(heap.valtree)  # internal, will be released in DataStructures v0.19
    fill!(parents, zero(T))
    fill!(edge_ids, zero(T))
    fill!(dists, typemax(W))
    return nothing
end

function dijkstra!(
    storage::DijkstraStorage,
    g::SimpleWeightedDiGraph,
    s::Integer;
    forbidden_intermediate_vertices,
)
    T, W = vertex_type(g), weight_type(g)
    reset!(storage)
    (; heap, parents, edge_ids, dists) = storage
    # Add source
    push!(heap, s => zero(W))
    # Main loop
    while !isempty(heap)
        u, Δu = pop!(heap)
        if u != s && u in forbidden_intermediate_vertices
            continue
        end
        if Δu <= dists[u]
            dists[u] = Δu
            for (v, w_uv, e_uv) in outneighbors_weights_edgeids(g, u)
                if Δu + w_uv < dists[v]
                    parents[v] = u
                    edge_ids[v] = e_uv
                    dists[v] = Δu + w_uv
                    push!(heap, v => Δu + w_uv)
                end
            end
        end
    end
end

function dijkstra(g::SimpleWeightedDiGraph, s::Integer; forbidden_intermediate_vertices=())
    storage = DijkstraStorage(g)
    dijkstra!(storage, g, s; forbidden_intermediate_vertices)
    return storage
end
