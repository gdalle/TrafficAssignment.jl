struct SimpleWeightedDiGraph{W,T<:Integer}
    adj_trans::SparseMatrixCSC{W,T}
    id_trans::SparseMatrixCSC{T,T}

    function SimpleWeightedDiGraph(
        adj::AbstractMatrix{W}, id::AbstractMatrix{T}
    ) where {W,T}
        @assert checksquare(adj) == checksquare(id)
        @assert nnz(adj) == nnz(id)
        return new{W,T}(sparse(transpose(adj)), sparse(transpose(id)))
    end
end

vertex_type(::SimpleWeightedDiGraph{W,T}) where {T,W} = T
weight_type(::SimpleWeightedDiGraph{W,T}) where {T,W} = W

nv(g::SimpleWeightedDiGraph) = size(g.adj_trans, 1)

function outneighbors_weights_edgeids(g::SimpleWeightedDiGraph, u::Integer)
    w, e = g.adj_trans, g.id_trans
    interval = w.colptr[u]:(w.colptr[u + 1] - 1)
    return zip(  #
        view(w.rowval, interval),
        view(w.nzval, interval),
        view(e.nzval, interval),
    )
end
