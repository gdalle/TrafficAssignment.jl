mynz(A::SparseMatrixCSC) = nonzeros(A)
mynz(a::Number) = a

function sparse_by_link(problem::TrafficAssignmentProblem, nzval::AbstractVector)
    A = problem.link_free_flow_time
    return SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, nzval)
end

@inline function _edge_obj(::Val{:equilibrium}, f, t0, α, β, c)
    return t0 * (f + (α * c / (β + 1)) * (f / c)^(β + 1))
end
@inline function _edge_grad(::Val{:equilibrium}, f, t0, α, β, c)
    return t0 * (1 + α * (f / c)^β)
end

@inline _edge_obj(::Val{:centralized}, f, t0, α, β, c) = t0 * (1 + α * (f / c)^β) * f
@inline _edge_grad(::Val{:centralized}, f, t0, α, β, c) = t0 * (1 + α * (β + 1) * (f / c)^β)

function _objective(
    problem::TrafficAssignmentProblem, control::Val, flow_vec::AbstractVector
)
    (; link_free_flow_time, link_bpr_mult, link_capacity, link_bpr_power) = problem
    t0, c = mynz(link_free_flow_time), mynz(link_capacity)
    α, β = mynz(link_bpr_mult), mynz(link_bpr_power)
    f = flow_vec
    s = zero(eltype(flow_vec))
    for i in eachindex(t0, c, f)
        s += _edge_obj(control, f[i], t0[i], α[i], β[i], c[i])
    end
    return s
end

function _gradient!(
    g::AbstractVector,
    problem::TrafficAssignmentProblem,
    control::Val,
    flow_vec::AbstractVector,
)
    (; link_free_flow_time, link_bpr_mult, link_capacity, link_bpr_power) = problem
    t0, c = mynz(link_free_flow_time), mynz(link_capacity)
    α, β = mynz(link_bpr_mult), mynz(link_bpr_power)
    f = flow_vec
    g .= _edge_grad.(control, f, t0, α, β, c)
    return nothing
end

struct ShortestPathOracle{P<:TrafficAssignmentProblem} <:
       FrankWolfe.LinearMinimizationOracle
    problem::P
end

function FrankWolfe.compute_extreme_point(
    spo::ShortestPathOracle, cost_vec::AbstractVector; kwargs...
)
    yield()
    (; problem) = spo
    (; link_id, demand, origins, destinations) = problem
    cost = sparse_by_link(problem, cost_vec)
    graph = SimpleWeightedDiGraph(cost, link_id)
    flow_vec = fill(one(float(valtype(demand))), length(cost_vec))
    @tasks for o in origins
        @local storage = DijkstraStorage(graph)
        dijkstra!(storage, graph, o)
        @one_by_one begin
            for d in destinations
                if haskey(demand, (o, d))
                    dem = demand[o, d]
                    v = d
                    while v != o
                        e = storage.edge_ids[v]
                        flow_vec[e] += dem
                        v = storage.parents[v]
                    end
                end
            end
        end
    end
    return flow_vec
end

"""
    solve_frank_wolfe(
        problem::TrafficAssignmentProblem,
        control::Val,
        frank_wolfe_alg=frank_wolfe;
        verbose, relative_gap, kwargs...
    )

Solve a traffic assignment `problem` with a specific `control` setting (either `:centralized` or `:equilibrium`), using an algorithm from the FrankWolfe library.

Remaining keyword arguments are passed to `frank_wolfe_alg`.
"""
function solve_frank_wolfe(
    problem::TrafficAssignmentProblem,
    control::Val=Val(:equilibrium),
    frank_wolfe_alg::A=frank_wolfe;
    verbose::Bool=true,
    relative_gap::Real=1e-4,
    kwargs...,
) where {A}
    lmo = ShortestPathOracle(problem)
    obj(flow_vec) = _objective(problem, control, flow_vec)
    grad!(g, flow_vec) = _gradient!(g, problem, control, flow_vec)
    direction_init = zeros(Float64, nb_links(problem))
    flow_init = FrankWolfe.compute_extreme_point(lmo, direction_init)
    callback = RelativeGapCallback(relative_gap)
    flow_opt, _ = frank_wolfe_alg(obj, grad!, lmo, flow_init; verbose, callback, kwargs...)
    return sparse_by_link(problem, flow_opt)
end

struct RelativeGapCallback
    relative_gap::Float64
end

function (callback::RelativeGapCallback)(state, args...)
    return state.dual_gap / state.primal > callback.relative_gap
end

"""
    social_cost(problem::TrafficAssignmentProblem, flow::SparseMatrixCSC)

Return the social cost (total travel time) of a given solution `flow` to `problem`.
"""
function social_cost(problem::TrafficAssignmentProblem, flow::SparseMatrixCSC)
    return _objective(problem, Val(:centralized), nonzeros(flow))
end

function wardrop_objective(problem::TrafficAssignmentProblem, flow::SparseMatrixCSC)
    return _objective(problem, Val(:equilibrium), nonzeros(flow))
end

"""
    price_of_anarchy(problem::TrafficAssignmentProblem, args...; kwargs...)

Compute the price of anarchy (cost of equilibrium flow over cost of centralized flow) for `problem`.

Positional and keyword arguments are forwarded to [`solve_frank_wolfe`](@ref).
"""
function price_of_anarchy(problem::TrafficAssignmentProblem, args...; kwargs...)
    flow_centralized = solve_frank_wolfe(problem, Val(:centralized), args...; kwargs...)
    flow_equilibrium = solve_frank_wolfe(problem, Val(:equilibrium), args...; kwargs...)
    cost_centralized = social_cost(problem, flow_centralized)
    cost_equilibrium = social_cost(problem, flow_equilibrium)
    return cost_equilibrium / cost_centralized
end
