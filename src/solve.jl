mynz(A::SparseMatrixCSC) = nonzeros(A)
mynz(a::Number) = a

function sparse_by_link(problem::TrafficAssignmentProblem, nzval::AbstractVector)
    A = problem.link_free_flow_time
    return SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, nzval)
end

function link_travel_time(problem::TrafficAssignmentProblem, flow::AbstractMatrix)
    (; link_free_flow_time, link_bpr_mult, link_capacity, link_bpr_power) = problem
    t0, c = mynz(link_free_flow_time), mynz(link_capacity)
    α, β = mynz(link_bpr_mult), mynz(link_bpr_power)
    f = mynz(flow)
    t = @. t0 * (1 + α * (f / c)^β)
    return sparse_by_link(problem, t)
end

function link_travel_time_integral(problem::TrafficAssignmentProblem, flow::AbstractMatrix)
    (; link_free_flow_time, link_bpr_mult, link_capacity, link_bpr_power) = problem
    t0, c = mynz(link_free_flow_time), mynz(link_capacity)
    α, β = mynz(link_bpr_mult), mynz(link_bpr_power)
    f = mynz(flow)
    t_integral = @. t0 * (f + (α * c / (β + 1)) * (f / c)^(β + 1))
    return sparse_by_link(problem, t_integral)
end

function objective(problem::TrafficAssignmentProblem, flow_vec::AbstractVector)
    flow = sparse_by_link(problem, flow_vec)
    return objective(problem, flow)
end

function objective(problem::TrafficAssignmentProblem, flow::AbstractMatrix)
    obj = sum(link_travel_time_integral(problem, flow))
    return obj
end

function objective_gradient(problem::TrafficAssignmentProblem, flow_vec::AbstractVector)
    flow = sparse_by_link(problem, flow_vec)
    return mynz(link_travel_time(problem, flow))
end

"""
    social_cost(problem::TrafficAssignmentProblem, flow::AbstractMatrix)

Compute the social cost induced by a matrix of link flows.
"""
function social_cost(problem::TrafficAssignmentProblem, flow::AbstractMatrix)
    return dot(link_travel_time(problem, flow), flow)
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
    flow_vec = zeros(float(valtype(demand)), length(cost_vec))
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
        frank_wolfe_alg;
        verbose, kwargs...
    )

Solve a traffic assignment problem using an algorithm from the FrankWolfe library.

Keyword arguments are passed to `frank_wolfe_alg`.
"""
function solve_frank_wolfe(
    problem::TrafficAssignmentProblem,
    frank_wolfe_alg::A=frank_wolfe;
    verbose::Bool=true,
    kwargs...,
) where {A}
    lmo = ShortestPathOracle(problem)
    f(cost_vec) = objective(problem, cost_vec)
    function grad!(storage, cost_vec)
        grad = objective_gradient(problem, cost_vec)
        return copyto!(storage, grad)
    end
    direction_init = zeros(Float64, nb_links(problem))
    flow_init = FrankWolfe.compute_extreme_point(lmo, direction_init)
    flow_opt, _ = frank_wolfe_alg(f, grad!, lmo, flow_init; verbose, kwargs...)
    return sparse_by_link(problem, flow_opt)
end
