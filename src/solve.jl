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
    (; demand, destination_free_flow_time) = problem
    cost = sparse_by_link(problem, cost_vec)
    graph = SimpleWeightedDiGraph(cost)
    I, J, C = findnz(cost)
    flow_vec = zero(C)
    edge_index = Dict((i, j) => e for (e, (i, j)) in enumerate(zip(I, J)))
    od_pairs = collect(keys(demand))
    origins = unique(first.(od_pairs))
    destinations = unique(last.(od_pairs))
    # Dijkstra version
    @tasks for o in origins
        @local (; heap, parents, dists) = init_dijkstra(graph)
        dijkstra!(heap, parents, dists, graph, o)
        for d in destinations
            if haskey(demand, (o, d))
                dem = demand[o, d]
                v = d
                while v != o
                    u = parents[v]
                    e = edge_index[u, v]
                    Atomix.@atomic flow_vec[e] += dem
                    v = parents[v]
                end
            end
        end
    end
    # A* version, much slower
    #=
    @tasks for (o, d) in od_pairs
        @local (; heap, parents, dists, path) = init_astar(graph)
        dem = demand[o, d]
        heuristic = destination_free_flow_time[d]
        astar!(heap, parents, dists, path, graph, o, d, heuristic)
        for k in eachindex(path)[1:(end - 1)]
            v, u = path[k], path[k + 1]
            if u == 0
                break
            else
                e = edge_index[u, v]
                Atomix.@atomic flow_vec[e] += dem
            end
        end
    end
    =#
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
