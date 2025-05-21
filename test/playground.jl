using TrafficAssignment
import TrafficAssignment as TA
using FrankWolfe
using SimpleWeightedGraphs, Graphs
using Test, LinearAlgebra, SparseArrays

dataset = "UnifiedTrafficDataset"
city = "New York"
@profview pb = TrafficAssignmentProblem(dataset, city)

@time solve_frank_wolfe(pb; max_iteration=5, print_iter=1, verbose=true);
@profview solve_frank_wolfe(pb; max_iteration=500, print_iter=10, verbose=true);

problem = TrafficAssignmentProblem("TransportationNetworks", "SiouxFalls")
(; optimal_flow) = problem
flow = solve_frank_wolfe(problem; verbose=true, max_iteration=1_000, print_iter=100)
@test reldist(optimal_flow, flow) < 1e-3
@test TA.objective(problem, flow) < 1.05 * TA.objective(problem, optimal_flow)
