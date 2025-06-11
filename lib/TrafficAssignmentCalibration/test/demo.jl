using Accessors
using DifferentiableFrankWolfe
using SparseArrays
using TrafficAssignment
import TrafficAssignment as TA
using Zygote

problem = TrafficAssignmentProblem("TransportationNetworks", "SiouxFalls")

flow = solve_frank_wolfe(problem, Val(:equilibrium))

function f(flow_vec, capa_vec)
    new_problem = @set problem.link_capacity = TA.sparse_by_link(problem, capa_vec)
    return TA._objective(new_problem, Val(:equilibrium), flow_vec)
end

function f_grad1(flow_vec, capa_vec)
    new_problem = @set problem.link_capacity = TA.sparse_by_link(problem, capa_vec)
    return TA._gradient(new_problem, Val(:equilibrium), flow_vec)
end

lmo = TA.ShortestPathOracle(problem)

dfw = DiffFW(f, f_grad1, lmo)
callback = TA.RelativeGapCallback(1e-4)

x0 = nonzeros(flow)
θ = nonzeros(problem.link_capacity)
function todiff(_θ)
    _f_vec = dfw(_θ, x0; callback)
    return TA._objective_vec(problem, Val(:centralized), _f_vec)
end

∇θ = Zygote.gradient(todiff, θ)[1]
@profview for _ in 1:1000
    Zygote.gradient(todiff, θ)[1]
end

capa_grad = TA.sparse_by_link(problem, ∇θ)
