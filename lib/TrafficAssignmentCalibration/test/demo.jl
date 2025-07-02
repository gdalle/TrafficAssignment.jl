using Accessors
using ImplicitDifferentiation
using DifferentiableFrankWolfe
using FrankWolfe
using SparseArrays
using TrafficAssignment
import TrafficAssignment as TA
using Zygote
using ForwardDiff

problem = TrafficAssignmentProblem(TransportationNetworks, "SiouxFalls")

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

dfw = DiffFW(
    f,
    f_grad1,
    lmo;
    representation=OperatorRepresentation(),
    linear_solver=IterativeLinearSolver(),
)
callback = TA.RelativeGapCallback(1e-4)

x0 = nonzeros(flow)
θ = nonzeros(problem.link_capacity)
θ_and_dθ = ForwardDiff.Dual.(θ, 1e3 .* rand(length(θ)))

function todiff(_θ; kwargs...)
    _f_vec = dfw(_θ, x0; callback, kwargs...)
    return TA._objective_vec(problem, Val(:centralized), _f_vec)
end

todiff(θ_and_dθ; away_steps=false)

∇θ = ForwardDiff.gradient(todiff, θ)
∇θ = Zygote.gradient(todiff, θ)

capa_grad = Matrix(TA.sparse_by_link(problem, ∇θ))
