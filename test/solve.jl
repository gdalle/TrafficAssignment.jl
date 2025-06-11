using TestItems

@testitem "Comparing results" begin
    import TrafficAssignment as TA
    using LinearAlgebra

    reldist(a, b) = norm(a - b) / norm(a)

    @testset "Sioux Falls" begin
        problem = TrafficAssignmentProblem("TransportationNetworks", "SiouxFalls")
        (; optimal_flow) = problem
        flow = solve_frank_wolfe(problem; verbose=false, max_iteration=1_000)
        @test reldist(optimal_flow, flow) < 1e-3
        @test TA.wardrop_objective(problem, flow) <
            1.05 * TA.wardrop_objective(problem, optimal_flow)
        @test price_of_anarchy(problem; verbose=false) > 1
    end

    @testset "Anaheim" begin
        problem = TrafficAssignmentProblem("TransportationNetworks", "Anaheim")
        (; optimal_flow) = problem
        flow = solve_frank_wolfe(problem; verbose=false, max_iteration=100)
        @test_broken reldist(optimal_flow, flow) < 1e-2
        @test TA.wardrop_objective(problem, flow) <
            1.05 * TA.wardrop_objective(problem, optimal_flow)
        @test price_of_anarchy(problem; verbose=false) > 1
    end
end
