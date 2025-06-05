using TestItems

@testitem "Plotting" begin
    using CairoMakie, Tyler
    @testset "Makie" begin
        pb = TrafficAssignmentProblem("TransportationNetworks", "SiouxFalls")
        plot_network(pb, pb.optimal_flow; nodes=true)
    end
    @testset "Tyler" begin
        pb = TrafficAssignmentProblem("UnifiedTrafficDataset", "San Francisco")
        plot_network(pb; tiles=true, zones=true)
    end
end
