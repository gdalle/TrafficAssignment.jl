using TestItems

@testset "Plotting" begin
    pb = TrafficAssignmentProblem("TransportationNetworks", "SiouxFalls")
    plot_network(pb, pb.optimal_flow; nodes=true)
    if VERSION >= v"1.11"
        using Pkg
        Pkg.add("Tyler")
        using Tyler
        pb = TrafficAssignmentProblem("UnifiedTrafficDataset", "San Francisco")
        plot_network(pb; tiles=true, zones=true)
    end
end
