using TrafficAssignment
using SimpleWeightedGraphs
using Graphs
using SparseArrays
using Test

@testset "Dijkstra" begin
    for p in 0.01:0.01:0.1
        w = sprand(100, 100, p)
        id = similar(w, Int)

        g1 = SimpleWeightedGraphs.SimpleWeightedDiGraph(w)
        d1 = Graphs.dijkstra_shortest_paths(g1, 1).dists

        g2 = TrafficAssignment.SimpleWeightedDiGraph(w, id)
        d2 = TrafficAssignment.dijkstra(g2, 1).dists
        @test d1 == d2
    end
end
