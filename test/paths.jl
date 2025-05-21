using TrafficAssignment
using SimpleWeightedGraphs
using Graphs
using SparseArrays
using Test

function path_from_edges(edge_list)
    if isempty(edge_list)
        return Int[]
    else
        path = [src(first(edge_list))]
        for e in edge_list
            push!(path, dst(e))
        end
        return path
    end
end

@testset "Dijkstra" begin
    for p in 0.01:0.01:0.1
        w = sprand(100, 100, p)
        g = SimpleWeightedDiGraph(w)

        d1 = Graphs.dijkstra_shortest_paths(g, 1).dists
        d2 = TrafficAssignment.dijkstra(g, 1)[2]
        @test d1 == d2
    end
end

@testset "A*" begin
    for p in 0.01:0.01:0.1
        w = sprand(100, 100, p)
        g = SimpleWeightedDiGraph(w)
        p1 = path_from_edges(Graphs.a_star(g, 1, nv(g)))
        p2 = try
            TrafficAssignment.astar(g, 1, nv(g), zeros(nv(g)))[1]
        catch e
            if e isa TrafficAssignment.PathNotFoundError
                Int[]
            else
                rethrow(e)
            end
        end
        @test p1 == reverse(p2)
    end
end
