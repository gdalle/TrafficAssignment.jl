using TestItems

@testitem "Parsing" begin
    import TrafficAssignment as TA

    problem = TrafficAssignmentProblem("TransportationNetworks", "Anaheim")
    @test TA.nb_nodes(problem) == 416
    @test TA.nb_links(problem) == 914
    @test TA.nb_zones(problem) == 38
    @test startswith(string(problem), "Traffic")

    problem = TrafficAssignmentProblem("UnifiedTrafficDataset", "San Francisco")
    @test TA.nb_nodes(problem) == 4986
    @test TA.nb_links(problem) == 18002
    @test TA.nb_zones(problem) == 194
    @test startswith(string(problem), "Traffic")

    pb1 = TrafficAssignmentProblem(
        "UnifiedTrafficDataset", "San Francisco"; solution="TransCAD"
    )
    pb2 = TrafficAssignmentProblem(
        "UnifiedTrafficDataset", "San Francisco"; solution="AequilibraE"
    )
end

@testitem "Read all instances" begin
    summary = TrafficAssignment.summarize_instances()
    @testset "$(row[:instance])" for row in eachrow(summary)
        if row[:instance] == "Munich"
            @test_broken row[:valid]
        else
            @test row[:valid]
        end
    end
end
