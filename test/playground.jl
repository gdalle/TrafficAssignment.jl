using TrafficAssignment
import TrafficAssignment as TA
using SimpleWeightedGraphs, Graphs

dataset = "UnifiedTrafficDataset"
city = "New York"
pb = TrafficAssignmentProblem(dataset, city)

@time solve_frank_wolfe(pb; max_iteration=1, verbose=false);
@profview solve_frank_wolfe(pb; max_iteration=5, verbose=false);

using JET, Cthulhu

g = SimpleWeightedDiGraph(pb.link_free_flow_time)

@test_opt TrafficAssignment.dijkstra(g, 1)
@profview solve_frank_wolfe(pb; max_iteration=1, verbose=false);
@test_opt solve_frank_wolfe(pb; max_iteration=1, verbose=false);
