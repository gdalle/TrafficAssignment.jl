using TrafficAssignment
import TrafficAssignment as TA
# using GLMakie, Tyler
using SimpleWeightedGraphs
using Graphs

# Debug cities

pb = TrafficAssignmentProblem("UnifiedTrafficDataset", "New York"; solution="TransCAD")
graph = SimpleWeightedDiGraph(pb.link_free_flow_time)
components = sort!(strongly_connected_components(graph); by=length, rev=true)
map(length, components)

s = 27986
t = 27521
components
s in components[1]
t in components[2]
dijkstra_shortest_paths(graph, s).parents[t]
s_original = s - minimum(pb.zone_nodes) + 10^7
t_original = t - minimum(pb.zone_nodes) + 10^7
s_original, t_original

dataset = "UnifiedTrafficDataset"
for (_, city) in list_instances(dataset)
    pb = TrafficAssignmentProblem(dataset, city)
    graph = SimpleWeightedDiGraph(pb.link_free_flow_time)
    components = strongly_connected_components(graph)
    if length(components) > 1
        @info "$city"
    end
end
