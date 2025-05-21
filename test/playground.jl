using TrafficAssignment
import TrafficAssignment as TA

dataset = "UnifiedTrafficDataset"
city = "Honolulu"
pb = TrafficAssignmentProblem(dataset, city)

@time solve_frank_wolfe(pb; max_iteration=5, verbose=false);
@profview solve_frank_wolfe(pb; max_iteration=5, verbose=false);
