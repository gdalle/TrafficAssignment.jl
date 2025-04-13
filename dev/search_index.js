var documenterSearchIndex = {"docs":
[{"location":"api/#API-reference","page":"API reference","title":"API reference","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"CollapsedDocStrings = true\nCurrentModule = TrafficAssignment","category":"page"},{"location":"api/#Public","page":"API reference","title":"Public","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [TrafficAssignment]\nPrivate = false","category":"page"},{"location":"api/#TrafficAssignment.TrafficAssignment","page":"API reference","title":"TrafficAssignment.TrafficAssignment","text":"TrafficAssignment\n\nA Julia package for studying traffic assignment problems, using data from https://github.com/bstabler/TransportationNetworks.\n\n\n\n\n\n","category":"module"},{"location":"api/#TrafficAssignment.TrafficAssignmentProblem","page":"API reference","title":"TrafficAssignment.TrafficAssignmentProblem","text":"TrafficAssignmentProblem(instance_name; ...)\nTrafficAssignmentProblem(\n    instance_name,\n    files;\n    toll_factor,\n    distance_factor\n)\n\n\nUser-friendly constructor for TrafficAssignmentProblem.\n\nThe provided instance_name must be one of the subfolders in https://github.com/bstabler/TransportationNetworks.\n\nWhen you run this function for the first time, the DataDeps package will ask you to confirm download. If you want to skip this check, for instance during CI, set the environment variable ENV[\"DATADEPS_ALWAYS_ACCEPT\"] = true.\n\n\n\n\n\n","category":"type"},{"location":"api/#TrafficAssignment.TrafficAssignmentProblem-2","page":"API reference","title":"TrafficAssignment.TrafficAssignmentProblem","text":"struct TrafficAssignmentProblem{C<:Union{Nothing, Vector{Float64}}, F<:Union{Nothing, SparseArrays.SparseMatrixCSC{Float64, Int64}}}\n\nInstance of the static traffic assignment problem.\n\nDetails\n\nThe link travel time is given by travel_time = free_flow_time * ( 1 + b * (flow/capacity)^power).\n\nThe generalized cost is cost = travel_time + toll_factor * toll + distance_factor * distance.\n\nFields\n\ninstance_name::String\nnumber_of_zones::Int64\nnumber_of_nodes::Int64\nfirst_thru_node::Int64\nnumber_of_links::Int64\ninit_node::Vector{Int64}\nterm_node::Vector{Int64}\ncapacity::SparseArrays.SparseMatrixCSC{Float64, Int64}\nlink_length::SparseArrays.SparseMatrixCSC{Float64, Int64}\nfree_flow_time::SparseArrays.SparseMatrixCSC{Float64, Int64}\nb::SparseArrays.SparseMatrixCSC{Float64, Int64}\npower::SparseArrays.SparseMatrixCSC{Float64, Int64}\nspeed_limit::SparseArrays.SparseMatrixCSC{Float64, Int64}\ntoll::SparseArrays.SparseMatrixCSC{Float64, Int64}\nlink_type::SparseArrays.SparseMatrixCSC{Int64, Int64}\ntotal_od_flow::Float64\ntravel_demand::Matrix{Float64}\nod_pairs::Vector{Tuple{Int64, Int64}}\nX::Union{Nothing, Vector{Float64}}\nY::Union{Nothing, Vector{Float64}}\noptimal_flow_volume::Union{Nothing, SparseArrays.SparseMatrixCSC{Float64, Int64}}\noptimal_flow_cost::Union{Nothing, SparseArrays.SparseMatrixCSC{Float64, Int64}}\ntoll_factor::Float64\ndistance_factor::Float64\n\n\n\n\n\n","category":"type"},{"location":"api/#TrafficAssignment.list_instances-Tuple{}","page":"API reference","title":"TrafficAssignment.list_instances","text":"list_instances()\n\n\nReturn a list of available instance names.\n\n\n\n\n\n","category":"method"},{"location":"api/#TrafficAssignment.plot_network-Tuple{TrafficAssignmentProblem}","page":"API reference","title":"TrafficAssignment.plot_network","text":"plot_network(problem)\n\n\nPlot a transportation network with Makie.\n\nThis function requires loading one of Makie's backends beforehand.\n\n\n\n\n\n","category":"method"},{"location":"api/#TrafficAssignment.social_cost-Tuple{TrafficAssignmentProblem, AbstractMatrix}","page":"API reference","title":"TrafficAssignment.social_cost","text":"social_cost(problem, flow)\n\n\nCompute the social cost induced by a matrix of link flows.\n\n\n\n\n\n","category":"method"},{"location":"api/#TrafficAssignment.solve_frank_wolfe-Union{Tuple{TrafficAssignmentProblem}, Tuple{A}, Tuple{TrafficAssignmentProblem, A}} where A","page":"API reference","title":"TrafficAssignment.solve_frank_wolfe","text":"solve_frank_wolfe(problem; ...)\nsolve_frank_wolfe(\n    problem,\n    frank_wolfe_alg;\n    verbose,\n    kwargs...\n)\n\n\nSolve a traffic assignment problem using an algorithm from the FrankWolfe library.\n\nKeyword arguments are passed to frank_wolfe_alg.\n\n\n\n\n\n","category":"method"},{"location":"api/#TrafficAssignment.summarize_instances-Tuple{}","page":"API reference","title":"TrafficAssignment.summarize_instances","text":"summarize_instances()\n\n\nReturn a DataFrame summarizing the dimensions of all available instances.\n\n\n\n\n\n","category":"method"},{"location":"api/#Private","page":"API reference","title":"Private","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [TrafficAssignment]\nPublic = false","category":"page"},{"location":"api/#TrafficAssignment.instance_files-Tuple{AbstractString}","page":"API reference","title":"TrafficAssignment.instance_files","text":"instance_files(instance_name)\n\n\nReturn a named tuple (; flow_file, net_file, node_file, trips_file) containing the absolute paths to the 4 data tables of an instance.\n\n\n\n\n\n","category":"method"},{"location":"#TrafficAssignment.jl","page":"Home","title":"TrafficAssignment.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Build Status) (Image: Coverage) (Image: Dev Documentation)","category":"page"},{"location":"","page":"Home","title":"Home","text":"This is a Julia package for studying traffic assignment on road networks.","category":"page"},{"location":"#Getting-started","page":"Home","title":"Getting started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install the latest development version, run this in your Julia Pkg REPL:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add https://github.com/gdalle/TrafficAssignment.jl","category":"page"},{"location":"","page":"Home","title":"Home","text":"You can easily load networks from the TransportationNetworks repository:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using TrafficAssignment\n\njulia> problem = TrafficAssignmentProblem(\"SiouxFalls\")\nTraffic assignment problem on the SiouxFalls network with 24 nodes and 76 links","category":"page"},{"location":"","page":"Home","title":"Home","text":"And then you can solve the equilibrium problem and compute the total system travel time:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> flow = solve_frank_wolfe(problem; max_iteration=1000, verbose=false)\n24×24 SparseArrays.SparseMatrixCSC{Float64, Int64} with 76 stored entries:\n⎡⠎⡡⡐⠀⠀⡠⠀⠀⠀⠀⠀⠀⎤\n⎢⠐⠈⢊⡰⡁⠀⠀⢀⠠⠀⠀⠀⎥\n⎢⠀⡠⠁⠈⠪⡢⡠⠒⠂⠀⠀⠀⎥\n⎢⠀⠀⠀⢀⢠⠊⠠⠂⣀⠄⠠⠊⎥\n⎢⠀⠀⠀⠂⠈⠀⠀⠜⢄⡱⣀⠀⎥\n⎣⠀⠀⠀⠀⠀⠀⡠⠂⠀⠘⡪⡪⎦\n\njulia> round(social_cost(problem, flow), sigdigits=4)\n7.481e6","category":"page"},{"location":"#Credits","page":"Home","title":"Credits","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package was originally written and maintained by Changhyun Kwon.","category":"page"}]
}
