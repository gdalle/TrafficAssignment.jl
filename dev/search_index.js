var documenterSearchIndex = {"docs":
[{"location":"api/#API-reference","page":"API reference","title":"API reference","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"CollapsedDocStrings = true\nCurrentModule = TrafficAssignment","category":"page"},{"location":"api/#Public","page":"API reference","title":"Public","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [TrafficAssignment]\nPrivate = false","category":"page"},{"location":"api/#Private","page":"API reference","title":"Private","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [TrafficAssignment]\nPublic = false","category":"page"},{"location":"#TrafficAssignment.jl","page":"Home","title":"TrafficAssignment.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Build Status) (Image: Coverage)","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package for the Julia Language does basically two tasks: (1) loading a network data and (2) finding a user equilibrium traffic pattern. See Traffic Assignment.","category":"page"},{"location":"#Install","page":"Home","title":"Install","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"julia> Pkg.add(\"TrafficAssignment\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"This will install Graphs.jl and Optim.jl, if you don't have them already.","category":"page"},{"location":"","page":"Home","title":"Home","text":"To check if works","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> Pkg.test(\"TrafficAssignment\")","category":"page"},{"location":"#load_ta_network","page":"Home","title":"load_ta_network","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This function loads a network data available in this TNTP github repository. The network name must match with the directory name in the TNTP repository.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Example:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using TrafficAssignment\nta_data = load_ta_network(\"SiouxFalls\")\n# ta_data = load_ta_network(\"Anaheim\")\n# ta_data = load_ta_network(\"Barcelona\")\n# ta_data = load_ta_network(\"Winnipeg\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"The return value is of the TA_Data type, which is defined as","category":"page"},{"location":"","page":"Home","title":"Home","text":"mutable struct TA_Data\n    network_name::String\n\n    number_of_zones::Int\n    number_of_nodes::Int\n    first_thru_node::Int\n    number_of_links::Int\n\n    init_node::Array{Int,1}\n    term_node::Array{Int,1}\n    capacity::Array{Float64,1}\n    link_length::Array{Float64,1}\n    free_flow_time::Array{Float64,1}\n    b::Array{Float64,1}\n    power::Array{Float64,1}\n    speed_limit::Array{Float64,1}\n    toll::Array{Float64,1}\n    link_type::Array{Int64,1}\n\n    total_od_flow::Float64\n\n    travel_demand::Array{Float64,2}\n    od_pairs::Array{Tuple{Int64,Int64},1}\n\n    toll_factor::Float64\n    distance_factor::Float64\n\n    best_objective::Float64\nend","category":"page"},{"location":"#ta_frank_wolfe","page":"Home","title":"ta_frank_wolfe","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This function implements methods to find traffic equilibrium flows: currently, Frank-Wolfe (FW) method, Conjugate FW (CFW) method, and Bi-conjugate FW (BFW) method.","category":"page"},{"location":"","page":"Home","title":"Home","text":"References:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Mitradjieva, M., & Lindberg, P. O. (2013). The Stiff Is Moving-Conjugate Direction Frank-Wolfe Methods with Applications to Traffic Assignment. Transportation Science, 47(2), 280-293.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Example:","category":"page"},{"location":"","page":"Home","title":"Home","text":"link_flow, link_travel_time, objective = ta_frank_wolfe(ta_data, log=\"off\", tol=1e-2)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Available optional arguments:","category":"page"},{"location":"","page":"Home","title":"Home","text":"method=:fw / :cfw / :bfw (default=:bfw)\nstep=\"exact\" / \"newton\" : exact line search using golden section / a simple Newton-type step (default=:exact)\nlog=:on / :off : displays information from each iteration or not (default=:off)\nmaxiterno=integer value : maximum number of iterations (default=2000)\ntol=numeric value : tolerance for the Average Excess Cost (AEC) (default=1e-3)","category":"page"},{"location":"","page":"Home","title":"Home","text":"For example, one may do:","category":"page"},{"location":"","page":"Home","title":"Home","text":"ta_data = load_ta_network(\"SiouxFalls\")\nlink_flow, link_travel_time, objective =\nta_frank_wolfe(ta_data, method=:cfw, max_iter_no=50000, step=:newton, log=:on, tol=1e-5)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The total system travel time can be simply computed as","category":"page"},{"location":"","page":"Home","title":"Home","text":"using LinearAlgebra\nsystem_travel_time = dot(link_travel_time, link_flow)","category":"page"},{"location":"#Contributor","page":"Home","title":"Contributor","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package is written and maintained by Changhyun Kwon.","category":"page"}]
}
