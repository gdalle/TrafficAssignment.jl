"""
$(TYPEDEF)

Instance of the static traffic assignment problem.

# Details

The link travel time is given by the formula of the Bureau of Public Roads (BPR):

    t = t0 * (1 + α * (f/c)^β)

where

  - `t` is the travel time
  - `t0` is the free flow time
  - `f` is the flow along the link
  - `c` is the link capacity
  - `α` is a multiplicative coefficient (often taken to be `0.15`)
  - `β` is an exponent (often taken to be `4`)

# Fields

$(TYPEDFIELDS)
"""
@kwdef struct TrafficAssignmentProblem{
    Coord<:Union{Missing,NTuple{2,Number}},
    Capa<:Number,
    Length<:Number,
    Free<:Number,
    Speed<:Number,
    BPRMult<:Number,
    BPRPow<:Number,
    Toll<:Union{Missing,Number},
    LinkT<:Any,
    Flow<:Union{Missing,Number},
    Dem<:Number,
    TF<:Union{Number,Missing},
    DF<:Union{Number,Missing},
}
    "name of the dataset"
    dataset::TrafficAssignmentDataset
    "name of the instance (subfolder inside the dataset)"
    instance_name::String

    # nodes
    "number of nodes in the network (nodes are numbered from `1` to `nb_nodes`)"
    nb_nodes::Int
    "number of directed links in the network"
    nb_links::Int
    "interval of nodes that correspond to real intersections"
    real_nodes::UnitRange{Int}
    "interval of nodes that correspond to artificial zones"
    zone_nodes::UnitRange{Int}
    "coordinates of the nodes for plotting"
    node_coord::Vector{Coord}
    "whether `node_coord` corresponds to the longitude and latitude"
    valid_longitude_latitude::Bool

    # links
    "matrix of link ids starting at one"
    link_id::SparseMatrixCSC{Int,Int}
    "matrix of link capacities (`c` in the BPR formula)"
    link_capacity::SparseMatrixCSC{Capa,Int}
    "matrix of link lengths"
    link_length::SparseMatrixCSC{Length,Int}
    "matrix of link free flow times (`t0` in the BPR formula)"
    link_free_flow_time::SparseMatrixCSC{Free,Int}
    "matrix of link speed limits"
    link_speed_limit::SparseMatrixCSC{Speed,Int}
    "link multiplicative factors `α` in the BPR formula, either a single scalar or a matrix"
    link_bpr_mult::SparseMatrixCSC{BPRMult,Int}
    "link exponents `β` in the BPR formula, either a single scalar or a matrix"
    link_bpr_power::SparseMatrixCSC{BPRPow,Int}
    "matrix of link tolls"
    link_toll::SparseMatrixCSC{Toll,Int}
    "matrix of link types"
    link_type::SparseMatrixCSC{LinkT,Int}

    # demand
    "demand by OD pair"
    demand::Dict{Tuple{Int,Int},Dem}
    "vector of unique origins for the OD pairs"
    origins::Vector{Int}
    "vector of unique destinations for the OD pairs"
    destinations::Vector{Int}
    "OD pairs removed because no path between them exists"
    removed_od_pairs::Vector{Tuple{Int,Int}}
    "dictionary mapping each destination to a vector of free flow path times for every possible origin"
    destination_free_flow_time::Dict{Int,Vector{Free}}

    # cost parameters
    "conversion factor turning toll costs into temporal costs, expressed in time/toll"
    toll_factor::TF
    "conversion factor turning distance costs into temporal costs, expressed in time/length"
    distance_factor::DF

    # solution
    "software used to compute the provided solution, if known"
    solution_software::Union{TrafficAssignmentSoftware,Missing}
    "provided matrix of optimal link flows"
    optimal_flow::SparseMatrixCSC{Flow,Int}
end

function Base.show(io::IO, problem::TrafficAssignmentProblem)
    (; instance_name) = problem
    return print(
        io,
        "Traffic assignment problem on the $instance_name network with $(nb_nodes(problem)) nodes, $(nb_links(problem)) links and $(nb_zones(problem)) zones",
    )
end

"""
    nb_nodes(problem::TrafficAssignmentProblem)

Return the number of nodes in the network (including zone nodes).
"""
nb_nodes(problem::TrafficAssignmentProblem) = problem.nb_nodes

"""
    nb_links(problem::TrafficAssignmentProblem)

Return the number of links in the network (including links to and from zone nodes).
"""
nb_links(problem::TrafficAssignmentProblem) = problem.nb_links

"""
    nb_zones(problem::TrafficAssignment)

Return the number of fake nodes in the network that represent zones.
"""
nb_zones(problem::TrafficAssignmentProblem) = length(problem.zone_nodes)
