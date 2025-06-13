"""
    TrafficAssignment

A Julia package for studying traffic assignment problems, using instaces from two datasets:

  - `"TransportationNetworks"`, available at [https://github.com/bstabler/TransportationNetworks](https://github.com/bstabler/TransportationNetworks)
  - `"Unified"`, available at [https://figshare.com/articles/dataset/A_unified_and_validated_traffic_dataset_for_20_U_S_cities/24235696](https://figshare.com/articles/dataset/A_unified_and_validated_traffic_dataset_for_20_U_S_cities/24235696)
"""
module TrafficAssignment

"""
    TrafficAssignmentDataset

Enum type listing the two possible sources of data: [`TransportationNetworks`](@ref) or [`Unified`](@ref).
"""
@enum TrafficAssignmentDataset begin
    TransportationNetworks
    Unified
end

@doc """
    TransportationNetworks

Dataset available at [https://github.com/bstabler/TransportationNetworks](https://github.com/bstabler/TransportationNetworks).
""" TransportationNetworks

@doc """
    Unified

Dataset available at [https://figshare.com/articles/dataset/A_unified_and_validated_traffic_dataset_for_20_U_S_cities/24235696](https://figshare.com/articles/dataset/A_unified_and_validated_traffic_dataset_for_20_U_S_cities/24235696).
""" Unified

"""
    TrafficAssignmentSoftware

Enum type listing the two possible software tools used to create solutions for instances from the [`Unified`](@ref) dataset: [`AequilibraE`](@ref) or [`TransCAD`](@ref).
"""
@enum TrafficAssignmentSoftware begin
    AequilibraE
    TransCAD
end

@doc """
    AequilibraE

Solutions computed with the Python package [AequilibraE](https://aequilibrae.com).
""" AequilibraE

@doc """
    TransCAD

Solutions computed with the commercial software [TransCAD](https://www.caliper.com/transcad/default.htm).
""" TransCAD

# outside packages
using BinDeps: unpack_cmd
using Colors
using CSV: CSV
using DataDeps: DataDeps, DataDep, @datadep_str
using DataFrames
using DataFramesMeta
using DataStructures
using DocStringExtensions
using FrankWolfe
using OhMyThreads
using OrderedCollections
using Proj
# standard libraries
using Distributed
using Printf
using LinearAlgebra
using LinearAlgebra: checksquare
using SparseArrays
using Statistics

include("download.jl")
include("types.jl")
include("graph.jl")
include("paths.jl")
include("read.jl")
include("solve.jl")
include("plot.jl")

export TrafficAssignmentDataset, TransportationNetworks, Unified
export TrafficAssignmentSoftware, AequilibraE, TransCAD
export TrafficAssignmentProblem
export nb_nodes, nb_links, nb_zones
export list_instances, summarize_instances
export solve_frank_wolfe, social_cost, price_of_anarchy
export plot_network

end # module
