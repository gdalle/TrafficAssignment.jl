"""
    plot_network(
        problem::TrafficAssignmentProblem, flow=nothing;
        nodes=false, zones=false, tiles=false
    )

Plot a transportation network, possibly on top of real-world map tiles.

If a `flow` is provided, network edges will be colored according to their congestion level.

!!! warning

    This function requires loading one of Makie.jl's backends beforehand, ideally GLMakie.jl.
    Using `tiles=true` requires loading Tyler.jl in addition.
"""
function plot_network(args...; kwargs...)
    return error("Please load a Makie backend")
end

function add_tiles!(args...; kwargs...)
    return error("Please load Tyler")
end
