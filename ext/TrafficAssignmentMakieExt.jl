module TrafficAssignmentMakieExt

using Colors
using Makie
using Proj
using LinearAlgebra
using SparseArrays
using TrafficAssignment
import TrafficAssignment as TA

const WebMercator = "EPSG:3857"

function TrafficAssignment.plot_network(
    problem::TrafficAssignmentProblem,
    edge_quantity::Union{Nothing,SparseMatrixCSC}=nothing;
    offset_ratio=0,
    edge_quantity_type=:flow,
    nodes=false,
    zones=false,
    tiles=false,
)
    (;
        instance_name,
        real_nodes,
        zone_nodes,
        node_coord,
        valid_longitude_latitude,
        link_capacity,
    ) = problem

    # coordinate conversion
    if valid_longitude_latitude
        trans = Proj.Transformation("WGS84", WebMercator; always_xy=true)
        XY = trans.(node_coord)
        X, Y = first.(XY), last.(XY)
    else
        X, Y = first.(node_coord), last.(node_coord)
    end

    # link filtering
    I, J, _ = findnz(link_capacity)

    IJ_real = [(i, j) for (i, j) in zip(I, J) if i in real_nodes && j in real_nodes]
    I_real, J_real = first.(IJ_real), last.(IJ_real)

    IJ_zone = [(i, j) for (i, j) in zip(I, J) if i in zone_nodes || j in zone_nodes]
    I_zone, J_zone = first.(IJ_zone), last.(IJ_zone)

    # point construction
    real_points = Point2f.(X[real_nodes], Y[real_nodes])
    real_start_points = Point2f.(X[I_real], Y[I_real])
    real_end_points = Point2f.(X[J_real], Y[J_real])

    X_q1 = (3 .* X[I_real] .+ X[J_real]) ./ 4
    Y_q1 = (3 .* Y[I_real] .+ Y[J_real]) ./ 4
    X_q2 = (X[I_real] .+ X[J_real]) ./ 2
    Y_q2 = (Y[I_real] .+ Y[J_real]) ./ 2
    X_q3 = (X[I_real] .+ 3 .* X[J_real]) ./ 4
    Y_q3 = (Y[I_real] .+ 3 .* Y[J_real]) ./ 4

    X_diff = (X[J_real] .- X[I_real]) ./ 2
    Y_diff = (Y[J_real] .- Y[I_real]) ./ 2

    X_q1 .+= Y_diff .* offset_ratio ./ 1.2
    X_q2 .+= Y_diff .* offset_ratio
    X_q3 .+= Y_diff .* offset_ratio ./ 1.2
    Y_q1 .-= X_diff .* offset_ratio ./ 1.2
    Y_q2 .-= X_diff .* offset_ratio
    Y_q3 .-= X_diff .* offset_ratio ./ 1.2

    real_q1_points = Point2f.(X_q1, Y_q1)
    real_q2_points = Point2f.(X_q2, Y_q2)
    real_q3_points = Point2f.(X_q3, Y_q3)
    real_point_couples = collect(zip(real_start_points, real_end_points))
    real_point_couples_q0_q1 = collect(zip(real_start_points, real_q1_points))
    real_point_couples_q1_q2 = collect(zip(real_q1_points, real_q2_points))
    real_point_couples_q2_q3 = collect(zip(real_q2_points, real_q3_points))
    real_point_couples_q3_q4 = collect(zip(real_q3_points, real_end_points))

    zone_points = Point2f.(X[zone_nodes], Y[zone_nodes])
    zone_start_points = Point2f.(X[I_zone], Y[I_zone])
    zone_end_points = Point2f.(X[J_zone], Y[J_zone])
    zone_point_couples = collect(zip(zone_start_points, zone_end_points))

    if !isnothing(edge_quantity)
        if edge_quantity_type == :flow
            flow = getindex.(Ref(edge_quantity), I_real, J_real)
            capacity = getindex.(Ref(link_capacity), I_real, J_real)
            # linewidth ∝ flow
            segment_linewidth_real = 2
            # TODO: toggle once error is fixed: "Cairo doesn't support two different line widths (1.3642795 and 2.4439764 at the endpoints of a line)."
            # segment_linewidth_real = 7 .* flow ./ maximum(flow)
            # color ∝ flow / capacity
            segment_color_real = flow ./ capacity
            segment_colormap = vcat(
                range(HSL(colorant"green"); stop=HSL(colorant"yellow"), length=4),
                range(HSL(colorant"yellow"); stop=HSL(colorant"red"), length=4)[2:end],
            )
            segment_colorrange = (0, 2)
        else
            throw(ArgumentError("Edge quantity of type $edge_quantity_type not supported"))
        end
    end

    # figure size
    Δx = maximum(X) - minimum(X)
    Δy = maximum(Y) - minimum(Y)

    # figure creation
    fig = Figure()
    ax = Axis(
        fig[1, 1:4];
        title=instance_name,
        subtitle="$(nb_nodes(problem)) nodes, $(nb_links(problem)) links",
        aspect=tiles ? nothing : DataAspect(),
    )
    hidedecorations!(ax)
    hidespines!(ax)
    Label(fig[2, 1], "Show nodes:"; tellwidth=false)
    Label(fig[2, 3], "Show zones:"; tellwidth=false)
    Label(fig[3, 1], "Show network:"; tellwidth=false)
    Label(fig[3, 3], "Show $edge_quantity_type:"; tellwidth=false)
    show_nodes = Toggle(fig[2, 2]; active=nodes, tellwidth=false)
    show_zones = Toggle(fig[2, 4]; active=zones, tellwidth=false)
    show_network = Toggle(fig[3, 2]; active=true, tellwidth=false)
    show_edge_quantity = Toggle(
        fig[3, 4]; active=(!isnothing(edge_quantity)), tellwidth=false
    )

    # actual plotting
    ls_real = linesegments!(
        ax,
        real_point_couples;
        linewidth=isnothing(edge_quantity) ? 2 : 0.5,
        color=:black,
        linecap=:round,
        visible=lift(identity, show_network.active),
    )

    ls_zone = linesegments!(
        ax,
        zone_point_couples;
        linewidth=isnothing(edge_quantity) ? 2 : 0.5,
        color=:gray,
        linecap=:round,
        visible=lift(&, show_zones.active, show_network.active),
    )

    ls_to_shift = []
    if !isnothing(edge_quantity)
        for couples in (
            real_point_couples_q0_q1, #
            real_point_couples_q1_q2, #
            real_point_couples_q2_q3, #
            real_point_couples_q3_q4, #
        )
            ls = linesegments!(
                ax,
                couples;
                linecap=:round,
                linewidth=segment_linewidth_real,
                color=segment_color_real,
                colorrange=segment_colorrange,
                colormap=segment_colormap,
                visible=lift(identity, show_edge_quantity.active),
            )
            push!(ls_to_shift, ls)
        end
        Colorbar(
            fig[1, 5],
            ls_to_shift[1];
            label=if edge_quantity_type == :flow
                "Link congestion (flow / capacity)"
            else
                nothing
            end,
            tellheight=false,
        )
    else
        nothing
    end

    sc_real = scatter!(
        ax, real_points; color=:black, visible=lift(identity, show_nodes.active)
    )
    sc_zone = scatter!(
        ax,
        zone_points;
        color=:black,
        marker=:diamond,
        visible=lift(identity, show_zones.active),
    )

    if tiles && valid_longitude_latitude
        TA.add_tiles!(fig, ax, node_coord)
        translate!(ls_real, 0, 0, 10)
        translate!(ls_zone, 0, 0, 10)
        translate!(sc_real, 0, 0, 10)
        translate!(sc_zone, 0, 0, 10)
        if !isnothing(edge_quantity)
            for ls in to_shift
                translate!(ls, 0, 0, 10)
            end
        end
    end

    if !tiles
        colsize!(fig.layout, 1, Aspect(1, 0.25 * Δx / Δy))
        colsize!(fig.layout, 2, Aspect(1, 0.25 * Δx / Δy))
        colsize!(fig.layout, 3, Aspect(1, 0.25 * Δx / Δy))
        colsize!(fig.layout, 4, Aspect(1, 0.25 * Δx / Δy))
    end
    resize_to_layout!(fig)

    return fig
end

end
