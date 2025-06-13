module TrafficAssignmentMakieExt

using Colors
using Makie
using Proj
using LinearAlgebra
using SparseArrays
using TrafficAssignment
import TrafficAssignment as TA

const WebMercator = "EPSG:3857"

zero_to_inf(x) = ifelse(iszero(x), Inf, x)

function aggregate_symmetric(A::SparseMatrixCSC, I, J, agg::F=max) where {F}
    nzvals = getindex.(Ref(A), I, J)
    nzvals_transpose = getindex.(Ref(A), J, I)
    return agg.(nzvals, nzvals_transpose)
end

function aggregate_symmetric_normalized(
    A::SparseMatrixCSC, B::SparseMatrixCSC, I, J, agg::F=max
) where {F}
    A_nzvals = getindex.(Ref(A), I, J)
    B_nzvals = zero_to_inf.(getindex.(Ref(B), I, J))
    A_nzvals_transpose = getindex.(Ref(A), J, I)
    B_nzvals_transpose = zero_to_inf.(getindex.(Ref(B), J, I))
    return agg.(A_nzvals ./ B_nzvals, A_nzvals_transpose ./ B_nzvals_transpose)
end

function TrafficAssignment.plot_network(
    problem::TrafficAssignmentProblem,
    edge_quantity::Union{Nothing,SparseMatrixCSC}=nothing;
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
    I, J, _ = findnz(tril(link_capacity + transpose(link_capacity)))

    IJ_real = [(i, j) for (i, j) in zip(I, J) if i in real_nodes && j in real_nodes]
    I_real, J_real = first.(IJ_real), last.(IJ_real)

    IJ_zone = [(i, j) for (i, j) in zip(I, J) if i in zone_nodes || j in zone_nodes]
    I_zone, J_zone = first.(IJ_zone), last.(IJ_zone)

    # point construction
    real_points = Point2f.(X[real_nodes], Y[real_nodes])
    real_start_points = Point2f.(X[I_real], Y[I_real])
    real_end_points = Point2f.(X[J_real], Y[J_real])
    real_point_couples = collect(zip(real_start_points, real_end_points))

    zone_points = Point2f.(X[zone_nodes], Y[zone_nodes])
    zone_start_points = Point2f.(X[I_zone], Y[I_zone])
    zone_end_points = Point2f.(X[J_zone], Y[J_zone])
    zone_point_couples = collect(zip(zone_start_points, zone_end_points))

    if !isnothing(edge_quantity)
        if edge_quantity_type == :flow
            flow = edge_quantity
            # linewidth ∝ flow
            segment_linewidth_real = aggregate_symmetric(flow, I_real, J_real, max)
            segment_linewidth_real .*= 7 ./ maximum(segment_linewidth_real)
            # color ∝ flow / capacity
            segment_color_real = aggregate_symmetric_normalized(
                flow, link_capacity, I_real, J_real, max
            )
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
        fig[3, 4]; active=!isnothing(edge_quantity), tellwidth=false
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

    ls_color_real = if !isnothing(edge_quantity)
        ls_color_real = linesegments!(
            ax,
            real_point_couples;
            linecap=:round,
            linewidth=segment_linewidth_real,
            color=segment_color_real,
            colorrange=segment_colorrange,
            colormap=segment_colormap,
            visible=lift(identity, show_edge_quantity.active),
        )
        Colorbar(
            fig[1, 5],
            ls_color_real;
            label=if edge_quantity_type == :flow
                "Link congestion (flow / capacity)"
            else
                nothing
            end,
            tellheight=false,
        )
        ls_color_real
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
        !isnothing(edge_quantity) && translate!(ls_color_real, 0, 0, 10)
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
