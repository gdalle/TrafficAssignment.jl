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

function TrafficAssignment.plot_network(
    problem::TrafficAssignmentProblem,
    flow::Union{Nothing,SparseMatrixCSC}=nothing;
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
    if !zones
        IJ_filtered = [(i, j) for (i, j) in zip(I, J) if i in real_nodes && j in real_nodes]
        I = first.(IJ_filtered)
        J = last.(IJ_filtered)
    end

    # point construction
    real_points = Point2f.(X[real_nodes], Y[real_nodes])
    zone_points = Point2f.(X[zone_nodes], Y[zone_nodes])
    start_points = Point2f.(X[I], Y[I])
    end_points = Point2f.(X[J], Y[J])
    point_couples = collect(zip(start_points, end_points))
    segment_linewidth = if isnothing(flow)
        1
    else
        segment_flow = getindex.(Ref(flow), I, J)
        segment_flow_reverse = getindex.(Ref(flow), J, I)
        segment_linewidth = max.(segment_flow, segment_flow_reverse)
        7 .* segment_linewidth ./ maximum(segment_linewidth)
    end
    segment_color = if isnothing(flow)
        :black
    else
        segment_congestion = (
            getindex.(Ref(flow), I, J) ./ zero_to_inf.(getindex.(Ref(link_capacity), I, J))
        )
        segment_congestion_reverse = (
            getindex.(Ref(flow), J, I) ./ zero_to_inf.(getindex.(Ref(link_capacity), J, I))
        )
        max.(segment_congestion, segment_congestion_reverse)
    end

    # figure size
    Δx = maximum(X) - minimum(X)
    Δy = maximum(Y) - minimum(Y)

    # figure creation
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title=instance_name,
        subtitle="$(nb_nodes(problem)) nodes, $(nb_links(problem)) links",
        aspect=tiles ? nothing : DataAspect(),
    )
    hidedecorations!(ax)
    hidespines!(ax)

    # actual plotting
    ls1 = linesegments!(
        ax,
        point_couples;
        linewidth=isnothing(flow) ? 2 : 0.5,
        color=isnothing(flow) ? :black : :gray,
        linecap=:round,
    )
    ls2 = if !isnothing(flow)
        red_green_yellow = vcat(
            range(HSL(colorant"green"); stop=HSL(colorant"yellow"), length=4),
            range(HSL(colorant"yellow"); stop=HSL(colorant"red"), length=4)[2:end],
        )
        ls2 = linesegments!(
            ax,
            point_couples;
            linecap=:round,
            linewidth=segment_linewidth,
            color=segment_color,
            colorrange=(0, 2),
            colormap=red_green_yellow,
        )
        Colorbar(fig[2, 1], ls2; vertical=false, label="Link flow (flow / capacity)")
        if !tiles
            colsize!(fig.layout, 1, Aspect(1, Δx / Δy))
        end
        ls2
    end
    sc1 = if nodes
        scatter!(ax, real_points; color=:black)
    end
    sc2 = if zones
        scatter!(ax, zone_points; color=:black, marker=:diamond)
    end
    if tiles && valid_longitude_latitude
        TA.add_tiles!(fig, ax, node_coord)
        translate!(ls1, 0, 0, 10)
        nodes && translate!(sc1, 0, 0, 10)
        zones && translate!(sc2, 0, 0, 10)
        !isnothing(flow) && translate!(ls2, 0, 0, 10)
    end
    resize_to_layout!(fig)
    return fig
end

end
