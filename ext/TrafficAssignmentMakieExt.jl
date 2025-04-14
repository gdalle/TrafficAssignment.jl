module TrafficAssignmentMakieExt

using Makie
using Proj
using SparseArrays
using TrafficAssignment
import TrafficAssignment as TA

const WebMercator = "EPSG:3857"

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
    I, J, _ = findnz(link_capacity)
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
    segment_color = if isnothing(flow)
        :black
    else
        congestion = flow ./ link_capacity
        getindex.(Ref(congestion), I, J)
    end
    segment_colorrange = if isnothing(flow)
        nothing
    else
        congestion_vec = nonzeros(flow) ./ nonzeros(link_capacity)
        (0, maximum(congestion_vec))
    end

    # figure size
    Δx = maximum(X) - minimum(X)
    Δy = maximum(Y) - minimum(Y)
    Δmax = max(Δx, Δy)
    figsize = (700 / Δmax) .* (Δx, Δy)

    # figure creation
    fig = Figure(; size=figsize)
    ax = Axis(
        fig[1, 1];
        title=instance_name,
        subtitle="$(nb_nodes(problem)) nodes, $(nb_links(problem)) links",
        aspect=DataAspect(),
    )
    hidedecorations!(ax)
    hidespines!(ax)

    # actual plotting
    ls = linesegments!(
        ax,
        point_couples;
        linewidth=2,
        color=segment_color,
        colorrange=segment_colorrange,
        colormap=:plasma,
    )
    if nodes
        sc1 = scatter!(ax, real_points; color=:black)
    end
    if zones
        sc2 = scatter!(ax, zone_points; color=:red, marker=:diamond)
    end
    if !isnothing(flow)
        Colorbar(fig[2, 1], ls; vertical=false, label="Link congestion (flow / capacity)")
    end
    if tiles && valid_longitude_latitude
        TA.add_tiles!(fig, ax, node_coord)
        nodes && translate!(sc1, 0, 0, 10)
        zones && translate!(sc2, 0, 0, 10)
        translate!(ls, 0, 0, 10)
    end
    return fig
end

end
