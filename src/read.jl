"""
    instance_files(dataset_name, instance_name)

Return a named tuple containing the absolute paths to the individual data tables of an instance.
"""
function instance_files(dataset_name::AbstractString, instance_name::AbstractString)
    return _instance_files(Val(Symbol(dataset_name)), instance_name)
end

function _instance_files(::Val{:TransportationNetworks}, instance_name)
    instance_dir = datapath("TransportationNetworks", instance_name)
    @assert ispath(instance_dir)

    flow_file = net_file = node_file = trips_file = nothing
    for f in readdir(instance_dir; join=true)
        if occursin("_flow", lowercase(f)) && endswith(lowercase(f), ".tntp")
            flow_file = f
        elseif occursin("_net", lowercase(f)) && endswith(lowercase(f), ".tntp")
            net_file = f
        elseif occursin("_node", lowercase(f)) && endswith(lowercase(f), ".tntp")
            node_file = f
        elseif occursin("_trips", lowercase(f)) && endswith(lowercase(f), ".tntp")
            trips_file = f
        end
    end

    if instance_name == "Munich"
        # https://github.com/bstabler/TransportationNetworks/issues/59
        net_file, trips_file = trips_file, net_file
    end

    @assert !isnothing(net_file)
    @assert !isnothing(trips_file)

    return (; flow_file, net_file, node_file, trips_file)
end

function _instance_files(::Val{:UnifiedTrafficDataset}, instance_name)
    dataset_name = "UnifiedTrafficDataset"
    underscored_instance_name = replace(instance_name, " " => "_")
    instance_dir = ""
    for potential_instance_dir in readdir(datapath(dataset_name); join=true)
        if endswith(potential_instance_dir, underscored_instance_name)
            instance_dir = potential_instance_dir
            break
        end
    end
    @assert ispath(instance_dir)
    node_file = joinpath(
        instance_dir, "01_Input_data", "$(underscored_instance_name)_node.csv"
    )
    link_file = joinpath(
        instance_dir, "01_Input_data", "$(underscored_instance_name)_link.csv"
    )
    od_file = joinpath(instance_dir, "01_Input_data", "$(underscored_instance_name)_od.csv")
    aequilibrae_network_file = joinpath(
        instance_dir, "03_AequilibraE_results", "network.csv"
    )
    aequilibrae_assignment_file = joinpath(
        instance_dir, "03_AequilibraE_results", "assignment_result.csv"
    )
    transcad_linkflows_file = joinpath(
        instance_dir,
        "02_TransCAD_results",
        "$(underscored_instance_name)_link_LinkFlows.csv",
    )
    return (;
        node_file,
        link_file,
        od_file,
        aequilibrae_network_file,
        aequilibrae_assignment_file,
        transcad_linkflows_file,
    )
end

"""
    TrafficAssignmentProblem(dataset_name, instance_name)

User-friendly constructor for [`TrafficAssignmentProblem`](@ref).

The dataset must be one of `$DATASET_NAMES`, the instance can be chosen from [`list_instances`](@ref).

!!! tip

    When you run this function for the first time, the DataDeps package will ask you to confirm download.
    If you want to skip this check, for instance during CI, set the environment variable `ENV["DATADEPS_ALWAYS_ACCEPT"] = true`.
"""
function TrafficAssignmentProblem(
    dataset_name::AbstractString, instance_name::AbstractString; kwargs...
)
    if !(dataset_name in DATASET_NAMES)
        throw(ArgumentError("The dataset name must be one of $DATASET_NAMES"))
    end
    pb = _TrafficAssignmentProblem(Val(Symbol(dataset_name)), instance_name; kwargs...)
    postprocess!(pb)
    return pb
end

function _TrafficAssignmentProblem(
    ::Val{:TransportationNetworks},
    instance_name,
    toll_factor::Real=0.0,
    distance_factor::Real=0.0,
)
    (; net_file, trips_file, node_file, flow_file) = instance_files(
        "TransportationNetworks", instance_name
    )
    @assert ispath(net_file)
    @assert ispath(trips_file)

    # network table

    nb_zones = 0
    nb_links = 0
    nb_nodes = 0
    first_thru_node = 0

    any_nb_of_spaces = r"[ \t]+"

    header_row = 0
    net_lines = readlines(net_file)
    for (k, line) in enumerate(net_lines)
        if startswith(line, "<NUMBER OF ZONES>")
            nb_zones = parse(Int, split(line, any_nb_of_spaces)[4])
        elseif startswith(line, "<NUMBER OF NODES>")
            nb_nodes = parse(Int, split(line, any_nb_of_spaces)[4])
        elseif startswith(line, "<NUMBER OF LINKS>")
            nb_links = parse(Int, split(line, any_nb_of_spaces)[4])
        elseif startswith(line, "<FIRST THRU NODE>")
            first_thru_node = parse(Int, split(line, any_nb_of_spaces)[4])
        elseif startswith(line, "~")
            header_row = k
            break
        end
    end

    zone_nodes = 1:(first_thru_node - 1)
    real_nodes = first_thru_node:nb_nodes

    net_col_names = string.(split(net_lines[header_row], any_nb_of_spaces))
    if first(net_col_names) == "~" &&
        !startswith(net_lines[header_row + 1], any_nb_of_spaces)
        # first column is not empty
        deleteat!(net_col_names, 1)
    end
    net_df = DataFrame(
        CSV.File(
            net_file;
            delim='\t',
            skipto=header_row + 1,
            header=net_col_names,
            maxwarnings=1,
            silencewarnings=true,
            drop=(i, name) -> i > length(net_col_names),
        ),
    )
    @assert size(net_df, 1) == nb_links

    n, m = nb_nodes, nb_links
    I = net_df[!, :init_node]
    J = net_df[!, :term_node]

    link_id = sparse(I, J, ones(Int, m), n, n)
    nonzeros(link_id) .= 1:nnz(link_id)  # doubles
    link_capacity = sparse(I, J, float.(net_df[!, :capacity]), n, n)
    link_length = sparse(I, J, float.(net_df[!, :length]), n, n)
    link_free_flow_time = sparse(I, J, float.(net_df[!, :free_flow_time]), n, n)
    link_speed_limit = sparse(I, J, float.(net_df[!, :speed]), n, n)
    if ==(extrema(net_df[!, :b])...)
        link_bpr_mult = float(first(net_df[!, :b]))  # single b value
    else
        link_bpr_mult = sparse(I, J, float.(net_df[!, :b]), n, n)
    end
    if ==(extrema(net_df[!, :power])...)
        link_bpr_power = float(first(net_df[!, :power]))  # single power value
    else
        link_bpr_power = sparse(I, J, float.(net_df[!, :power]), n, n)
    end
    if "toll" in names(net_df)
        link_toll = sparse(I, J, float.(net_df[!, :toll]), n, n)
    else
        link_toll = missing
    end
    if "link_type" in names(net_df)
        link_type_nzval = if eltype(net_df[!, :link_type]) <: AbstractString
            # in some instances, the semicolon is stuck at the end
            parse.(Int, strip.(net_df[!, :link_type], ';'))
        else
            net_df[!, :link_type]
        end
        link_type = sparse(I, J, link_type_nzval, n, n)
    else
        link_type = missing
    end

    # trips table
    nb_zones_trip = 0
    total_od_flow = 0

    demand = Dict{Tuple{Int,Int},Float64}()

    trips_lines = readlines(trips_file)
    origin = -1
    for line in trips_lines
        if startswith(line, "<NUMBER OF ZONES>")
            nb_zones_trip = parse(Int, split(line, any_nb_of_spaces)[4])
        elseif startswith(line, "<TOTAL OD FLOW>")
            total_od_flow = parse(Float64, split(line, any_nb_of_spaces)[4])
        elseif line == ""
            origin = -1
        elseif occursin("Origin", line)
            origin = parse(Int, split(line)[2])
        elseif occursin(";", line)
            pairs = split(line, ";")
            for i in 1:size(pairs)[1]
                if occursin(":", pairs[i])
                    pair = split(pairs[i], ":")
                    destination = parse(Int, strip(pair[1]))
                    od_flow = parse(Float64, strip(pair[2]))
                    demand[origin, destination] = od_flow
                end
            end
        end
    end
    origins_destinations = collect(keys(demand))
    origins = unique(map(first, origins_destinations))
    destinations = unique(map(last, origins_destinations))
    removed_od_pairs = Tuple{Int,Int}[]

    @assert nb_zones_trip == nb_zones # Check if nb_zone is same in both txt files
    @assert total_od_flow > 0

    destination_free_flow_time = Dict{Int,Vector{Float64}}()

    # node table

    if !isnothing(node_file)
        coord_lines = readlines(node_file)
        if startswith(lowercase(coord_lines[1]), "node")
            coord_lines = @view(coord_lines[2:end])
        end
        coord_lines_split = split.(coord_lines, Ref(r"[\t ]+"))
        inds = parse.(Int, getindex.(coord_lines_split, 1))
        @assert inds == 1:nb_nodes
        x = parse.(Float64, getindex.(coord_lines_split, 2))
        y = parse.(Float64, getindex.(coord_lines_split, 3))
        source_crs = if occursin("Birmingham", instance_name)
            "EPSG:27700"  # right
        elseif occursin("chicago", lowercase(instance_name))
            "EPSG:26771"  # slightly off
        elseif instance_name in ("GoldCoast", "SiouxFalls", "Sydney")
            "WGS84"
        else
            nothing
        end
        if source_crs !== nothing
            trans = Proj.Transformation(source_crs, "WGS84"; always_xy=true)
            node_coord = trans.(collect(zip(x, y)))
            valid_longitude_latitude = true
        else
            node_coord = collect(zip(x, y))
            valid_longitude_latitude = false
        end
    else
        node_coord = missing
        valid_longitude_latitude = false
    end

    if !isnothing(flow_file) && instance_name != "chicago-regional"
        optimal_flow_df = DataFrame(CSV.File(flow_file))
        I, J = optimal_flow_df[!, "From "], optimal_flow_df[!, "To "]
        optimal_flow = sparse(I, J, optimal_flow_df[!, "Volume "], n, n)
        optimal_flow_cost = sparse(I, J, optimal_flow_df[!, "Cost "], n, n)
    else
        optimal_flow = missing
        optimal_flow_cost = missing
    end

    return TrafficAssignmentProblem(;
        dataset_name="TransportationNetworks",
        instance_name,
        # nodes
        nb_nodes,
        nb_links,
        real_nodes,
        zone_nodes,
        node_coord,
        valid_longitude_latitude,
        # links
        link_id,
        link_capacity,
        link_length,
        link_free_flow_time,
        link_speed_limit,
        link_bpr_mult,
        link_bpr_power,
        link_toll,
        link_type,
        # demand
        demand,
        origins,
        destinations,
        removed_od_pairs,
        destination_free_flow_time,
        # cost
        toll_factor,
        distance_factor,
        # solution
        optimal_flow,
    )
end

function _TrafficAssignmentProblem(
    ::Val{:UnifiedTrafficDataset}, instance_name; solution=nothing
)
    dataset_name = "UnifiedTrafficDataset"
    underscored_instance_name = replace(instance_name, " " => "_")
    files = TrafficAssignment.instance_files(dataset_name, instance_name)

    # nodes

    node_df = DataFrame(CSV.File(files.node_file::String))
    node_df[!, :New_Node_ID] = 1:size(node_df, 1)

    nb_nodes = size(node_df, 1)
    node_coord = collect(zip(node_df[!, :Lon], node_df[!, :Lat]))

    min_zone_node, max_zone_node = extrema(
        @rsubset(node_df, :Tract_Node == true)[!, :New_Node_ID]
    )
    zone_nodes = min_zone_node:max_zone_node
    real_nodes = 1:(min_zone_node - 1)

    # links

    link_df = DataFrame(CSV.File(files.link_file))
    link_df[!, :Free_Flow_Time] =
        (60 / 1000) .* link_df[!, :Length] ./ link_df[!, :Free_Speed]
    link_df = leftjoin(
        link_df,
        @select(node_df, :Node_ID, :New_Node_ID);
        on=:From_Node_ID => :Node_ID,
        renamecols="" => "_From",
    )
    link_df = leftjoin(
        link_df,
        @select(node_df, :Node_ID, :New_Node_ID);
        on=:To_Node_ID => :Node_ID,
        renamecols="" => "_To",
    )

    nb_links = size(link_df, 1)

    n, m = nb_nodes, nb_links
    I, J = link_df[!, :New_Node_ID_From], link_df[!, :New_Node_ID_To]
    link_id = sparse(I, J, ones(Int, m), n, n)
    nonzeros(link_id) .= 1:m
    link_capacity = sparse(I, J, link_df[!, :Capacity], n, n)
    link_length = sparse(I, J, link_df[!, :Length], n, n)
    link_free_flow_time = sparse(I, J, link_df[!, :Free_Flow_Time], n, n)
    link_speed_limit = sparse(I, J, link_df[!, :Free_Speed], n, n)
    link_type = sparse(I, J, link_df[!, :Link_Type], n, n)

    bpr_df = DataFrame(
        CSV.File(joinpath(dirname(@__DIR__), "data", "UnifiedTrafficDataset", "bpr.csv"))
    )

    link_bpr_mult = only(
        @rsubset(bpr_df, :city == underscored_instance_name)[!, :bpr_alpha]
    )
    link_bpr_power = only(
        @rsubset(bpr_df, :city == underscored_instance_name)[!, :bpr_beta]
    )

    # demand

    od_df = DataFrame(CSV.File(files.od_file))

    od_df = leftjoin(
        od_df,
        @select(node_df, :Node_ID, :New_Node_ID);
        on=:O_ID => :Node_ID,
        renamecols="" => "_O",
    )
    od_df = leftjoin(
        od_df,
        @select(node_df, :Node_ID, :New_Node_ID);
        on=:D_ID => :Node_ID,
        renamecols="" => "_D",
    )

    demand = Dict(
        collect(zip(od_df[!, :New_Node_ID_O], od_df[!, :New_Node_ID_D])) .=>
            od_df[!, :OD_Number],
    )
    origins_destinations = collect(keys(demand))
    origins = unique(map(first, origins_destinations))
    destinations = unique(map(last, origins_destinations))
    removed_od_pairs = Tuple{Int,Int}[]
    destination_free_flow_time = Dict{Int,Vector{Float64}}()

    # solution

    if solution == "AequilibraE"
        sol_network_df = DataFrame(CSV.File(files.aequilibrae_network_file))
        sol_assignment_df = DataFrame(CSV.File(files.aequilibrae_assignment_file))
        sol_df = leftjoin(sol_assignment_df, sol_network_df; on=:link_id)
        sol_df = leftjoin(
            sol_df,
            @select(node_df, :Node_ID, :New_Node_ID);
            on=:a_node => :Node_ID,
            renamecols="" => "_a_node",
        )
        sol_df = leftjoin(
            sol_df,
            @select(node_df, :Node_ID, :New_Node_ID);
            on=:b_node => :Node_ID,
            renamecols="" => "_b_node",
        )
        A = sol_df[!, :New_Node_ID_a_node]
        B = sol_df[!, :New_Node_ID_b_node]
        VAB = coalesce.(sol_df[!, :matrix_ab], 0.0)
        VBA = coalesce.(sol_df[!, :matrix_ba], 0.0)
        optimal_flow =
            dropzeros(sparse(A, B, VAB, n, n)) + dropzeros(sparse(B, A, VBA, n, n))
    elseif solution == "TransCAD"
        linkflows_df = DataFrame(CSV.File(files.transcad_linkflows_file))
        linkflows_df = @select(
            linkflows_df,
            :ID,
            :Dir,
            :AB_From_Node_,
            :AB_To_Node_ID,
            :BA_From_Node_,
            :BA_To_Node_ID,
            :AB_Flow,
            :BA_Flow,
            :AB_VOC,
            :BA_VOC,
            :Max_VOC,
        )
        linkflows_df = leftjoin(
            linkflows_df,
            @select(node_df, :Node_ID, :New_Node_ID);
            on=:AB_From_Node_ => :Node_ID,
            renamecols="" => "_A",
        )
        linkflows_df = leftjoin(
            linkflows_df,
            @select(node_df, :Node_ID, :New_Node_ID);
            on=:AB_To_Node_ID => :Node_ID,
            renamecols="" => "_B",
        )
        A = linkflows_df[!, :New_Node_ID_A]
        B = linkflows_df[!, :New_Node_ID_B]
        VAB = linkflows_df[!, :AB_Flow]
        VBA = linkflows_df[!, :BA_Flow]
        optimal_flow_AB = sparse(A, B, VAB, n, n)
        reverse_exists = (!ismissing).(VBA)
        optimal_flow_BA = sparse(
            B[reverse_exists], A[reverse_exists], identity.(VBA[reverse_exists]), n, n
        )
        optimal_flow = dropzeros(optimal_flow_AB) + dropzeros(optimal_flow_BA)
    else
        optimal_flow = missing
    end

    return TrafficAssignmentProblem(;
        dataset_name="UnifiedTrafficDataset",
        instance_name,
        # nodes
        nb_nodes,
        nb_links,
        real_nodes,
        zone_nodes,
        node_coord,
        valid_longitude_latitude=true,
        # links
        link_id,
        link_capacity,
        link_length,
        link_free_flow_time,
        link_speed_limit,
        link_bpr_mult,
        link_bpr_power,
        link_toll=missing,
        link_type,
        # demand
        demand,
        origins,
        destinations,
        removed_od_pairs,
        destination_free_flow_time,
        # cost
        toll_factor=missing,
        distance_factor=missing,
        # solution
        optimal_flow,
    )
end

"""
    postprocess!(pb)

Perform some data cleaning on a `TrafficAssignmentProblem`:

  - precompute free flow times to every destination present in the demand
  - remove OD pairs without a path between them (the free flow time is infinite)
"""
function postprocess!(pb::TrafficAssignmentProblem)
    (;
        link_id,
        link_free_flow_time,
        demand,
        destinations,
        destination_free_flow_time,
        removed_od_pairs,
    ) = pb
    W = eltype(link_free_flow_time)
    g_rev = SimpleWeightedDiGraph(transpose(link_free_flow_time), transpose(link_id))
    @tasks for d in destinations
        @set collect = true
        @local storage = DijkstraStorage(g_rev)
        dijkstra!(storage, g_rev, d)
        @one_by_one destination_free_flow_time[d] = copy(storage.dists)
    end
    for (o, d) in keys(demand)
        if destination_free_flow_time[d][o] == typemax(W)
            push!(removed_od_pairs, (o, d))
            delete!(demand, (o, d))  # modifies the instance
        end
    end
end

"""
    list_instances()
    list_instances(dataset_name)

Return a list of the available instances, given as a tuple with their dataset.
"""
function list_instances(dataset_name::AbstractString)
    data_dir = datapath(dataset_name)
    if dataset_name == "TransportationNetworks"
        instance_names = String[]
        for potential_name in readdir(data_dir)
            isdir(joinpath(data_dir, potential_name)) || continue
            for instance_file in readdir(joinpath(data_dir, potential_name))
                if endswith(instance_file, ".tntp")
                    push!(instance_names, potential_name)
                    break
                end
            end
        end
        return [(dataset_name, instance_name) for instance_name in instance_names]
    else
        instance_names = String[]
        for underscored_name in readdir(data_dir)
            instance_name = replace(underscored_name[4:end], "_" => " ")
            push!(instance_names, instance_name)
        end
        return [(dataset_name, instance_name) for instance_name in instance_names]
    end
end

list_instances() = mapreduce(list_instances, vcat, DATASET_NAMES)

"""
    summarize_instances()
    summarize_instances(dataset_name)

Return a `DataFrame` summarizing the dimensions of the available instances inside a datase.
"""
function summarize_instances(dataset_name::AbstractString)
    df = DataFrame(;
        dataset=String[],
        instance=String[],
        valid=Bool[],
        nodes=Int[],
        links=Int[],
        zones=Int[],
    )
    for (_, instance_name) in list_instances(dataset_name)
        yield()
        valid = false
        nn, nl, nz = (-1, -1, -1)
        try
            problem = TrafficAssignmentProblem(dataset_name, instance_name)
            valid = true
            nn = nb_nodes(problem)
            nl = nb_links(problem)
            nz = nb_zones(problem)
        catch exception
            @warn "Loading $instance_name from $dataset_name failed" exception
            # nothing
        end
        push!(
            df,
            (;
                dataset=dataset_name,
                instance=instance_name,
                valid,
                nodes=nn,
                links=nl,
                zones=nz,
            ),
        )
    end
    return df
end

summarize_instances() = mapreduce(summarize_instances, vcat, DATASET_NAMES)
