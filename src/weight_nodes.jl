"""
    create_graph_weights(graph, geodata, weightcols, distance)

Using a graph and a spatial data frame, assign weights to graph nodes.

Each row in the spatial data frame should have a weight associated with it
(if they are all equivalent, the weight can be set to a column of ones). For each
row, all graph edges within `distance` of the geometry of that row will be identified
(the spatial data frame should be in the same projection as the network data used to
create the graph). An equal amount of weight will be assigned to the nodes at each end of
all identified edges. (In a corner, twice as much weight will be assigned to the node at the
corner, as it will receive weight from each of the links adjacent to the property. The effect
may be even more pronounced if the streets on the other side of the intersection from the row
are also within `distance`).

`weightcols` specifies which columns to use as weights. It must be a vector. Multiple weight columns
can be specified. The function will return an n x w matrix, where n is the number of nodes in the graph,
and w is the number of weightcols specified. This will contain the weight for each node for each weightcol.
"""
function create_graph_weights(G, gdf, weightcols, distance)
    geomcol = first(metadata(gdf, "geometrycolumns"))
    weights = zeros(Float64, (nv(G), length(weightcols)))
    @info "Indexing graph"
    edgeidx, edges = index_graph_edges(G)

    for (i, row) in enumerate(eachrow(gdf))
        if i % 1000 == 0
            @info "Processed $i / $(nrow(gdf)) locations"
        end

        geom = row[geomcol]
        env = ArchGDAL.envelope(geom)
        nodes = Int64[] # node codes, not ids
        # expand the envelope to account for buffer
        for candidate_idx in intersects(edgeidx, [env.MinX - distance, env.MinY - distance], [env.MaxX + distance, env.MaxY + distance])
            edge = edges[candidate_idx]
            edge_geom = G[edge[1], edge[2]].geom
            if ArchGDAL.distance(geom, edge_geom) < distance
                # this edge is close by
                # some nodes will be added more than once, and that's okay - they are attached
                # to two nearby edges, so get twice the weight
                push!(nodes, code_for(G, edge[1]))
                push!(nodes, code_for(G, edge[2]))
            end
        end

        for (i, weightcol) in enumerate(weightcols)
            # can't use .+= here as there may (often) be duplicate nodes when
            # there are multiple close-by edges. x[[1, 1]] .+= 2 only adds two once.
            # https://discourse.julialang.org/t/unexpected-behavior-of-vectorized-with-duplicate-indices/112537/4
            for node in nodes
                weights[node, i] += row[weightcol] / length(nodes)
            end
        end
    end

    return weights
end

"""
    partition_weights(G, Gsub, weights)

Given weights `weights` that go with graph `G`, extract a weight array that goes with GraphPartition `Gsub`.
"""
function partition_weights(G::T, Gsub::GraphPartition{T}, weights::AbstractVector{W}) where {T, W}
    output = zeros(W, nv(Gsub.G))

    for (idx, wt) ∈ pairs(weights)
        label = label_for(G, idx)
        if haskey(Gsub.G, label)
            output[code_for(Gsub.G, label)] = wt
        end
    end

    return output
end