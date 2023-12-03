"""
Using a graph and a spatial data frame, assign weights to graph nodes. This works as follows:

1. For each entry in the data frame
    - Find all graph edges that are within a given distance
    - Divide each of the weights by the number of graph edges times two
    - Assign that divided weight to each of the start and end nodes of those edges
"""
function create_graph_weights(G, gdf, weightcols, distance; geomcol=:geometry)
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
            weights[nodes, i] .+= row[weightcol] / length(nodes)
        end
    end

    return weights
end