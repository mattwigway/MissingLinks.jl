# This contains code used to troubleshoot the graph and detect common errors

"""
    find_disconnected_crossings(G, dmat; tolerance=10)

Find locations where edges of graph G cross without intersecting, and where the network distance between
the intersecting points is more than `tolerance`. 

These are often graph errors, but may also be overpasses,
tunnels, etc. Returns a GeoDataFrames-compatible point DataFrame with the locations, for further examination
in GIS.
"""
function find_disconnected_crossings(G, dmat; tolerance=10)
    # index graph edges
    eidx, edges = index_graph_edges(G)

    # why two different types of geometry??
    result = Union{ArchGDAL.IGeometry{ArchGDAL.wkbPoint}, ArchGDAL.IGeometry{ArchGDAL.wkbPoint25D}}[]

    for (i, sedge) in enumerate(edge_labels(G))
        if i % 10000 == 0
            @info "Processed $i / $(ne(G)) edges"
        end
        sgeom = G[sedge...].geom
        source_edge_envelope = ArchGDAL.envelope(sgeom)
        candidates = intersects(eidx,
            [source_edge_envelope.MinX, source_edge_envelope.MinY],
            [source_edge_envelope.MaxX, source_edge_envelope.MaxY]
        )

        for candidate in candidates
            tedge = edges[candidate]

            # don't duplicate points in both directions
            if tedge != sedge && tedge[1] ≥ sedge[1]
                tgeom = G[tedge...].geom

                if ArchGDAL.intersects(sgeom, tgeom)
                    for_each_geom(ArchGDAL.intersection(sgeom, tgeom)) do pt
                        if pt isa ArchGDAL.IGeometry{ArchGDAL.wkbLineString} || pt isa ArchGDAL.IGeometry{ArchGDAL.wkbLineString25D}
                            # mark start and end of intersection
                            #push!(result, ArchGDAL.createpoint(get_first_point(pt)))
                            #push!(result, ArchGDAL.createpoint(get_last_point(pt)))
                        else
                            sdist = LibGEOS.project(sgeom, pt)
                            tdist = LibGEOS.project(tgeom, pt)
                            dist = compute_net_distance(G, dmat,
                                code_for(G, sedge[1]), code_for(G, sedge[2]), sdist, G[sedge...].length_m - sdist,
                                code_for(G, tedge[1]), code_for(G, tedge[2]), tdist, G[tedge...].length_m - tdist
                            )

                            if dist > tolerance
                                push!(result, pt)
                            end
                        end
                    end
                end
            end
        end
    end

    df = DataFrame(:geom=>result)
    metadata!(df, "geometrycolumns", (:geom,))

    return df
end

"""
    find_dead_ends(G)

Find locations in graph `G` where there is a dead end.

Returns a GeoDataFrame. This is useful for network data creation as by opening
this layer in GIS you can easily see where there are existing
dead ends and determine if they are correct or not.
"""
function find_dead_ends(G)
    result = ArchGDAL.IGeometry{ArchGDAL.wkbPoint}[]

    for node in 1:nv(G)
        if Graphs.degree(G, node) == 1
            position = G[label_for(G, node)]
            push!(result, ArchGDAL.createpoint(collect(position)))
        end
    end

    @info "Found $(length(result)) degree-one nodes"

    df = DataFrame(:geom => result)
    metadata!(df, "geometrycolumns", (:geom,))

    return df
end