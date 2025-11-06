"""
    regional_access(decay_func, G, dmat, link; resolution=50)

Perform a regional access analysis (change in access for each point) as a result of the link. Return
a Rasters.jl raster with the values at each point in the area.
"""
function regional_access_oneway(decay_func, G, dmat, link::CandidateLink, weights, decay_cutoff; resolution=25, index=nothing, max_snap_dist=100)
    # first, figure out extent of analysis
    fredge = G[link.fr_edge_src, link.fr_edge_tgt].geom
    frpt = LibGEOS.interpolate(gdal_to_geos(fredge), link.fr_dist_from_start)
    toedge = G[link.to_edge_src, link.to_edge_tgt].geom
    topt = LibGEOS.interpolate(gdal_to_geos(toedge), link.to_dist_from_start)

    # the extents where this could have an effect
    minx = min(LibGEOS.getGeomX(frpt), LibGEOS.getGeomX(topt)) - decay_cutoff
    maxx = max(LibGEOS.getGeomX(frpt), LibGEOS.getGeomX(topt)) + decay_cutoff
    miny = min(LibGEOS.getGeomY(frpt), LibGEOS.getGeomY(topt)) - decay_cutoff
    maxy = max(LibGEOS.getGeomY(frpt), LibGEOS.getGeomY(topt)) + decay_cutoff

    # y iterates negatively for north-up in most projected coordinate systems
    xdim, ydim = Rasters.X(minx:resolution:(maxx + resolution)), Rasters.Y(maxy:-resolution:(miny - resolution))
    rast = Rasters.Raster(zeros(xdim, ydim))

    # Index the graph edges
    index = isnothing(index) ? index_graph_edges(G) : index

    @showprogress for (xi, yi) in Iterators.product(1:size(rast, 1), 1:size(rast, 2))
        x = Rasters.lookup(xdim)[xi]
        y = Rasters.lookup(ydim)[yi]

        edge, dist, dist_from_start, dist_to_end = closest_edge(G, [x, y], index)
        if dist > max_snap_dist
            # not close to edge
            continue
        end

        # only looking at end of link. whole function gets called again for traversing link in other direction.
        distance_from_end_of_link_to_pixel = add_unless_typemax(
            compute_net_distance(G, dmat,
                link.to_edge_src, link.to_edge_tgt, link.fr_dist_from_start, link.fr_dist_to_end,
                edge..., dist_from_start, dist_to_end),
            dist
        )

        if distance_from_end_of_link_to_pixel > decay_cutoff
            continue # too far from link
        end

        # other vertex - could be origin or destination depending on whether origin or dest weights
        # were supplied.
        # TODO this won't work on digraphs
        for other_v in eachindex(weights)
            if weights[other_v] == 0
                continue # can't affect access 
            end

            distance_from_vertex_to_start_of_link = min(
                add_unless_typemax(dmat[other_v, code_for(G, link.fr_edge_src)], link.fr_dist_from_start),
                add_unless_typemax(dmat[other_v, code_for(G, link.fr_edge_tgt)], link.fr_dist_to_end)
            )

            if distance_from_vertex_to_start_of_link > decay_cutoff
                continue # can't have any effect
            end

            old_distance = min(
                add_unless_typemax(dmat[other_v, code_for(G, edge[1])], dist_from_start + dist),
                add_unless_typemax(dmat[other_v, code_for(G, edge[2])], dist_to_end + dist)
            )

            # this should not overflow, b/c both to start and from end lead to short circuits above if typemax
            new_distance = OverflowContexts.@checked distance_from_vertex_to_start_of_link + link.geographic_length_m + distance_from_end_of_link_to_pixel

            if new_distance < old_distance
                Δdecayed = decay_func(new_distance)
                if old_distance < typemax(eltype(dmat))
                    Δdecayed -= decay_func(old_distance)
                end

                rast[xi, yi] += Δdecayed * weights[other_v]
            end
        end
    end

    # TODO this won't work with digraphs
    return rast
end

function regional_access(decay_func, G, dmat, link::CandidateLink, weights, decay_cutoff; resolution=50, index=nothing, max_snap_dist=100)
    index = isnothing(index) ? index_graph_edges(G) : index

    regional_access_oneway(decay_func, G, dmat, link::CandidateLink, weights, decay_cutoff; resolution=resolution, index=index, max_snap_dist=max_snap_dist) +
        regional_access_oneway(decay_func, G, dmat, reverse(link::CandidateLink), weights, decay_cutoff; resolution=resolution, index=index, max_snap_dist=max_snap_dist)
end