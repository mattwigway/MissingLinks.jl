"""
Snap the point origin to the nearest edge in G, returning a namedtuple with fields:

- src: Source vertex label
- tgt: Target vertex label
- dist_from_src: distance from the source to the closest point on the edge
- dist_from_tgt: distance from the target to the closest point on the edge
- perpendicular_dist: distance from the requested point to the closest point on the edge

or nothing if there is no edge within snap_tolerance
"""
function find_point_on_edge(origin::ArchGDAL.IGeometry{ArchGDAL.wkbPoint}, G, edge_index, snap_tolerance)
    x, y, _ = ArchGDAL.getpoint(origin, 0)

    candidates = map(e -> edge_index[2][e], intersects(edge_index[1], [x, y] .- snap_tolerance, [x, y] .+ snap_tolerance))

    best_candidate = nothing
    best_dist = Inf64

    for candidate in candidates
        candidate_geom = G[candidate...].geom
        dist = ArchGDAL.distance(origin, candidate_geom)
        # gotta check distance, spatial index overselects
        if dist < snap_tolerance && dist < best_dist
            best_candidate = candidate
            best_dist = dist
        end
    end

    if !isnothing(best_candidate)
        # find the closest point
        egeom = G[best_candidate...].geom
        _, point_on_edge = LibGEOS.nearestPoints(origin, egeom)
        dist_from_start = LibGEOS.project(egeom, point_on_edge)
        dist_from_end = LibGEOS.geomLength(egeom) - dist_from_start

        return (
            src=best_candidate[1],
            tgt=best_candidate[2],
            dist_from_src=dist_from_start,
            dist_from_tgt=dist_from_end,
            perpendicular_dist=best_dist
        )
    else
        # nothing within snap_tolerance
        return nothing
    end
end

"""
Compute a service area starting at a particular location, based on a precomputed distance matrix, graph,
and an optional list of links to include.

Supply edge_index=<result of index_graph_edges(G)> if you precomputed the edge spatial index to avoid recomputing
each time the function is run.
"""
function service_area(origin::ArchGDAL.IGeometry{ArchGDAL.wkbPoint}, G, dmat, threshold;
        links=nothing, scores=nothing, edge_index=nothing, snap_tolerance=100)
    eidx = if isnothing(edge_index)
        index_graph_edges(G)
    else
        edge_index
    end

    # snap the origin to an edge
    snapped = find_point_on_edge(origin, G, eidx, snap_tolerance)

    if isnothing(snapped)
        return nothing
    end

    src_from_origin = round(eltype(dmat), snapped.dist_from_src + snapped.perpendicular_dist)
    tgt_from_origin = round(eltype(dmat), snapped.dist_from_tgt + snapped.perpendicular_dist)

    src = code_for(G, snapped.src)
    tgt = code_for(G, snapped.tgt)

    # find distances to all nodes
    dists_from_src = add_unless_typemax.(dmat[:, src], src_from_origin)
    dists_from_tgt = add_unless_typemax.(dmat[:, tgt], tgt_from_origin)

    # account for new links
    if !isnothing(links)
        for origlink in links
            for dists in (dists_from_src, dists_from_tgt)
                for dest in eachindex(dists)
                    for link in (origlink, reverse(origlink))
                        # shortest distance to start of link
                        dist_to_start_of_link = min(
                            add_unless_typemax(dists[link.fr_edge_src], link.fr_dist_from_start), # already includes dist from origin
                            add_unless_typemax(dists[link.fr_edge_tgt], link.fr_dist_to_end) # already includes dist from origin
                        )

                        dist_from_end_of_link = min(
                            add_unless_typemax(dmat[link.to_edge_src, dest], link.to_dist_from_start),
                            add_unless_typemax(dmat[link.to_edge_tgt, dest], link.to_dist_to_end),
                        )

                        if dist_to_start_of_link < threshold && dist_from_end_of_link < threshold
                            dist_to_dest_via_link = dist_to_start_of_link + dist_from_end_of_link + link.geographic_length_m
                            if dist_to_dest_via_link < dists[dest]
                                dists[dest] = dist_to_dest_via_link
                            end
                        end
                    end
                end
            end
        end
    end

    dists_from_point = min.(dists_from_src, dists_from_tgt)

    geoms = ArchGDAL.IGeometry{ArchGDAL.wkbLineString}[]

    # make the service area geography
    for (esrc, etgt) in edge_labels(G)
        esrccode = code_for(G, esrc)
        etgtcode = code_for(G, etgt)

        if dists_from_point[esrccode] < threshold && dists_from_point[etgtcode] < threshold
            # the whole edge is in the service area
            push!(geoms, G[esrc, etgt].geom)
        elseif dists_from_point[esrccode] < threshold
            # only the start is in the service area
            dist_to_threshold = threshold - dists_from_point[esrccode]
            push!(geoms, geom_between(G[esrc, etgt].geom, 0, dist_to_threshold))
        elseif dists_from_point[etgtcode] < threshold
            # only the end is in the service area
            dist_to_threshold = threshold - dists_from_point[etgtcode]
            len = ArchGDAL.geomlength(G[esrc, etgt].geom)
            push!(geoms, geom_between(G[esrc, etgt].geom, len - dist_to_threshold, len))
        end
    end

    # union everything to a single feature so multiple service areas can easily be exported in a single file
    outgeom = ArchGDAL.createmultilinestring()
    for g in geoms
        ArchGDAL.addgeom!(outgeom, g)
    end

    score = if !isnothing(scores)
        sum(scores[dists_from_point .≤ threshold])
    else
        0.0
    end

    return (geom=outgeom, score=score)
end