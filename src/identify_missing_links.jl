# The algorithm works like this:

# First, we identify pairs of nodes that are geographically within a short distance, but are at least
# a minimum distance away on the network

# Then, for each pair, we identify how much making this connection would increase the aggregate number
# of origin-weighted destinations (i.e. the sum across origins of the number of destinations accessible within
# a threshold distance). We do this in several steps:
#  1. Identify all origins within the threshold of the first node of the proposed new connection
#  2. For each origin:
#    a. generate distances to each destination via the new link by summing
#      - Distance from origin to link
#      - Length of link
#      - Distance from link to each destination
#    b. generate the shortest distance to each destination, by taking min(existing distances, distances via new link)
#    c. calculate the (weighted) access implied these distances
#    d. subtract the (weighted) access for this origin originally
#  3. Sum the weighted access deltas
#  4. Repeat with the second node of the new link
#  5. Sum the results. This is the contribution to access of this link.

"""
Add y to x, but only if x is not already typemax (used to indicate unreachable). This assumes that
if x < typemax(typeof(T)), x + y < typemax(typeof(T)), which is generally reasonable as x and y are in meters
and anything we're adding is well less than 65km, unless it's unreachable, which should only be the case for
x (coming from the distance matrix). However, we do use a checked_add to throw an error if that assumption is violated.

If x == typemax(T), return x
If x < typemax(T) and x + y <= typemax(T), return x + y
If x < typemax(T) and x + y > typemax(T), throw OverflowError
"""
add_unless_typemax(x::T, y) where T = x == typemax(T) ? x : Base.Checked.checked_add(x, round(T, y))

"""
    identify_potential_missing_links(graph, distance_matrix, max_link_dist, min_net_dist)

Identify locations in `graph` that are less than `max_link_dist` apart geographically, but at least `min_net_dist` apart
in network distance (or disconnected entirely).

`dmat` should be a distance matrix for all nodes in the graph, generally created by [`fill_distance_matrix!`](@ref fill_distance_matrix!)
"""
function identify_potential_missing_links(G, dmat::Matrix{T}, max_link_dist, min_net_dist) where T
    @info "Indexing graph"
    sidx, edges = index_graph_edges(G)

    @info "Searching for candidate links"

    links = Vector{CandidateLink{T}}()

    # Find edges that are nearby geographically, but far in network space
    # This is kinda slow, multithreading might help, but need to look into thread safety
    # order_vertices ensures that v1 < v2 so geometry is coded correctly
    # TODO directed graph
    for (i, source_edge) in enumerate(order_vertices.(edge_labels(G)))
        if i % 10000 == 0
            @info "Processed $i / $(ne(G)) edges"
        end

        source_edge_geom = G[source_edge...].geom
        source_edge_length_m = G[source_edge...].length_m
        source_edge_envelope = ArchGDAL.envelope(source_edge_geom)
        source_edge_fr = code_for(G, source_edge[1])
        source_edge_to = code_for(G, source_edge[2])

        # find other edges whose bounding boxes intersect the source_edge edge bounding box
        # expanded by the max_link_dist
        candidates = map(
            # optimization - could discard edges that have already been used as a source
            x -> edges[x],
            intersects(sidx,
                [source_edge_envelope.MinX, source_edge_envelope.MinY] .- max_link_dist,
                [source_edge_envelope.MaxX, source_edge_envelope.MaxY] .+ max_link_dist
            )
        )

        # all geometries are coded with lower-numbered node first
        # TODO directed graph
        for candidate in order_vertices.(candidates)
            # first, figure out if this is even worth doing - often they are going to be connected
            # to one another fairly directly
            candidate_edge_fr = code_for(G, candidate[1])
            candidate_edge_to = code_for(G, candidate[2])
            candidate_edge_length_m = G[candidate...].length_m

            # an upper bound on the net dist is the shortest distance from either end to either other
            # end plus the total length of both edges (if the closest points geographically
            # happened to be opposite the closest points topologically).

            # To avoid overflow, instead of adding the lengths here, we subtract them from the max
            # in the next line.
            upper_bound_net_dist = min(
                dmat[source_edge_fr, candidate_edge_to],
                dmat[source_edge_to, candidate_edge_to],
                dmat[source_edge_fr, candidate_edge_fr],
                dmat[source_edge_to, candidate_edge_fr],
                dmat[candidate_edge_to, source_edge_to],
                dmat[candidate_edge_fr, source_edge_to],
                dmat[candidate_edge_to, source_edge_fr],
                dmat[candidate_edge_fr, source_edge_fr]
            )
            
            # we do the subtraction here instead of adding above to avoid overflow
            # if upper_bound_net_dist is typemax(UInt16)
            # the edge lengths are Float64, so this should not overflow as promotion will happen
            # Furthermore, even if types were not promoted, this cannot underflow with the default parameters,
            # as min_net_dist is 1000m and edges can be no longer than 250m.
            # doing a checked subtract here would not be the world's worst idea
            if upper_bound_net_dist > min_net_dist - source_edge_length_m - candidate_edge_length_m
                candidate_edge_geom = G[candidate...].geom
                # we might have a missing link. Calculate geographic distance.
                geo_distance = LibGEOS.distance(source_edge_geom, candidate_edge_geom)

                if geo_distance ≤ max_link_dist
                    # Now, find the places that are closest
                    point_on_source, point_on_candidate = LibGEOS.nearestPoints(source_edge_geom, candidate_edge_geom)
                    length_m = LibGEOS.distance(point_on_source, point_on_candidate)
                    isapprox(length_m, geo_distance, atol=1e-6) ||
                        error("Nearest points do not have nearest distance: expected $geo_distance, got $length_m, linking edge $source_edge to $candidate")

                    source_dist = LibGEOS.project(source_edge_geom, point_on_source)
                    candidate_dist = LibGEOS.project(candidate_edge_geom, point_on_candidate)

                    # calculate the distances in UInt16 units
                    source_dist_from_start = round(T, source_dist)
                    source_dist_to_end = round(T, source_edge_length_m) - source_dist_from_start

                    candidate_dist_from_start = round(T, candidate_dist)
                    candidate_dist_to_end = round(T, candidate_edge_length_m) - candidate_dist_from_start

                    # calculate the actual network distance from these points
                    net_dist_m = compute_net_distance(G, dmat, source_edge..., source_dist_from_start, source_dist_to_end,
                        candidate..., candidate_dist_from_start, candidate_dist_to_end)

                    # re-check net dist now that we have an actual value, not an upper bound
                    if net_dist_m > min_net_dist
                        push!(links, CandidateLink{T}(
                            source_edge...,
                            source_dist_from_start,
                            source_dist_to_end,
                            candidate...,
                            candidate_dist_from_start,
                            candidate_dist_to_end,
                            round(T, length_m),
                            net_dist_m
                        ))
                    end
                end
            end            
        end
    end

    return links
end

"""
    links_to_gis(graph, links, pairs...)

Convenience function that converts a vector of links (and associated graph) into a GeoDataFrame, suitable
for export to a GIS file. 

`pairs`` should be pairs of :col_name=>values that will be added to the output. For example, you will often
want to add scores to your output.
"""
function links_to_gis(G, links, pairs...)
    gdf = DataFrame(links)

    gdf.network_length_m = ifelse.(
        gdf.network_length_m .≠ typemax(eltype(gdf.network_length_m)),
        convert.(Float64, gdf.network_length_m),
        Inf64
    )

    # add geometry
    gdf.geometry = map(links) do link
        startg = GeoInterface.convert(ArchGDAL, LibGEOS.interpolate(gdal_to_geos(G[label_for(G, link.fr_edge_src), label_for(G, link.fr_edge_tgt)].geom), link.fr_dist_from_start))
        endg = GeoInterface.convert(ArchGDAL, LibGEOS.interpolate(gdal_to_geos(G[label_for(G, link.to_edge_src), label_for(G, link.to_edge_tgt)].geom), link.to_dist_from_start))
        ArchGDAL.createlinestring([
            # this can't possibly be the best way to do this

            [ArchGDAL.getx(startg, 0), ArchGDAL.gety(startg, 0)],
            [ArchGDAL.getx(endg, 0), ArchGDAL.gety(endg, 0)]
        ])
    end

    for pair in pairs
        gdf[!, pair[1]] = pair[2]
    end

    return gdf
end
