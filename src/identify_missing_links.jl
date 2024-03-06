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

struct CandidateLink{T}
    fr_edge_src::VertexID
    fr_edge_tgt::VertexID
    fr_dist_from_start::T
    fr_dist_to_end::T
    to_edge_src::VertexID
    to_edge_tgt::VertexID
    to_dist_from_start::T
    to_dist_to_end::T
    geographic_length_m::Float64
    network_length_m::T
end

saturated_add(x::T, y) where T = first(x) == typemax(T) ? x : x + round(T, y)

"""
This identifies possible missing links
"""
function identify_potential_missing_links(G, dmat::Matrix{T}, max_link_dist, min_net_dist) where T
    @info "Indexing graph"
    sidx, edges = index_graph_edges(G)

    @info "Searching for candidate links"

    links = Vector{CandidateLink}()

    # Find edges that are nearby geographically, but far in network space
    # This is kinda slow, multithreading might help, but need to look into thread safety
    for (i, source_edge) in enumerate(edge_labels(G))
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
            x -> edges[x],
            intersects(sidx,
                [source_edge_envelope.MinX, source_edge_envelope.MinY] .- max_link_dist,
                [source_edge_envelope.MaxX, source_edge_envelope.MaxY] .+ max_link_dist
            )
        )

        for candidate in candidates
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
            if upper_bound_net_dist > min_net_dist - source_edge_length_m - candidate_edge_length_m
                candidate_edge_geom = G[candidate...].geom
                # we might have a missing link. Calculate geographic distance.
                geo_distance = LibGEOS.distance(source_edge_geom, candidate_edge_geom)

                if geo_distance ≤ max_link_dist
                    # Now, find the places that are closest
                    # possible optimization: move outside loop. But in many cases this will never
                    # get hit as the conditionals will be false.
                    point_on_source, point_on_candidate = LibGEOS.nearestPoints(source_edge_geom, candidate_edge_geom)
                    length_m = LibGEOS.distance(point_on_source, point_on_candidate)
                    isapprox(length_m, geo_distance, atol=1e-6) ||
                        error("Nearest points do not have nearest distance: expected $geo_distance, got $length_m, linking edge $source_edge to $candidate")

                    source_dist = LibGEOS.project(source_edge_geom, point_on_source)
                    candidate_dist = LibGEOS.project(candidate_edge_geom, point_on_candidate)

                    # calculate the actual network distance from these points
                    net_dist_m = compute_net_distance(dmat, source_edge_fr, source_edge_to, source_dist, source_edge_length_m, candidate_edge_fr, candidate_edge_to, candidate_dist, candidate_edge_length_m)

                    # calculate the distances in UInt16 units
                    source_dist_from_start = round(T, source_dist)
                    source_dist_to_end = round(T, source_edge_length_m) - source_dist_from_start

                    candidate_dist_from_start = round(T, candidate_dist)
                    candidate_dist_to_end = round(T, candidate_edge_length_m) - candidate_dist_from_start

                    # re-check net dist now that we have an actual value, not an upper bound
                    if net_dist_m > min_net_dist
                        push!(links, CandidateLink{T}(
                            source_edge[1],
                            source_edge[2],
                            source_dist_from_start,
                            source_dist_to_end,
                            candidate[1],
                            candidate[2],
                            candidate_dist_from_start,
                            candidate_dist_to_end,
                            length_m,
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
Compute the network distance between the two points on a candidate link, by computing
between from ends and adding fractions of the edge.
"""
function compute_net_distance(dmat::Matrix{T}, sfr, sto, sdist, slength, dfr, dto, ddist, dlength) where T
    sdist .≤ slength + 1e-12 || error("Source dist greater than overall source length")
    ddist .≤ dlength + 1e-12 || error("Source dist greater than overall source length")

    senddist = slength - sdist
    denddist = dlength - ddist
    min(
        saturated_add(dmat[sfr, dfr], sdist + ddist),
        saturated_add(dmat[sto, dfr], senddist + ddist),
        saturated_add(dmat[sfr, dto], sdist + denddist),
        saturated_add(dmat[sto, dto], senddist + denddist),
        saturated_add(dmat[dfr, sfr], ddist + sdist),
        saturated_add(dmat[dto, sfr], denddist + sdist),
        saturated_add(dmat[dfr, sto], ddist + senddist),
        saturated_add(dmat[dto, sto], denddist + senddist)
    )
end

"""
Convenience function that converts links to a GeoDataFrame
"""
function links_to_gdf(G, links, scores=ones(length(links)))
    gdf = DataFrame(links)

    gdf.fr_edge_src = [x.id for x in gdf.fr_edge_src]
    gdf.to_edge_src = [x.id for x in gdf.to_edge_src]
    gdf.fr_edge_tgt = [x.id for x in gdf.fr_edge_tgt]
    gdf.to_edge_tgt = [x.id for x in gdf.to_edge_tgt]

    gdf.network_length_m = ifelse.(
        gdf.network_length_m .≠ typemax(eltype(gdf.network_length_m)),
        convert.(Float64, gdf.network_length_m),
        Inf64
    )

    gdf.score = scores

    # add geometry
    gdf.geometry = map(links) do link
        startg = GeoInterface.convert(ArchGDAL, LibGEOS.interpolate(gdal_to_geos(G[link.fr_edge_src, link.fr_edge_tgt].geom), link.fr_dist_from_start))
        endg = GeoInterface.convert(ArchGDAL, LibGEOS.interpolate(gdal_to_geos(G[link.to_edge_src, link.to_edge_tgt].geom), link.to_dist_from_start))
        ArchGDAL.createlinestring([
            # this can't possibly be the best way to do this

            [ArchGDAL.getx(startg, 0), ArchGDAL.gety(startg, 0)],
            [ArchGDAL.getx(endg, 0), ArchGDAL.gety(endg, 0)]
        ])
    end

    return gdf
end

mutable struct SphereOfInfluence
    nodes::Vector{Int64}
    link::CandidateLink
end

"""
This is a heuristic function to deduplicate links. Consider the situation below.
Double lines are existing roads. Here, you have the ends of two blocks in different subdivisions,
that do not connect.

    ====a   c===
       ||   ||
    ====b   d===

There are four ways to connect the streets above (assuming we only make connections at nodes):
AC, AD, BC, and BD. Clearly, they all provide essentially the same level of access, and there's no
reason to consider all of them.

    ====a---c===
       || X ||
    ====b---d===

This function greedily groups connections using the following algorithm. For each connection,
it defines a "sphere of influence" which is all nodes within 100m network distance of the from node,
and all nodes within 100m network distance of the to node. If any subsequent link connects two nodes
in this sphere of influence, whichever link has the shortest geographic distance is retained.

Note that, with an undirected graph and as long as the min network distance for proposing a link is more
than twice the radius of the sphere of influence, any link connecting two nodes in the sphere of influence is
perforce connecting the area around the from node to the area around the to node or vice versa. Consider the
counterfactual: for any two nodes that are within x meters of one of the nodes of the link, they are at most
2x meters apart (because there is a path through the node).
"""
function deduplicate_links(links, dmat, sphere_of_influence_radius)
    spheres_of_influence = SphereOfInfluence[]
    
    for link in links
        # check if this link is in a sphere of influence
        in_soi = false
        for soi in spheres_of_influence
            if link.fr ∈ soi.nodes && link.to ∈ soi.nodes
                # note: order matters here. it is possible for a link to be in two spheres of influence,
                # so ordering might change results.
                in_soi = true

                # if it's better than the best link found for this soi, retain it
                if link.geographic_length_m < soi.link.geographic_length_m
                    soi.link = link
                end

                # don't keep considering spheres of influence
                break
            end
        end

        if !in_soi
            # not in a sphere of influence, create one
            nodes = vcat(
                findall(dmat[:, link.fr] .< sphere_of_influence_radius),
                findall(dmat[:, link.to] .< sphere_of_influence_radius)
            )

            push!(spheres_of_influence, SphereOfInfluence(nodes, link))
        end
    end

    return collect(getfield.(spheres_of_influence, :link))
end