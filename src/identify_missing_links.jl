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
    fr::Int64
    to::Int64
    geographic_length_m::Float64
    network_length_m::T
end

"""
This identifies possible missing links
"""
function identify_potential_missing_links(G, dmat::Matrix{T}, max_link_dist, min_net_dist) where T
    @info "Indexing graph"
    sidx = RTree(2)

    for node in labels(G)
        loc = collect(G[node])
        insert!(sidx, code_for(G, node), loc, loc)
    end

    @info "Searching for candidate links"

    links = Vector{CandidateLink}()

    # now, find nodes that are nearby geographically, but far away in network space
    for node in labels(G)
        # find nearby nodes
        loc = collect(G[node])
        nidx = code_for(G, node)
        map(intersects(sidx, loc .- max_link_dist, loc .+ max_link_dist)) do candidate
            # from is always less than to, to avoid duplicate links a->b and b->a
            if candidate > nidx
                geo_dist = norm2(loc .- G[label_for(G, candidate)])
                # check distance in both directions. should be no-op with walking, but with
                # biking we may ultimately use a directed graph.
                net_dist = min(dmat[nidx, candidate], dmat[candidate, nidx])
                if geo_dist ≤ max_link_dist && net_dist ≥ min_net_dist
                    push!(links, CandidateLink{T}(nidx, candidate, geo_dist, net_dist))
                end
            end
        end
    end

    return links
end

"""
Convenience function that converts links to a GeoDataFrame
"""
function links_to_gdf(G, links)
    gdf = DataFrame(links)

    gdf.network_length_m = ifelse.(
        gdf.network_length_m .≠ typemax(eltype(gdf.network_length_m)),
        convert.(Float64, gdf.network_length_m),
        Inf64
    )

    # add geometry
    gdf.geometry = map(zip(gdf.fr, gdf.to)) do (fr, to)
        orig = G[label_for(G, fr)]
        dest = G[label_for(G, to)]
        ArchGDAL.createlinestring([orig, dest])
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