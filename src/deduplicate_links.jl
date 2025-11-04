
mutable struct SphereOfInfluence
    nodes::Vector{VertexID}
    original_link::CandidateLink
    links::Vector{CandidateLink}
end

"""
    deduplicate_links(list, distance_matrix, sphere_of_influence_radius)

Deduplicate links using a heuristic.

Arguments are a list of links (e.g. as produced by [`identify_potential_missing_links`](@ref identify_potential_missing_links)),
a distance matrix between nodes for the graph used to identify those links, and the radius of the sphere of influence (see below).
    
Consider the situation below. Lines are existing
roads. Here, you have the ends of two blocks in different subdivisions, that do not connect.

```
--a--o    o--d--
     |    |
     c    e
     |    |
--b--o    o--f--     
```

Since there are three edges in each subdivision in the graph above, all of which pass within (say)
100m of each other, there are nine possible candidate links: ad, ae, af, bd, be, bf, cd, ce, cf.
Clearly, all of these provide essentially the same level of access, so we want to deduplicate them.

We consider two candidate links to be duplicates if they connect pairs of locations within a
`sphere_of_influence_radius` of one another (we use 100m in the paper). To identify these links, we use the
following greedy algorithm. Note that since it is a greedy algorithm, the order the links are
identified in matters. Currently, link identification is not multithreaded, so the order of the
links should be deterministic, and therefore you should get the same results from this deduplication
whenever it is run.

For the first candidate link, we identify the "sphere of influence" of that link; the sphere of
influence is all nodes that are within 100m of either end of the candidate link. We calculate this
using the node-to-node distance matrix and the distances of the link from the start and end of the
edges it connects. We do not need to define separate spheres of influence for the start and end of
the link—as long as the radius of the sphere of influence is less than half the minimum network
distance for proposing a link, any candidate link that connects two nodes in the sphere of influence
must have one end near each end of the original link that defined the sphere of influence.

For subsequent candidate links, we determine if the link is within an already-identified sphere of
influence. We check to see if one of the ends of each of the edges the candidate link connects is in
an already identified sphere of influence. If it is, we further check to see if the network
distances between the ends of this candidate link and the corresponding ends of candidate link that
    defined the sphere of influence are less than the threshold for both ends of the links (note
that links are undirected, so we explicitly consider links duplicates even if the start of one is
near the end of the other and vice versa). If these network distances are both less than the radius
of the sphere of influence, we add the link to this sphere of influence and continue to the next
link. If they are not, we continue to the next sphere of influence. If a link is not within any
existing sphere of influence, we define a new sphere of influence based on the link.

We then return the link with the shortest geographic distance from each sphere of influence. Note
that there may still be links that are close to one another; if two links define adjacent or even
overlapping spheres of influence, but with neither link in either sphere of influence, the best
links in the two spheres of influence may be very similar ones in between the original two links.
You could iterate the algorithm if you wanted to to prevent this, but this could result in links
being supplanted by ones more than 100m away.
"""
function deduplicate_links(G, links::AbstractVector{<:CandidateLink{<:Any}}, dmat, sphere_of_influence_radius)
    spheres_of_influence = SphereOfInfluence[]
    
    for (i, link) in enumerate(links)
        if i % 10000 == 0
            @info "Processed $i / $(length(links)) links"
        end

        # check if this link is in a sphere of influence
        in_soi = false
        for soi in spheres_of_influence
            # one end of each of the edges must be in the sphere of influence
            if (link.fr_edge_src ∈ soi.nodes || link.fr_edge_tgt ∈ soi.nodes) &&
                (link.to_edge_src ∈ soi.nodes || link.to_edge_tgt ∈ soi.nodes) &&
                # That still doesn't mean it's in the sphere of influence - an end node could be but the actual point
                # where the link connects might still not be
                (
                    # we don't know whether the links face the same direction or not
                    (
                        # links face same direction: start near start, end near end
                        compute_net_distance(G, dmat,
                            link.fr_edge_src, link.fr_edge_tgt, link.fr_dist_from_start, link.fr_dist_to_end,
                            soi.original_link.fr_edge_src, soi.original_link.fr_edge_tgt, soi.original_link.fr_dist_from_start, soi.original_link.fr_dist_to_end) ≤ sphere_of_influence_radius &&
                        compute_net_distance(G, dmat,
                            link.to_edge_src, link.to_edge_tgt, link.to_dist_from_start, link.to_dist_to_end,
                            soi.original_link.to_edge_src, soi.original_link.to_edge_tgt, soi.original_link.to_dist_from_start, soi.original_link.to_dist_to_end) ≤ sphere_of_influence_radius
                    ) ||
                    (
                        # links face opposite directions: start near end, end near start
                        compute_net_distance(G, dmat,
                            link.fr_edge_src, link.fr_edge_tgt, link.fr_dist_from_start, link.fr_dist_to_end,
                            soi.original_link.to_edge_src, soi.original_link.to_edge_tgt, soi.original_link.to_dist_from_start, soi.original_link.to_dist_to_end) ≤ sphere_of_influence_radius &&
                        compute_net_distance(G, dmat,
                            link.to_edge_src, link.to_edge_tgt, link.to_dist_from_start, link.to_dist_to_end,
                            soi.original_link.fr_edge_src, soi.original_link.fr_edge_tgt, soi.original_link.fr_dist_from_start, soi.original_link.fr_dist_to_end) ≤ sphere_of_influence_radius
                    )
                )
                # note: order matters here. it is possible for a link to be in two spheres of influence,
                # so ordering might change results.
                in_soi = true
                push!(soi.links, link)

                # don't keep considering spheres of influence
                break
            end
        end

        if !in_soi
            # not in a sphere of influence, create one - everything near either end
            nodes = unique(vcat(
                findall(dmat[:, code_for(G, link.fr_edge_src)] .< sphere_of_influence_radius),
                findall(dmat[:, code_for(G, link.fr_edge_tgt)] .< sphere_of_influence_radius),
                findall(dmat[:, code_for(G, link.to_edge_src)] .< sphere_of_influence_radius),
                findall(dmat[:, code_for(G, link.to_edge_tgt)] .< sphere_of_influence_radius)
            ))

            push!(spheres_of_influence, SphereOfInfluence(label_for.(Ref(G), nodes), link, [link]))
        end
    end

    return map(spheres_of_influence) do soi
        best, rest = Iterators.peel(soi.links)

        for link in rest
            if link.geographic_length_m < best.geographic_length_m
                best = link
            end
        end

        return best
    end
end