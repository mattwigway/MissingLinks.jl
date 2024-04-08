# This file contains code to score the missing links

"""
Score potential missing links

For each link, we identify how much making this connection would increase the aggregate number of
origin-weighted destinations (i.e. the sum across origins of the number of destinations accessible
within a threshold distance). We do this in several steps:
 1. Identify all origins and destinations within the threshold of the start of the proposed new link
 2. For each origin and destination:
    a. generate distances from the origin to the destination via the new link by summing
        - Distance from origin to link
        - Length of link
        - Distance from link to each destination
    b. if the distance via the new link is longer than the existing distance, continue to the next pair
    c. otherwise, compute the distance weight using the weighting function in the accessibility metric for the distance via the new link and the existing distance
    d. take the difference of these weights
    e. multiply by the origin and destination weights
 3. Sum the results
 4. Reverse the start and end of the link, and repeat
 5. Sum the results
"""
function compute_link_score(link::CandidateLink, dmat, origin_weights, dest_weights, decay_function, decay_cutoff_m)
    # Performance optimization
    origin_fr_distances = @view dmat[:, link.fr_edge_src]
    origin_to_distances = @view dmat[:, link.fr_edge_tgt]
    dest_fr_distances = @view dmat[link.to_edge_src, :]
    dest_to_distances = @view dmat[link.to_edge_tgt, :]

    new_access = 0.0
    # could flip this loop so we loop over destinations first, and thus access the GBMatrix in row-order
    for origin in eachindex(origin_fr_distances)
        origin_distance = min(
            add_unless_typemax(origin_fr_distances[origin], link.fr_dist_from_start),
            add_unless_typemax(origin_to_distances[origin], link.fr_dist_to_end)
        )
        if origin_distance ≤ (decay_cutoff_m - link.geographic_length_m)
            for dest in eachindex(dest_fr_distances)
                dest_distance = min(
                    add_unless_typemax(dest_fr_distances[dest], link.to_dist_from_start),
                    add_unless_typemax(dest_to_distances[dest], link.to_dist_to_end)
                )

                if dest_distance ≤ (decay_cutoff_m - link.geographic_length_m)
                    new_dist = origin_distance + link.geographic_length_m + dest_distance
                    # For construction performance, the array is i
                    old_dist = dmat[dest, origin]
                    if (new_dist < old_dist)
                        # only affects access if it makes the trip shorter
                        Δdecayed = decay_function(new_dist) - decay_function(old_dist)
                        new_access += Δdecayed * origin_weights[origin] * dest_weights[dest]
                    end
                end
            end
        end
    end

    return new_access
end

function score_links(decay_function, links, dmat, origin_weights, dest_weights, decay_cutoff_m)
    # First, score all origins using the base network
    @info "Processing links"
    link_scores = ThreadsX.mapi(enumerate(links)) do (i, link)
        if i % 100 == 0
            @info "Scored $i / $(length(links)) links"
        end
        # compute in both directions
        compute_link_score(link, dmat, origin_weights, dest_weights, decay_function, decay_cutoff_m) +
            compute_link_score(reverse(link), dmat, origin_weights, dest_weights, decay_function, decay_cutoff_m)
    end

    return link_scores
end