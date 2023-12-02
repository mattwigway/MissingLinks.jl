# This file contains code to score the missing links

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
#
# This algorithm is O(nm^2), where n is the number of links, and m is the average number of nodes within decay_cutoff_m
# of each link

function compute_link_score(from, to, link_dist_m, dmat, origin_weights, dest_weights, decay_function, decay_cutoff_m)
    origin_distances = @view dmat[:, from]
    dest_distances = @view dmat[:, to]

    new_access = 0.0
    for origin in eachindex(origin_distances)
        if origin_distances[origin] ≤ (decay_cutoff_m - link_dist_m)
            for dest in eachindex(dest_distances)
                if dest_distances[dest] ≤ (decay_cutoff_m - link_dist_m)
                    new_dist = origin_distances[origin] + link_dist_m + dest_distances[dest]
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
        compute_link_score(link.fr, link.to, link.geographic_length_m, dmat, origin_weights, dest_weights, decay_function, decay_cutoff_m) +
            compute_link_score(link.to, link.fr, link.geographic_length_m, dmat, origin_weights, dest_weights, decay_function, decay_cutoff_m)
    end

    return link_scores
end