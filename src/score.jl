# This file contains code to score the missing links

"""
Score a potential missing link
"""
function compute_link_score(G, link::CandidateLink, dmat, origin_weights, dest_weights, decay_function, decay_cutoff_m)
    # note: assumes undirected graph
    origin_fr_distances = @view dmat[:, code_for(G, link.fr_edge_src)]
    origin_to_distances = @view dmat[:, code_for(G, link.fr_edge_tgt)]
    dest_fr_distances = @view dmat[:, code_for(G, link.to_edge_src)]
    dest_to_distances = @view dmat[:, code_for(G, link.to_edge_tgt)]

    new_access = 0.0
    for origin in eachindex(origin_fr_distances)
        if origin_weights[origin] == 0
            continue
        end

        # how far is this origin from the start of this link?
        origin_distance = min(
            add_unless_typemax(origin_fr_distances[origin], link.fr_dist_from_start),
            add_unless_typemax(origin_to_distances[origin], link.fr_dist_to_end)
        )

        # the subtraction should never overflow, because decay_cutoff will generally be more than the
        # maximum length of the link, but make sure
        if origin_distance ≤ Base.Checked.checked_sub(decay_cutoff_m, link.geographic_length_m)
            # It is close enough it could possibly provide access to something (distance to start of link
            # plus length of link itself)
            for dest in eachindex(dest_fr_distances)
                if dest_weights[dest] == 0
                    continue
                end

                dest_distance = min(
                    add_unless_typemax(dest_fr_distances[dest], link.to_dist_from_start),
                    add_unless_typemax(dest_to_distances[dest], link.to_dist_to_end)
                )

                if dest_distance ≤ Base.Checked.checked_sub(decay_cutoff_m, link.geographic_length_m)
                    # also should not overflow as the 2x the decay cutoff + the link length should not be more than typemax,
                    # but check
                    new_dist = Base.Checked.checked_add(Base.Checked.checked_add(origin_distance, link.geographic_length_m), dest_distance)
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

"""
    score_links(decay_function, graph, links, distance_matrix, origin_weights, dest_weights, decay_cutoff_m)

Score the contribution of each link in `links` to aggregate accessibility, using the decay function and weights specified.

`decay_function` is a Julia function that takes a single argument, a distance, and returns a weight; smaller weights
mean less influence. For a two-mile cumulative opportunities function, this could just be e.g. `x -> x < 3218`. It is the
first argument to allow use of Julia do-block syntax.

`links` is the vector of links (e.g. from [`deduplicate_links`](@ref deduplicate_links)). The distance matrix is the same
one used throughout the analysis. `origin_weights` and `dest_weights` are vectors of weights associated with each node,
e.g. computed by [`create_graph_weight`](@ref create_graph_weights).

`decay_cutoff_m` is the distance at which `decay_function` goes to zero (for a smooth function, we recommend
reimplementing as piecewise with a dropoff to zero at some point where additional access is largely immaterial).
It is in meters if the input data were in meters (which we recommend as we have not tested with other units to
ensure units are not hardcoded anywhere. Just use meters. Be a scientist.)

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
function score_links(decay_function, G, links, dmat, origin_weights, dest_weights, decay_cutoff_m)
    # First, score all origins using the base network
    @info "Processing links"
    link_scores = ThreadsX.mapi(enumerate(links)) do (i, link)
        if i % 100 == 0
            @info "Scored $i / $(length(links)) links"
        end
        # compute in both directions
        compute_link_score(G, link, dmat, origin_weights, dest_weights, decay_function, decay_cutoff_m) +
            compute_link_score(G, reverse(link), dmat, origin_weights, dest_weights, decay_function, decay_cutoff_m)
    end

    return link_scores
end