# This file contains code to score the missing links

"""
Score a potential missing link
"""
function compute_link_score(link::CandidateLink, pool, origin_weights, dest_weights, decay_function, decay_cutoff_m)
    new_access::Float64 = 0.0

    with(pool) do dmat
        # distances _to_ each end of the from link
        fr_distances = distances_to(dmat, link.fr_edge_src, link.fr_edge_tgt)
        # distances _from_ each end of the target link
        to_distances = distances_from(dmat, link.to_edge_src, link.to_edge_tgt)

        # touching the database to get the old distance of every O-D pair is expensive
        # instead, touch the database once per origin to get distances to every destination, and store them here
        # we put this in a vector so we can write and read without allocations and indirection.
        # This is effectively one row of the non-sparse distance matrix. Even in very large networks this should not
        # take much memory (e.g. with 1 million vertices, this is 8mb).
        distances = Vector{Union{Int64, Missing}}(undef, nv(dmat.graph))

        for (origin, fr_dist_start, fr_dist_end) in fr_distances
            # how far is this origin from the start of this link?
            origin_distance = minimum_nonmissing(fr_dist_start + link.fr_dist_from_start, fr_dist_end + link.fr_dist_to_end)

            # not checking for missing... it should not be missing (one of the distances should be nonmissing or it would not be
            # returned). Let it error below if it is as that's a bug.
            if origin_distance + link.geographic_length_m ≤ decay_cutoff_m
                # It is close enough it could possibly provide access to something (distance to start of link
                # plus length of link itself)

                # cache distances through existing graph to everything nearby
                fill!(distances, missing)
                for (dest, dist) in distances_from(dmat, origin)
                    distances[code_for(dmat.graph, dest)] = dist
                end

                for (dest, to_dist_start, to_dist_end) in to_distances
                    dest_distance = minimum_nonmissing(to_dist_start + link.to_dist_from_start, to_dist_end + link.to_dist_to_end)

                    if origin_distance + link.geographic_length_m + dest_distance ≤ decay_cutoff_m
                        new_dist = origin_distance + link.geographic_length_m + dest_distance
                        old_dist = distances[code_for(dmat.graph, dest)]
                        if (ismissing(old_dist) || new_dist < old_dist)
                            # only affects access if it makes the trip shorter
                            Δdecayed = decay_function(new_dist)
                            if !ismissing(old_dist) # if it was unreachable, weight is zero
                                Δdecayed -= decay_function(old_dist)
                            end
                            new_access += Δdecayed * origin_weights[code_for(dmat.graph, origin)] * dest_weights[code_for(dmat.graph, dest)]
                        end
                    end
                end
            end
        end
    end

    return new_access
end

"""
    score_links(decay_function, links, distance_matrix, origin_weights, dest_weights, decay_cutoff_m)

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
function score_links(decay_function, links, dmat, origin_weights, dest_weights, decay_cutoff_m)
    # First, score all origins using the base network
    @info "Processing links"

    p = pool(dmat)

    # threading doesn't seem to make this faster in my toy example, and Julia peaks at 190% CPU usage.
    # Probably SQLite is in serial mode.
    link_scores = ThreadsX.mapi(enumerate(links)) do (i, link)
        if i % 1 == 0
            @info "Scored $i / $(length(links)) links"
        end
        # compute in both directions
        compute_link_score(link, p, origin_weights, dest_weights, decay_function, decay_cutoff_m) +
            compute_link_score(reverse(link), p, origin_weights, dest_weights, decay_function, decay_cutoff_m)
    end

    return link_scores
end