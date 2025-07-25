"""
    partition(G, nrow, ncol, overlap)

Partition graph G into nrow x ncol subgraphs. This is useful with graphs that are too large
to process at once. `buffer` is the amount the graphs should overlap, in meters/graph units
(though using units other than meters is untested and unscientific).

Specifically, the bounding box of each subgraph will be extended by `buffer` meters, so the actual
overlap will be 2 * `overlap`.

Returns an nrow x ncol matrix of graphs.
"""
function partition(G::T, nrow, ncol, buffer) where T
    result = Matrix{T}(undef, nrow, ncol)

    # Get the envelope of the entire graph
    bbox = extent(G)

    col_x = (bbox.max_x - bbox.min_x) / ncol
    row_y = (bbox.max_y - bbox.min_y) / nrow

    for r in 1:nrow, c in 1:ncol
        Gsub = new_graph()

        xmin = bbox.min_x + col_x * (c - 1)
        xmax = xmin + col_x
        ymin = bbox.min_y + row_y * (r - 1)
        ymax = ymin + row_y

        xmin -= buffer
        ymin -= buffer
        xmax += buffer
        ymax += buffer

        # copy over all vertices in area
        for v in labels(G)
            x, y = G[v]
            if x ≥ xmin && x ≤ xmax && y ≥ ymin && y ≤ ymax
                Gsub[v] = G[v]
            end
        end

        # copy all edges
        for v in labels(Gsub)
            nbrs = neighbor_labels(G, v)
            
            for nbr in nbrs
                if !haskey(Gsub, nbr)
                    # include edges that start within and go outside area
                    Gsub[nbr] = G[nbr]
                end

                v1, v2 = order_vertices(v, nbr)
                Gsub[v1, v2] = G[v1, v2]
            end
        end

        result[r, c] = Gsub
    end

    return result
end

"""
    extent(G)

Get the spatial extent of graph G. Returns a NamedTuple with members `min_x`, `max_x`, `min_y`, `max_y`.
The bounding box is based only on the vertices.
"""
function extent(G)
    min_x = Inf64
    max_x = -Inf64
    min_y = Inf64
    max_y = -Inf64

    for v in labels(G)
        x, y = G[v]
        min_x = min(x, min_x)
        max_x = max(x, max_x)
        min_y = min(y, min_y)
        max_y = max(y, max_y)
    end

    @assert all(isfinite.((min_x, max_x, min_y, max_y)))

    return (min_x=min_x, min_y=min_y, max_x=max_x, max_y=max_y)
end

"""
    merge_links(G, Gs, links, scores)

Merge links from partitioned graphs Gs back into a single list referencing graph G. Return
a tuple of (links, scores).

Note: currently this returns the maximum score for each link. That is usually but not always correct. Consider a portion
of a graph that looks like this:
  
    ✂
A --✂--- B
    ✂    |
    ✂    |
    ✂    |
D---✂--- C
    ✂

A is relatively directly connected to D, but suppose where the scissor are is the edge of a partition that extends
to the left. In the partition that extends to the left, AB and CD will be present, but BC will not as it is completely
outside the partition. This will result in a candidate link AB->CD being identified and may score highly as it closes
a major gap in the graph. This subgraph will be completely within the partition to the right, and we will get the correct
link and score there, but the score may be lower or AB->CD may not be identified at all.

I think what we need to do is only include links where at least one end of them is in the "main" part of the partition,
and then increase the buffer to maxdist + max_link_dist.

"""
function merge_links(G::T, Gs::Matrix{T}, links::Matrix{Vector{CandidateLink}}, scores::Matrix{Vector{S}}) where {S <: Real, T}
    size(Gs) == size(links) || error("links must have same size as Gs")
    size(Gs) == size(scores) || error("links must have same size as Gs")

    # map from (edge, edge) to (candidatelink, score)
    out_links = Dict{NTuple{2, NTuple{2, VertexID}}, @NamedTuple{link::CandidateLink, score::S}}()

    for r in 1:size(links, 1), c in 1:size(links, 2)
        Gsub = Gs[r, c]
        lsub = links[r, c]
        ssub = scores[r, c]

        for (link, score) in zip(lsub, ssub)
            fr = order_vertices(label_for(Gsub, link.fr_edge_src), label_for(Gsub, link.fr_edge_tgt))
            to = order_vertices(label_for(Gsub, link.to_edge_src), label_for(Gsub, link.to_edge_tgt))

            # translate link to original graph vertex codes

            # near an edge, the network distance of a link may be unreliable, but there will be another
            # graph where the distance is reliable because the edges overlap. use the minimum network distance.
            net_dist = link.network_length_m

            if haskey(out_links, (fr, to))
                net_dist = min(net_dist, out_links[(fr, to)].link.network_length_m)
            elseif haskey(out_links, (to, fr))
                net_dist = min(net_dist, out_links[(to, fr)].link.network_length_m)
            end

            link_tr = CandidateLink(
                code_for(G, fr[1]),
                code_for(G, fr[2]),
                link.fr_dist_from_start,
                link.fr_dist_to_end,
                code_for(G, to[1]),
                code_for(G, to[2]),
                link.to_dist_from_start,
                link.to_dist_to_end,
                link.geographic_length_m,
                net_dist
            )

            # treat links in both directions as substitutes, and don't retain both,
            # but keep them coded in the right order so the distances match up.
            # the order within fr and to is constant thanks to order_vertices above.

            # should only have one version
            @assert !haskey(out_links, (fr, to)) || !haskey(out_links, (to, fr))

            if haskey(out_links, (fr, to))
                o = out_links[(fr, to)].link

                # they should be identical
                @assert o.fr_dist_from_start == link_tr.fr_dist_from_start
                @assert o.fr_dist_to_end == link_tr.fr_dist_to_end
                @assert o.to_dist_from_start == link_tr.to_dist_from_start
                @assert o.to_dist_to_end == link_tr.to_dist_to_end
                @assert o.geographic_length_m == link_tr.geographic_length_m
                # network length might not be identical, see above

                if out_links[(fr, to)].score < score
                    out_links[(fr, to)] = (link=link_tr, score=score)
                end
            elseif haskey(out_links, (to, fr))
                o = out_links[(to, fr)].link

                # they should be reverses of one another
                @assert o.fr_dist_from_start == link_tr.to_dist_from_start
                @assert o.fr_dist_to_end == link_tr.to_dist_to_end
                @assert o.to_dist_from_start == link_tr.fr_dist_from_start
                @assert o.to_dist_to_end == link_tr.fr_dist_to_end
                @assert o.geographic_length_m == link_tr.geographic_length_m
                # network length might not be identical, see above

                if out_links[(to, fr)].score < score
                    out_links[(fr, to)] = (link=link_tr, score=score)
                    delete!(out_links, (to, fr))
                end
            else
                out_links[(fr, to)] = (link=link_tr, score=score)
            end
        end
    end

    merged_links = collect(values(out_links))

    return (
        collect(map(x -> x.link, merged_links)),
        collect(map(x -> x.score, merged_links))
    )
end
