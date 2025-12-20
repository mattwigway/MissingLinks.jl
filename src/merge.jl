
Base.contains(p::GraphPartition, x::Real, y::Real) = x ≥ p.min_x && x ≤ p.max_x && y ≥ p.min_y && y ≤ p.max_y
Base.contains(p::GraphPartition, pt::LibGEOS.Point) = Base.contains(p, LibGEOS.getGeomX(pt), LibGEOS.getGeomY(pt))

"""
    contains(p, l)

Is link l in the non-buffer area of graph partition p? Link l is assumed to be identified from this graph partition.
"""
function Base.contains(p::GraphPartition, l::CandidateLink)
    frgeom = gdal_to_geos(p.G[l.fr_edge_src, l.fr_edge_tgt].geom)
    frpt = LibGEOS.interpolate(frgeom, l.fr_dist_from_start)
    if Base.contains(p, frpt)
        # it only has to contain one point, short circuit if it does
        return true
    end

    togeom = gdal_to_geos(p.G[l.to_edge_src, l.to_edge_tgt].geom)
    topt = LibGEOS.interpolate(togeom, l.to_dist_from_start)
    return Base.contains(p, topt)
end

"""
    merge_links(G, Gs, links, scores)

Merge links from partitioned graphs Gs back into a single list referencing graph G. Return
a tuple of (links, scores).

"""
function merge_links(Gs::Matrix{GraphPartition{T}}, links::Matrix{Vector{CandidateLink}}, scores::Matrix{Vector{S}}) where {S <: Real, T}
    size(Gs) == size(links) || error("links must have same size as Gs")
    size(Gs) == size(scores) || error("links must have same size as Gs")

    # map from (edge, edge) to (candidatelink, score)
    out_links = Dict{NTuple{2, NTuple{2, VertexID}}, @NamedTuple{link::CandidateLink, score::S}}()

    for r in 1:size(links, 1), c in 1:size(links, 2)
        Gsub = Gs[r, c]
        lsub = links[r, c]
        ssub = scores[r, c]

        for (link, score) in zip(lsub, ssub)
            # figure out if this link is actually in this partition as opposed to in the buffer around it
            if contains(Gsub, link)
                fr = order_vertices(link.fr_edge_src, link.fr_edge_tgt)
                to = order_vertices(link.to_edge_src, link.to_edge_tgt)

                # translate link to original graph vertex codes
                link_tr = CandidateLink(
                    fr[1],
                    fr[2],
                    link.fr_dist_from_start,
                    link.fr_dist_to_end,
                    to[1],
                    to[2],
                    link.to_dist_from_start,
                    link.to_dist_to_end,
                    link.geographic_length_m,
                    # network distances are not reliable in partitioned graph, as the old network distance may have been
                    # much larger than a partition
                    zero(typeof(link.network_length_m))
                )

                # treat links in both directions as substitutes, and don't retain both,
                # but keep them coded in the right order so the distances match up.
                # the order within fr and to is constant thanks to order_vertices above.

                # should only have one version
                @assert !haskey(out_links, (fr, to)) || !haskey(out_links, (to, fr))

                # Links that cross a graph partition will be found in both partitions (because
                # only one end has to be in the partition for contains() to return true). But
                # they should be identical.
                if haskey(out_links, (fr, to))
                    o = out_links[(fr, to)].link

                    # they should be identical
                    @assert o.fr_dist_from_start == link_tr.fr_dist_from_start
                    @assert o.fr_dist_to_end == link_tr.fr_dist_to_end
                    @assert o.to_dist_from_start == link_tr.to_dist_from_start
                    @assert o.to_dist_to_end == link_tr.to_dist_to_end
                    @assert o.geographic_length_m == link_tr.geographic_length_m
                    @assert o.network_length_m == link_tr.network_length_m
                    @assert out_links[(fr, to)].score ≈ score
                elseif haskey(out_links, (to, fr))
                    o = out_links[(to, fr)].link

                    # they should be reverses of one another
                    @assert o.fr_dist_from_start == link_tr.to_dist_from_start
                    @assert o.fr_dist_to_end == link_tr.to_dist_to_end
                    @assert o.to_dist_from_start == link_tr.fr_dist_from_start
                    @assert o.to_dist_to_end == link_tr.fr_dist_to_end
                    @assert o.geographic_length_m == link_tr.geographic_length_m
                    @assert o.network_length_m == link_tr.network_length_m
                    @assert out_links[(to, fr)].score ≈ score
                else
                    out_links[(fr, to)] = (link=link_tr, score=score)
                end
            end
        end
    end

    merged_links = collect(values(out_links))

    return (
        collect(map(x -> x.link, merged_links)),
        collect(map(x -> x.score, merged_links))
    )
end
