"""
    collapse_realized_graph!(G)

Remove degree-2 nodes and combine the associated edges. _Please note_ that this graph
has not been tested with link identification code and unexpected things may happen because of
vertex removal changing the order of nodes.
"""
function collapse_realized_graph!(G)
    @warn "This function should only be applied to a realized graph for export, it should not be used prior to link identification."
    # collect so we're not looping over something tied to graph as we're modifying
    # the graph
    for vertex ∈ collect(labels(G))
        vidx = code_for(G, vertex)
        if degree(G, vidx) == 2
            # merge the edges
            nbr1, nbr2 = neighbor_labels(G, vertex)

            @assert nbr1 != vertex
            @assert nbr2 != vertex

            # TODO need to fix for directed graphs
            if code_for(G, nbr2) < code_for(G, nbr1)
                nbr1, nbr2 = nbr2, nbr1
            end

            e1 = G[nbr1, vertex]
            e2 = G[vertex, nbr2]

            if nbr1 == nbr2 || has_edge(G, code_for(G, nbr1), code_for(G, nbr2)) ||
                    e1.link_type == "candidate" ||
                    e2.link_type == "candidate"
                # this vertex is necessary to preserve topology
                continue
            end

            # the geometry should go nbr1 -> vertex -> nbr2
            # confirmed above that nbr1 < nbr2
            g1 = code_for(G, nbr1) < code_for(G, vertex) ? e1.geom : reverse_geom(e1.geom)
            g2 = code_for(G, vertex) < code_for(G, nbr2) ? e2.geom : reverse_geom(e2.geom)

            combined_geom = ArchGDAL.createlinestring([
                get_xy(g1)[begin:end-1]..., # remove the end point, replace with vertex geom
                collect(G[vertex]),
                # remove the start point, replaced with vertex geom
                get_xy(g2)[begin+1:end]...
            ])

            # when multiple edges get merged, we want sidewalk|crosswalk, not sidewalk|crosswalk|sidewalk|sidewalk
            etypes = Set([split(e1.link_type, "|")..., split(e2.link_type, "|")...])
            etype = join(etypes, "|")

            G[nbr1, nbr2] = EdgeData((
                e1.length_m + e2.length_m,
                etype,
                combined_geom
            ))

            rem_edge!(G, code_for(G, nbr1), vidx)
            rem_edge!(G, vidx, code_for(G, nbr2))
            rem_vertex!(G, vidx)
        end
    end
end