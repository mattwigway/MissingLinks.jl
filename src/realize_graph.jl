mutable struct RealizedCandidateLink{T}
    link::CandidateLink{T}
    srcnode::Union{VertexID, Nothing}
    dstnode::Union{VertexID, Nothing}
end

"""
    index_candidate_links(G, links)

Return a NamedTuple with three members:
- links: Vector of RealizedCandidateLinks, one for each link in `links`
- src: Dict mapping (VertexID, VertexID) (i.e. an edge ID) to a vector of RealizedCandidateLinks
    starting from that edge
- dest: Similar dict mapping edge to links ending on that edge

`G` should be the original graph the links were built with (although since the indices are the same,
the realized graph should also work).

The RealizedCandidateLinks in src and dest are references to those in links, so modifying
a link in any of these will modify it in all of them.
"""
function index_candidate_links(G, links::Vector{CandidateLink{T}}) where T
    realized = [RealizedCandidateLink{T}(link, nothing, nothing) for link ∈ links]

    src = DefaultDict{NTuple{2, VertexID}, Vector{RealizedCandidateLink{T}}}(Vector{RealizedCandidateLink})
    dest = DefaultDict{NTuple{2, VertexID}, Vector{RealizedCandidateLink{T}}}(Vector{RealizedCandidateLink})

    for link in realized
        v1 = label_for(G, link.link.fr_edge_src)
        v2 = label_for(G, link.link.fr_edge_tgt)
        push!(src[(v1, v2)], link)

        v3 = label_for(G, link.link.to_edge_src)
        v4 = label_for(G, link.link.to_edge_tgt)
        push!(dest[(v3, v4)], link)
    end

    return (realized=realized, src=src, dest=dest)
end

"""
    realize_graph(G, links)

"Realize" the graph `G` and candidate links `links`.

In identifying, deduplicating, and scoring missing links, we never actually put those
links into the graph. Instead, each link stores information about where it connects to
the network, and scoring is based only on the distance matrix, not on the graph itself.

However, for some applications, we need to "realize" the graph—i.e. create a graph that
has edges for all the candidate links, labeled as such. This function does that.
"""
function realize_graph(G, links::Vector{CandidateLink{T}}) where T
    # Create the new graph
    G2 = new_graph()

    # Pass 1: store all current vertices
    for v in vertices(G)
        G2[label_for(G, v)] = G[label_for(G, v)]
    end

    # Pass 2: index candidate links by edge
    realized_links = index_candidate_links(G, links)

    # Step 3: loop over edges
    for (v1, v2) in edge_labels(G)
        # locations where we will break this edge
        breakpoints = T[]

        # links that start on this edge
        for link in realized_links.src[(v1, v2)]
            push!(breakpoints, link.link.fr_dist_from_start)
        end

        # links that end on this edge
        for link in realized_links.dest[(v1, v2)]
            push!(breakpoints, link.link.to_dist_from_start)
        end

        unique!(breakpoints)
        sort!(breakpoints)

        # create nodes
        nodes = map(breakpoints) do breakpoint
            if breakpoint == zero(T)
                v1
            elseif breakpoint == round(T, G[v1, v2].length_m)
                v2
            else
                new_id = VertexID(nv(G2) + 1)
                @assert !haskey(G, new_id)

                # figure out where it goes
                new_coord = LibGEOS.interpolate(gdal_to_geos(G[v1, v2].geom), breakpoint)

                G2[new_id] = (LibGEOS.getGeomX(new_coord), LibGEOS.getGeomY(new_coord))

                new_id
            end
        end

        nodes_by_breakpoint = Dict(breakpoints .=> nodes)

        # mark links
        for link in realized_links.src[(v1, v2)]
            link.srcnode = nodes_by_breakpoint[link.link.fr_dist_from_start]
        end

        for link in realized_links.dest[(v1, v2)]
            link.dstnode = nodes_by_breakpoint[link.link.to_dist_from_start]
        end

        # add original links to graph
        if isempty(nodes) || first(nodes) ≠ v1
            pushfirst!(nodes, v1)
            pushfirst!(breakpoints, zero(T))
        end

        if last(nodes) ≠ v2
            push!(nodes, v2)
            push!(breakpoints, round(T, G[v1, v2].length_m))
        end

        node_and_dist = collect(zip(nodes, breakpoints))

        for ((n1, d1), (n2, d2)) ∈ zip(node_and_dist[begin:end-1], node_and_dist[begin+1:end])
            # if it's the end, use the actual edge length for end
            d2f = n2 == v2 ? G[v1, v2].length_m : convert(Float64, d2)

            G2[n1, n2] = EdgeData((
                d2f - d1,
                G[v1, v2].link_type,
                geom_between(G[v1, v2].geom, d1, d2f)
            ))
        end
    end

    # add candidate links to graph
    for link in realized_links.realized
        @assert !isnothing(link.srcnode)
        @assert !isnothing(link.dstnode)
        @assert !haskey(G2, (link.srcnode, link.dstnode))

        G2[link.srcnode, link.dstnode] = EdgeData((
            link.link.geographic_length_m,
            "candidate",
            ArchGDAL.createlinestring([G2[link.srcnode], G2[link.dstnode]])
        ))
    end

    return G2
end