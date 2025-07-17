"""
    split_link!(G, fr, to, dist)

Split the link from fr to to at dist, and return the new vertex number.
"""
function split_link!(G, fr, to, dist)
    @assert fr < to

    edge = G[fr, to]
    g1 = geom_between(edge.geom, 0, dist)
    g2 = geom_between(edge.geom, dist, edge.length_m)

    rem_edge!(G, code_for(G, fr), code_for(G, to))

    g1_len = ArchGDAL.geomlength(g1)

    # split edges will not be measured sequentially, but this avoids scanning the graph to
    # find the next unused split id
    split = VertexID(nv(G) + 1, :split)
    G[split] = (ArchGDAL.getx(g2, 0), ArchGDAL.gety(g2, 0))

    if fr < split
        G[fr, split] = EdgeData((
            g1_len,
            edge.link_type,
            g1
        ))
    else
        G[split, fr] = EdgeData((
            g1_len,
            edge.link_type,
            reverse_geom(g1)
        ))
    end

    if split < to
        G[split, to] = EdgeData((
            edge.length_m - g1_len,
            edge.link_type,
            g2
        ))
    else
        G[to, split] = EdgeData((
            edge.length_m - g1_len,
            edge.link_type,
            reverse_geom(g2)
        ))
    end

    return split
end

"""
    link_point!(G, pts, idx; tol=20, min_split_length=1, create=false)

Link points pts (Vector of ArchGDAL objects) into graph destructively by splitting edges.

`tol` represents the tolerance (in graph units) for linking, `min_split_length` represents the tolerance
for splitting (if the closest point is closer to the end of the edge than min_split_length, we just connect
it to the end), and if create is true we will create an orphaned edge to connect the point to if there is no
nearby edge. Since the algorithm only connects edges to edges, this allows missing links to the point to be
identified even if there is no nearby edge.
"""
function link_points!(G, pts; tol=20, min_split_length=1, create=false)
    @info "indexing graph edges"
    sidx, edges = index_graph_edges(G)

    return map(pts) do pt
        coords = [ArchGDAL.getx(pt, 0), ArchGDAL.gety(pt, 0)]
        candidates = map(eid -> edges[eid], intersects(sidx, coords .- tol, coords .+ tol))
        best_dist = Inf64
        best_candidate = nothing
        for candidate in candidates
            if !haskey(G, candidate...)
                # this edge has already been split, we will come across the split edges at the end
                continue
            end

            dist = ArchGDAL.distance(G[candidate...].geom, pt)
            if dist < best_dist
                best_candidate = candidate
                best_dist = dist
            end
        end

        if best_dist < tol
            # there is a link within tolerance
            dist = LibGEOS.project(gdal_to_geos(G[best_candidate...].geom), gdal_to_geos(pt))
            if dist < min_split_length
                return min(best_candidate...)
            elseif dist > G[best_candidate...].length_m - min_split_length
                return max(best_candidate...)
            else
                split = split_link!(G, min(best_candidate...), max(best_candidate...), dist)
                # update index
                e1 = tuple(sort([best_candidate[1], split])...)
                push!(edges, e1)
                env = ArchGDAL.envelope(G[e1...].geom)
                insert!(sidx, length(edges), [env.MinX, env.MinY], [env.MaxX, env.MaxY])
                
                e2 = tuple(sort([best_candidate[2], split])...)
                push!(edges, e2)
                env = ArchGDAL.envelope(G[e2...].geom)
                insert!(sidx, length(edges), [env.MinX, env.MinY], [env.MaxX, env.MaxY])

                return split
            end
        elseif create
            v1 = VertexID(nv(G) + 1, :island)
            v2 = VertexID(nv(G) + 2, :island)
            G[v1] = tuple(coords...)
            G[v2] = tuple((coords .+ 1e-6)...)
            G[v1, v2] = EdgeData((
                1e-6 * sqrt(2),
                "island",
                ArchGDAL.createlinestring([coords, coords .+ 1e-6])
            ))

            return v1
        else
            return missing
        end
    end
end
