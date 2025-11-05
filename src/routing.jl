@kwdef struct OneToOneRouteResult
    distance::Float64
    geom::ArchGDAL.IGeometry{ArchGDAL.wkbLineString}
end

"Get the closest edge to a location (tuple/vector of coordinates)"
function closest_edge(G, loc, index)
    eidx, edges = index

    candidates = map(x -> edges[x], knn(eidx, loc, loc, 1))
    pt = ArchGDAL.createpoint(loc)
    dists = map(c -> ArchGDAL.distance(pt, G[c...].geom), candidates)
    best = argmin(dists)

    egeom = G[candidates[best]...].geom
    _, point_on_edge = LibGEOS.nearestPoints(loc, egeom)
    dist_from_start = LibGEOS.project(egeom, point_on_edge)
    dist_from_end = LibGEOS.geomLength(egeom) - dist_from_start

    return candidates[best], dists[best], dist_from_start, dist_from_end
end

"""
    route_one_to_one(G, origin, dest; crs=nothing, index=nothing)

Do a one-to-one route. If CRS is specified, origin and dest can be in lat-lon (not lon-lat) coordinates
and will be reprojected to the specified CRS before routing.

If you're doing many routes, setting `index` to the result of [index_graph_edges](@ref) will speed up results.
"""
function route_one_to_one(G, origin, dest; crs=nothing, index=nothing, max_snap_dist=100)
    if isnothing(index)
        index = index_graph_edges(G)
    end

    if !isnothing(crs)
        origin = ArchGDAL.reproject(origin, GFT.EPSG(4326), crs)
        dest = ArchGDAL.reproject(dest, GFT.EPSG(4326), crs)
    end

    origin_edge, odist, odist_from_start, odist_from_end = closest_edge(G, origin, index)
    dest_edge, ddist, ddist_from_start, ddist_from_end = closest_edge(G, dest, index)

    odist <= max_snap_dist || error("No edges near origin")
    ddist <= max_snap_dist || error("No edges near destination")

    # can only do dijkstra from nodes (multi-origin dijkstra theoretically possible,
    # but not implemented in Graphs.jl)
    paths_by_origin = [
        dijkstra_shortest_paths(G, [code_for(G, origin_edge[1])])
        dijkstra_shortest_paths(G, [code_for(G, origin_edge[2])])
    ]

    best_distance = typemax(Int64)
    best_ostartend = nothing
    best_dstartend = nothing
    for ostartend in 1:2, dstartend in 1:2
        paths = paths_by_origin[ostartend]
        dist = paths.dists[code_for(G, dest_edge[dstartend])]

        # add distance from origin, destination
        dist += ostartend == 1 ? odist_from_start : odist_from_end
        dist += dstartend == 1 ? ddist_from_start : ddist_from_end

        if dist < best_distance
            best_distance = dist
            best_ostartend = ostartend
            best_dstartend = dstartend
        end
    end
        
    # trace the path
    # what's your vector victor? ... and stop calling me shirley
    coords = Vector{Vector{Float64}}()

    # end partial edge
    if best_dstartend == 1 && ddist_from_start != 0
        # start of the last edge
        push!.(Ref(coords), reverse(get_xy(geom_between(G[dest_edge...].geom, 0, ddist_from_start)))[begin+1:end])
    elseif best_dstartend == 2 && ddist_from_end != 0
        push!.(Ref(coords), reverse(get_xy(geom_between(reverse_geom(G[dest_edge...].geom), 0, ddist_from_end))[begin+1:end]))
    end

    # main edges
    vx = dest_edge[best_dstartend]
    
    while true
        prev_vxidx = paths_by_origin[best_ostartend].parents[code_for(G, vx)]

        if prev_vxidx == 0
            break
        end

        prev_vx = label_for(G, prev_vxidx)

        geom = G[order_vertices(prev_vx, vx)...].geom
        # this seems backward, b/c we are walking the tree backwards; the whole coords array
        # gets reversed at the end
        if prev_vx < vx
            geom = reverse_geom(geom)
        end

        push!.(Ref(coords), get_xy(geom)[begin+1:end])

        vx = prev_vx
    end

    # begin partial edge
    if best_ostartend == 1 && odist_from_start != 0
        push!.(Ref(coords), get_xy(geom_between(G[origin_edge...].geom, 0, odist_from_start))[begin+1:end])
    elseif best_dstartend == 2 && odist_from_end != 0
        push!.(Ref(coords), get_xy(geom_between(reverse_geom(G[origin_edge...].geom), 0, odist_from_end))[begin+1:end])
    else
        # no partial, just the vertex
        push!(coords, G[vx...])
    end

    return OneToOneRouteResult(
        distance = best_distance,
        geom = ArchGDAL.createlinestring(collect(reverse(coords)))
    )
end