"""
    distance_surface(G, origin; maxdist=5000, resolution=25, crs=nothing)

Get a Raster.jl raster of the travel times from origin.
"""
function distance_surface(G, dmat, origin; maxdist=5000, resolution=25, index=nothing, crs=nothing, max_snap_dist=100)
    if isnothing(index)
        index = index_graph_edges(G)
    end

    if !isnothing(crs)
        origin = ArchGDAL.reproject(origin, GFT.EPSG(4326), crs)
    end

    origin_edge, dist_to_edge, dist_from_start, dist_from_end = closest_edge(G, collect(origin), index)

    dist_to_edge ≤ maxdist || error("No edges near origin")

    minx, miny = origin .- maxdist .- resolution
    maxx, maxy = origin .+ maxdist .+ resolution

    # y iterates negatively for north-up in most projected coordinate systems
    xdim, ydim = Rasters.X(minx:resolution:maxx), Rasters.Y(maxy:-resolution:miny)
    # rasters.jl doesn't write the UInt16 to geotiff correctly
    rast = Rasters.Raster(fill(convert(Float64, typemax(UInt16)), xdim, ydim))

    @showprogress for (xi, yi) in Iterators.product(1:size(rast, 1), 1:size(rast, 2))
        x = Rasters.lookup(xdim)[xi]
        y = Rasters.lookup(ydim)[yi]

        dest_edge, dest_dist, dest_fr_start, dest_to_end = closest_edge(G, [x, y], index)

        if dest_dist ≤ max_snap_dist
            rast[xi, yi] = round(UInt16,
            add_unless_typemax(
                compute_net_distance(
                    G, dmat,
                    origin_edge..., dist_from_start, dist_from_end,
                    dest_edge..., dest_fr_start, dest_to_end
                ),
                dist_to_edge + dest_dist
            ))
        end
    end

    return rast
end
