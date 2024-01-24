"""
This converts GDAL data sources where the sources are "semi-noded" - that is,
coincident ends, or an end coincident with a line segment, count as intersections - into a fully
noded graph with nodes at each intersection, ready to be built into a graph.

The snap tolerance is how close an end has to be to a line to be considered connected. The split
tolerance is how close two splits must be to each other to be considered the same split,
and should generally be much smaller, to avoid a situation where two ends snap to a line at
points near each other, the line is split at one of them, and the split is more than snap_tolerance
distance from the other one. I further recommend setting the tolerance when building the graph to
snap_tolerance + split_tolerance.
"""
function semi_to_fully_noded(data...; snap_tolerance=1e-6, split_tolerance=1e-6)
    # flatten all the data files
    geoms = map(data) do file
        map(file.geom) do geom
            result = ArchGDAL.IGeometry{ArchGDAL.wkbLineString}[]
            MissingLinks.for_each_geom(geom) do part
                push!(result, part)
            end
            return result
        end
    end |>
        Iterators.flatten |>
        Iterators.flatten |>
        collect
    @info "$(length(geoms)) geometry records"
    # next, we build a spatial index
    index = RTree(2)

    @info "Building spatial index"
    for (i, geom) in enumerate(geoms)
        if i % 10000 == 0
            @info "$i / $(length(geoms))"
        end
        envl = ArchGDAL.envelope(geom)
        insert!(index, i, [envl.MinX, envl.MinY], [envl.MaxX, envl.MaxY])
    end

    @info "Linking ends to middles"
    # okay - new approach. We record which end is closest, and the offset for the closest point on the line
    # then we scan over that and quantize those offsets
    # Then we create nodes for each geom and offset, and link them up
    geos_geoms = map(geoms) do gdalgeom
        coords = map(0:(ArchGDAL.ngeom(gdalgeom) - 1)) do i
            [ArchGDAL.getx(gdalgeom, i), ArchGDAL.gety(gdalgeom, i)]
        end
        LibGEOS.LineString(coords)
    end

    connections = []

    for (i, geom) in enumerate(geoms)
        if i % 10000 == 0
            @info "$i / $(length(geoms))"
        end

        # check both ends
        startpt = ArchGDAL.getpoint(geom, 0)
        endpt = ArchGDAL.getpoint(geom, ArchGDAL.ngeom(geom) - 1)

        for (pt, which) in zip([startpt, endpt], [:start, :end])
            x, y, _ = pt
            xy = [x, y]
            xygeom = ArchGDAL.createpoint(xy)
            candidates = intersects(index, xy .- snap_tolerance, xy .+ snap_tolerance)
            for candidate in candidates
                if ArchGDAL.distance(geoms[candidate], xygeom) ≤ snap_tolerance
                    # It is within the distance, find the closest point on the candidate
                    geosgeom = geos_geoms[candidate]
                    closest_offset = LibGEOS.project(geosgeom, LibGEOS.Point(xy))
                    push!(connections, (from=i, to=candidate, from_end=which, to_pos=closest_offset))
                end
            end
        end
    end

    # calculate breaks for each geom
    geom_breaks = Dict()

    for (i, geom) in enumerate(geos_geoms)
        geom_breaks[i] = [zero(Float64), LibGEOS.geomLength(geom)]
    end

    for connection in connections
        push!(geom_breaks[connection.to], connection.to_pos)
    end

    result = ArchGDAL.IGeometry{ArchGDAL.wkbLineString}[]

    for (i, geom) in enumerate(geos_geoms)
        breaks = sort(geom_breaks[i])
        prevbreak, restbreaks = Iterators.peel(breaks)
        for brk in restbreaks
            if brk - prevbreak < split_tolerance
                continue # don't make a new split
                # TODO this makes some lines slightly shorter if the end of the line is snapped back to an earlier break
            end

            push!(result, geom_between(geom, prevbreak, brk))
            prevbreak = brk
        end
    end

    return DataFrame(:geom=>result)
end

function geom_between(geom::LibGEOS.LineString, pos1, pos2)
    pos1 < pos2 || error("pos1 must be less than pos2")
    startpt = LibGEOS.interpolate(geom, pos1)

    coords = [[LibGEOS.getGeomX(startpt), LibGEOS.getGeomY(startpt)]]

    cumulative_dist = zero(Float64)
    # skip first and last as they will be handled by startpt and endpt
    for ptidx in 2:(LibGEOS.numPoints(geom) - 1)
        lastpt = LibGEOS.getPoint(geom, ptidx - 1)
        thispt = LibGEOS.getPoint(geom, ptidx)
        cumulative_dist += LibGEOS.distance(lastpt, thispt)
        
        if cumulative_dist > pos2
            break
        end

        if cumulative_dist > pos1
            newcoord = [LibGEOS.getGeomX(thispt), LibGEOS.getGeomY(thispt)]
            if !(newcoord ≈ coords[end])
                push!(coords, newcoord)
            end
        end
    end

    endpt = LibGEOS.interpolate(geom, pos2)
    newcoord = [LibGEOS.getGeomX(endpt), LibGEOS.getGeomY(endpt)]
    if !(newcoord ≈ coords[end])
        push!(coords, newcoord)
    end

    return ArchGDAL.createlinestring(coords)
end