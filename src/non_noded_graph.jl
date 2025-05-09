"""
    semi_to_fully_noded(data...; snap_tolerance=1e-6, split_tolerance=1e-6)

Convert GDAL data sources where the sources are "semi-noded" - that is,
coincident ends, or an end coincident with a line segment, count as intersections - into a fully
noded graph with nodes at each intersection, ready to be built into a graph.

Data is a vararg, so arbitrarily many layers may be supplied. They will be merged into a
single output layer.

The snap tolerance is how close an end has to be to a line to be considered connected. The split
tolerance is how close two splits must be to each other to be considered the same split,
and should generally be much smaller, to avoid a situation where two ends snap to a line at
points near each other, the line is split at one of them, and the split is more than `snap_tolerance`
distance from the other one. I further recommend setting the tolerance when building the graph to
`snap_tolerance` + `split_tolerance`. These should both be quite small; larger gaps in the network
should be addressed with [`add_short_edges!`](@ref add_short_edges!).
"""
function semi_to_fully_noded(data...; snap_tolerance=1e-6, split_tolerance=1e-6)
    # flatten all the data files
    geoms_and_types = map(data) do file
        geomcol = if haskey(metadata(file), "geometrycolumns")
            first(metadata(file, "geometrycolumns"))
        elseif "geom" ∈ names(file)
            :geom
        elseif "geometry" ∈ names(file)
            :geometry
        else
            error("No geometry column detected. Columns: $(join(names(file), ", "))")
        end

        link_types = "link_type" ∈ names(file) ? file.link_type : fill(missing, nrow(file))
        
        map(zip(link_types, file[!, geomcol])) do (link_type, geom)
            result = Tuple{Union{String, Missing}, ArchGDAL.IGeometry{ArchGDAL.wkbLineString}}[]
            MissingLinks.for_each_geom(geom) do part
                push!(result, (link_type, part))
            end
            return result
        end
    end |>
        Iterators.flatten |>
        Iterators.flatten |>
        collect

    geoms = [g[2] for g in geoms_and_types]
    types = [g[1] for g in geoms_and_types]

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
    geos_geoms = map(gdal_to_geos, geoms)

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
                # TODO skip snaps to the same link as the start and end point are from
                # This should not change outputs as they'll be snapped to the ends already
                # We aren't currently handling p-shaped streets
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
        # we manually insert "breaks" at the start and end to ensure we don't lose the ends of lines
        # though this should be a no-op at the moment because the ends of lines get snapped to themselves
        geom_breaks[i] = [zero(Float64), LibGEOS.geomLength(geom)]
    end

    for connection in connections
        push!(geom_breaks[connection.to], connection.to_pos)
    end

    result = Vector{@NamedTuple{geom::ArchGDAL.IGeometry{ArchGDAL.wkbLineString}, link_type::Union{String, Missing}}}()

    for (i, link_type, geom) in zip(1:length(geoms), types, geos_geoms)
        breaks = sort(geom_breaks[i])
        prevbreak, restbreaks = Iterators.peel(breaks)
        for brk in restbreaks
            if brk - prevbreak < split_tolerance
                continue # don't make a new split
                # TODO this makes some lines slightly shorter if the end of the line is snapped back to an earlier break
                # This is probably actually desirable though, as it will clean up overshoots automatically, and should not affect connectivity
                # if add_short_edges is used with a tolerance of split_tolerance + snap_tolerance
            end

            push!(result, (geom=geom_between(geom, prevbreak, brk), link_type=link_type))
            prevbreak = brk
        end
    end

    return DataFrame(result)
end