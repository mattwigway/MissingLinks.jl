function gdal_to_geos(gdalgeom::ArchGDAL.IGeometry{ArchGDAL.wkbLineString})
    coords = map(0:(ArchGDAL.ngeom(gdalgeom) - 1)) do i
        [ArchGDAL.getx(gdalgeom, i), ArchGDAL.gety(gdalgeom, i)]
    end
    LibGEOS.LineString(coords)
end

geom_between(geom::ArchGDAL.IGeometry{ArchGDAL.wkbLineString}, pos1, pos2) = geom_between(gdal_to_geos(geom), pos1, pos2)

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
    if !(newcoord ≈ coords[end]) || length(coords) == 1
        push!(coords, newcoord)
    end

    return ArchGDAL.createlinestring(coords)
end

"""
Break long lines into pieces of at most `length`
"""
function break_long_line(line::ArchGDAL.IGeometry{ArchGDAL.wkbLineString}, max_length)::Vector{ArchGDAL.IGeometry{ArchGDAL.wkbLineString}}
    gline = GeoInterface.convert(LibGEOS, line)

    len = LibGEOS.geomLength(gline)
    if len ≤ max_length
        return [line]
    else
        geoms = ArchGDAL.IGeometry{ArchGDAL.wkbLineString}[]
        n_segs = len ÷ max_length + 1
        for seg in 1:n_segs
            start = (seg - 1) / n_segs * len
            endd = seg / n_segs * len
            push!(geoms, GeoInterface.convert(ArchGDAL, geom_between(gline, start, endd)))
        end
        return geoms
    end
end