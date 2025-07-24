function gdal_to_geos(gdalgeom::ArchGDAL.IGeometry{ArchGDAL.wkbLineString})
    coords = map(0:(ArchGDAL.ngeom(gdalgeom) - 1)) do i
        [ArchGDAL.getx(gdalgeom, i), ArchGDAL.gety(gdalgeom, i)]
    end
    LibGEOS.LineString(coords)
end

function gdal_to_geos(gdalgeom::ArchGDAL.IGeometry{ArchGDAL.wkbPoint})
    coords = [ArchGDAL.getx(gdalgeom, 0), ArchGDAL.gety(gdalgeom, 0)]
    LibGEOS.Point(coords)
end

geom_between(geom::ArchGDAL.IGeometry{ArchGDAL.wkbLineString}, pos1, pos2) = geom_between(gdal_to_geos(geom), pos1, pos2)

function geom_between(geom::LibGEOS.LineString, pos1, pos2)
    pos1 < pos2 || error("pos1 ($pos1) must be less than pos2 ($pos2), total length $(LibGEOS.geomLength(geom))")
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

reverse_geom(g::ArchGDAL.IGeometry{ArchGDAL.wkbLineString}) =
    ArchGDAL.createlinestring([ArchGDAL.getpoint(g, i) for i in (ArchGDAL.ngeom(g) - 1):-1:0])

"Run the provided function on each constituent geometry of a multigeometry"
function for_each_geom(f, g::Union{
        ArchGDAL.IGeometry{ArchGDAL.wkbMultiLineString},
        ArchGDAL.IGeometry{ArchGDAL.wkbMultiLineString25D},
        ArchGDAL.IGeometry{ArchGDAL.wkbMultiPoint},
        ArchGDAL.IGeometry{ArchGDAL.wkbMultiPoint25D},
        ArchGDAL.IGeometry{ArchGDAL.wkbGeometryCollection},
        ArchGDAL.IGeometry{ArchGDAL.wkbGeometryCollection25D}
    })
    for i in 1:ArchGDAL.ngeom(g)
        # archgdal has zero based indexing
        f(ArchGDAL.getgeom(g, i - 1))
    end
end

for_each_geom(f, g::Union{ArchGDAL.IGeometry{ArchGDAL.wkbLineString}, ArchGDAL.IGeometry{ArchGDAL.wkbLineString25D}, ArchGDAL.IGeometry{ArchGDAL.wkbPoint}, ArchGDAL.IGeometry{ArchGDAL.wkbPoint25D}}) = f(g)

"get the first point of a line"
get_first_point(g::Union{ArchGDAL.IGeometry{ArchGDAL.wkbLineString}, ArchGDAL.IGeometry{ArchGDAL.wkbLineString25D}}) = ArchGDAL.getpoint(g, 0)
# 0 based indexing in GDAL
get_last_point(g::Union{ArchGDAL.IGeometry{ArchGDAL.wkbLineString}, ArchGDAL.IGeometry{ArchGDAL.wkbLineString25D}}) = ArchGDAL.getpoint(g, ArchGDAL.ngeom(g) - 1)


"Get the geometry column from a GeoDataFrame"
function get_geometry(gdf)
    if haskey(metadata(gdf), "geometrycolumns")
        gdf[!, first(metadata(gdf, "geometrycolumns"))]
    elseif "geom" ∈ names(gdf)
        gdf.geom
    elseif "geometry" ∈ names(gdf)
        gdf.geometry
    else
        error("Could not autodetect geometry column; geometrycolumns metadata not specified and neither :geom nor :geometry columns present.")
    end
end

"""
    remove_elevation!(geom)

Remove the Z component of a geometry, and change the ArchGDAL wrapper type to match.

This does mutate its argument, but in order to re-wrap in a new ArchGDAL type a new wrapper must be returned. So you should
use the return value of the function.

Workaround for https://github.com/yeesian/ArchGDAL.jl/issues/333
"""
function remove_elevation!(geom)
    # keep geom in scope so it doesn't get destroyed before we're ready
    # atm this shouldn't happen as flattento2d! returns its argument, but if #333 is fixed
    # we don't want to rely on implementation details.
    new_geom = ArchGDAL.flattento2d!(geom)
    return ArchGDAL.IGeometry(GDAL.ogr_g_clone(new_geom.ptr))
end

"""
    get_xy(geom)

Convert a geometry into a Vector{Vector{Float64}} with the xy coordinates (used in testing).
"""
function get_xy(geom)
    map(0:(ArchGDAL.ngeom(geom) - 1)) do i
        [ArchGDAL.getx(geom, i), ArchGDAL.gety(geom, i)]
    end
end

geomapprox(a::ArchGDAL.IGeometry, b::ArchGDAL.IGeometry; kwargs...) = all(isapprox.(get_xy(a), get_xy(b); kwargs...))