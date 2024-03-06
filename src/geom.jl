function gdal_to_geos(gdalgeom::ArchGDAL.IGeometry{ArchGDAL.wkbLineString})
    coords = map(0:(ArchGDAL.ngeom(gdalgeom) - 1)) do i
        [ArchGDAL.getx(gdalgeom, i), ArchGDAL.gety(gdalgeom, i)]
    end
    LibGEOS.LineString(coords)
end