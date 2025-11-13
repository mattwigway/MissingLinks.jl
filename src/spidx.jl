# wrap LibSpatialIndex functions to work around https://github.com/JuliaGeo/LibSpatialIndex.jl/pull/36

insert!(sidx::RTree, idx, min, max) = GC.@preserve min max LibSpatialIndex.insert!(sidx, idx, min, max)
intersects(sidx::RTree, min, max) = GC.@preserve min max LibSpatialIndex.intersects(sidx, min, max)
knn(sidx::RTree, min, max, n) = GC.@preserve min max LibSpatialIndex.knn(sidx, min, max, n)