module MissingLinks
import MetaGraphsNext: MetaGraph, labels, edge_labels
import Graphs: DiGraph, dijkstra_shortest_paths, nv
import LibSpatialIndex: RTree, insert!, intersects
import GeoInterface, ArchGDAL
import DataFrames: nrow
import LinearAlgebra: norm2
import Logging: @info, @warn, @error
import EzXML: XMLDocument, ElementNode, TextNode, link!, setroot!, prettyprint

include("graph.jl")

export graph_from_gdal

end
