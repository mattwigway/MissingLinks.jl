module MissingLinks
import MetaGraphsNext: MetaGraph, labels, edge_labels, code_for, label_for
import Graphs: Graph, dijkstra_shortest_paths, nv, is_directed, connected_components, strongly_connected_components, rem_vertex!,
    outneighbors
import LibSpatialIndex: RTree, insert!, intersects
import GeoInterface, ArchGDAL
import DataFrames: DataFrame, nrow
import LinearAlgebra: norm2
import Logging: @info, @warn, @error
import EzXML: XMLDocument, ElementNode, TextNode, link!, setroot!, prettyprint
import ThreadsX
import LibGEOS

include("graph.jl")
include("node_to_node.jl")
include("identify_missing_links.jl")
include("score.jl")
include("weight_nodes.jl")
include("non_noded_graph.jl")

export graph_from_gdal, identify_potential_missing_links, links_to_gdf, remove_tiny_islands, deduplicate_links,
    score_links, create_graph_weights, semi_to_fully_noded

end
