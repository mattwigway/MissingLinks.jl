module MissingLinks
import MetaGraphsNext: MetaGraph, labels, edge_labels, code_for, label_for
import Graphs: Graph, dijkstra_shortest_paths, nv, ne, is_directed, connected_components, strongly_connected_components, rem_vertex!,
    outneighbors, has_edge, has_vertex, vertices, edges
import LibSpatialIndex: RTree, insert!, intersects
import DataFrames: DataFrame, nrow, metadata, metadata!
import LinearAlgebra: norm2
import Logging: @info, @warn, @error
import EzXML: XMLDocument, ElementNode, TextNode, link!, setroot!, prettyprint
import Compat: @compat
import Artifacts: @artifact_str
import GeoInterface, ArchGDAL, ThreadsX, LibGEOS, Graphs, GeoDataFrames, GDAL


include("graph.jl")
include("candidate_link.jl")
include("dist.jl")
include("identify_missing_links.jl")
include("deduplicate_links.jl")
include("score.jl")
include("weight_nodes.jl")
include("non_noded_graph.jl")
include("geom.jl")
include("service_area.jl")
include("graph_troubleshooting.jl")
include("example_data.jl")


export graph_from_gdal, identify_potential_missing_links, links_to_gdf, remove_tiny_islands, deduplicate_links,
    score_links, create_graph_weights, semi_to_fully_noded, add_short_edges!, index_graph_edges, service_area,
    graph_to_gis, graph_to_graphml, links_to_gis, find_dead_ends, find_disconnected_crossings, fill_distance_matrix!,
    nodes_to_gis

@compat public remove_elevation!, get_example_data
end
