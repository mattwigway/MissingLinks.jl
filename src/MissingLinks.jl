module MissingLinks
import MetaGraphsNext: MetaGraph, labels, edge_labels, code_for, label_for, neighbor_labels
import Graphs: Graph, dijkstra_shortest_paths, nv, ne, is_directed, connected_components, strongly_connected_components,
    outneighbors, has_edge, has_vertex, vertices, edges, degree, rem_edge!, rem_vertex!
import LibSpatialIndex: RTree, insert!, intersects
import DataFrames: DataFrame, nrow, metadata, metadata!
import LinearAlgebra: norm2
import Logging: @info, @warn, @error
import EzXML: XMLDocument, ElementNode, TextNode, link!, setroot!, prettyprint
import Compat: @compat
import Artifacts: @artifact_str
import GeoInterface, ArchGDAL, ThreadsX, LibGEOS, Graphs, GeoDataFrames, GDAL
import DataStructures: DefaultDict
import CSV, Dates
import OpenStreetMapPBF: Way, scan_nodes, scan_ways
import GeoFormatTypes as GFT

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
include("realize_graph.jl")
include("simplify_graph.jl")
include("tntp.jl")
include("osm.jl")


export graph_from_gdal, identify_potential_missing_links, remove_tiny_islands, deduplicate_links,
    score_links, create_graph_weights, semi_to_fully_noded, add_short_edges!, index_graph_edges, service_area,
    graph_to_gis, graph_to_graphml, links_to_gis, find_dead_ends, find_disconnected_crossings, fill_distance_matrix!,
    nodes_to_gis, realize_graph, collapse_realized_graph!

@compat public remove_elevation!, get_example_data, write_tntp
end
