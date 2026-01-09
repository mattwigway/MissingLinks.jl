# API Documentation

## Graph building and troubleshooting

```@docs
MissingLinks.get_example_data
semi_to_fully_noded
graph_from_gdal
graph_to_graphml
graph_to_gis
nodes_to_gis
remove_tiny_islands
add_short_edges!
find_disconnected_crossings
find_dead_ends
MissingLinks.remove_elevation!
TraversalPermissionSettings
graph_from_osm
MissingLinks.is_traversable
MissingLinks.extent
```
## Opportunity data

```@docs
create_graph_weights
```

## Link identification, deduplication, and scoring

```@docs
fill_distance_matrix!
identify_potential_missing_links
deduplicate_links
score_links
links_to_gis
```

## Routing and mapping of impacts of individual links

```@docs
route_one_to_one
distance_surface
regional_access
```

## Graph partitioning

Large graphs are intractable; these tools can be used to split up a graph and analyze it.

```@docs
MissingLinks.GraphPartition
MissingLinks.partition
MissingLinks.identify_and_score
MissingLinks.partition_weights
MissingLinks.merge_links
Base.contains
```

## Internal

```@docs
MissingLinks.new_graph
MissingLinks.for_each_geom
MissingLinks.get_first_point
MissingLinks.get_geometry
MissingLinks.CandidateLink
MissingLinks.add_geom_to_graph!
MissingLinks.compute_net_distance
MissingLinks.add_unless_typemax
MissingLinks.break_long_line
MissingLinks.index_graph_nodes
MissingLinks.compute_link_score
MissingLinks.split_link!
MissingLinks.collapse_realized_graph!
MissingLinks.get_xy
MissingLinks.link_points!
MissingLinks.write_tntp
MissingLinks.realize_graph
MissingLinks.index_candidate_links
MissingLinks.index_graph_edges
MissingLinks.closest_edge
Base.reverse
```