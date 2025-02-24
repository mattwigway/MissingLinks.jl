# API Documentation

## Graph building and troubleshooting

```@docs
semi_to_fully_noded
graph_from_gdal
graph_to_graphml
graph_to_gis
remove_tiny_islands
add_short_edges!
find_disconnected_crossings
find_dead_ends
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

## Service area analysis

```@docs
service_area
```

## Internal

```@docs
MissingLinks.new_graph
MissingLinks.for_each_geom
MissingLinks.get_first_point
MissingLinks.CandidateLink
MissingLinks.add_geom_to_graph!
MissingLinks.compute_net_distance
MissingLinks.add_unless_typemax
MissingLinks.break_long_line
MissingLinks.index_graph_nodes
MissingLinks.find_point_on_edge
MissingLinks.compute_link_score
Base.reverse
```