EdgeData = @NamedTuple{length_m::Float64, link_type::String, geom::ArchGDAL.IGeometry{ArchGDAL.wkbLineString}}

# if the geographic distance from two nodes is less than the snapping tolerance, we connect them
# if the network distance is more than NODE_CONNECT_FACTOR * the geographic distance. Note that
# this is for automated fixing of graph errors, not for missing link identification, which uses a
# different set of thresholds
const NODE_CONNECT_FACTOR = 2

# wrap int64, avoid confusion with Int64 Graphs.jl numbers
struct VertexID
    id::Int64
end

Base.isless(a::VertexID, b::VertexID) = a.id < b.id
Base.isequal(a::VertexID, b::VertexID) = a.id == b.id

function find_or_create_vertex!(G, end_node_idx, location::NTuple{2, Float64}, tolerance)
    loc = collect(location) # tuple to vector
    existing = filter(map(x -> label_for(G, x), intersects(end_node_idx, loc .- tolerance, loc .+ tolerance))) do candidate
        # is it closer than the tolerance (l2-norm is euclidean distance)
        candidate_loc = G[candidate]
        norm2(loc .- candidate_loc) ≤ tolerance
    end

    sort!(existing, by=candidate -> norm2(loc .- G[candidate]))

    if isempty(existing)
        # create a new vertex and add it to the index
        vidx = nv(G) + 1
        vid = VertexID(vidx)
        G[vid] = location
        insert!(end_node_idx, vidx, loc, loc)
        return vid
    else
        if length(existing) > 1
            # This could happen e.g. if there are three nearby nodes each tolerance units apart, in a line
            # if the end ones are found first, the middle one will have multiple candidates. We connect all the candidates
            # to the first (closest) candidate, then return that one. This may create a bit of a rat's nest of short links,
            # but that should not affect routing much.
            @warn "At location $location, found $(length(existing)) candidate nodes within tolerance $tolerance; choosing one" map(x -> G[x], existing)
            # closest, rest = Iterators.peel(existing)
            # frloc = G[closest]
            # for node in rest
            #     if !has_edge(G, code_for(G, closest), code_for(G, node))
            #         toloc = G[node]
            #         # NB will not work for digraph
            #         G[closest, node] = (
            #             length_m = norm2(frloc .- toloc),
            #             geom = ArchGDAL.createlinestring([frloc, toloc])
            #         )
            #     end
            # end
        end
        # TODO could use label_forhere
        return first(existing)
    end
end

"Create a new graph with appropriately set types"
function new_graph()
    MetaGraph(
        Graph(), # for walking - undirected
        label_type=VertexID,
        vertex_data_type=NTuple{2, Float64},
        edge_data_type=EdgeData,
        graph_data=nothing,
        weight_function=ed -> ed.length_m,
        default_weight=Inf64
    )
end

"""
This adds edges between nodes that are close to one another.

This used to be done by snapping nodes during graph build, but that can be prone to errors
because it results in nodes being actually moved. Consider this situation:

---a (0.3 m) b--3.1m---c (0.2m) d----

Assuming the vertices are read in a b c d order, b and c will get snapped to a, and d will not get snapped to anything.

Instead, we only snap a very small distance (1e-6 meters). We then add edges between nodes close together.
"""
function add_short_edges!(G, max_edge_length)
    @info "Indexing graph nodes"
    idx = index_graph_nodes(G)
    for (i, vlabel) in enumerate(labels(G))
        if i % 10000 == 0
            @info "Processed $i / $(nv(G)) vertices"
        end

        # find nearby vertices
        loc = G[vlabel]

        # find all candidates ≤ max_edge_length away
        candidates = filter(
                map(candidate -> (label_for(G, candidate), norm2(loc .- G[label_for(G, candidate)])),
                intersects(idx, collect(loc .- max_edge_length), collect(loc .+ max_edge_length)))
            ) do (_, dist)
            dist .≤ max_edge_length            
        end

        for (candidate, geographic_dist) in candidates
            # doing this inside the loop because each added edge may affect other distances. This is, of course, non-ideal for performance.
            dists = dijkstra_shortest_paths(G, [code_for(G, vlabel)], maxdist=geographic_dist * NODE_CONNECT_FACTOR).dists

            if dists[code_for(G, candidate)] > geographic_dist * NODE_CONNECT_FACTOR
                G[vlabel, candidate] = (
                    length_m=geographic_dist,
                    link_type="short",
                    geom=ArchGDAL.createlinestring([loc, G[candidate]])
                )
            end
        end
    end
end

"""
This adds edge(s) representing geom to the graph. It may add multiple edges in cases where
two edges would otherwise be duplicates because they connect the same intersections (see #2).

Returns a tuple of tuples (frv, tov) - usually just one, but some edges may be split before adding to graph.
"""
function add_geom_to_graph!(G, geom, link_type, end_node_idx, tolerance)
    # add to graph
    startpt = get_first_point(geom)[1:2]
    endpt = get_last_point(geom)[1:2]
    frv = find_or_create_vertex!(G, end_node_idx, startpt, tolerance)
    tov = find_or_create_vertex!(G, end_node_idx, endpt, tolerance)

    # geometry always goes from lower-numbered to higher-numbered node
    if tov < frv
        frv, tov = tov, frv
        geom = reverse(geom)
    end

    if !haskey(G, frv, tov)
        G[frv, tov] = EdgeData((ArchGDAL.geomlength(geom), link_type, geom))
        return ((frv, tov),)
    else
        # we already have an edge between these vertices. This can happen when the network looks like this:
        # 
        #   ________
        #  /        \
        # *          *
        #  \________/
        # 
        # break the edge and try again.
        # Note that this is not 100% correct, because we're introducing a node where
        # there wasn't one before. If that happens to be very close to a node on a disconnected
        # segment of the network, those sections of the network will be fused (unlikely). We
        # attempt to ameliorate the situation by placing the extra node very close to the original node,
        # but if it is just far enough to be within the tolerance of another node this could be an issue
        length_m = ArchGDAL.geomlength(geom)
        offset = min(1e-2, length_m * 0.1)
        return (
            add_geom_to_graph!(G, geom_between(geom, 0.0, offset), link_type, end_node_idx, tolerance)...,
            add_geom_to_graph!(G, geom_between(geom, offset, length_m), link_type, end_node_idx, tolerance)...
        )
    end
end

"""
This builds a graph from one or more GeoDataFrame layers.
"""
function graph_from_gdal(layers...; tolerance=1e-6, max_edge_length=250)
    G = new_graph()

    # pass 1: index all end nodes (connections), and add vertices to graph
    end_node_idx = RTree(2)

    total = sum(nrow.(layers))

    for (i, type_geom) in enumerate(Iterators.flatten((zip(l.link_type, l.geom) for l in layers)))
        if i % 1000 == 0
            @info "Processed $i / $total geometries"
        end

        link_type, multigeom = type_geom

        for_each_geom(multigeom) do biggeom
            # break_long_line is a no-op on lines shorter than max_edge_length
            for geom in break_long_line(biggeom, max_edge_length)
                add_geom_to_graph!(G, geom, link_type, end_node_idx, tolerance)
            end
        end
    end

    return G
end

function graph_to_graphml(out, G; pretty=false)
    doc = XMLDocument()
    root = ElementNode("graphml")
    root["xmlns"] = "http://graphml.graphdrawing.org/xmlns"  
    root["xmlns:xsi"] ="http://www.w3.org/2001/XMLSchema-instance"
    root["xsi:schemaLocation"]="http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd"

    setroot!(doc, root)

    length_m = ElementNode("key")
    length_m["id"] = "length_m"
    length_m["for"] = "edge"
    length_m["attr.name"] = "length_m"
    length_m["attr.type"] = "double"
    link!(root, length_m)
    
    x = ElementNode("key")
    x["id"] = "x"
    x["for"] = "node"
    x["attr.name"] = "x"
    x["attr.type"] = "double"
    link!(root, x)

    y = ElementNode("key")
    y["id"] = "y"
    y["for"] = "node"
    y["attr.name"] = "y"
    y["attr.type"] = "double"
    link!(root, y)

    gph = ElementNode("graph")
    gph["id"] = "G"
    gph["edgedefault"] = "directed"
    link!(root, gph)

    for node in labels(G)
        data = G[node]
        n = ElementNode("node")
        n["id"] = node.id
        link!(gph, n)

        x = ElementNode("data")
        x["key"] = "x"
        link!(x, TextNode("$(data[1])"))
        link!(n, x)

        y = ElementNode("data")
        y["key"] = "y"
        link!(y, TextNode("$(data[2])"))
        link!(n, y)
    end

    for edge in edge_labels(G)
        data = G[edge[1], edge[2]]::EdgeData

        e = ElementNode("edge")
        e["source"] = edge[1].id
        e["target"] = edge[2].id
        link!(gph, e)

        length_m = ElementNode("data")
        length_m["key"] = "length_m"
        link!(length_m, TextNode("$(data.length_m)"))
        link!(e, length_m)
    end

    if pretty
        open(out, "w") do os
            prettyprint(os, doc)
        end
    else
        write(out, doc)
    end
end

"""
    graph_to_gis(fn, G)

Write graph G in GIS format to file fn (format determined by extension)
"""
function graph_to_gis(fn, G; crs=nothing)
    gdf = DataFrame((G[t...] for t in edge_labels(G)))
    metadata!(gdf, "geometrycolumns", (:geom,))
    GeoDataFrames.write(fn, gdf, crs=crs)
end

"""
Data cleaning: remove islands that are smaller than a certain size
"""
function remove_tiny_islands(G, min_vertices_to_retain)
    components = if is_directed(G)
        strongly_connected_components(G)
    else
        connected_components(G)
    end

    component_sizes = length.(components)
    @info "Removing $(sum(component_sizes .< min_vertices_to_retain)) components (out of $(length(components)) total) with $(min_vertices_to_retain - 1) or fewer vertices. " *
        "($(sum(component_sizes[component_sizes .< min_vertices_to_retain])) / $(nv(G)) vertices)"

    # Workaround for https://github.com/JuliaGraphs/MetaGraphsNext.jl/issues/74: create
    # a new graph retaining only what we want
    vertices_to_keep = Set(label_for.(Ref(G), Iterators.flatten(components[component_sizes .≥ min_vertices_to_retain])))

    G2 = new_graph()

    for v in labels(G)
        if v ∈ vertices_to_keep
            G2[v] = G[v]
        end
    end

    for (v1, v2) in edge_labels(G)
        if v1 ∈ vertices_to_keep && v2 ∈ vertices_to_keep
            G2[v1, v2] = G[v1, v2]
        end
    end

    return G2
end

function index_graph_edges(G)
    edgeidx = RTree(2)
    edges = Vector{NTuple{2, VertexID}}()
    for (v1, v2) in edge_labels(G)
        env = ArchGDAL.envelope(G[v1, v2].geom)
        push!(edges, (v1, v2))
        insert!(edgeidx, length(edges), [env.MinX, env.MinY], [env.MaxX, env.MaxY])
    end
    return (edgeidx, edges)
end

"""
Create a spatial index from Graphs.jl vertex codes.
"""
function index_graph_nodes(G)
    nodeidx = RTree(2)
    for v in vertices(G)
        loc = G[label_for(G, v)]
        insert!(nodeidx, v, collect(loc), collect(loc))
    end

    return nodeidx
end