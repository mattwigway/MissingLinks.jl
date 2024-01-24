EdgeData = @NamedTuple{length_m::Float64, geom::ArchGDAL.IGeometry{ArchGDAL.wkbLineString}}

"Run the provided function on each constituent geometry of a multigeometry"
function for_each_geom(f, g::ArchGDAL.IGeometry{ArchGDAL.wkbMultiLineString})
    for i in 1:ArchGDAL.ngeom(g)
        # archgdal has zero based indexing
        f(ArchGDAL.getgeom(g, i - 1))
    end
end

for_each_geom(f, g::ArchGDAL.IGeometry{ArchGDAL.wkbLineString}) = f(g)

"get the first point of a line"
get_first_point(g::ArchGDAL.IGeometry{ArchGDAL.wkbLineString}) = ArchGDAL.getpoint(g, 0)
# 0 based indexing in GDAL
get_last_point(g::ArchGDAL.IGeometry{ArchGDAL.wkbLineString}) = ArchGDAL.getpoint(g, ArchGDAL.ngeom(g) - 1)

# wrap int64, avoid confusion with Int64 Graphs.jl numbers
struct VertexID
    id::Int64
end

function find_or_create_vertex!(G, end_node_idx, location, tolerance)
    loc = collect(location) # tuple to vector
    existing = filter(intersects(end_node_idx, loc .- tolerance, loc .+ tolerance)) do candidate
        # is it closer than the tolerance (l2-norm is euclidean distance)
        candidate_loc = G[VertexID(candidate)]
        norm2(loc .- candidate_loc) ≤ tolerance
    end

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
            # if the end ones are found first, the middle one will have multiple candidates
            @warn "At location $location, found $(length(existing)) candidate nodes within tolerance $tolerance:" map(x -> G[VertexID(x)], existing)
        end

        # TODO could use label_for here
        return VertexID(first(existing))
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
This builds a graph from one or more GeoDataFrame layers.
"""
function graph_from_gdal(layers...; tolerance=1e-6)
    G = new_graph()

    # pass 1: index all end nodes (connections), and add vertices to graph
    end_node_idx = RTree(2)

    total = sum(nrow.(layers))

    for (i, multigeom) in enumerate(Iterators.flatten((l.geom for l in layers)))
        if i % 1000 == 0
            @info "Processed $i / $total geometries"
        end

        for_each_geom(multigeom) do geom
            # add to graph
            startpt = get_first_point(geom)[1:2]
            endpt = get_last_point(geom)[1:2]
            frv = find_or_create_vertex!(G, end_node_idx, startpt, tolerance)
            tov = find_or_create_vertex!(G, end_node_idx, endpt, tolerance)

            G[frv, tov] = EdgeData((ArchGDAL.geomlength(geom), geom))
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