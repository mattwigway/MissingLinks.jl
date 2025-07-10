const Tag = Pair{String, String}

"""
Settings about what streets are considered walkable
"""
struct TraversalPermissionSettings
    walkable_tags::Set{Tag}
    not_walkable_tags::Set{Tag}
    override_walkable_tags::Set{Tag}
end

# default settings
function TraversalPermissionSettings()
    TraversalPermissionSettings(
        Set([
            ("highway" .=> ["footway", "pedestrian", "track", "sidewalk", "service", "road", "steps", "path", "crossing", "residential"])...,
            ("sidewalk" .=> ["yes", "both", "left", "right"])...,
            "sidewalk:left" .=> "yes",
            "sidewalk:right" .=> "yes",
            "sidewalk:both" .=> "yes"
        ]),
        Set([
            "foot" => "no",
            "access" => "no"
        ])
    )
end

function is_traversable(settings, way::Way)
    traversable = false
    override_no = false
    override_yes = false
    for tag in pairs(way.tags)
        if tag ∈ settings.walkable_tags
            traversable = true
        elseif tag ∈ settings.not_walkable_tags
            override_no = true
        elseif tag ∈ settings.override_walkable_tags
            override_yes = true
        end
    end

    return override_yes || traversable && !override_no
end

function graph_from_osm(pbf, settings, projection)
    nodes = Dict{Int64, NTuple{2, Float64}}()
    nodes_encountered = Set{Int64}()
    intersection_nodes = Set{Int64}()

    # pass 1: nodes
    scan_nodes(pbf) do node
        nodes[node.id] = (node.lon, node.lat)
    end

    # pass 2: intersection nodes
    scan_ways(pbf) do way
        if is_traversable(settings, way)
            for node in way.nodes
                if node ∈ nodes_encountered
                    push!(intersection_nodes, node)
                else
                    push!(nodes_encountered, node)
                end
            end
        end
    end

    # pass 3: build graph
    G = new_graph()

    for node ∈ intersection_nodes
        G[VertexID(node)] = nodes[node]
    end

    scan_ways(pbf) do way
        if is_traversable(settings, way) && length(way.nodes) > 1
            first_nid = first(way.nodes)
            coords = []

            for (idx, nids) in enumerate(zip(way.nodes, [way.nodes[begin + 1:end]..., -1]))
                nid, next_nid = nids
                push!(coords, collect(nodes[nid]))
                # split if there is already an edge from first_nid to the next node
                if (idx > 1 && nid ∈ intersection_nodes) || idx == length(way.nodes) ||
                        has_key(G, (first_nid, next_nid))
                    # create an edge
                    geom = ArchGDAL.createlinestring(coords)

                    # reproject from WGS 84 to local projection
                    geom = ArchGDAL.reproject(geom, GFT.EPSG(4326), projection)

                    G[first_nid, nid] = EdgeData(
                        ArchGDAL.length(geom),
                        has_key(way.tags, "highway") ? way.tags["highway"] : missing,
                        geom
                    )

                    # get ready for next edge
                    first_nid = nid
                    coords = [nodes[nid]]
                end
            end
        end
    end

    return G
end