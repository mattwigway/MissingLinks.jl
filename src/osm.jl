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
            ("highway" .=> ["footway", "cycleway", "pedestrian", "track", "sidewalk", "service", "road", "steps", "path", "crossing", "residential"])...,
            ("sidewalk" .=> ["yes", "both", "left", "right"])...,
            "sidewalk:left" .=> "yes",
            "sidewalk:right" .=> "yes",
            "sidewalk:both" .=> "yes"
        ]),
        Set([
            "foot" => "no",
            "access" => "no"
        ]),
        Set([
            # TODO currently foot=yes on e.g. an arterial that doesn't have a sidewalk overrides
            # that we don't want to walk there generally? 
            "foot" => "yes"
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
        if is_traversable(settings, way) && length(way.nodes) > 1
            for (node, next_node) in zip(way.nodes[begin:end-1], way.nodes[begin+1:end])
                if node ∈ nodes_encountered || next_node == first(way.nodes)
                    push!(intersection_nodes, node)
                else
                    push!(nodes_encountered, node)
                end
            end

            push!(intersection_nodes, first(way.nodes))
            push!(intersection_nodes, last(way.nodes))
        end
    end

    # pass 3: build graph
    G = new_graph()

    scan_ways(pbf) do way
        split_nodeid = -1

        if is_traversable(settings, way) && length(way.nodes) > 1
            first_nid = first(way.nodes)
            coords = NTuple{2, Float64}[]

            for (idx, nids) in enumerate(zip(
                way.nodes,
                [way.nodes[begin + 1:end]..., -1],
                ))
                nid, next_nid = nids
                push!(coords, nodes[nid])

                if (idx > 1 && nid ∈ intersection_nodes) || idx == length(way.nodes) ||
                        first_nid == next_nid

                    nid1 = first_nid
                    nid2 = nid

                    # make sure vertices are in graph
                    for osmid in (nid1, nid2)
                        label = VertexID(osmid)

                        if !haskey(G, label)
                            pt = ArchGDAL.createpoint(nodes[osmid])
                            # Project from OSM WGS 84 to desired projection
                            ptpr = ArchGDAL.reproject(pt, GFT.EPSG(4326), projection, order=:trad)
                            G[label] = (ArchGDAL.getx(ptpr, 0), ArchGDAL.gety(ptpr, 0))
                        end
                    end

                    if nid1 > nid2
                        nid1, nid2 = nid2, nid1
                        reverse!(coords)
                    end
                    
                    # create an edge
                    geom = ArchGDAL.createlinestring(coords)

                    # reproject
                    geom = ArchGDAL.reproject(geom, GFT.EPSG(4326), projection, order=:trad)

                    if haskey(G, VertexID(nid1), VertexID(nid2))
                        # can't have a duplicate edge, split it in the middle
                        # we have to do this after reprojection
                        len = ArchGDAL.geomlength(geom)
                        g1 = geom_between(geom, 0, len / 2)
                        g2 = geom_between(geom, len / 2, len)
                        pt = (ArchGDAL.getx(g2, 0), ArchGDAL.gety(g2, 0))
                        
                        G[VertexID(split_nodeid)] = pt

                        G[VertexID(nid1), VertexID(split_nodeid)] = EdgeData((
                            len / 2,
                            haskey(way.tags, "highway") ? way.tags["highway"] : "missing",
                            reverse_geom(g1) # should always go from lower numbered vertex ID, which is always the split since it is negative (not lower numbered Graphs.jl code)
                        ))

                        G[VertexID(nid2), VertexID(split_nodeid)] = EdgeData((
                            len / 2,
                            haskey(way.tags, "highway") ? way.tags["highway"] : "missing",
                            g2
                        ))

                        split_nodeid -= 1
                    else
                        G[VertexID.((nid1, nid2))...] = EdgeData((
                            ArchGDAL.geomlength(geom),
                            haskey(way.tags, "highway") ? way.tags["highway"] : "missing",
                            geom
                        ))
                    end

                    # get ready for next edge
                    first_nid = nid
                    coords = [nodes[nid]]
                end
            end
        end
    end

    return G
end