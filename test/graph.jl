@testitem "find_or_create_vertex!" begin
    import ArchGDAL as AG
    import LibSpatialIndex: RTree
    import MissingLinks: VertexID, find_or_create_vertex!
    import Graphs: ne, nv

    G = MissingLinks.new_graph()

    end_node_idx = RTree(2)

    vx_a = find_or_create_vertex!(G, end_node_idx, (1.0, 1.0), 0.011)
    @test vx_a == VertexID(1)

    vx_b = find_or_create_vertex!(G, end_node_idx, (1.02, 1.0), 0.011)
    @test vx_b == VertexID(2)

    vx_c = find_or_create_vertex!(G, end_node_idx, (0.99, 1.0), 0.011)
    @test vx_c == VertexID(1)

    # near vx_c, but that wasn't actually added to graph
    vx_d = find_or_create_vertex!(G, end_node_idx, (0.98, 1.0), 0.011)
    @test vx_d == VertexID(3)

    # between vid1 and vid2, slightly closer to 2
    vx_e = find_or_create_vertex!(G, end_node_idx, (1.0101, 1.0), 0.011)
    @test vx_e == VertexID(2)

    # and diagonal distance
    doff = 0.01 / √2

    # make sure we're using Euclidean distance
    @test doff * 2 > 0.011
    vx_f = find_or_create_vertex!(G, end_node_idx, (1.0 - doff, 1.0 - doff), 0.011)
    @test vx_f == VertexID(1)
end

@testitem "add_geom_to_graph!, add_short_edges!, and remove_islands!" begin
    #=
    The graph we're building looks like this:

    1------2
    |
    |    /\
    |   /  \  
    |  /    \
    | /      \
    3 ------- 4
    5----6
        8  \
         \  \
         |  7
          \_/
 
    9 ---- 10

    3 and 5 are 1m apart but are not connected, add_short_edges will connect them
    6, 7, and 8 are <2m apart as the horse flies and <1m apart as the horse flies, add_short_edges should not connect them
    3 and 4 have two edges between them
    9 and 10 are several meters away and should be eliminated by remove_islands!
    =#

    import ArchGDAL as AG
    import ArchGDAL
    import LibSpatialIndex: RTree
    import MissingLinks: VertexID, add_geom_to_graph!, add_short_edges!, remove_tiny_islands
    import Graphs: ne, nv
    import Base: +
    include("geom/util.jl")

    +(x::VertexID, y::Integer) = VertexID(x.id + y)

    G = MissingLinks.new_graph()
    end_node_idx = RTree(2)

    tol = 1e-6

    # start by building up the simple graph
    v1v2 = add_geom_to_graph!(G, AG.createlinestring([[0.0, 0.0], [10.0, 0.0]]), missing, end_node_idx, tol)
    @test length(v1v2) == 1
    v1, v2 = first(v1v2)
    @test get_coords_xy(G[v1, v2].geom) ≈ [[0.0, 0.0], [10.0, 0.0]]
    @test G[v1, v2].length_m ≈ 10.0

    v3v4 = add_geom_to_graph!(G, AG.createlinestring([[0.0, 10.0], [10.0, 10.0]]), missing, end_node_idx, tol)
    @test v3v4 == ((v2 + 1, v2 + 2),)
    v3, v4 = first(v3v4)

    # v3v1 is backwards; lower vertex ID always is start of edge and geometry should get reversed
    v3v1 = add_geom_to_graph!(G, AG.createlinestring([[0.0, 10.0], [0.0, 0.0]]), missing, end_node_idx, tol)
    @test v3v1 == ((v1, v3),)
    @test get_coords_xy(G[v1, v3].geom) ≈ [[0.0, 0.0], [0.0, 10.0]]

    # there is already v3v4, so a new node should get inserted
    v4v3 = add_geom_to_graph!(G, AG.createlinestring([[10.0, 10.0], [5.0, 5.0], [0.0, 10.0]]), missing, end_node_idx, tol)
    @test v4v3 == (
        (v3, VertexID(5)),
        (v4, VertexID(5))
    )

    # v3' is the added vertex
    v3b = v4v3[1][2]

    @test G[v3, v3b].length_m + G[v4, v3b].length_m ≈ AG.geomlength(AG.createlinestring([[10.0, 10.0], [5.0, 5.0], [0.0, 10.0]]))
    # the breakpoint is offset 1e-2 units down the line. Since the line is 45 degrees, this should be 1e-2 / √2 in the X and Y dimensions
    expected_breakpoint = [1e-2 / √2, 10.0 - 1e-2 / √2]
    @test get_coords_xy(G[v3, v3b].geom) ≈ [[0.0, 10.0], expected_breakpoint]
    @test get_coords_xy(G[v4, v3b].geom) ≈ [[10.0, 10.0], [5.0, 5.0], expected_breakpoint]

    # now the second component, which will get linked to the first by add_short_edges!
    v5v6 = add_geom_to_graph!(G, AG.createlinestring([[0.0, 11.0], [5.0, 11.0]]), missing, end_node_idx, tol)
    @test length(v5v6) == 1
    v5, v6 = v5v6[1]
    
    v6v7 = add_geom_to_graph!(G, AG.createlinestring([[5.0, 11.0], [5.1, 11.1], [5.11, 11.0]]), missing, end_node_idx, tol)
    @test v6v7 == ((v6, v6 + 1),)
    v7 = v6v7[1][2]

    v7v8 = add_geom_to_graph!(G, AG.createlinestring([[5.11, 11.0], [5.1, 11.0]]), missing, end_node_idx, tol)
    @test v7v8 == ((v7, v7 + 1),)
    v8 = v7v8[1][2]

    # will be removed by island pruning
    v9v10 = add_geom_to_graph!(G, AG.createlinestring([[0.0, 13.1], [13.1, 15.0]]), missing, end_node_idx, tol)
    @test v9v10 == ((v8 + 1, v8 + 2),)
    v9, v10 = v9v10[1]


    # one more than in the picture above because v4v3 is split
    @test nv(G) == 11
    @test ne(G) == 9

    add_short_edges!(G, 2)
    @test nv(G) == 11
    # one more edge, v3v5
    @test ne(G) == 10
    @test haskey(G, v3, v5)
    @test G[v3, v5].length_m ≈ 1.0
    @test get_coords_xy(G[v3, v5].geom) ≈ [[0.0, 10.0], [0.0, 11.0]]

    # v6v8 should not be in the graph even though they're close because they're also close in terms of network distance
    @test !haskey(G, v6, v8)

    # v9v10 should disappear when removing islands
    @test haskey(G, v9, v10)
    G = remove_tiny_islands(G, 3)
    # v9, v10 should be gone
    @test nv(G) == 9
    @test ne(G) == 9
    @test !haskey(G, v9, v10)
end