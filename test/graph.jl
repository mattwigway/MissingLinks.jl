@testitem "find_or_create_vertex!" begin
    import ArchGDAL as AG
    import LibSpatialIndex: RTree
    import MissingLinks: VertexID, find_or_create_vertex!

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
