@testitem "link points" begin
    import MissingLinks: VertexID, link_points!, add_geom_to_graph!, EdgeData
    import Graphs: nv, ne
    import ArchGDAL as AG
    import LibSpatialIndex: RTree
    import LinearAlgebra: norm2

    # Graph looks like this
    #   a  f
    # 1 --- 2
    # |
    # | b
    # | c
    # | d     g
    # | e
    # 3
    # The lowercase letters are the points that get linked in
    # a should get snapped directly to 1
    # b should split 13, and then c should be linked to b. d should split b3
    # e should get snapped directly to 3
    # f should get snapped to 2
    # g is isolated, it should have an edge created when create = true and not when create = false

    for create in [true, false]
        G = MissingLinks.new_graph()
        end_node_idx = RTree(2)
        tol = 1e-6
        a!(g...) = add_geom_to_graph!(G, AG.createlinestring([g...]), "existing", end_node_idx, tol)

        a!([0., 0.], [100., 0.])
        a!([0., 0.], [0., 100.])

        @test nv(G) == 3
        @test ne(G) == 2

        points = AG.createpoint.([
            [0.1, -1.], # a
            [0.2, 25.], # b
            [0.25, 25.5], # c
            [5.0, 75.], # d
            [19.0, 100.0], # e
            [99.5, 0.2], # f
            [100., 100.,] # g
        ])

        vxs = link_points!(G, points; create=create)

        if create
            @test ne(G) == 5
            @test nv(G) == 7
        else
            @test ne(G) == 4
            @test nv(G) == 5
        end

        @test length(vxs) == length(points)
        @test vxs[1] == VertexID(1)

        # split at B
        # we don't need to test for absence as that is covered by the ne() call above
        @test haskey(G, VertexID(1), vxs[2]) 
        @test G[VertexID(1), vxs[2]] ≈ EdgeData((
            25.,
            "existing",
            AG.createlinestring([[0., 0.], [0., 25.]])
        ))
        @test vxs[2] == VertexID(4, :split)
        @test all(G[vxs[2]] .≈ (0., 25.)) # will be snapped to line
        # TODO test geometry
        
        # no split at C
        @test vxs[3] == vxs[2]

        # split at D
        @test vxs[4] == VertexID(5, :split)
        @test all(G[vxs[4]] .≈ (0., 75.))
        @test haskey(G, vxs[2], vxs[4])
        @test G[vxs[2], vxs[4]] ≈ EdgeData((
            50.,
            "existing",
            AG.createlinestring([[0., 25.], [0., 75.]])
        ))

        # no split at E
        @test vxs[5] == VertexID(3)
        @test haskey(G, vxs[4], VertexID(3))
        @test G[vxs[4], VertexID(3)] ≈ EdgeData((
            25.,
            "existing",
            # vxs[4] > VertexID(3)
            AG.createlinestring([[0., 100.], [0., 75.]])
        ))

        # f snapped directly to 2
        @test vxs[6] == VertexID(2)

        if create
            # g has its own edge
            @test vxs[7] == VertexID(6, :island)
            @test all(G[vxs[7]] .≈ (100., 100.))
            @test all(abs.(G[VertexID(7, :island)] .- [100., 100.]) .< 1e-5)
            @test haskey(G, vxs[7], VertexID(7, :island))
            @test G[vxs[7], VertexID(7, :island)] ≈ EdgeData((
                norm2(G[vxs[7]] .- G[VertexID(7, :island)]),
                "island",
                AG.createlinestring([[100., 100.], collect(G[VertexID(7, :island)])])
            ))

        else
            @test ismissing(vxs[7])
        end
    end
end