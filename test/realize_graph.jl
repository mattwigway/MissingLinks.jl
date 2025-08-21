@testitem "Realize graph" begin
    # Graph looks like this, single lines are existing and double lines are candidates
    # each row or column is 10m.
    # top: read down
    # m             111111111122222222223333
    #      123456789012345678901234567890123
    #      000000000000000000000000000000000
    # 10   A       C== 8 ==E                  10
    # 20   |       |       |       J          20
    # 30   |== 1 ==|       |       |          30
    # 40   |       |== 2 ==F== 3 ==|          40
    # 50   B       |               |          50
    # 60           D== 4 ==G== 5 ==K          60
    # 70                   |       |          70
    # 80                   H== 6 ==|==10= M   80
    # 90                   |       |      |   90
    # 100                  |== 9 ==L==7 ==N   100
    # 110                  |              |   110
    # 120                  I              O   120
    # 130                                 |   130
    # 140                                 P   140
    #
    # m             111111111122222222223333
    #      123456789012345678901234567890123
    #      000000000000000000000000000000000
    #
    # Candidate link 1 connects a middle to a middle
    # Candidate link 2 connects a middle to a start
    # Candidate link 3 connects a start to a middle
    # Candidate link 4 connects an end to a start
    # Candidate link 5 connects a start to an end
    # Candidate link 6 connects and end to a middle
    # Candidate link 7 connects an end to an end
    # Candidate link 8 connects a start to a start
    # Candidate link 9 connects a middle to an end
    # Candidate link 10 connects a middle to a start, but should share start with end of candidate link 6

    # edge AB has one link in middle
    # Edge CD has one link at start, two in middle, one at end
    # Edges EF, GH have links at start and end
    # Edge HN has link at start and middle
    # Edge IJ has link at middle and end
    # Edge JK has links at start, middle, end
    # Edge LM has links at end
    import MissingLinks: VertexID, add_geom_to_graph!, CandidateLink
    import ArchGDAL as AG
    import LibSpatialIndex: RTree
    import Graphs: nv, ne
    import MetaGraphsNext: code_for, label_for, edge_labels
    import LinearAlgebra: norm

    G = MissingLinks.new_graph()
    end_node_idx = RTree(2)
    tol = 1e-6

    a!(g...) = add_geom_to_graph!(G, AG.createlinestring([g...]), "existing", end_node_idx, tol)

    # coordinates
    c = (
        A = [10., 10.],
        B = [10., 50.],
        C = [90., 10.],
        D = [90., 60.],
        E = [170., 10.],
        F = [170., 40.],
        G = [170., 60.],
        H = [170., 80.],
        I = [170., 120.],
        J = [250., 20.],
        K = [250., 60.],
        L = [250., 100.],
        M = [320., 80.],
        N = [320., 100.],
        O = [320., 120.],
        P = [320., 140.]
    )

    a!(c.A, c.B)
    a!(c.C, c.D)
    a!(c.E, c.F)
    a!(c.G, c.H)
    a!(c.H, c.I)
    a!(c.J, c.K)
    a!(c.K, c.L)
    a!(c.M, c.N)
    a!(c.N, c.O)
    a!(c.O, c.P)

    @test nv(G) == 16
    @test ne(G) == 10

    v = (
        A = VertexID(1),
        B = VertexID(2),
        C = VertexID(3),
        D = VertexID(4),
        E = VertexID(5),
        F = VertexID(6),
        G = VertexID(7),
        H = VertexID(8),
        I = VertexID(9),
        J = VertexID(10),
        K = VertexID(11),
        L = VertexID(12),
        M = VertexID(13),
        N = VertexID(14),
        O = VertexID(15),
        P = VertexID(16)
    )


    for k in keys(v)
        # collect converts tuple in graph to vector in test
        @test collect(G[v[k]]) == c[k]
    end

    # Add a link between v1v2 -> v3v4 at y
    function l(v1, v2, v3, v4, y)
        x1 = c[v1][1]
        x2 = c[v3][1]
        CandidateLink(
            code_for(G, v[v1]),
            code_for(G, v[v2]),
            abs(y - c[v1][2]),
            abs(c[v2][2] - y),
            code_for(G, v[v3]),
            code_for(G, v[v4]),
            abs(y - c[v3][2]),
            abs(c[v4][2] - y),
            abs(x2 - x1),
            missing # network length: not relevant to test
        )
    end

    links = [
        l(:A, :B, :C, :D, 30),
        l(:C, :D, :E, :F, 40),
        l(:E, :F, :J, :K, 40),
        l(:C, :D, :G, :H, 60),
        l(:G, :H, :K, :L, 60),
        l(:G, :H, :K, :L, 80),
        l(:K, :L, :M, :N, 100),
        l(:C, :D, :E, :F, 10),
        l(:H, :I, :K, :L, 100),
        l(:K, :L, :M, :N, 80)
    ]

    G2 = realize_graph(G, links)

    # should add six vertices (all the middles)
    @test nv(G2) == nv(G) + 6

    # should add 10 edges (links) + 6 edges (each time an edge is split for a middle)
    @test ne(G2) == ne(G) + 10 + 6

    # all original nodes should be preserved
    for n in v
        @test all(G[n] .≈ G2[n])
    end

    function test_edge(v1, v2, kind)
        if v1 > v2
            v2, v1 = v1, v2
        end

        @test (v1, v2) ∈ collect(edge_labels(G2))
        e = G2[v1, v2]
        @test e.link_type == kind
        @test e.length_m ≈ norm(G2[v2] .- G2[v1], 2)

        # can't use ≈ here, not defined for ArchGDAL geoms. but they should
        # be exact anyways.
        # TOD for some reason some graph geometries have elevation and some don't - why?
        @test all(MissingLinks.get_xy(e.geom) .≈ [
            collect(G2[v1]),
            collect(G2[v2])
        ])
    end

    newv(i) = VertexID(16 + i)

    # A - B should be split at 30
    @test all(G2[newv(1)] .≈ [10., 30.])
    test_edge(v.A, newv(1), "existing")
    test_edge(newv(1), v.B, "existing")

    # C-D should be split at 30 and 40
    @test all(G2[newv(2)] .≈ [90., 30.])
    @test all(G2[newv(3)] .≈ [90., 40.])
    test_edge(v.C, newv(2), "existing")
    test_edge(newv(2), newv(3), "existing")
    test_edge(newv(3), v.D, "existing")

    # Candidate link 1 should exist
    test_edge(newv(1), newv(2), "candidate")

    # E-F, G-H should not be split
    test_edge(v.E, v.F, "existing")
    test_edge(v.G, v.H, "existing")

    # Candidate link 2 should exist and connects to existing vx
    test_edge(newv(3), v.F, "candidate")

    # Candidate link 4 should exist
    test_edge(v.D, v.G, "candidate")

    # Candidate link 8 should exist
    test_edge(v.C, v.E, "candidate")

    # HI should be split by at 100
    # splits happen in order of edges, so HI will be split before JK
    # even though JK is split by a lower numbered candidate
    @test all(collect(G2[newv(4)]) .≈ [170., 100.])
    test_edge(v.H, newv(4), "existing")
    test_edge(newv(4), v.H, "existing")

    # JK should be split at 40 by link 3
    @test all(collect(G2[newv(5)]) .≈ [250., 40.])
    test_edge(v.J, newv(5), "existing")
    test_edge(newv(5), v.K, "existing")

    # KL should be split at 80 (this split shared between link 6 and 10)
    @test all(collect(G2[newv(6)]) .≈ [250., 80.])
    test_edge(v.K, newv(6), "existing")
    test_edge(newv(6), v.L, "existing")

    # Candidate link 3 should exist
    test_edge(v.F, newv(5), "candidate")

    # Candidate link 5 should exist
    test_edge(v.G, v.K, "candidate")

    # Candidate link 6 should exist
    test_edge(v.H, newv(6), "candidate")

    # Candidate link 9 should exist
    test_edge(newv(4), v.L, "candidate")

    # MN, NO, OP not split
    test_edge(v.M, v.N, "existing")
    test_edge(v.N, v.O, "existing")
    test_edge(v.O, v.P, "existing")

    # candidate link 7 should exist
    test_edge(v.L, v.N, "candidate")

    # candidate link 10 should exist
    test_edge(newv(6), v.M, "candidate")
end