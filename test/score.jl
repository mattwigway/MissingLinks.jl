@testitem "Scoring" begin
    import MissingLinks: graph_from_gdal, score_links, CandidateLink, fill_matrix!
    import DataFrames: DataFrame
    import ArchGDAL as AG
    import Graphs: nv
    
    @testset "Missing link across disconnected components" begin
        # This tests a missing link that connects two previously disconnected components
        G = graph_from_gdal(DataFrame(:geom=>[
            # Component 1
            AG.createlinestring([[0.0, 0.0], [1.0, 0.0]]),
            AG.createlinestring([[1.0, 0.0], [1.0, 1.1]]),

            # Component 2
            AG.createlinestring([[2.0, 0.0], [3.0, 0.0]]),
            AG.createlinestring([[2.0, 0.0], [2.0, 1.0]])
        ]))

        dmat = zeros(UInt16, nv(G), nv(G))
        fill_matrix!(G, dmat)

        links = [
            CandidateLink(
                fr_edge_src=1,
                fr_edge_tgt=2,
                fr_dist_from_start=UInt16(1),
                fr_dist_to_end=UInt16(0),
                to_edge_src=4,
                to_edge_tgt=5,
                to_dist_from_start=UInt16(0),
                to_dist_to_end=UInt16(1),
                geographic_length_m=UInt16(1),
                network_length_m=typemax(UInt16)
            )
        ]

        # Every combination of these is unique so we can be sure we got the right combination
        # Though I suppose that is not true since some nodes are accessible more than once 🤔 - still better than
        # just using ones everywhere
        weights = convert.(Float64, [
            0b1,
            0b10,
            0b100,
            0b1000,
            0b10000,
            0b100000
        ])

        # each of six origin nodes gets access to three new destinations
        @test score_links(x -> x < 2000, links, dmat, weights, weights, 2000) == [
            sum(weights[1:3] .* sum(weights[4:6])) + # Nodes 1-3 get access to 4-6
            sum(weights[4:6] .* sum(weights[1:3]))
        ]
    end

    @testset "Missing links across connected component" begin
        # in this case, the locations at each end of the link are already connected, but by a roundabout route
        # that means some nodes are too far from others
        # network looks like this, and we use a 1900m cutoff
        #              5---500m-----6
        #              |            |
        #             500m        500m
        #              |            |
        #  1 --500m -- 2 ***500m*** 3 --500m-- 4
        #
        # In the base network, the following access is possible:
        # 1: 1, 2, 5, 6
        # 2: 1, 3, 5, 6
        # 3: 2, 3, 4, 5, 6
        # 4: 3, 4, 5, 6
        # 5: 1, 2, 3, 4, 5, 6
        # 6: 1, 2, 3, 4, 5, 6

        # When we add the starred "missing link", we additionally can access:
        # 1: 3, 4
        # 2: 4
        # 3: 1
        # 4: 1, 2


        G = graph_from_gdal(DataFrame(:geom=>[
            AG.createlinestring([[0.0, 500.0], [500.0, 500.0]]), # 12
            AG.createlinestring([[1000.0, 500.0], [1500.0, 500.0]]), # 34
            AG.createlinestring([[500.0, 500.0], [500.0, 0.0]]), # 25
            AG.createlinestring([[500.0, 0.0], [1000.0, 0.0]]), # 56
            AG.createlinestring([[1000.0, 0.0], [1000.0, 500.0]]) # 36
        ]); max_edge_length=50_000)

        dmat = zeros(UInt16, nv(G), nv(G))
        fill_matrix!(G, dmat)

        links = [
            CandidateLink(
                fr_edge_src=1,
                fr_edge_tgt=2,
                fr_dist_from_start=UInt16(0),
                fr_dist_to_end=UInt16(500),
                to_edge_src=3,
                to_edge_tgt=4,
                to_dist_from_start=UInt16(0),
                to_dist_to_end=UInt16(500),
                geographic_length_m=UInt16(500),
                network_length_m=UInt16(500)
            )
        ]

        weights = convert.(Float64, [
            0b1,
            0b10,
            0b100,
            0b1000,
            0b10000,
            0b100000
        ])

        @test score_links(x -> x < 1900, links, dmat, weights, weights, 1900) ≈ [
            weights[1] * sum(weights[[3, 4]]) +
            weights[2] * weights[4] +
            weights[3] * weights[1] +
            weights[4] * sum(weights[[1, 2]])
        ]
    end

    @testset "Middles" begin
        # This tests that links in the middle of edges are scored correctly
        # Network looks like this
        #
        # 1------2        5-----6
        #         \      /
        #          \ ** /
        #          /    \
        #  3------4      7------8
        #
        # 56 is coded backwards to test that as well
        # After the link is added, all nodes in each component gain access to all nodes in the other component,
        # _except_:
        #
        # 1 cannot reach 6 or 8
        # 3 cannot reach 6
        # 6 cannot reach 1 or 3
        # 8 cannot reach 1
        #
        # Cutoff distance is 2000m. All of the horizontal links are 900m. The missing link
        # is 100m, and the link is 45m from 4/7 and 72m from 2/5.

        G = graph_from_gdal(DataFrame(:geom=>[
            AG.createlinestring([[0.0, 0.0], [900.0, 0.0]]), # 12
            AG.createlinestring([[0.0, 96.0], [900.0, 96.0]]), # 34
            AG.createlinestring([[1064.0, 0.0], [1964.0, 0.0]]), # 56
            AG.createlinestring([[1064.0, 96.0], [1964.0, 96.0]]), # 78
            AG.createlinestring([[900.0, 0.0], [932.0, 64.0], [900.0, 96.0]]), # 24
            AG.createlinestring([[1064.0, 96.0], [1032.0, 64.0], [1064.0, 0.0]]) # 75 coded backwards
        ]), max_edge_length=50_000)

        dmat = zeros(UInt16, nv(G), nv(G))
        fill_matrix!(G, dmat)

        links = [
            CandidateLink(
                fr_edge_src=2,
                fr_edge_tgt=4,
                fr_dist_from_start=UInt16(72),
                fr_dist_to_end=UInt16(45),
                to_edge_src=5,
                to_edge_tgt=7,
                to_dist_from_start=UInt16(72),
                to_dist_to_end=UInt16(45),
                geographic_length_m=UInt16(100),
                network_length_m=typemax(UInt16)
            )
        ]

        weights = convert.(Float64, [
            0b1,
            0b10,
            0b100,
            0b1000,
            0b10000,
            0b100000,
            0b1000000,
            0b10000000
        ])

        @test score_links(x -> x < 2000, links, dmat, weights, weights, 2000) == [
            weights[1] * sum(weights[[5, 7]]) +
            weights[2] * sum(weights[5:8]) +
            weights[3] * sum(weights[[5, 7, 8]]) +
            weights[4] * sum(weights[5:8]) +
            weights[5] * sum(weights[1:4]) +
            weights[6] * sum(weights[[2, 4]]) +
            weights[7] * sum(weights[1:4]) +
            weights[8] * sum(weights[[2, 3, 4]])
        ]

    end

    # TODO test with gravity metric
end