@testitem "compute_net_distance and fill_matrix" begin
    import MissingLinks: graph_from_gdal, compute_net_distance, fill_distance_matrix!
    import DataFrames: DataFrame
    import ArchGDAL as AG

    #= Graph looks like this:
    1234567

1    1----2
2    |
3    |
4    |
5    3----4
6
7    5----6

    Each row/column is one unit
    =#

    G = graph_from_gdal(DataFrame(:geom => [
        AG.createlinestring([[1.0, 1.0], [7.0, 1.0]]),
        AG.createlinestring([[1.0, 1.0], [1.0, 5.0]]),
        AG.createlinestring([[1.0, 5.0], [7.0, 5.0]]),
        AG.createlinestring([[1.0, 7.0], [7.0, 7.0]])
    ]))

    dmat = zeros(UInt16, (6, 6))
    fill_distance_matrix!(G, dmat; maxdist=12) # 1, 4 should be too far

    U = typemax(UInt16)
    @test dmat ≈ [
        0  6  4  10 U U;
        6  0  10 U  U U;
        4  10 0  6  U U;
        10 U  6  0  U U;
        U  U  U  U  0 6;
        U  U  U  U  6 0
    ]

    # from 2 meters down 1->2 to two meters down 3->4 - 2 meters each end to access node, plus d[2, 3]
    @test compute_net_distance(G, dmat, 1, 2, UInt16(2), UInt16(4), 3, 4, UInt16(2), UInt16(4)) == UInt16(2 + 4 + 2)

    # adjacent edges - two meters down 1->2 to three meters down 1->3 = 2 + 3
    @test compute_net_distance(G, dmat, 1, 2, UInt16(2), UInt16(4), 1, 3, UInt16(3), UInt16(1)) == UInt16(2 + 3)

    # same edge
    @test compute_net_distance(G, dmat, 1, 2, UInt16(2), UInt16(4), 1, 2, UInt16(3), UInt16(3)) == UInt16(1)

    # and backwards
    @test compute_net_distance(G, dmat, 1, 2, UInt16(3), UInt16(3), 1, 2, UInt16(2), UInt16(4)) == UInt16(1)

    # disconnected
    @test compute_net_distance(G, dmat, 1, 2, UInt16(2), UInt16(4), 5, 6, UInt16(2), UInt16(4)) == U
end

@testitem "all possible combinations of compute_net_dist" begin
    # Test all possible combinations of which ends are closer.
    # The graph looks like this:
    #  |\
    #  | \ 
    #  |
    #  | /
    #  |/
    #  We build it with the two angled edges in all possible orientations to make sure the correct distance is calculated
    import MissingLinks: graph_from_gdal, compute_net_distance, fill_distance_matrix!
    import DataFrames: DataFrame
    import ArchGDAL as AG

    for reverse1 in [false, true]
        for reverse2 in [false, true]
            for reversemid in [false, true]
                line1 = [[0.0, 0.0], [0.0, 1.0]]
                line2 = [[3.0, 0.0], [3.0, 1.0]]
                mid = [[0.0, 0.0], [3.0, 0.0]]

                G = graph_from_gdal(DataFrame(:geom=>[
                    AG.createlinestring(reverse1 ? reverse(line1) : line1),
                    AG.createlinestring(reverse2 ? reverse(line2) : line2),
                    AG.createlinestring(reversemid ? reverse(mid) : mid)
                ]))

                # true distance from endpoints should be 1 + 1 + 3
                dmat = zeros(UInt16, (4, 4))
                fill_distance_matrix!(G, dmat)

                @test compute_net_distance(
                    G,
                    dmat,
                    1,
                    2,
                    # if one is reversed, route from the start of one (the disconnected part)
                    # otherwise, the end
                    reverse1 ? zero(UInt16) : one(UInt16),
                    reverse1 ? one(UInt16) : zero(UInt16),
                    3,
                    4,
                    reverse2 ? zero(UInt16) : one(UInt16),
                    reverse2 ? one(UInt16) : zero(UInt16),
                ) ≈ 5
            end
        end
    end
end