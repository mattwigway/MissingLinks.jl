@testitem "Distance matrix functions" begin
    import MissingLinks: graph_from_gdal, compute_net_distance, create_matrix, margin
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

    dmat = create_matrix(G, UInt16; maxdist=12) # 1, 4 should be too far

    U = typemax(UInt16)
    @test dmat ≈ [
        0  6  4  10 U U;
        6  0  10 U  U U;
        4  10 0  6  U U;
        10 U  6  0  U U;
        U  U  U  U  0 6;
        U  U  U  U  6 0
    ]

    # make sure margins work
    # TODO better tests when we use digraphs
    @test collect(margin(dmat, origin=1)) == [
        (1, convert(UInt16, 0)),
        (2, convert(UInt16, 6)),
        (3, convert(UInt16, 4)),
        (4, convert(UInt16, 10))
    ]

    @test collect(margin(dmat, destination=1)) == [
        (1, convert(UInt16, 0)),
        (2, convert(UInt16, 6)),
        (3, convert(UInt16, 4)),
        (4, convert(UInt16, 10))
    ]

    @test collect(margin(dmat, origin=6)) == [
        (5, convert(UInt16, 6)),
        (6, convert(UInt16, 0))
    ]

    @test collect(margin(dmat, destination=6)) == [
        (5, convert(UInt16, 6)),
        (6, convert(UInt16, 0))
    ]

    @test_throws BoundsError margin(dmat, origin=7)
    @test_throws BoundsError margin(dmat, origin=-1)
    @test_throws BoundsError margin(dmat, destination=7)
    @test_throws BoundsError margin(dmat, destination=-1)

    # These are syntactically allowed but don't make sense, so we should get an error
    @test_throws ErrorException margin(dmat)
    @test_throws ErrorException margin(dmat, origin=2, destination=3)

    @test_throws TypeError margin(dmat, origin="UNC")
    @test_throws TypeError margin(dmat, destination=(1, 2))


    # from 2 meters down 1->2 to two meters down 3->4 - 2 meters each end to access node, plus d[2, 3]
    @test compute_net_distance(dmat, 1, 2, UInt16(2), UInt16(4), 3, 4, UInt16(2), UInt16(4)) == UInt16(2 + 4 + 2)

    # adjacent edges - two meters down 1->2 to three meters down 1->3 = 2 + 3
    @test compute_net_distance(dmat, 1, 2, UInt16(2), UInt16(4), 1, 3, UInt16(3), UInt16(1)) == UInt16(2 + 3)

    # same edge
    @test compute_net_distance(dmat, 1, 2, UInt16(2), UInt16(4), 1, 2, UInt16(3), UInt16(3)) == UInt16(1)

    # and backwards
    @test compute_net_distance(dmat, 1, 2, UInt16(3), UInt16(3), 1, 2, UInt16(2), UInt16(4)) == UInt16(1)

    # disconnected
    @test compute_net_distance(dmat, 1, 2, UInt16(2), UInt16(4), 5, 6, UInt16(2), UInt16(4)) == U
end