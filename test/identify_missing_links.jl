@testitem "add_unless_typemax" begin
    import MissingLinks: add_unless_typemax
    for T in [Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64]
        # No overflow
        @test add_unless_typemax(one(T), one(T)) == convert(T, 2)
        @test add_unless_typemax(typemax(T), one(T)) == typemax(T)
        @test_throws OverflowError println(T, add_unless_typemax(typemax(T) - one(T), convert(T, 42)))
    end
end

@testitem "compute_net_distance and fill_matrix" begin
    import MissingLinks: graph_from_gdal, compute_net_distance, fill_matrix!
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
    fill_matrix!(G, dmat; maxdist=12) # 1, 4 should be too far

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
    @test compute_net_distance(dmat, 1, 2, 2, 4, 3, 4, 2, 4) == 2 + 4 + 2

    # adjacent edges - two meters down 1->2 to three meters down 1->3 = 2 + 3
    @test compute_net_distance(dmat, 1, 2, 2, 4, 1, 3, 3, 1) == 2 + 3

    # same edge
    @test compute_net_distance(dmat, 1, 2, 2, 4, 1, 2, 3, 3) == 1

    # and backwards
    @test compute_net_distance(dmat, 1, 2, 3, 3, 1, 2, 2, 4) == 1

    # disconnected
    @test compute_net_distance(dmat, 1, 2, 2, 4, 5, 6, 2, 4) == U
end