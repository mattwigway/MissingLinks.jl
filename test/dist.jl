@testitem "compute_net_distance and fill_matrix" begin
    import MissingLinks: graph_from_gdal, compute_net_distance, calculate_distances, VertexID
    import DataFrames: DataFrame
    import ArchGDAL as AG
    import Graphs: nv
    import MetaGraphsNext: label_for

    function to_mtx(dmat)
        mtx = Matrix{Int64}(undef, nv(G), nv(G))
        for fr in 1:nv(G)
            for to in 1:nv(G)
                val = dmat[label_for(dmat.graph, fr), label_for(dmat.graph, to)]
                # recode missing to -1 so we don't have to deal with missing != missing
                mtx[fr, to] = ismissing(val) ? -1 : val
            end
        end
        return mtx
    end

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

    dmat = calculate_distances(G; maxdist=12) # 1, 4 should be too far

    U = -1
    @test to_mtx(dmat) == [
        0  6  4  10 U U;
        6  0  10 U  U U;
        4  10 0  6  U U;
        10 U  6  0  U U;
        U  U  U  U  0 6;
        U  U  U  U  6 0
    ]

    # from 2 meters down 1->2 to two meters down 3->4 - 2 meters each end to access node, plus d[2, 3]
    @test compute_net_distance(dmat, VertexID(1), VertexID(2), 2, 4, VertexID(3), VertexID(4), 2, 4) == 2 + 4 + 2

    # adjacent edges - two meters down 1->2 to three meters down 1->3 = 2 + 3
    @test compute_net_distance(dmat, VertexID(1), VertexID(2), 2, 4, VertexID(1), VertexID(3), 3, 1) == 2 + 3

    # same edge
    @test compute_net_distance(dmat, VertexID(1), VertexID(2), 2, 4, VertexID(1), VertexID(2), 3, 3) == 1

    # and backwards
    @test compute_net_distance(dmat, VertexID(1), VertexID(2), 3, 3, VertexID(1), VertexID(2), 2, 4) == 1

    # disconnected
    @test ismissing(compute_net_distance(dmat, VertexID(1), VertexID(2), 2, 4, VertexID(5), VertexID(6), 2, 4))
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
    import MissingLinks: graph_from_gdal, compute_net_distance, calculate_distances, VertexID
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
                dmat = calculate_distances(G)

                @test compute_net_distance(
                    dmat,
                    VertexID(1),
                    VertexID(2),
                    # if one is reversed, route from the start of one (the disconnected part)
                    # otherwise, the end
                    reverse1 ? 0 : 1,
                    reverse1 ? 1 : 0,
                    VertexID(3),
                    VertexID(4),
                    reverse2 ? 0 : 1,
                    reverse2 ? 1 : 0,
                ) ≈ 5
            end
        end
    end
end