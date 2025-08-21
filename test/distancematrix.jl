@testitem "Can create multiple distance matrices" setup=[InitializeIdentify] begin
    # When I first implemented the DistanceMatrix struct based on SQLite, I was so confused
    # because everything seemed to work but most of the tests failed. I eventually figured out that
    # using DBInterface.@prepare meant that all SQL statements were going to whatever database was
    # opened first, so subsequent distance matrices were empty. Apparently because DBInterface.@prepare
    # always prepares a statement for the first database it is used with.
        G = graph_from_gdal(DataFrame(:geom=>[
        AG.createlinestring([[0.0, 0.0], [10.0, 0.0]]),
        AG.createlinestring([[0.0, 0.0], [5.0, 2.0], [6.0, 1.0]])
    ]))

    dmat = calculate_distances(G)
    dmat2 = calculate_distances(G)

    @test !ismissing(dmat[VertexID(1), VertexID(1)])
    @test !ismissing(dmat[VertexID(1), VertexID(2)])
    @test !ismissing(dmat[VertexID(1), VertexID(3)])
    @test !ismissing(dmat[VertexID(2), VertexID(1)])
    @test !ismissing(dmat[VertexID(2), VertexID(2)])
    @test !ismissing(dmat[VertexID(2), VertexID(3)])
    @test !ismissing(dmat[VertexID(3), VertexID(1)])
    @test !ismissing(dmat[VertexID(3), VertexID(2)])
    @test !ismissing(dmat[VertexID(3), VertexID(3)])

    @test !ismissing(dmat2[VertexID(1), VertexID(1)])
    @test !ismissing(dmat2[VertexID(1), VertexID(2)])
    @test !ismissing(dmat2[VertexID(1), VertexID(3)])
    @test !ismissing(dmat2[VertexID(2), VertexID(1)])
    @test !ismissing(dmat2[VertexID(2), VertexID(2)])
    @test !ismissing(dmat2[VertexID(2), VertexID(3)])
    @test !ismissing(dmat2[VertexID(3), VertexID(1)])
    @test !ismissing(dmat2[VertexID(3), VertexID(2)])
    @test !ismissing(dmat2[VertexID(3), VertexID(3)])

    @test dmat[VertexID(1), VertexID(1)] == dmat2[VertexID(1), VertexID(1)]
    @test dmat[VertexID(1), VertexID(2)] == dmat2[VertexID(1), VertexID(2)]
    @test dmat[VertexID(1), VertexID(3)] == dmat2[VertexID(1), VertexID(3)]
    @test dmat[VertexID(2), VertexID(1)] == dmat2[VertexID(2), VertexID(1)]
    @test dmat[VertexID(2), VertexID(2)] == dmat2[VertexID(2), VertexID(2)]
    @test dmat[VertexID(2), VertexID(3)] == dmat2[VertexID(2), VertexID(3)]
    @test dmat[VertexID(3), VertexID(1)] == dmat2[VertexID(3), VertexID(1)]
    @test dmat[VertexID(3), VertexID(2)] == dmat2[VertexID(3), VertexID(2)]
    @test dmat[VertexID(3), VertexID(3)] == dmat2[VertexID(3), VertexID(3)]


end