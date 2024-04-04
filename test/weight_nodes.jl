@testitem "Weight nodes" begin
    import ArchGDAL as AG
    import MissingLinks: graph_from_gdal, create_graph_weights
    import DataFrames: DataFrame, metadata!

    @testset "Points" begin
        # Graph looks like this, *s are points + are nodes
        # ------+------+ *       *
        # |  *  |
        # +-----+---+
        # Point in the middle should get assigned to nodes 1-3. Point on right should get assigned
        # to node 2/4. Point on far right should be unassigned. Nothing should be assigned to node 5
        G = graph_from_gdal(DataFrame(:geom=>[
            AG.createlinestring([[0.0, 38.0], [0.0, 0.0], [38.0, 0.0]]),
            AG.createlinestring([[0.0, 38.0], [38.0, 38.0]]),
            AG.createlinestring([[38.0, 38.0], [38.0, 0.0]]),
            AG.createlinestring([[38.0, 0.0], [80.0, 0.0]]),
            AG.createlinestring([[38.0, 38.0], [40.0, 38.0]])
        ]))

        points = DataFrame([
            (weight1=1.0, weight2=3.0, geom=AG.createpoint([19.0, 19.0])),
            (weight1=3.0, weight2=5.0, geom=AG.createpoint([81.0, 0.0])),
            (weight1=42.0, weight2=7.0, geom=AG.createpoint([101.0, 0.0]))
        ])

        metadata!(points, "geometrycolumns", (:geom,))

        @test create_graph_weights(G, points, [:weight1, :weight2], 20) == [
            1.0 / 3          1.0;
            1.5 + 1.0 / 3    3.5;
            1.0 / 3          1.0;
            1.5              2.5;
            0.0              0.0
        ]
    end

    @testset "Polygons" begin
        # Graph is the same as before, but scaled up a bit so that the center of the polygon is not within 20 meters of the edges, but the edges are.
        G = graph_from_gdal(DataFrame(:geom=>[
            AG.createlinestring([[0.0, 42.0], [0.0, 0.0], [42.0, 0.0]]),
            AG.createlinestring([[0.0, 42.0], [42.0, 42.0]]),
            AG.createlinestring([[42.0, 42.0], [42.0, 0.0]]),
            AG.createlinestring([[42.0, 0.0], [84.0, 0.0]]),
            AG.createlinestring([[42.0, 42.0], [44.0, 42.0]])
        ]))

        polys = DataFrame([
            (weight1=1.0, weight2=3.0,  geom=AG.createpolygon([[[19.0, 19.0], [19.0, 23.0], [23.0, 23.0], [23.0, 19.0], [19.0, 19.0]]])),
            (weight1=3.0, weight2=5.0,  geom=AG.createpolygon([[[86.0, 0.0], [87.0, 0.0], [87.0, 1.0], [86.0, 0.0]]])),
            (weight1=42.0, weight2=7.0, geom=AG.createpolygon([[[105.0, 0.0], [106.0, 0.0], [106.0, 1.0], [105.0, 0.0]]])),
        ])

        metadata!(polys, "geometrycolumns", (:geom,))

        @test create_graph_weights(G, polys, [:weight1, :weight2], 20) == [
            1.0 / 3          1.0;
            1.5 + 1.0 / 3    3.5;
            1.0 / 3          1.0;
            1.5              2.5;
            0.0              0.0
        ]
    end

    @testset "Multipolygons" begin
        # Graph is same as Polygons, but the center polygon is now split in two and one part is near the top, one near the bottom, to make sure it is still
        # assigned to both edges.
        G = graph_from_gdal(DataFrame(:geom=>[
            AG.createlinestring([[0.0, 42.0], [0.0, 0.0], [42.0, 0.0]]),
            AG.createlinestring([[0.0, 42.0], [42.0, 42.0]]),
            AG.createlinestring([[42.0, 42.0], [42.0, 0.0]]),
            AG.createlinestring([[42.0, 0.0], [84.0, 0.0]]),
            AG.createlinestring([[42.0, 42.0], [44.0, 42.0]])
        ]))

        polys = DataFrame([
            (weight1=1.0, weight2=3.0,  geom=AG.createmultipolygon([
                [[[19.0, 19.0], [19.0, 21.0], [21.0, 21.0], [21.0, 19.0], [19.0, 19.0]]],
                [[[21.0, 21.0], [21.0, 23.0], [23.0, 23.0], [23.0, 21.0], [21.0, 21.0]]]
            ])),
            (weight1=3.0, weight2=5.0,  geom=AG.createmultipolygon([[[[86.0, 0.0], [87.0, 0.0], [87.0, 1.0], [86.0, 0.0]]]])),
            (weight1=42.0, weight2=7.0, geom=AG.createmultipolygon([[[[105.0, 0.0], [106.0, 0.0], [106.0, 1.0], [105.0, 0.0]]]])),
        ])

        metadata!(polys, "geometrycolumns", (:geom,))

        @test create_graph_weights(G, polys, [:weight1, :weight2], 20) == [
            1.0 / 3          1.0;
            1.5 + 1.0 / 3    3.5;
            1.0 / 3          1.0;
            1.5              2.5;
            0.0              0.0
        ]
    end
end