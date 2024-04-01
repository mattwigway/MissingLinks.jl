@testitem "Non noded graph" begin
    import ArchGDAL as AG
    import ArchGDAL
    import DataFrames: DataFrame, metadata!

    include("util.jl")

    # test that the simple case of a T intersection works
    @testset "T intersection" begin
        tdf = DataFrame(:geom=>[
            AG.createlinestring([[-1.0, 0.0], [1.0, 0.0]]),
            AG.createlinestring([[0.0, 0.0], [0.0, 1.0]])
        ])
        metadata!(tdf, "geometrycolumns", (:geom,))

        tsplit = MissingLinks.semi_to_fully_noded(tdf)

        # make sure the geometries are in order, sorted by flattened coordinates
        sort!(tsplit, :geom, by=collect ∘ Iterators.flatten ∘ get_coords_xy)

        @test get_coords_xy.(tsplit.geom) ≈ [
                [[-1.0, 0.0], [0.0, 0.0]],
                [[0.0, 0.0], [0.0, 1.0]],
                [[0.0, 0.0], [1.0, 0.0]]
            ]

        # Reverse the intersecting edge and put it before the through edge
        tdf = DataFrame(:geom=>[
            AG.createlinestring([[0.0, 1.0], [0.0, 0.0]]),
            AG.createlinestring([[-1.0, 0.0], [1.0, 0.0]])
        ])
        metadata!(tdf, "geometrycolumns", (:geom,))

        tsplit = MissingLinks.semi_to_fully_noded(tdf)

        # make sure the geometries are in order, sorted by flattened coordinates
        sort!(tsplit, :geom, by=collect ∘ Iterators.flatten ∘ get_coords_xy)

        @test get_coords_xy.(tsplit.geom) ≈ [
                [[-1.0, 0.0], [0.0, 0.0]],
                [[0.0, 0.0], [1.0, 0.0]],
                [[0.0, 1.0], [0.0, 0.0]]
            ]
    end

    @testset "X crossing" begin
        # A crossing should not trigger a split
        xdf = DataFrame(:geom=>[
            AG.createlinestring([[0.0, -1.0], [0.0, 1.0]]),
            AG.createlinestring([[-1.0, 0.0], [1.0, 0.0]])
        ])


        xsplit = MissingLinks.semi_to_fully_noded(xdf)

        # make sure the geometries are in order, sorted by flattened coordinates
        sort!(xsplit, :geom, by=collect ∘ Iterators.flatten ∘ get_coords_xy)

        @test get_coords_xy.(xsplit.geom) ≈ [
            [[-1.0, 0.0], [1.0, 0.0]],
            [[0.0, -1.0], [0.0, 1.0]]
        ]
    end

    # snapping tests
    @testset "Snapping" begin
        @testset "Undershoot T intersection" begin
            # here we have a T intersection where the non-through leg undershoots the through leg
            undershoot = DataFrame(:geom => [
                AG.createlinestring([[-1.0, 0.0], [1.0, 0.0]]),
                AG.createlinestring([[0.0, 0.5], [0.0, 1.5]])
            ])

            udf = MissingLinks.semi_to_fully_noded(undershoot, snap_tolerance = 0.75)

            @test get_coords_xy.(udf.geom) ≈ [
                [[-1.0, 0.0], [0.0, 0.0]],
                [[0.0, 0.0], [1.0, 0.0]],
                [[0.0, 0.5], [0.0, 1.5]]
            ]

            # and again with reverse ordering
            undershoot = DataFrame(:geom => [
                AG.createlinestring([[0.0, 1.5], [0.0, 0.5]]),
                AG.createlinestring([[-1.0, 0.0], [1.0, 0.0]])
            ])

            udf = MissingLinks.semi_to_fully_noded(undershoot, snap_tolerance = 0.75)
            sort!(udf, :geom, by=collect ∘ Iterators.flatten ∘ get_coords_xy)

            @test get_coords_xy.(udf.geom) ≈ [
                [[-1.0, 0.0], [0.0, 0.0]],
                [[0.0, 0.0], [1.0, 0.0]],
                [[0.0, 1.5], [0.0, 0.5]]
            ]
        end

        @testset "Overshoot T intersection" begin
            # here we have a T intersection where the non-through leg undershoots the through leg
            overshoot = DataFrame(:geom => [
                AG.createlinestring([[-1.0, 0.0], [1.0, 0.0]]),
                AG.createlinestring([[0.0, -0.5], [0.0, 1.5]])
            ])

            odf = MissingLinks.semi_to_fully_noded(overshoot, snap_tolerance = 0.75)
            sort!(odf, :geom, by=collect ∘ Iterators.flatten ∘ get_coords_xy)


            @test get_coords_xy.(odf.geom) ≈ [
                [[-1.0, 0.0], [0.0, 0.0]],
                [[0.0, -0.5], [0.0, 1.5]],
                [[0.0, 0.0], [1.0, 0.0]]
            ]

            overshoot = DataFrame(:geom => [
                AG.createlinestring([[0.0, 1.5], [0.0, -0.5]]),
                AG.createlinestring([[-1.0, 0.0], [1.0, 0.0]]),
            ])

            odf = MissingLinks.semi_to_fully_noded(overshoot, snap_tolerance = 0.75)
            sort!(odf, :geom, by=collect ∘ Iterators.flatten ∘ get_coords_xy)

            @test get_coords_xy.(odf.geom) ≈ [
                [[-1.0, 0.0], [0.0, 0.0]],
                [[0.0, 0.0], [1.0, 0.0]],
                [[0.0, 1.5], [0.0, -0.5]]
            ]
        end
    end

    @testset "Offset X intersection" begin
        # This tests an intersection like this:
        #     |
        #  ---|
        #     |---
        #     |
        # where the offset between the roads is very small, so they get snapped to the same point

        offset = DataFrame(:geom => [
            AG.createlinestring([[0.0, -2.0], [0.0, 2.0]]),
            # NB slight undershoot as well
            AG.createlinestring([[-1.0, 0.1], [-0.1, 0.0]]),
            AG.createlinestring([[0.1, 0.1], [1.0, 0.1]])
        ])

        offdf = MissingLinks.semi_to_fully_noded(offset, snap_tolerance=0.5, split_tolerance=0.5)
        sort!(offdf, :geom, by=collect ∘ Iterators.flatten ∘ get_coords_xy)

        @test get_coords_xy.(offdf.geom) ≈ [
            [[-1.0, 0.1], [-0.1, 0.0]],
            [[0.0, -2.0], [0.0, 0.0]], # combined split happens at earliest point in line
            [[0.0, 0.0], [0.0, 2.0]],
            [[0.1, 0.1], [1.0, 0.1]]
        ]
    end
end