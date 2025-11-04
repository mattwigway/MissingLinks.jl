@testitem "add_unless_typemax" begin
    import MissingLinks: add_unless_typemax
    for T in [Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64]
        # No overflow
        @test add_unless_typemax(one(T), one(T)) == convert(T, 2)
        @test add_unless_typemax(typemax(T), one(T)) == typemax(T)
        @test_throws OverflowError println(T, add_unless_typemax(typemax(T) - one(T), convert(T, 42)))
    end
end

@testitem "Identify missing links" begin
    import MissingLinks: graph_from_gdal, identify_potential_missing_links, fill_distance_matrix!, CandidateLink, VertexID
    import DataFrames: DataFrame
    import ArchGDAL as AG
    import StructEquality: @struct_hash_equal
    import Graphs: nv

    # make sure we can compare two links by value rather than identity
    @struct_hash_equal CandidateLink

    # Sort the links into a stable order. There can only be one link between an edge pair, so there will
    # be no ties.
    sortlinks!(links) = sort!(links, by=l -> (l.fr_edge_src, l.fr_edge_tgt, l.to_edge_src, l.to_edge_tgt))

    # This first graph is very simple, it looks like > | 
    # We should find a link from the middle of one edge to the middle of the other
    @testset "Snapping to middles\nGraph > |" begin
        gebar = graph_from_gdal(DataFrame(:geom =>[
            AG.createlinestring([[0.0, 0.0], [1.0, 1.0], [0.0, 2.0]]),
            AG.createlinestring([[2.0, 0.0], [2.0, 2.0]])
        ]))

        dmat = zeros(UInt16, (4, 4))
        fill_distance_matrix!(gebar, dmat)
        links = identify_potential_missing_links(gebar, dmat, 2, 5)
        sortlinks!(links)

        @test length(links) == 2
        # only one link, but we found it in both directions
        @test links[1] == reverse(links[2])

        @test links[1] == CandidateLink(
            fr_edge_src=VertexID(1),
            fr_edge_tgt=VertexID(2),
            fr_dist_from_start=one(UInt16), # actually √2 but is rounded
            fr_dist_to_end=UInt16(2),
            to_edge_src=VertexID(3),
            to_edge_tgt=VertexID(4),
            to_dist_from_start=one(UInt16),
            to_dist_to_end=one(UInt16),
            geographic_length_m=one(UInt16),
            network_length_m=typemax(UInt16)
        )
    end

    # Now we check one where the snapping is to the ends
    # Network looks like / \
    @testset "Snapping to ends\nGraph / \\" begin
        fbslash = graph_from_gdal(DataFrame(:geom=>[
            AG.createlinestring([[0.0, 1.0], [1.0, 0.0]]),
            AG.createlinestring([[2.0, 0.0], [3.0, 1.0]])
        ]))

        dmat = zeros(UInt16, (4, 4))
        fill_distance_matrix!(fbslash, dmat)
        links = identify_potential_missing_links(fbslash, dmat, 2, 5)
        sortlinks!(links)

        @test length(links) == 2
        # only one link, but we found it in both directions
        @test links[1] == reverse(links[2])

        @test links[1] == CandidateLink(
            fr_edge_src=VertexID(1),
            fr_edge_tgt=VertexID(2),
            fr_dist_from_start=one(UInt16),
            fr_dist_to_end=zero(UInt16),
            to_edge_src=VertexID(3),
            to_edge_tgt=VertexID(4),
            to_dist_from_start=zero(UInt16),
            to_dist_to_end=one(UInt16),
            geographic_length_m=one(UInt16),
            network_length_m=typemax(UInt16)
        )
    end

    @testset """One end and one middle
    Graph ----------------
          ______/
    """ begin

        # One end and one middle
        # Network looks like this:
        # ----------------
        # ______/

        G = graph_from_gdal(DataFrame(:geom=>[
            AG.createlinestring([[0.0, 0.0], [10.0, 0.0]]),
            AG.createlinestring([[0.0, 2.0], [5.0, 2.0], [6.0, 1.0]])
        ]))

        dmat = zeros(UInt16, (4, 4))
        fill_distance_matrix!(G, dmat)
        links = identify_potential_missing_links(G, dmat, 2, 5)
        sortlinks!(links)

        @test length(links) == 2
        # only one link, but we found it in both directions
        @test links[1] == reverse(links[2])

        @test links[1] == CandidateLink(
            fr_edge_src=VertexID(1),
            fr_edge_tgt=VertexID(2),
            fr_dist_from_start=UInt16(6), 
            fr_dist_to_end=UInt16(4),
            to_edge_src=VertexID(3),
            to_edge_tgt=VertexID(4),
            to_dist_from_start=UInt16(6), # Actually 6.414
            to_dist_to_end=UInt16(0),
            geographic_length_m=one(UInt16),
            network_length_m=typemax(UInt16)
        )
    end

    @testset "Multiple links - Graph = <" begin
        G = graph_from_gdal(DataFrame(:geom=>[
            AG.createlinestring([[0.0, 0.0], [5.0, 0.0]]),
            AG.createlinestring([[5.0, 1.0], [0.0, 1.0]]),
            AG.createlinestring([[7.0, 0.0], [6.0, 1.0], [7.0, 2.0]])
        ]))

        dmat = zeros(UInt16, (6, 6))
        fill_distance_matrix!(G, dmat)
        links = identify_potential_missing_links(G, dmat, 2, 5)
        sortlinks!(links)

        @test length(links) == 6
        @test links[2] == reverse(links[5])
        @test links[4] == reverse(links[6])

        @test links[2] == CandidateLink(
            fr_edge_src=VertexID(1),
            fr_edge_tgt=VertexID(2),
            fr_dist_from_start=UInt16(5), 
            fr_dist_to_end=UInt16(0),
            to_edge_src=VertexID(5),
            to_edge_tgt=VertexID(6),
            to_dist_from_start=UInt16(1), # Actually 1.414
            to_dist_to_end=UInt16(2),
            geographic_length_m=one(UInt16),
            network_length_m=typemax(UInt16)
        )

        @test links[4] == CandidateLink(
            fr_edge_src=VertexID(3),
            fr_edge_tgt=VertexID(4),
            fr_dist_from_start=UInt16(0), 
            fr_dist_to_end=UInt16(5),
            to_edge_src=VertexID(5),
            to_edge_tgt=VertexID(6),
            to_dist_from_start=UInt16(1), # Actually 1.414
            to_dist_to_end=UInt16(2),
            geographic_length_m=one(UInt16),
            network_length_m=typemax(UInt16)
        )

        # there is also a link 1-2 -> 3-4, but offsets are implementation-defined since the lines are parallel
    end

    @testset "Connections between already connected edges that are far apart; backwards edges" begin
        # in this test we look at this network:
        #        B------------------------------C
        #        |                              |
        #        |                              |
        #        |                              |
        #        |                              |
        #        |                              |
        #        |                              |
        #        |                              |
        #        |                              |
        #        |_________A   D________________|

        # A and D are 50m geographic distance, 3950m network distance, so should get a link
        # DC also has reversed geometry from the way it is specified, so this checks that distances are correct

        G = graph_from_gdal(DataFrame(:geom=>[
            AG.createlinestring([[475.0, 1000.0], [0.0, 1000.0], [0.0, 0.0]]), # AB
            AG.createlinestring([[0.0, 0.0], [1000.0, 0.0]]), # BC
            AG.createlinestring([[525.0, 1000.0], [1000.0, 1000.0], [1000.0, 0.0]]) # DC - will be reversed
        ]); max_edge_length=50_000)
        dmat = zeros(UInt16, (nv(G), nv(G)))
        fill_distance_matrix!(G, dmat)
        links = identify_potential_missing_links(G, dmat, 100, 1000)

        @test length(links) == 2
        @test links[1] == reverse(links[2])
        @test links[1] == CandidateLink(
            fr_edge_src=VertexID(1),
            fr_edge_tgt=VertexID(2),
            fr_dist_from_start=UInt16(0), 
            fr_dist_to_end=UInt16(1475),
            to_edge_src=VertexID(3),
            to_edge_tgt=VertexID(4),
            to_dist_from_start=UInt16(1475), # Geometry will be reversed
            to_dist_to_end=UInt16(0),
            geographic_length_m=UInt16(50),
            network_length_m=UInt16(3950)
        )
    end

    @testset "Connected links should not be identified" begin
        # Network looks like this
        #   ___B
        #  /   |
        #  |   |
        #  A   C
        # AC is 40m, but ABC is

        G = graph_from_gdal(DataFrame(:geom=>[
            AG.createlinestring([[0.0, 100.0], [0.0, 20.0], [20.0, 0.0], [40.0, 0.0]]), # AB
            AG.createlinestring([[40.0, 0.0], [40.0, 100.0]]) # BC
        ]); max_edge_length=50_000)
        dmat = zeros(UInt16, (nv(G), nv(G)))
        fill_distance_matrix!(G, dmat)
        links = identify_potential_missing_links(G, dmat, 100, 1000)
        @test isempty(links)
    end

    @testset "Edge breaking" begin
        # Network looks like this
        #       __B_____C__
        #       |         |
        #       |         |
        #       |         |
        #       |         |
        #       |         |
        #       |         |
        #       |         |
        #       |         |
        #       |_A     D_|
        # The algorithm would ordinarily not find AD, because the closest point between them
        # is at BC, but AB and CD get broken into two links each because of their length.
        G = graph_from_gdal(DataFrame(:geom=>[
            # NB AB and CD are each 400m
            AG.createlinestring([[10.0, 370.0], [0.0, 370.0], [0.0, 0.0], [20.0, 0.0]]), # AB
            AG.createlinestring([[20.0, 0.0], [30.0, 0.0]]), # BC
            AG.createlinestring([[30.0, 0.0], [50.0, 0.0], [50.0, 370.0], [40.0, 370.0]])
        ]); max_edge_length=250)
        dmat = zeros(UInt16, (nv(G), nv(G)))
        fill_distance_matrix!(G, dmat)
        links = identify_potential_missing_links(G, dmat, 100, 800)

        @test length(links) == 2
        @test links[1] == reverse(links[2])

        @test links[1] == CandidateLink(
            fr_edge_src=VertexID(1),
            fr_edge_tgt=VertexID(2),
            fr_dist_from_start=UInt16(0), 
            fr_dist_to_end=UInt16(200), # 200m split edge
            to_edge_src=VertexID(5),
            to_edge_tgt=VertexID(6),
            to_dist_from_start=UInt16(200), # Geometry will be reversed
            to_dist_to_end=UInt16(0),
            geographic_length_m=UInt16(30),
            network_length_m=UInt16(810)
        )
    end

    @testset "Unbroken long edges" begin
        # this tests a known issue in the algorithm - we have a network that looks like this:
        # 
        #
        #      B----------C
        #     /            \
        #    /              \
        #   /                \
        #  /                  \
        # /                    \
        # A                     D
        # 
        # A link between A and D would be optimal, as this is less than 100m and ABCD is more than 1000 m. However,
        # when AB and CD are compared, the candidate link found will be BC, which is already connected.
        # So we don't find the link in this corner case. However, we do split long edges in the graph build to 
        # minimize the impact of this issue.

        G = graph_from_gdal(DataFrame(:geom=>[
            AG.createlinestring([[0.0, 1000.0], [2.0, 0.0]]),
            AG.createlinestring([[2.0, 0.0], [3.0, 0.0]]),
            AG.createlinestring([[3.0, 0.0], [4.0, 1000.0]])
        ]); max_edge_length=50_000)
        dmat = zeros(UInt16, (nv(G), nv(G)))
        fill_distance_matrix!(G, dmat)
        links = identify_potential_missing_links(G, dmat, 100, 1000)
        @test isempty(links)
    end
end