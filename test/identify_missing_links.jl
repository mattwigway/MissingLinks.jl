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
    import MissingLinks: graph_from_gdal, identify_potential_missing_links, fill_matrix!, CandidateLink
    import DataFrames: DataFrame
    import ArchGDAL as AG
    import StructEquality: @struct_hash_equal

    # make sure we can compare two links by value rather than identity
    @struct_hash_equal CandidateLink

    # Sort the links into a stable order. There can only be one link between an edge pair, so there will
    # be no ties.
    sortlinks!(links) = sort!(links, by=l -> (l.fr_edge_src, l.fr_edge_tgt, l.to_edge_src, l.to_edge_tgt))

    # This first graph is very simple, it looks like > | 
    # We should find a link from the middle of one edge to the middle of the other
    gebar = graph_from_gdal(DataFrame(:geom =>[
        AG.createlinestring([[0.0, 0.0], [1.0, 1.0], [0.0, 2.0]]),
        AG.createlinestring([[2.0, 0.0], [2.0, 2.0]])
    ]))

    dmat = zeros(UInt16, (4, 4))
    fill_matrix!(gebar, dmat)
    links = identify_potential_missing_links(gebar, dmat, 2, 5)
    sortlinks!(links)

    @test length(links) == 2
    # only one link, but we found it in both directions
    @test links[1] == reverse(links[2])

    @test links[1] == CandidateLink(
        fr_edge_src=1,
        fr_edge_tgt=2,
        fr_dist_from_start=one(UInt16), # actually √2 but is rounded
        fr_dist_to_end=UInt16(2),
        to_edge_src=3,
        to_edge_tgt=4,
        to_dist_from_start=one(UInt16),
        to_dist_to_end=one(UInt16),
        geographic_length_m=one(UInt16),
        network_length_m=typemax(UInt16)
    )

    # Now we check one where the snapping is to the ends
    # Network looks like / \
    fbslash = graph_from_gdal(DataFrame(:geom=>[
        AG.createlinestring([[0.0, 1.0], [1.0, 0.0]]),
        AG.createlinestring([[2.0, 0.0], [3.0, 1.0]])
    ]))

    dmat = zeros(UInt16, (4, 4))
    fill_matrix!(fbslash, dmat)
    links = identify_potential_missing_links(fbslash, dmat, 2, 5)
    sortlinks!(links)

    @test length(links) == 2
    # only one link, but we found it in both directions
    @test links[1] == reverse(links[2])

    @test links[1] == CandidateLink(
        fr_edge_src=1,
        fr_edge_tgt=2,
        fr_dist_from_start=one(UInt16),
        fr_dist_to_end=zero(UInt16),
        to_edge_src=3,
        to_edge_tgt=4,
        to_dist_from_start=zero(UInt16),
        to_dist_to_end=one(UInt16),
        geographic_length_m=one(UInt16),
        network_length_m=typemax(UInt16)
    )

    # One end and one middle
    # Network looks like this:
    # ----------------
    # ______/

    G = graph_from_gdal(DataFrame(:geom=>[
        AG.createlinestring([[0.0, 0.0], [10.0, 0.0]]),
        AG.createlinestring([[0.0, 2.0], [5.0, 2.0], [6.0, 1.0]])
    ]))

    dmat = zeros(UInt16, (4, 4))
    fill_matrix!(G, dmat)
    links = identify_potential_missing_links(G, dmat, 2, 5)
    sortlinks!(links)

    @test length(links) == 2
    # only one link, but we found it in both directions
    @test links[1] == reverse(links[2])

    @test links[1] == CandidateLink(
        fr_edge_src=1,
        fr_edge_tgt=2,
        fr_dist_from_start=UInt16(6), 
        fr_dist_to_end=UInt16(4),
        to_edge_src=3,
        to_edge_tgt=4,
        to_dist_from_start=UInt16(6), # Actually 6.414
        to_dist_to_end=UInt16(0),
        geographic_length_m=one(UInt16),
        network_length_m=typemax(UInt16)
    )

    # Multiple links - network looks like = <
    G = graph_from_gdal(DataFrame(:geom=>[
        AG.createlinestring([[0.0, 0.0], [5.0, 0.0]]),
        AG.createlinestring([[5.0, 1.0], [0.0, 1.0]]),
        AG.createlinestring([[7.0, 0.0], [6.0, 1.0], [7.0, 2.0]])
    ]))

    dmat = zeros(UInt16, (6, 6))
    fill_matrix!(G, dmat)
    links = identify_potential_missing_links(G, dmat, 2, 5)
    sortlinks!(links)

    @test length(links) == 6
    @test links[2] == reverse(links[5])
    @test links[4] == reverse(links[6])

    @test links[2] == CandidateLink(
        fr_edge_src=1,
        fr_edge_tgt=2,
        fr_dist_from_start=UInt16(5), 
        fr_dist_to_end=UInt16(0),
        to_edge_src=5,
        to_edge_tgt=6,
        to_dist_from_start=UInt16(1), # Actually 1.414
        to_dist_to_end=UInt16(2),
        geographic_length_m=one(UInt16),
        network_length_m=typemax(UInt16)
    )

    @test links[4] == CandidateLink(
        fr_edge_src=3,
        fr_edge_tgt=4,
        fr_dist_from_start=UInt16(0), 
        fr_dist_to_end=UInt16(5),
        to_edge_src=5,
        to_edge_tgt=6,
        to_dist_from_start=UInt16(1), # Actually 1.414
        to_dist_to_end=UInt16(2),
        geographic_length_m=one(UInt16),
        network_length_m=typemax(UInt16)
    )

    # there is also a link 1-2 -> 3-4, but offsets are implementation-defined since the lines are parallel

    # add'l tests: network where the links are connected distantly and missing link should be identified, network
    # where links are connected and missing link should not be identified.
end