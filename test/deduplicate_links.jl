# Test deduplication of links
@testsnippet InitializeDedupe begin
    import ArchGDAL as AG
    import MissingLinks: deduplicate_links, graph_from_gdal, calculate_distances, CandidateLink, VertexID
    import DataFrames: DataFrame
    import Graphs: nv


    # the network we have is just two parallel lines split into a number of edges
    # We then have a number of candidate links between them to test deduplication
    # 1 --50m-- 2 --50m-- 3 --200m-- 4
    # 5 --50m-- 6 --50m-- 8 --200m-- 7 Note that 78 is backwards

    G = graph_from_gdal(DataFrame(:geom => [
        AG.createlinestring([[0.0, 0.0], [50.0, 0.0]]), # 12
        AG.createlinestring([[50.0, 0.0], [100.0, 0.0]]), # 23
        AG.createlinestring([[100.0, 0.0], [300.0, 0.0]]), # 34

        AG.createlinestring([[0.0, 10.0], [50.0, 10.0]]), # 56
        AG.createlinestring([[300.0, 10.0], [100.0, 10.0]]), # 78
        AG.createlinestring([[50.0, 10.0], [100.0, 10.0]]) # 67
    ]))

    dmat = calculate_distances(G)
end

@testitem "Basic deduplication" setup=[InitializeDedupe] begin
    # Here we have three links, two on the left two edges that are duplicates, and one further right that is not
    # Since the link chosen is based on geographic length, we fudge the geographic length a bit, since otherwise
    # they'd all be equivalent.
    links = [
        CandidateLink(
            fr_edge_src=VertexID(1),
            fr_edge_tgt=VertexID(2),
            fr_dist_from_start=25,
            fr_dist_to_end=25,
            to_edge_src=VertexID(5),
            to_edge_tgt=VertexID(6),
            to_dist_from_start=25,
            to_dist_to_end=25,
            geographic_length_m=9, # shortest of this SOI
            network_length_m=missing
        ),

        CandidateLink(
            fr_edge_src=VertexID(2),
            fr_edge_tgt=VertexID(3),
            fr_dist_from_start=25,
            fr_dist_to_end=25,
            to_edge_src=VertexID(6),
            to_edge_tgt=VertexID(8),
            to_dist_from_start=25,
            to_dist_to_end=25,
            geographic_length_m=10, #longer
            network_length_m=missing
        ),

        # this one is not a dupe
        CandidateLink(
            fr_edge_src=VertexID(3),
            fr_edge_tgt=VertexID(4),
            fr_dist_from_start=0, # start is in previous SOI
            fr_dist_to_end=200,
            to_edge_src=VertexID(7),
            to_edge_tgt=VertexID(8),
            to_dist_from_start=0, # end is not, because link is backwards
            to_dist_to_end=200,
            geographic_length_m=10,
            network_length_m=missing
        )
    ]

    # we will always find the backwards versions as well
    links = [links..., reverse.(links)...]

    dedupe = deduplicate_links(links, dmat, 100)

    @test length(dedupe) == 2
    @test dedupe[1] === links[1]
    @test dedupe[2] === links[3]
end

@testitem "Does not always return first link" setup=[InitializeDedupe] begin
    # we add another link that's even shorter to the first sphere of influence, and make sure it is the
    # one that gets retained.
    links = [
        CandidateLink(
            fr_edge_src=VertexID(1),
            fr_edge_tgt=VertexID(2),
            fr_dist_from_start=25,
            fr_dist_to_end=25,
            to_edge_src=VertexID(5),
            to_edge_tgt=VertexID(6),
            to_dist_from_start=25,
            to_dist_to_end=25,
            geographic_length_m=9,
            network_length_m=missing
        ),

        CandidateLink(
            fr_edge_src=VertexID(2),
            fr_edge_tgt=VertexID(3),
            fr_dist_from_start=25,
            fr_dist_to_end=25,
            to_edge_src=VertexID(6),
            to_edge_tgt=VertexID(8),
            to_dist_from_start=25,
            to_dist_to_end=25,
            geographic_length_m=10, #longer
            network_length_m=missing
        ),

        CandidateLink(
            fr_edge_src=VertexID(1),
            fr_edge_tgt=VertexID(2),
            fr_dist_from_start=25,
            fr_dist_to_end=25,
            to_edge_src=VertexID(6),
            to_edge_tgt=VertexID(8),
            to_dist_from_start=25,
            to_dist_to_end=25,
            geographic_length_m=8, # shortest of this SOI
            network_length_m=missing
        ),

        # this one is not a dupe
        CandidateLink(
            fr_edge_src=VertexID(3),
            fr_edge_tgt=VertexID(4),
            fr_dist_from_start=0, # start is in previous SOI
            fr_dist_to_end=200,
            to_edge_src=VertexID(7),
            to_edge_tgt=VertexID(8),
            to_dist_from_start=0, # end is not, because link is backwards
            to_dist_to_end=200,
            geographic_length_m=10,
            network_length_m=missing
        )
    ]

    # we will always find the backwards versions as well
    links = [links..., reverse.(links)...]

    dedupe = deduplicate_links(links, dmat, 100)

    @test length(dedupe) == 2
    @test dedupe[1] === links[3]
    @test dedupe[2] === links[4]
end

@testitem "Connections to same edges" setup=[InitializeDedupe] begin
    links = [
        # both of these originate on 34, one goes to 56, the other to 68,
        # they should be treated as duplicates
        CandidateLink(
            fr_edge_src=VertexID(3),
            fr_edge_tgt=VertexID(4),
            fr_dist_from_start=100,
            fr_dist_to_end=100,
            to_edge_src=VertexID(5),
            to_edge_tgt=VertexID(6),
            to_dist_from_start=50,
            to_dist_to_end=0,
            geographic_length_m=9, # shortest of this SOI
            network_length_m=missing
        ),
        CandidateLink(
            fr_edge_src=VertexID(3),
            fr_edge_tgt=VertexID(4),
            fr_dist_from_start=100,
            fr_dist_to_end=100,
            to_edge_src=VertexID(6),
            to_edge_tgt=VertexID(8),
            to_dist_from_start=0,
            to_dist_to_end=50,
            geographic_length_m=10,
            network_length_m=missing
        )
    ]

    links = [links..., reverse.(links)...]

    # threshold set at 75 b/c we want to make sure that the distance measurement for the start of
    # the links is not requiring a walk to the end of the street and back.
    dedupe = deduplicate_links(links, dmat, 75)

    @test length(dedupe) == 1
    @test dedupe[1] === links[1]
end
