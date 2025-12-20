@testitem "Partition real world graph" begin
    import GeoFormatTypes as GFT
    import MetaGraphsNext: labels
    import Graphs: nv

    Base.isless(a::MissingLinks.CandidateLink, b::MissingLinks.CandidateLink) =
        (a.fr_edge_src, a.fr_edge_tgt, a.to_edge_src, a.to_edge_tgt) <
        (b.fr_edge_src, b.fr_edge_tgt, b.to_edge_src, b.to_edge_tgt)

    remove_distance(l::MissingLinks.CandidateLink{T}) where T = MissingLinks.CandidateLink{T}(
        fr_edge_src=l.fr_edge_src,
        fr_edge_tgt=l.fr_edge_tgt,
        fr_dist_from_start=l.fr_dist_from_start,
        fr_dist_to_end=l.fr_dist_to_end,
        to_edge_src=l.to_edge_src,
        to_edge_tgt=l.to_edge_tgt,
        to_dist_from_start=l.to_dist_from_start,
        to_dist_to_end=l.to_dist_to_end,
        geographic_length_m=l.geographic_length_m,
        network_length_m=zero(T)
    )

    # this is an area of Huntersville NC with lots of disconnected locations
    G = graph_from_osm(joinpath(MissingLinks.TestData.huntersville_osm(), "huntersville.osm.pbf"), TraversalPermissionSettings(), GFT.EPSG(32119))
    
    # we want the weights to be "deterministically random" i.e. not show any pattern but be the same every time
    # the test is run, so we use the LSB of the OSM ID
    weights = [v.id & 0xff for v in labels(G)]

    # identify and score links with monolithic graph
    dmat = Matrix{UInt16}(undef, nv(G), nv(G))
    fill_distance_matrix!(G, dmat)
    links = identify_potential_missing_links(G, dmat, 100, 1000)
    # We are intentionally not deduplicating. Since the deduplication is a heuristic that depends on other links
    # in the graph, partitioning will mean that some links are not exactly perfectly deduplicated identically.
    scores = score_links(x -> x < 1000, G, links, dmat, fill(1, nv(G)), weights, 1000)

    # network distances are not reliable in partitioned graph, as the old network distance may have been
    # much larger than a partition
    links = collect(map(remove_distance, links))

    # ensure they are the same up to ordering
    # sort links instead of scores in case of ties
    sorter = sortperm(links)
    links = links[sorter]
    scores = scores[sorter]

    # partition graph
    Gs = MissingLinks.partition(G, 3, 3, 1100)

    plinks, pscores = MissingLinks.identify_and_score(
        x -> x < 1000,
        G,
        Gs,
        cutoff_distance = 1000,
        origin_weights = fill(1, nv(G)),
        dest_weights = weights,
        deduplicate = false
    )

    # since we've disabled deduplication, we'll find all links in both directions above
    # but merge_links will only return links in one direction
    # so make the reverse of every link

    plinks = [plinks..., reverse.(plinks)...]
    pscores = [pscores..., pscores...]

    sorter = sortperm(plinks)
    plinks = plinks[sorter]
    pscores = pscores[sorter]

    @test length(links) == 7430
    @test length(plinks) == 7430
    @test all(links .== plinks)
    @test all(scores .≈ pscores)
end

