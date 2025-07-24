@testitem "Aqua" begin
    import Aqua, MissingLinks
    Aqua.test_all(MissingLinks;
        # used in notebook code or scripts
        stale_deps=(ignore=[:ArgParse, :Gadfly, :Humanize, :StatsBase, :Cairo, :Fontconfig],),
        
        # EdgeData is an alias of a specific NamedTuple type. Now technically defining Base.isapprox
        # is type piracy, because we don't own NamedTuple, but the chance that there is another
        # NamedTuple with the same members and types is basically nil.

        # I tried putting MissingLinks.EdgeData as well as the full definition of the NamedTuple type
        # here and it did not work, the test still claimed type piracy.
        piracies=(treat_as_own=[NamedTuple],)
    )
end
