@testitem "Aqua" begin
    import Aqua
    Aqua.test_all(MissingLinks;
        # used in notebook code or scripts
        stale_deps=(ignore=[:ArgParse, :Gadfly, :Humanize, :StatsBase, :Cairo, :Fontconfig],)
    )
end
