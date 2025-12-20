# datasets used in testing, but which need to be referred to from within MissingLinks itself
module TestData
    import Pkg.Artifacts: @artifact_str

    huntersville_osm() = artifact"huntersville_osm"
end