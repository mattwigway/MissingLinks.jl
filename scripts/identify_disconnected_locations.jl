# This script identifies locations where ends of streets are disconnected
# i.e. vertices that connect to only one edge

using ArgParse, TOML, ArchGDAL, GeoDataFrames, MissingLinks, DataFrames, Logging, Graphs, MetaGraphsNext
import GeoFormatTypes as GFT

config = TOML.parsefile(joinpath(dirname(@__FILE__), "..", "config.toml"))

function main(raw_args)
    s = ArgParseSettings()

    @add_arg_table! s begin
        "output"
            help = "output geopackage"
        "input"
            help = "input files"
            nargs='+'
    end

    args = parse_args(raw_args, s)

    @info "Reading files:" args["input"]
    data = map(GeoDataFrames.read.(args["input"])) do file
        if "geometry" ∈ names(file)
            rename!(file, :geometry=>:geom) # ?
        end

        file.geom = ArchGDAL.lineargeom.(file.geom)

        file
    end

    # the input files are semi-noded, convert to fully noded
    @info "Converting semi-noded data to fully noded"
    noded_inputs = semi_to_fully_noded(data..., snap_tolerance=config["snap_tolerance"], split_tolerance=config["split_tolerance"])

    @info "building graph"
    G = MissingLinks.graph_from_gdal(noded_inputs)

    add_short_edges!(G, config["snap_tolerance"] + config["split_tolerance"])

    @info "Identifying degree-1 nodes"
    result = ArchGDAL.IGeometry{ArchGDAL.wkbPoint}[]

    for node in 1:nv(G)
        if Graphs.degree(G, node) == 1
            position = G[label_for(G, node)]
            push!(result, ArchGDAL.createpoint(collect(position)))
        end
    end

    @info "Found $(length(result)) degree-one nodes"

    resultdf = DataFrame(:geometry => result)
    GeoDataFrames.write(args["output"], resultdf, crs=GFT.EPSG(config["coord_system"]))
end

main(ARGS)