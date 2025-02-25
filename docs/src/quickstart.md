# Quickstart

In this page, we will quickly work through the Julia code necessary to load a network dataset, assign weights to it, and identify and score missing links. I recommend running this code in an interactive, notebook style environment; I use [Quarto](https://quarto.org) and [Visual Studio Code](https://code.visualstudio.com/). More details about all of the functions from the MissingLinks package that are used herein are available in the [API Documentation](@ref).

## Load packages

First, we need to load some packages. The Missing Links code is distributed as a Julia package, so we need to load MissingLinks. We also need to load the GeoDataFrames package for
spatial data manipulation, the Plots library for visualization, the Graphs library for working with graph data structures, the StatsBase library for statistics, and the Mmap package for working with memory-mapped files. Most of these packages are available from the Julia general registry, and can be installed [with the Pkg-mode add command](https://docs.julialang.org/en/v1/stdlib/Pkg/). Mmap is built into Julia and does not need to be installed. MissingLinks, however, is not in the general registry; to install it you should enter Pkg-mode (press `]`) and then type `add https://github.com/mattwigway/MissingLinks.jl`. In a nutshell, you will run this code, after entering package mode by pressing `]`. You'll run this code in the REPL (in VSCode, choose View -> Command Palette and search for "Julia: Start REPL"). All of the other code you'll run should be in your notebook interface so you can save it. While it is optional, I do recommend [creating a Julia environment](https://pkgdocs.julialang.org/v1/environments/) to keep MissingLinks code separate from any other Julia projects you may work on.

```
add GeoDataFrames Plots Graphs StatsBase
add https://github.com/mattwigway/MissingLinks.jl
```

Press backspace to exit package mode.

Once packages are installed, we are ready to load them.

```@example main
using MissingLinks, GeoDataFrames, Plots, Mmap, Graphs, StatsBase
```

If you're using Quarto and VSCode, place this code between triple-backticks to create a Julia cell, then hover over it and click "run cell":

````
```{julia}
using MissingLinks, GeoDataFrames, Plots, Mmap, Graphs, StatsBase
```
````

## Loading data

Next, we will read our data. To keep things fast for this example, we are working with a tiny segment of the Charlotte pedestrian network, in the [Northlake area](https://www.openstreetmap.org/#map=15/35.34618/-80.85734). We need data on the pedestrian network, and origins and destinations. In this case we are using residential parcels as origins, and a selection of commercial destinations in the Northlake area digitized from OpenStreetMap as destinations. The data are included with the MissingLinks package; running

```@example main
MissingLinks.get_example_data()
```

will create a folder `example_data` in the current working directory, containing the data.

All of our layers are already projected to EPSG:32119 (State Plane North Carolina, meters). I recommend using a meter-based coordinate system. The system will not work correctly with a geographic coordinate system, and has not been tested with feet.

```@example main
sidewalks = GeoDataFrames.read("example_data/Northlake.gpkg")
parcels = GeoDataFrames.read("example_data/Northlake_Parcels.gpkg")
destinations = GeoDataFrames.read("example_data/Northlake_Destinations.gpkg")
```

```@setup main
rm("example_data", recursive=true)
```

We can plot what we read in:

```@example main

# sidewalks have an elevation component that we need to remove to plot
sidewalks.geom = MissingLinks.remove_elevation!.(sidewalks.geom)


plot(parcels.geom, color="#0f0", legend=false, aspect_ratio=:equal)
plot!(sidewalks.geom, color="#f00")
plot!(destinations.geometry, color="#f00")
```

## Building the network

Next, we need to convert our network GIS layer into a mathematical "graph"—a data structure of nodes and edges representing intersections and streets. This network dataset is already "fully noded"—i.e. all intersections occur at the ends of line segments. The tool can also work with "semi-noded" data—where all intersections are at the end of one line segment, but may be in the middle of another. To use that kind of data, run the preprocessing function [`semi_to_fully_noded`](@ref) before building the graph.

To build the graph, we will use the [`graph_from_gdal`](@ref) function, which will convert GeoDataFrames into a graph.

```@example main
graph = graph_from_gdal(sidewalks)
```

It is possible to specify more than one layer (for example, if you have sidewalks and crosswalks in separate layers) by separating the layers with commas. It is also possible to set a tolerance for when nodes are considered coincident. However, I recommend leaving this at the default small value and using [`add_short_edges!`](@ref) to close such gaps, to avoid issues that are described in the linked documentation. Similarly, [`remove_tiny_islands`](@ref) can be used to remove small "islands" in the graph that may be data errors. There are also tools for troubleshooting the graph, namely [`find_dead_ends`](@ref) to find locations that are dead ends, and [`find_disconnected_crossings`](@ref) to find places where sidewalks cross without intersecting. While both of these things will be present in a correct graph, they may also represent data errors.

## Assigning weights to nodes

Next, we need to assign weights to the nodes. We will create two sets of weights, for origins and destinations. For origins, the weights are just the number of units on the parcel. The weight of each parcel will get divided among the graph edges within 20 meters, and then the weight for each edge will be evenly divided between the nodes it connects. `[:,1]` retrieves the first column of the weights (if multiple weight columns are specified instead of just units, there will be one column in the output for each column specified. This avoids re-doing computationally-intensive geographic calculations if there are more weight columns—for instance, if both our origins and our destinations were coded in the parcel layer, or if we wanted to measure accessibility from multiple types of origin).

```@example main
origin_weights = create_graph_weights(graph, parcels, [:units], 20)[:,1]
```

We can look at what proportion of the units are served by sidewalks and therefore were assigned to the network. Some units will not be assigned to any edge and therefore will not count towards origin_weights.

```@example main
sum(origin_weights) / sum(parcels.units)
```

We can similarly create our destination weights. There is no weight column here—all of the destinations are equally weighted. However, we need a weight column so we just create a column of all ones. Since destinations are further from the sidewalk network due to parking lots, we use a 100 meter distance threshold. In the paper we use parcels and a 20 meter threshold, on the assumption that the edges of parcels will be close to the road, but here we have point data on destinations.

```@example main
destinations.one .= 1
destination_weights = create_graph_weights(graph, destinations, [:one], 100)[:, 1]
sum(destination_weights) / sum(destinations.one)
```

## Creating the distance matrix

A point-to-point distance matrix is used extensively throughout the algorithm for identify and scoring links. We store distances in this matrix as UInt16 integer meters, to save memory (each UInt16 can represent values from 0-65535, and takes only two bytes of memory, as opposed to four or eight for traditional integers and longs). For a network as small as the one we're working with here, it would be possible to store this distance matrix in memory. However, in real-world applications, often the matrix will be larger than available memory, so storing it on disk and using memory-mapping of the disk file is recommended. Memory-mapping lets the computer treat a portion of the disk as if it were memory, and the CPU and the operating system work to make sure that the data are available quickly when needed.

```@example main
matrix_file = open("quickstart_distances.bin", "w+")
matrix = Mmap.mmap(matrix_file, Matrix{UInt16}, (nv(graph), nv(graph)))
# this just makes sure the matrix is filled with 0s before we begin. This should not make any real difference as all values should be overwritten but is just a good practice.
fill!(matrix, 0)

fill_distance_matrix!(graph, matrix)
```

## Link identification

We are now ready to identify missing links. The [`identify_potential_missing_links`](@ref) will do this and return a vector of places in the network that meet the specified criteria (in this case, no longer than 100 meters, and with a network distance of 1000 meters or more; you can change these values in the code below).

```@example main
all_links = identify_potential_missing_links(graph, matrix, 100, 1000)
```

We can plot the links identified (shown here in blue):

```@example main
all_links_geo = links_to_gis(graph, all_links)
plot(sidewalks.geom, color="#000", legend=false, aspect_ratio=:equal)
plot!(all_links_geo.geometry, color="#4B9CD3")
```

## Link deduplication

Clearly, most of these links are essentially duplicates. The next step is to deduplicate them, by merging links where both ends are within 100m of one another (configurable by changing the number below).

```@example main
links = deduplicate_links(all_links, matrix, 100)
```

We can now plot these deduplicated links:

```@example main
links_geo = links_to_gis(graph, links)
plot(sidewalks.geom, color="#000", legend=false, aspect_ratio=:equal)
# set line width wider to make the links easier to see
plot!(links_geo.geometry, color="#4B9CD3", lw=2)
```

## Link scoring

The final step is to score the identified links, using the [`score_links`](@ref) function. For this we need to specify what our origin and destination weights are, as well as our distance decay function. Here, I use a 1-mile (1609 meter) cutoff. This will return a vector with the score for each link. We have to specify both the distance decay function (first argument) and the point at which that function goes to zero (last argument); for a smooth decay function, I recommend instead specifying a piecewise function that goes to zero when the result is small enough that additional destinations do not materially affect accessibility.

```@example main
link_scores = score_links(d -> d < 1609, links, matrix, origin_weights, destination_weights, 1609)
```

## Saving links to GIS

Generally, you will want to export your results to GIS for further analysis, including the scores. You can create a GeoDataFrame including the scores like this:

```@example main
links_geo = links_to_gis(graph, links, :score=>link_scores)
```

You can then write this to GIS file like this:

```julia
GeoDataFrames.write("links.gpkg", links_geo)
```

We can also use it in further plotting. For example, this plot shows the top five links (the `competerank` function just ranks values largest to smallest).

```@example main
plot(sidewalks.geom, color="#000", legend=false, aspect_ratio=:equal)
# set line width wider to make the links easier to see
plot!(links_geo[competerank(links_geo.score) .<= 5, :geometry], color="#4B9CD3", lw=2)
```

As we might expect, the new links that connect the completely disconnected neighborhood on the west side are most valuable. Those that shorten trips in already connected areas are less valuable.
