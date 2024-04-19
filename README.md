# MissingLinks

This code supports finding and scoring missing links in transportation networks (primarily pedestrian networks; some changes would be needed for other types of networks, because it assumes the network is undirected, and currently only supports distances up to 65 km). This repository contains a [Julia](https://julialang.org) package to support finding missing links, as well as code for the initial application in Charlotte, North Carolina, USA.

If you want to use this with your own network, the easiest way is to install it as a Julia package. It's currently not registered in the Julia general repository; you can install it with `]add https://github.com/mattwigway/MissingLinks.jl`.

There are several steps; to see how to run them, see `notebooks/run_analysis.qmd`. Generally, they are:

1. Build a network using MetaGraphs.jl. The `graph_from_gdal` function will build a graph from a fully-noded GIS network. There is also `semi_to_fully_noded` to convert a semi-noded (intersections occur where the end of one line intersects any part of another line) GIS network to a fully-noded GIS network that can then be used by `graph_to_gdal`
1. Clean up the network—the two main functions are `add_short_edges!` which will add edges any place where there are very short disconnections (configurable)–e.g. a gap of a meter or two is probably not an actual missing link but a data error. `remove_tiny_islands` will remove small components of the network which often represent data errors. The parameters for both of these methods need to be carefully considered to handle data errors without connecting or removing places that are actual missing links.
1. Create a distance matrix - see the notebook for how to create an mmapped matrix if you don't have enough memory for the full matrix. `UInt16` is the recommended type. Then use `fill_matrix!` to fill in the distance matrix.
1. Identify disconnected locations - the `identify_disconnected_locations` function does this.
1. Deduplicate the found links—many links are basically identical, we do not want to consider all of them. The deduplication algorithm clusters nearby links and retains only the shortest one. `deduplicate_links` does this.
1. Load origin and destination weights—ultimately weights need to be assigned to edges in the network. `create_graph_weights` assigns weights from GIS files to nearby edges.
1. Score the links - `score_links` does this.

## Running the code

The code is written in [Julia](https://julialang.org) for performance, and uses a [Quarto](https://quarto.org) notebook format. To use it, [install Julia](https://julialang.org/downloads/); I also recommend installing [Visual Studio Code](https://code.visualstudio.com/) to edit Quarto files across many languages, as well as the [Quarto](https://marketplace.visualstudio.com/items?itemName=quarto.quarto) and [Julia](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia) VSCode extensions.

It will use multiple cores if you have them to speed up routing and scoring. To enable this, in the VSCode settings, search for "julia num threads", click "edit in settings.json", and edit the selected line to read `"julia.NumThreads": "auto"`.

Open the `MissingLinks.jl` folder in VSCode, and press Cmd/Ctrl-shift-P to get the "command palette" and search for and select "Julia: Start REPL". A Julia REPL (console) should appear. Type `]` to enter the package manager. The prompt should change to say `(MissingLinks.jl) pkg>`. (If it says something like `(@1.9) pkg>` type `activate .` to activate the "environment" for the MissingLinks.jl folder, which contains the dependencies needed to run the model.) Type `instantiate` to download all the necessary packages, and then press backspace to exit the package manager. Type `Threads.nthreads()` and press enter; confirm you get a number greater than 1 to indicate that multithreaded processing is enabled.

You now have Julia and all dependencies installed, and should not need to repeat the above steps unless you get a new computer, upgrade Julia, or change the dependencies of the project.

The datasets used in the Charlotte project are too large to store in Git, so they are stored in the OneDrive folder shared with all project collaborators. This will potentially be in a different location on each person's computer. To tell Julia where to find the data, copy the `config.toml.template` file to `config.toml`, and enter the path to the `Missing links project/data` folder.

To run the analysis, open the `notebooks/run_analysis.qmd file`. You can run cells by clicking "Run cell" or pressing "Cmd/Ctrl-shift-enter". The model requires about 26 GB of free disk space to run (for the mmapped distance matrix).