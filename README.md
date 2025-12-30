# MissingLinks

![Automated test status](https://github.com/mattwigway/MissingLinks.jl/actions/workflows/ci.yml/badge.svg) [![Documentation](https://img.shields.io/badge/docs-latest-4b9cd3)](https://docs.missinglinks.city)

This code supports finding and scoring missing links in transportation networks (primarily pedestrian networks; some changes would be needed for other types of networks, because it assumes the network is undirected, and currently only supports distances up to 65 km). This repository contains a [Julia](https://julialang.org) package to support finding missing links, as well as code for the initial application in Charlotte, North Carolina, USA. The algorithm is described in [Bhagat-Conway et al. 2025](https://files.indicatrix.org/publications/bhagat-conway-compiano-ivie-missing-links.pdf), and [code documentation is available here](https://docs.missinglinks.city). The version of the code used in the paper is accessible under the git tag `ceus`.

If you want to use this with your own network, the easiest way is to install it as a Julia package. You can install it with the Julia package manager by typing `]add MissingLinks`. Full installation instructions are in [the documentation](https://projects.indicatrix.org/MissingLinks.jl).

## Running the code

The code is written in [Julia](https://julialang.org) for performance, and uses a [Quarto](https://quarto.org) notebook format. To use it, [install Julia](https://julialang.org/downloads/); I also recommend installing [Visual Studio Code](https://code.visualstudio.com/) to edit Quarto files across many languages, as well as the [Quarto](https://marketplace.visualstudio.com/items?itemName=quarto.quarto) and [Julia](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia) VSCode extensions.

It will use multiple cores if you have them to speed up routing and scoring. To enable this, in the VSCode settings, search for "julia num threads", click "edit in settings.json", and edit the selected line to read `"julia.NumThreads": "auto"`.

There is a quickstart for getting started with the package [available in the documentation](https://projects.indicatrix.org/MissingLinks.jl/quickstart).

## Reproducing the Charlotte results

The Charlotte results are produced by `run_analysis.qmd` and `rank_comparisons.qmd` in the notebook folder; the other `qmd` files in that folder prepare the data. Unfortunately, the Data Axle data used in the paper is not publicly available.

Open the `MissingLinks.jl` folder in VSCode, and press Cmd/Ctrl-shift-P to get the "command palette" and search for and select "Julia: Start REPL". A Julia REPL (console) should appear. Type `]` to enter the package manager. The prompt should change to say `(MissingLinks.jl) pkg>`. (If it says something like `(@1.9) pkg>` type `activate .` to activate the "environment" for the MissingLinks.jl folder, which contains the dependencies needed to run the model.) Type `instantiate` to download all the necessary packages, and then press backspace to exit the package manager. Type `Threads.nthreads()` and press enter; confirm you get a number greater than 1 to indicate that multithreaded processing is enabled.

You now have Julia and all dependencies installed, and should not need to repeat the above steps unless you get a new computer, upgrade Julia, or change the dependencies of the project.

The datasets used in the Charlotte project are too large to store in Git, so they are stored in the OneDrive folder shared with all project collaborators. This will potentially be in a different location on each person's computer. To tell Julia where to find the data, copy the `config.toml.template` file to `config.toml`, and enter the path to the `Missing links project/data` folder.

To run the analysis, open the `notebooks/run_analysis.qmd file`. You can run cells by clicking "Run cell" or pressing "Cmd/Ctrl-shift-enter". The model requires about 26 GB of free disk space to run (for the mmapped distance matrix).
