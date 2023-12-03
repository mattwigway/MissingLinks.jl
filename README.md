# MissingLinks

This code supports the bicycle/pedestrian missing links project; details of the algorithm are in the paper.

## Running the code

The code is written in [Julia](https://julialang.org) for performance, and uses a [Quarto](https://quarto.org) notebook format. To use it, [install Julia](https://julialang.org/downloads/); I also recommend installing [Visual Studio Code](https://code.visualstudio.com/) to edit Quarto files across many languages, as well as the [Quarto](https://marketplace.visualstudio.com/items?itemName=quarto.quarto) and [Julia](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia) VSCode extensions.

It will use multiple cores if you have them to speed up routing and scoring. To enable this, in the VSCode settings, search for "julia num threads", click "edit in settings.json", and edit the selected line to read `"julia.NumThreads": "auto"`.

Open the `MissingLinks.jl` folder in VSCode, and press Cmd/Ctrl-shift-P to get the "command palette" and search for and select "Julia: Start REPL". A Julia REPL (console) should appear. Type `]` to enter the package manager. The prompt should change to say `(MissingLinks.jl) pkg>`. (If it says something like `(@1.9) pkg>` type `activate .` to activate the "environment" for the MissingLinks.jl folder, which contains the dependencies needed to run the model.) Type `instantiate` to download all the necessary packages, and then press backspace to exit the package manager. Type `Threads.nthreads()` and press enter; confirm you get a number greater than 1 to indicate that multithreaded processing is enabled.

You now have Julia and all dependencies installed, and should not need to repeat the above steps unless you get a new computer, upgrade Julia, or change the dependencies of the project.

The datasets used are too large to store in Git, so they are stored in the OneDrive folder shared with all project collaborators. This will potentially be in a different location on each person's computer. To tell Julia where to find the data, copy the `config.toml.template` file to `config.toml`, and enter the path to the `Missing links project/data` folder.

To run the analysis, open the `notebooks/run_analysis.qmd file`. You can run cells by clicking "Run cell" or pressing "Cmd/Ctrl-shift-enter". The model requires about 26 GB of free disk space to run (for the mmapped distance matrix).