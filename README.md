# MissingLinks

![Automated test status](https://github.com/mattwigway/MissingLinks.jl/actions/workflows/ci.yml/badge.svg) [![Documentation](https://img.shields.io/badge/docs-latest-4b9cd3)](https://docs.missinglinks.city)

This code supports finding and scoring missing links in transportation networks (primarily pedestrian networks; some changes would be needed for other types of networks, because it assumes the network is undirected, and currently only supports distances up to 65 km). This repository contains a [Julia](https://julialang.org) package to support finding missing links, as well as code for the initial application in Charlotte, North Carolina, USA. The algorithm is described in [Bhagat-Conway et al. 2025](https://files.indicatrix.org/publications/bhagat-conway-compiano-ivie-missing-links.pdf), and [code documentation is available here](https://docs.missinglinks.city). The version of the code used in the paper is accessible under the git tag `ceus`.

If you want to use this with your own network, the easiest way is to install it as a Julia package. You can install it with the Julia package manager by typing `]add MissingLinks`. Full installation instructions are in [the documentation](https://projects.indicatrix.org/MissingLinks.jl).

## Running the code

The code is written in [Julia](https://julialang.org) for performance, and uses a [Quarto](https://quarto.org) notebook format. To use it, [install Julia](https://julialang.org/downloads/); I also recommend installing [Visual Studio Code](https://code.visualstudio.com/) to edit Quarto files across many languages, as well as the [Quarto](https://marketplace.visualstudio.com/items?itemName=quarto.quarto) and [Julia](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia) VSCode extensions.

It will use multiple cores if you have them to speed up routing and scoring. To enable this, in the VSCode settings, search for "julia num threads", click "edit in settings.json", and edit the selected line to read `"julia.NumThreads": "auto"`.

There is a quickstart for getting started with the package [available in the documentation](https://projects.indicatrix.org/MissingLinks.jl/quickstart).

## Reproducing the Charlotte results

As of 2026-01-09, code to run the Charlotte results used in the original paper have been moved to [mattwigway/missing-links-clt](https://github.com/mattwigway/missing-links-clt).
