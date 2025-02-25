using Documenter
using MissingLinks

makedocs(
    sitename = "MissingLinks",
    format = Documenter.HTML(),
    modules = [MissingLinks],
    pages = [
        "index.md",
        "quickstart.md",
        "api.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
