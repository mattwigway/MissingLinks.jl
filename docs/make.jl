using Documenter
using MissingLinks

makedocs(
    sitename = "MissingLinks",
    format = Documenter.HTML(),
    modules = [MissingLinks]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
