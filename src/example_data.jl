"""
    get_example_data()

Get the example dataset used in the [Quickstart](@ref) and save it into the current working directory,
in a folder called "example_data".
"""
function get_example_data()
    if ispath("example_data")
        @warn "Example data folder already exists, not extracting new data."
    else
        cp(artifact"example_data", "example_data")
    end
end