# This file is separate from partition.jl to avoid a circular dependency

"""
A GraphPartition represents a subgraph of a larger graph, defined by a bounding box. The subgraph
is actually `buffer` units larger than the bounding box, to avoid edge effects, but only candidate
links with an end inside the bounding box should be considered. Consider a portion
of a graph that looks like this:
  
    ✂
A --✂--- B
    ✂    |
    ✂    |
    ✂    |
D---✂--- C
    ✂

A is relatively directly connected to D, but suppose where the scissor are is the edge of a partition that extends
to the left. In the partition that extends to the left, AB and CD will be present, but BC will not as it is completely
outside the partition. This will result in a candidate link AB->CD being identified and may score highly as it closes
a major gap in the graph. This subgraph will be completely within the partition to the right, and we will get the correct
link and score there, but the score may be lower or AB->CD may not be identified at all.

I think what we need to do is only include links where at least one end of them is in the "main" part of the partition,
and then increase the buffer to maxdist + max_link_dist.
"""
struct GraphPartition{T}
    G::T
    min_x::Float64
    max_x::Float64
    min_y::Float64
    max_y::Float64
    buffer::Float64
end

Graphs.nv(G::GraphPartition) = Graphs.nv(G.G)
Graphs.ne(G::GraphPartition) = Graphs.ne(G.G)


"""
    partition(G, nrow, ncol, overlap)

Partition graph G into nrow x ncol subgraphs. This is useful with graphs that are too large
to process at once. `buffer` is the amount the graphs should overlap, in meters/graph units
(though using units other than meters is untested and unscientific).

Specifically, the bounding box of each subgraph will be extended by `buffer` meters, so the actual
overlap will be 2 * `overlap`.

Returns an nrow x ncol matrix of graphs.
"""
function partition(G::T, nrow, ncol, buffer) where T
    result = Matrix{GraphPartition{T}}(undef, nrow, ncol)

    # Get the envelope of the entire graph
    bbox = extent(G)

    col_x = (bbox.max_x - bbox.min_x) / ncol
    row_y = (bbox.max_y - bbox.min_y) / nrow

    for r in 1:nrow, c in 1:ncol
        Gsub = new_graph()

        xmin = bbox.min_x + col_x * (c - 1)
        xmax = xmin + col_x
        ymin = bbox.min_y + row_y * (r - 1)
        ymax = ymin + row_y

        xminb = xmin - buffer
        yminb = ymin - buffer
        xmaxb = xmax + buffer
        ymaxb = ymax + buffer

        # copy over all vertices in area
        for v in labels(G)
            x, y = G[v]
            if x ≥ xminb && x ≤ xmaxb && y ≥ yminb && y ≤ ymaxb
                Gsub[v] = G[v]
            end
        end

        # copy all edges
        for v in labels(Gsub)
            nbrs = neighbor_labels(G, v)
            
            for nbr in nbrs
                if !haskey(Gsub, nbr)
                    # include edges that start within and go outside area
                    Gsub[nbr] = G[nbr]
                end

                v1, v2 = order_vertices(v, nbr)
                Gsub[v1, v2] = G[v1, v2]
            end
        end

        result[r, c] = GraphPartition{T}(
            Gsub,
            xmin,
            xmax,
            ymin,
            ymax,
            buffer
        )
    end

    return result
end

"""
    extent(G)

Get the spatial extent of graph G. Returns a NamedTuple with members `min_x`, `max_x`, `min_y`, `max_y`.
The bounding box is based only on the vertices.
"""
function extent(G)
    min_x = Inf64
    max_x = -Inf64
    min_y = Inf64
    max_y = -Inf64

    for v in labels(G)
        x, y = G[v]
        min_x = min(x, min_x)
        max_x = max(x, max_x)
        min_y = min(y, min_y)
        max_y = max(y, max_y)
    end

    @assert all(isfinite.((min_x, max_x, min_y, max_y)))

    return (min_x=min_x, min_y=min_y, max_x=max_x, max_y=max_y)
end
