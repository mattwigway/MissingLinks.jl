# Architecture thoughts:
# Store everything in UInt16 meters. with 120k vertices, will "only" take 26 GB RAM - mmap if needed
# matrix is travel time from origin i to destination j

# Tested in identify_missing_links tests

# NB this codes -Inf as typemax as well. that should not matter.
function round_clamp(T, x)
    if isfinite(x)
        round(T, x)
    elseif x == Inf
        typemax(T)
    elseif x == -Inf
        T <: Signed ? typemin(T) : error("-Inf cannot be converted to unsigned")
    end
end

function one_origin_distances(G, origin, queue, maxdist)
    paths = dijkstra_shortest_paths(G, [origin], maxdist=maxdist)
    dists = round_clamp.(Int64, paths.dists) # replace Inf with typemax(T)
    insert_distances!(queue, origin, dists, maxdist)
end

"""
    calculate_distances(G; maxdist=5000)

Calculate distances  with shortest-path distances based on graph `G`.

The units are the same as the underlying data, and will be rounded to the integers.

To make this faster, you can set a maximum distance `maxdist`; for destinations beyond this distance (or destinations
that are unreachable altogether) the matrix will contain `typemax(T)`.

Will use multiple threads if Julia is started with multiple threads.
"""
function calculate_distances(G; maxdist=5000, origins=1:nv(G))
    @info "Routing with $(Threads.nthreads()) threads"

    # initialize the matrix
    mtx = DistanceMatrix(G)
    queue = DistanceMatrixQueue(mtx)

    writer_task = run(queue)

    # wait for all the routing to finish
    @sync begin
        for i in origins
            Threads.@spawn one_origin_distances(G, $i, queue, maxdist)
        end
    end

    # signal that we're done adding items, let the queue drain
    drain(queue)

    # wait for the writer task to finish
    wait(writer_task)

    index!(mtx)

    return mtx
end

"""
Compute the network distance between the two points on links, by computing
between from ends and adding fractions of the edge.
"""
function compute_net_distance(dmat::DistanceMatrix, sfr, sto, sdist, senddist, dfr, dto, ddist, denddist)
    if sfr == dfr && sto == dto
        # same edge
        # NB: we could remove the base.checked_sub here if perf is poor, as this should never overflow, this is a
        # belt-and-suspenders solution for #7
        sdist > ddist ? Base.checked_sub(sdist, ddist) : Base.checked_sub(ddist, sdist)
    elseif sfr == dto && sto == dfr
        # Links are always supposed to go from lower-numbered vertex to higher-numbered, so
        # this should not be possible
        error("One link is reverse of another!")
    else
        # get the minimum distance, all possible combinations of from and to
        # there are 8 comparisons because this is not assuming undirectedness (though
        # other places in the code still do).
        distances = (
            dmat[sfr, dfr] + sdist + ddist,
            dmat[sto, dfr] + senddist + ddist,
            dmat[sfr, dto] + sdist + denddist,
            dmat[sto, dto] + senddist + denddist,
            dmat[dfr, sfr] + ddist + sdist,
            dmat[dto, sfr] + denddist + sdist,
            dmat[dfr, sto] + ddist + senddist,
            dmat[dto, sto] + denddist + senddist
        )

        all(ismissing.(distances)) ? missing : minimum(skipmissing(distances))
    end
end