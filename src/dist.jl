# Architecture thoughts:
# Store everything in UInt16 meters. with 120k vertices, will "only" take 26 GB RAM - mmap if needed
# matrix is travel time from origin i to destination j

# Tested in identify_missing_links tests

function fill_matrix!(G, mtx::AbstractMatrix{T}; maxdist=5000, origins=1:nv(G)) where T
    size(mtx) == (nv(G), nv(G)) || error("Matrix must be $(nv(G))x$(nv(G))")

    @info "Routing with $(Threads.nthreads()) threads"

    ThreadsX.mapi(origins) do origin
        if origin % 1000 == 0
            @info "Processed $origin / $(nv(G)) trips"
        end
        
        paths = dijkstra_shortest_paths(G, [origin], maxdist=maxdist)
        mtx[:, origin] = round.(T, min.(paths.dists, typemax(T)))
    end
end

"""
Compute the network distance between the two points on links, by computing
between from ends and adding fractions of the edge.
"""
function compute_net_distance(dmat::Matrix{T}, sfr, sto, sdist, senddist, dfr, dto, ddist, denddist) where T
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
        min(
            add_unless_typemax(dmat[sfr, dfr], sdist + ddist),
            add_unless_typemax(dmat[sto, dfr], senddist + ddist),
            add_unless_typemax(dmat[sfr, dto], sdist + denddist),
            add_unless_typemax(dmat[sto, dto], senddist + denddist),
            add_unless_typemax(dmat[dfr, sfr], ddist + sdist),
            add_unless_typemax(dmat[dto, sfr], denddist + sdist),
            add_unless_typemax(dmat[dfr, sto], ddist + senddist),
            add_unless_typemax(dmat[dto, sto], denddist + senddist)
        )
    end
end
