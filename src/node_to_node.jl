# Architecture thoughts:
# Store everything in UInt16 meters. with 120k vertices, will "only" take 26 GB RAM - mmap if needed
# matrix is travel time from origin i to destination j

function fill_matrix!(G, mtx::Matrix{T}; maxdist=5000, origins=1:nv(G)) where T
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
