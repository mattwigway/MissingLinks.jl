"""
Code to manage inserts into a distance matrix in a threadsafe way.
"""

@kwdef struct DistanceMessage
    src_code::Int64
    dists::Vector{Int64}
    maxdist::Int64
end

struct DistanceMatrixQueue
    mtx::DistanceMatrix
    channel::Channel{DistanceMessage}

    DistanceMatrixQueue(mtx::DistanceMatrix) = new(mtx, Channel{DistanceMessage}())
end

"""
    run(q)

This task handles receiving distances from generating threads and writing them to SQLite. It starts the writer in a new thread.
"""
function run(q::DistanceMatrixQueue)
    t = Threads.@spawn _run_task(q)
    # fail fast if writer fails
    # Note this is not working
    bind(q.channel, t)
    return t
end

function _run_task(q::DistanceMatrixQueue)
    count = 0
    total = nv(q.mtx.graph)

    SQLite.transaction(q.mtx.db) do
        for msg ∈ q.channel
            count += 1
            insert_distances!(q.mtx, msg.src_code, msg.dists, msg.maxdist)

            if count % 1000 == 0
                @info "Processed $count / $total origins"
            end
        end
    end
end

"""
    insert_distances!(mtx::DistanceMatrix, src_code, distances, max_dist)

Insert the distances from origin src_code (a vertex _code_, not [VertexID](@ref)) to all other nodes, represented by vector `distances`.
Distances greater than max_dist will not be inserted.

This function is thread safe, and can be safely used from multiple writing threads.
"""
function insert_distances!(q::DistanceMatrixQueue, src_code::Int64, distances::Vector{Int64}, max_dist::Int64)
    put!(q.channel, DistanceMessage(src_code=src_code, dists=distances, maxdist=max_dist))
end

"""
    done(queue)

Indicate that there are no more values to be placed in the queue.
"""
drain(queue::DistanceMatrixQueue) = close(queue.channel)

