# should this subclass matrix?
struct DistanceMatrix{G <: MetaGraph}
    graph::MetaGraph
    db::SQLite.DB

    function DistanceMatrix(graph::G, db::SQLite.DB) where G
        mtx = new{G}(graph, db)
        initialize(mtx)
        return mtx
    end
end

"""
    DistanceMatrix(graph; memory=false)

Initialize a new distance matrix. If you have a lot of memory, set memory=true for best performance.
Otherwise set memory=false for lowest memory consumption.
"""
function DistanceMatrix(graph; memory=false)
    db = if memory
        SQLite.DB()
    else
        dbfile, io = mktemp()
        close(io)
        SQLite.DB(dbfile)
    end

    return DistanceMatrix(graph, db)
end

"""
    initialize(mtx)

Initializes the database with the distance matrix table.
"""
function initialize(mtx::DistanceMatrix)
    # check if table exists already and give a useful error if it does
    # https://stackoverflow.com/questions/1601151
    DBInterface.execute(mtx.db, "SELECT name FROM sqlite_master WHERE type = 'table' AND name = 'distances'") do result
        isnothing(iterate(result)) || error("Cache file is stale; table $TABLE_NAME already exists!")
    end

    # In benchmarks these don't seem to make it faster, so don't bother.
    # performance, don't require double writes. matrices are ephemeral
    #DBInterface.execute(identity, mtx.db, "PRAGMA journal_mode = OFF")
    # TODO is this safe? Is it possible that we read it later and it hasn't been written to disk
    #DBInterface.execute(identity, mtx.db, "PRAGMA synchronous = OFF")

    # or perhaps this should use labels? but that's a lot of strings
    DBInterface.execute(identity, mtx.db, "CREATE TABLE distances (src_code INTEGER, dst_code INTEGER, dist INTEGER) STRICT")
end

"""
    insert_distances!(mtx::DistanceMatrix, src_code, distances, max_dist)

Insert the distances from origin src_code (a vertex _code_, not [VertexID](@ref)) to all other nodes, represented by vector `distances`.
Distances greater than max_dist will not be inserted.

Note that this function is NOT THREAD SAFE. Use DistanceMatrixQueue for thread-safe insertion.
"""
function insert_distances!(mtx::DistanceMatrix, src_code::Int64, distances::AbstractVector{<:Integer}, max_dist)
    Base.require_one_based_indexing(distances)
    for (tgt, dist) in enumerate(distances)
        if dist ≤ max_dist
            # why does this use a function?
            DBInterface.execute(
                identity,
                DBInterface.@prepare(() -> mtx.db, "INSERT INTO distances (src_code, dst_code, dist) VALUES (?, ?, ?)"),
                (src_code, tgt, dist)
            )
        end
    end
end

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

function Base.filesize(d::DistanceMatrix)
    if d.db.file == ":memory:"
        0
    else
        filesize(d.db.file) + filesize(d.db.file * "-shm") + filesize(d.db.file * "-wal")
    end
end
