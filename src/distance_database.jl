const FROM_QUERY = "SELECT dst_code, dist FROM distances WHERE src_code = :src_code"
const TO_QUERY = "SELECT src_code, dist FROM distances WHERE dst_code = :dst_code"
const BOTH_QUERY = "SELECT dist FROM distances WHERE src_code = :src_code AND dst_code = :dst_code"

# should this subclass matrix?
@kwdef struct DistanceMatrix{G <: MetaGraph}
    graph::MetaGraph
    file::String
    db::SQLite.DB
    from_query::SQLite.Stmt
    to_query::SQLite.Stmt
    both_query::SQLite.Stmt
end

"""
    DistanceMatrix(graph)

Initialize a new distance matrix.

Note that DistanceMatrices are not thread safe _for either reading or writing_. For writing, see
[DistanceMatrixQueue](@ref) for a thread-safe writing method. For reading, use [fork()](@ref) to create
multiple independent connections, and then use a RevolvingPool to share connections safely.
"""
function DistanceMatrix(graph::G) where G <: MetaGraph
    dbfile, io = mktemp()
    close(io)

    mtx = DistanceMatrix(graph, dbfile; initialize=true)
    return mtx
end

function DistanceMatrix(graph::G, file; initialize=false) where G <: MetaGraph
    db = SQLite.DB(file)

    if initialize
        initialize!(db)
    end

    DistanceMatrix{G}(;
        graph=graph,
        file=file,
        db=db,
        from_query=DBInterface.prepare(db, FROM_QUERY),
        to_query=DBInterface.prepare(db, TO_QUERY),
        both_query=DBInterface.prepare(db, BOTH_QUERY)
    )
end

"""
    fork(mtx)

Create another connection to the same distance matrix, used for safe multithreaded access
in conjunction with RevolvingPool.
"""
function fork(mtx::DistanceMatrix)
    DistanceMatrix(mtx.graph, mtx.file; initialize=false)
end

"""
    pool(mtx, size=Threads.nthreads() * 2)

Create a [RevolvingPool](@ref) of this DistanceMatrix, for use in multithreaded applications.
Each matrix in the pool shares the same data, and can be used by one thread at a time. DistanceMatrices
are cheap; to avoid blocking, I recommend making the pool twice the number of concurrent tasks that may be
running.
"""
function pool(mtx::DistanceMatrix, size=Threads.nthreads() * 2)
    p = RevolvingPool{DistanceMatrix}(size)
    initialize!(() -> fork(mtx), p)
    return p
end

"""
    initialize!(mtx)

Initializes the database with the distance matrix table.
"""
function initialize!(db)
    # check if table exists already and give a useful error if it does
    # https://stackoverflow.com/questions/1601151
    DBInterface.execute(db, "SELECT name FROM sqlite_master WHERE type = 'table' AND name = 'distances'") do result
        isnothing(iterate(result)) || error("Cache file is stale; table $TABLE_NAME already exists!")
    end

    # In benchmarks these don't seem to make it faster, so don't bother.
    # performance, don't require double writes. matrices are ephemeral
    #DBInterface.execute(identity, mtx.db, "PRAGMA journal_mode = OFF")
    # TODO is this safe? Is it possible that we read it later and it hasn't been written to disk
    #DBInterface.execute(identity, mtx.db, "PRAGMA synchronous = OFF")

    # or perhaps this should use labels? but that's a lot of strings
    DBInterface.execute(identity, db, "CREATE TABLE distances (src_code INTEGER, dst_code INTEGER, dist INTEGER) STRICT")
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
            # why does @prepare use a function?
            DBInterface.execute(
                identity,
                DBInterface.@prepare(() -> mtx.db, "INSERT INTO distances (src_code, dst_code, dist) VALUES (?, ?, ?)"),
                (src_code, tgt, dist)
            )
        end
    end
end

function Base.filesize(d::DistanceMatrix)
    if d.db.file == ":memory:"
        0
    else
        filesize(d.db.file) + filesize(d.db.file * "-shm") + filesize(d.db.file * "-wal")
    end
end

"""
    index!(d::DistanceMatrix)

Index this distance matrix for fast access. Modifying the distance matrix after running index! will still result
in correct results, but will be slower, so you you should build the distance matrix and then index it before querying it.

This creates two indices in the SQLite file, one on (src_code, dst_code, dist) and one on (dst_code, src_code, dist).
We have two so that finding all the distances to _or_ from a node is fast. We include dist in the indices so these can be
index-only queries which are very fast.
"""
function index!(d::DistanceMatrix)
    SQLite.transaction(d.db) do
        @info "creating forward index"
        DBInterface.execute(identity, d.db, "CREATE INDEX fwd_idx ON distances (src_code, dst_code, dist)")
        @info "creating reverse index"
        DBInterface.execute(identity, d.db, "CREATE INDEX rev_idx ON distances (dst_code, src_code, dist)")
    end

    # Not sure this is necessary (i mean for correctness it is not, not sure if it helps performance enough
    # to take the hit of building it)
    # SQLite.transaction(d.db) do
    #     @info "VACUUM ANALYZE"
    #     DBInterface.execute(identity, d.db, "VACUUM")
    #     DBInterface.execute(identity, d.db, "ANALYZE")
    # end
end

"""
    distances_from(matrix, v::VertexID)

Return an iterator of VertexID, Int64 of all vertices reachable from v.
"""
function distances_from(m::DistanceMatrix, v::VertexID)::Tuple{VertexID, Int64}
    # get code
    code = code_for(m.graph, v)
    DBInterface.execute(m.from_query, (src_code=code,)) do cursor
        collect(map(row -> (label_for(m.graph, row.dst_code), row.dist), cursor))
    end
end

"""
    distances_to(matrix, v::VertexID)

Return an iterator of VertexID, Int64 of all vertices that can reach v.
"""
function distances_to(m::DistanceMatrix, v::VertexID)::Tuple{VertexID, Int64}
    code = code_for(m.graph, v)
    DBInterface.execute(m.to_query, (dst_code=code,)) do cursor
        collect(map(row -> (label_for(m.graph, row.src_code), row.dist), cursor))
    end
end

function Base.getindex(m::DistanceMatrix, from::VertexID, to::VertexID)
    DBInterface.execute(m.both_query, (src_code=code_for(m.graph, from), dst_code=code_for(m.graph, to))) do cursor
        itr = iterate(cursor)

        if isnothing(itr)
            # not found
            return missing
        end

        row, next = itr
        # have to do this before next iterate, row will get updated to missing
        dist = row.dist
        isnothing(iterate(cursor, next)) || error("Multiple rows returned")
        dist
    end
end