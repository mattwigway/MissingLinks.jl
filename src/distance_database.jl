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

Note that DistanceMatrices are not thread safe _for either reading or writing_. For writing, see
[DistanceMatrixQueue](@ref) for a thread-safe writing method. For reading, use [fork()](@ref) to create
multiple independent connections, and then use a RevolvingPool to share connections safely.
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
function distances_from(m::DistanceMatrix, v::VertexID)
    # get code
    code = code_for(m.graph, v)
    DBInterface.execute(
        DBInterface.@prepare(() -> m.db, "SELECT dst_code, dist FROM distances WHERE src_code = ?"),
        (code,)
    ) do cursor
        collect(map(row -> (label_for(m.graph, row.dst_code), row.dist), cursor))
    end
end

"""
    distances_to(matrix, v::VertexID)

Return an iterator of VertexID, Int64 of all vertices that can reach v.
"""
function distances_to(m::DistanceMatrix, v::VertexID)
    # get code
    code = code_for(m.graph, v)
    DBInterface.execute(
        DBInterface.@prepare(() -> m.db, "SELECT src_code, dist FROM distances WHERE dst_code = ?"),
        (code,)
    ) do cursor
        collect(map(row -> (label_for(m.graph, row.src_code), row.dist), cursor))
    end
end