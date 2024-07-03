# Architecture thoughts:
# Store everything in UInt16 meters. with 120k vertices, will "only" take 26 GB RAM - mmap if needed
# matrix is travel time from origin i to destination j

# Tested in identify_missing_links tests

function create_matrix(G, T; maxdist=5000, origins=1:nv(G))
    @info "Routing with $(Threads.nthreads()) threads"

    maxdist < typemax(T) || error("maxdist must be less than typemax(T) = $(typemax(T)); typemax is used to indicate unreachable")

    return Distances{T}(Iterators.flatten(ThreadsX.mapi(origins) do origin
        if origin % 1000 == 0
            @info "Processed $origin / $(nv(G)) trips"
        end
        
        paths = dijkstra_shortest_paths(G, [origin], maxdist=maxdist)
        return filter(x -> x[3] ≤ maxdist, tuple.(origin, 1:nv(G), round.(T, min.(paths.dists, typemax(T)))))
    end), nv(G))
end

"""
Store distances between all nodes in the graph. The distances are in two sorted lists, one of
(origin, destination, distance) and one of (destination, origin, distance) so that they can be
easily iterated over.
"""
struct Distances{T <: Real} <: AbstractMatrix{T}
    by_origin::Vector{Tuple{Int64, Int64, T}}
    by_destination::Vector{Tuple{Int64, Int64, T}}
    n::Int64
end

"""
    Distances([(origin, destination, dist)])
"""
function Distances{T}(iter, n::Int64) where T
    by_origin = collect(iter)
    sort!(by_origin) # since routing is multithreaded, it may not be sorted
    by_dest = [(d, o, dist) for (o, d, dist) in by_origin]
    sort!(by_dest)
    Distances{T}(by_origin, by_dest, n)
end

# Matrix interface - used in deduplication where we are not using the optimized
# iteration method used in score_links
Base.eltype(::Distances{T}) where T = T
Base.size(d::Distances) = (d.n, d.n)

"""
    getindex(d::Distances{T}, origin, destination)

Get the distance from origin to destination stored in d, or typemax(T) if unreachable.
"""
function Base.getindex(d::Distances{T}, origin, destination) where T
    idx = searchsortedfirst(d.by_origin, (origin, destination, typemin(T)))
    if idx ≤ lastindex(d.by_origin) && d.by_origin[idx][1] == origin && d.by_origin[idx][2] == destination
        return d.by_origin[idx][3]
    else
        return typemax(T)
    end
end



"""
An iterator over all destinations from a given origin or destination.
"""
struct DistanceIterator{T}
    origin_or_destination::Int64
    first_index::Int64
    sorted::Vector{Tuple{Int64, Int64, T}}
end

function DistanceIterator{T}(sorted, origin_or_destination) where T
    first_index = searchsortedfirst(sorted, (origin_or_destination, 0, typemin(T)))
    DistanceIterator{T}(origin_or_destination, first_index, sorted)
end

Base.iterate(d::DistanceIterator) = Base.iterate(d, d.first_index)

Base.IteratorSize(d::DistanceIterator) = Base.SizeUnknown()

function Base.eltype(d::DistanceIterator{T}) where T
    Tuple{Int64, T}
end

function Base.iterate(d::DistanceIterator, idx)
    if idx > lastindex(d.sorted) || d.sorted[idx][1] != d.origin_or_destination
        # nothing left to iterate over
        return nothing
    else
        # 2:3 will always return the origin/destination (whichever we are iterating over) and distance
        return (d.sorted[idx][2:3], idx + 1)
    end
end

"""
    CombinedDistanceIterator(iter1, dist1, iter2, dist2)

This combines two distance iterators into a single distance iterator, used when iterating
over distances to both ends of a source/target edge. It efficiently yields tuples of
(node, distance) where node is a node reachable in either iter1 or iter2, and distance
is the minimum distance to that node, i.e. (min(dist from iter1 + dist1, dist from iter2 + dist2)).
"""
struct CombinedDistanceIterator{T}
    iter1::DistanceIterator{T}
    dist1::T
    iter2::DistanceIterator{T}
    dist2::T
end

# Start out at the start of both 
Base.iterate(d::CombinedDistanceIterator) = Base.iterate(d, (d.iter1.first_index, d.iter2.first_index))

function Base.iterate(d::CombinedDistanceIterator, state)
    local newloc1, node1, dist1, res1, newloc2, node2, dist2, res2
    loc1, loc2 = state
    
    res1 = iterate(d.iter1, loc1)
    if !isnothing(res1)
        newloc1 = res1[2]
        node1, dist1 = res1[1]
        dist1 = add_unless_typemax(dist1, d.dist1)
    end

    res2 = iterate(d.iter2, loc2)
    if !isnothing(res2)
        newloc2 = res2[2]
        node2, dist2 = res2[1]
        dist2 = add_unless_typemax(dist2, d.dist2)
    end

    if isnothing(res1) && isnothing(res2)
        # both iterators have run out
        return nothing
    elseif isnothing(res1)
        # iterator 1 has run out, increment iter2
        return ((node2, dist2), (loc1, newloc2))
    elseif isnothing(res2)
        # iterator 2 has run out, increment iter1
        return ((node1, dist1), (newloc1, loc2))
    elseif node1 == node2
        # The node is in both iterators
        # increment both
        return ((node1, min(dist1, dist2)), (newloc1, newloc2))
    elseif node1 < node2
        # the next node should come from iterator 1, iterator 2 is already past this point and didn't have it
        # increment iter1, don't increment iter2
        return ((node1, dist1), (newloc1, loc2))
    elseif node2 < node1
        # the next node should come from iterator 2, iterator 1 is already past this point and didn't have it
        # increment iter2, don't increment iter1
        return ((node2, dist2), (loc1, newloc2))
    else
        error("Inconsistent state in CombinedDistanceIterator; this is bug. Report at https://github.com/mattwigway/MissingLinks.jl/issues")
    end
end

"""
    margin(distances, origin=n)

or

    margin(distances, destination=n)

Get a margin of a distance matrix - in the first form, all destinations from a particular origin, and in the second,
all origins to a particular destination. Returns a stateless iterator of tuples of (origin/destination, distance)—the first
element is whichever dimension was not fixed.
"""
function margin(distances::Distances{T}; origin::Union{Nothing, Int64}=nothing, destination::Union{Nothing, Int64}=nothing) where T
    ((isnothing(origin) && isnothing(destination)) || (!isnothing(origin) && !isnothing(destination))) &&
        error("Exactly one of origin or destination must be specified")

    return if !isnothing(origin)
        if origin < 1 || origin > distances.n
            throw(BoundsError(distances, origin))
        end
        DistanceIterator{T}(distances.by_origin, origin)
    else
        if destination < 1 || destination > distances.n
            throw(BoundsError(distances, destination))
        end
        DistanceIterator{T}(distances.by_destination, destination)
    end
end

"""
    sizeof(d::Distances{T})

Return the size in bytes of the stored distance matrix. Not 100% accurate as it does not account for struct padding,
but that is negligible in real-world distance matrices.
"""
sizeof(d::Distances) = sizeof(d.by_origin) + sizeof(d.by_destination) + sizeof(d.n)

"""
Compute the network distance between the two points on links, by computing
between from ends and adding fractions of the edge.
"""
function compute_net_distance(dmat::Distances{T}, sfr, sto, sdist, senddist, dfr, dto, ddist, denddist) where T
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
