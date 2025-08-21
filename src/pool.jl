"""
A RevolvingPool contains a number of interchangeable objects, that can be used in a threadsafe way.
"""
struct RevolvingPool{T}
    size::Int64
    channel::Channel{T}
end

RevolvingPool{T}(size::Int64) where T = RevolvingPool{T}(size, Channel{T}(size))

function with(f, pool::RevolvingPool)
    # grab something from the pool
    obj = take!(pool.channel)
    val = f(obj)
    # return it to the pool
    put!(pool.channel, obj)
    return val
end

function initialize!(f, pool::RevolvingPool)
    for _ in 1:pool.size
        put!(pool.channel, f())
    end
end