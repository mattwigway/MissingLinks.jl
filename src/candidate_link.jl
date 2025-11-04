"""
A CandidateLink represents a link between two existing edges. It stores the vertex codes (not labels)
from each end of each edge, as well as the distances from each vertex at the point where the link is.
"""
@kwdef struct CandidateLink{T<:Real}
    fr_edge_src::VertexID
    fr_edge_tgt::VertexID
    fr_dist_from_start::T
    fr_dist_to_end::T
    to_edge_src::VertexID
    to_edge_tgt::VertexID
    to_dist_from_start::T
    to_dist_to_end::T
    geographic_length_m::T
    network_length_m::T
end

Base.hash(l::CandidateLink, h::UInt64) = hash(
    l.fr_edge_src, hash(
        l.fr_edge_tgt, hash(
            l.fr_dist_from_start, hash(
                l.fr_dist_to_end, hash(
                    l.to_edge_src, hash(
                        l.to_edge_tgt, hash(
                            l.to_dist_from_start, hash(
                                l.to_dist_to_end, hash(
                                    l.geographic_length_m, hash(
                                        l.network_length_m, hash(
                                            :CandidateLink,
                                            h
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )
    )
)

Base.isequal(a::CandidateLink, b::CandidateLink) =
    a.fr_edge_src == fr_edge_src &&
    a.fr_edge_tgt == fr_edge_tgt &&
    a.fr_dist_from_start == fr_dist_from_start &&
    a.fr_dist_to_end == fr_dist_to_end &&
    a.to_edge_src == to_edge_src &&
    a.to_edge_tgt == to_edge_tgt &&
    a.to_dist_from_start == to_dist_from_start &&
    a.to_dist_to_end == to_dist_to_end &&
    a.geographic_length_m == geographic_length_m &&
    # this is why we can't use @struct_hash_equals - need to handle missings differently
    ((ismissing(a.network_length_m) && ismissing(b.network_length_m)) || a.network_length_m == network_length_m)

"""
Create a reversed version of a candidate link, used in evaluating accessibility to calculate
accessibility in both directions.
"""
function Base.reverse(link::CandidateLink)
    CandidateLink(
        link.to_edge_src,
        link.to_edge_tgt,
        link.to_dist_from_start,
        link.to_dist_to_end,
        link.fr_edge_src,
        link.fr_edge_tgt,
        link.fr_dist_from_start,
        link.fr_dist_to_end,
        link.geographic_length_m,
        link.network_length_m # in undirected graph this should be identical
    )
end

Base.show(io::IO, x::CandidateLink) =
    println(io, "$(round(Int64, x.geographic_length_m))m CandidateLink connecting edge $(x.fr_edge_src)-$(x.fr_edge_tgt)@$((x.fr_dist_from_start))m to $(x.to_edge_src)-$(x.to_edge_tgt)@$((x.to_dist_from_start))m; network distance $((x.network_length_m))m")