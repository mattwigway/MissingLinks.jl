"""
A CandidateLink represents a link between two existing edges. It stores the vertex codes (not labels)
from each end of each edge, as well as the distances from each vertex at the point where the link is.
"""
struct CandidateLink{T<:Real}
    fr_edge_src::Int64
    fr_edge_tgt::Int64
    fr_dist_from_start::T
    fr_dist_to_end::T
    to_edge_src::Int64
    to_edge_tgt::Int64
    to_dist_from_start::T
    to_dist_to_end::T
    geographic_length_m::T
    network_length_m::T
end

# constructor with keyword args for clarity in tests. All keywords must be passed, or will TypeError (intentionally)
CandidateLink(;
    fr_edge_src=nothing,
    fr_edge_tgt=nothing,
    fr_dist_from_start=nothing,
    fr_dist_to_end=nothing,
    to_edge_src=nothing,
    to_edge_tgt=nothing,
    to_dist_from_start=nothing,
    to_dist_to_end=nothing,
    geographic_length_m=nothing,
    network_length_m=nothing
) = CandidateLink(
        fr_edge_src,
        fr_edge_tgt,
        fr_dist_from_start,
        fr_dist_to_end,
        to_edge_src,
        to_edge_tgt,
        to_dist_from_start,
        to_dist_to_end,
        geographic_length_m,
        network_length_m
    )


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

Base.show(io::IO, x::CandidateLink{<:Any}) =
    println(io, "$(round(Int64, x.geographic_length_m))m CandidateLink connecting edge $(x.fr_edge_src)-$(x.fr_edge_tgt)@$(round(Int64, x.fr_dist_from_start))m to $(x.to_edge_src)-$(x.to_edge_tgt)@$(round(Int64, x.to_dist_from_start))m; network distance $(round(Int64, x.network_length_m))m")