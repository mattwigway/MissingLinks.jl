"""
A CandidateLink represents a link between two existing edges. It stores the vertex codes (not labels)
from each end of each edge, as well as the distances from each vertex at the point where the link is.
"""
@kwdef struct CandidateLink
    fr_edge_src::VertexID
    fr_edge_tgt::VertexID
    fr_dist_from_start::Int64
    fr_dist_to_end::Int64
    to_edge_src::VertexID
    to_edge_tgt::VertexID
    to_dist_from_start::Int64
    to_dist_to_end::Int64
    geographic_length_m::Int64
    network_length_m::Union{Int64, Missing}
end

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