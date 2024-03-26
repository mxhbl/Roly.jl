using LinearAlgebra
using StaticArrays
using Graphs

DefInt, DefFloat = Int16, Float32
SideLoc{T} = Tuple{T,T}                            # Location of a binding site within a structure, in the form (particle, side)
EdgeIdx{T} = Tuple{T,Symbol} # change to +T = src, -T = dst

struct InteractionEdge{T<:Integer}
    src::Tuple{T,T}       # Species, Side
    dst::Tuple{T,T}       # Species, Side
    index::T
end

InteractionEdge{T}(p::Pair, index) where {T} = InteractionEdge(p.first, p.second, index)
InteractionEdge(p::Pair, index) = InteractionEdge{DefInt}(p, index)
function InteractionEdge(v::Vector{T}, index) where {T}
    @assert length(v) == 4
    return InteractionEdge{T}((v[1], v[2]), (v[3], v[4]), index)
end
Base.reverse(e::InteractionEdge) = InteractionEdge(e.dst, e.src, e.index)

function connected_to(edge::InteractionEdge, n::Integer, s::Integer)
    return edge.src == (n, s) || edge.dst == (n, s)
end
function connected_to(edge::InteractionEdge, n::Integer)
    return edge.src[1] == n || edge.dst[1] == n
end
Base.:(==)(e_i::InteractionEdge, e_j::InteractionEdge) = (e_i.src == e_j.src && e_i.dst == e_j.dst) || (e_i.src == e_j.dst && e_i.dst == e_j.src)

Base.convert(::Type{InteractionEdge{T}}, int_edge::InteractionEdge) where {T} = InteractionEdge{T}(int_edge.src, int_edge.dst, int_edge.index)
struct InteractionDiagram{T<:Integer}
    edges::Vector{InteractionEdge{T}}

    function InteractionDiagram{T}(edges::Vector{<:InteractionEdge}) where {T} 
        filtered_edges = []
        k = 1
        for i in 1:length(edges)
            uqe = true
            
            for j in i+1:length(edges)
                if (edges[i] == edges[j])
                    uqe = false
                    break
                end
            end

            if uqe
                push!(filtered_edges, InteractionEdge{T}(edges[i].src, edges[i].dst, k))
                k += 1
            end
        end

        return new{T}(filtered_edges)
    end
end

function InteractionDiagram(A::Matrix{T}) where {T}
    return InteractionDiagram{T}([InteractionEdge(A[i, :], i) for i in 1:size(A, 1)])
end

function Base.push!(int_diag::InteractionDiagram, edge::InteractionEdge)
    return push!(int_diag.edges, edge)
end
Base.pop!(int_diag::InteractionDiagram) = pop!(int_diag.edges)

Base.getindex(int_diag::InteractionDiagram, i::Integer) = int_diag.edges[i]
function Base.size(int_diag::InteractionDiagram)
    return maximum(edge -> max(edge.src[1], edge.dst[1]), int_diag.edges),
           length(int_diag.edges)
end
function Base.iterate(int_diag::InteractionDiagram, state=1)
    return state > length(int_diag.edges) ? nothing : (int_diag[state], state + 1)
end
function Base.show(io::IO, int_diag::InteractionDiagram{T}) where {T}
    let (n, m) = size(int_diag)
        println(io, "InteractionDiagram{$T}[m=$n, k=$m]")
    end
end

function edges(int_diag::InteractionDiagram, n::Integer)
    return [edge for edge in int_diag if connected_to(edge, n)]
end

Base.convert(::Type{InteractionDiagram{T}}, int_diag::InteractionDiagram) where {T} = InteractionDiagram{T}(convert.(InteractionEdge{T}, int_diag.edges))
