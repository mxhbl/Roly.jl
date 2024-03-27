using LinearAlgebra
using StaticArrays
using Graphs
using NautyGraphs

mutable struct Polyform{T<:Integer,F<:AbstractFloat}
    anatomy::DirectedDenseNautyGraph{Cint}
    translator::EncTranslator{Int}
    species::Vector{T}
    positions::Coordinates2D{F}
    σ::T

    function Polyform{T,F}(anatomy, translator, species, positions) where {T,F}
        s = new{T,F}(anatomy, translator, species, positions, 0)
        canonize!(s)
        return s
    end

    function Polyform{T,F}(anatomy, translator, species, positions, σ) where {T,F}
        return new{T,F}(anatomy, translator, species, positions, σ)
    end
end
Polyform(args...) = Polyform{DefInt,DefFloat}(args...)
Polyform{T,F}() where {T,F} = Polyform{T,F}(DirectedDenseNautyGraph(0), EncTranslator{Int}(), zeros(T, 0), Coordinates2D{F}([], []), 0)

function canonize!(p::Polyform)
    canon_perm, n = NautyGraphs.canonize!(p.anatomy)
    permute!(p.translator, canon_perm)
    p.σ = n
    return
end

faces(::Type{T}, p::Polyform) where {T} = convert(Vector{T}, p.anatomy.labels)
faces(p::Polyform) = faces(Int, p)
species(p::Polyform) = p.species

function Base.size(p::Polyform)
    e = ne(p.anatomy)
    v = nv(p.anatomy)
    return length(p.positions), (e - v) ÷ 2, v
end #  TODO: nbonds only valid in 2D
function Base.show(io::IO, p::Polyform{T,F}) where {T,F}
    let (n, m, _) = size(p)
        println(io, "Polyform{$T, $F}[n=$n, k=$m]")
    end
end

function Base.hash(p::Polyform)
    return hash(p.anatomy)
end
function Base.hash(p::Polyform, h::UInt64)
    return hash(p.anatomy, h)
end
function Base.:(==)(si::Polyform, sj::Polyform)
    return hash(si.anatomy) == hash(sj.anatomy)
end

function Base.copy(p::Polyform{T,F}) where {T,F}
    return Polyform{T,F}(copy(p.anatomy), copy(p.translator), copy(p.species), deepcopy(p.positions), p.σ)
end
function Base.copy!(dest::Polyform{T,F}, src::Polyform{T,F}) where {T,F}
    copy!(dest.anatomy, src.anatomy)
    copy!(dest.translator, src.translator)
    copy!(dest.species, src.species)
    copy!(dest.positions, src.positions)
    dest.σ = src.σ
    return dest
end

function create_monomer(geometry::AbstractGeometry{T,F},
                        species_idx::T,
                        face_labels::AbstractVector{<:Integer}) where {T<:Integer,F<:Real}

    return Polyform{T,F}(DirectedDenseNautyGraph(geometry.anatomy, face_labels),
                          EncTranslator([collect(1:length(geometry))]),
                          [species_idx],
                        #   Coordinates2D{F}([SA[0.0, 0.0, 0.0]], [SA[1.0, 0.0, 0.0, 0.0]]))
                          Coordinates2D{F}([SA[0.0, 0.0]], [0.0]))

end