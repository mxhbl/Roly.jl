using LinearAlgebra
using StaticArrays
using Graphs
using NautyGraphs

mutable struct Polyform{D,T<:Integer,F<:AbstractFloat,R<:RotationOperator{F}}
    anatomy::DirectedDenseNautyGraph{Cint}
    encoder::PolyEncoder{T}
    bond_partners::Vector{T}
    species::Vector{T}
    xs::Vector{Point{D,F}}
    ψs::Vector{R}
    σ::T
    function Polyform{D,T,F}(anatomy, encoder, bond_status, species, xs, ψs, σ) where {D,T,F}
        if D == 2
            return new{D,T,F,Angle{F}}(anatomy, encoder, bond_status, species, xs, ψs, σ)
        elseif D == 3
            return new{D,T,F,Quaternion{F}}(anatomy, encoder, bond_status, species, xs, ψs, σ)
        else
            error("Polyforms are only defined in dimension 2 and 3.")
        end
    end
end
function Polyform{D,T,F}(anatomy, encoder, bond_status, species, xs, ψs) where {D,T,F}
    p = Polyform{D,T,F}(anatomy, encoder, bond_status, species, xs, ψs, 0)
    canonize!(p)
    return p
end 
Polyform(anatomy, encoder, bond_status, species::AbstractVector{T}, xs::AbstractVector{<:Point{D,F}}, ψs, σ) where {D,T,F} = 
Polyform{D,T,F}(anatomy, encoder, bond_status, species, xs, ψs, σ)
Polyform(anatomy, encoder, bond_status, species::AbstractVector{T}, xs::AbstractVector{<:Point{D,F}}, ψs) where {D,T,F} = 
Polyform{D,T,F}(anatomy, encoder, bond_status, species, xs, ψs)
function Polyform{D,T,F}() where {D,T,F}
    if D == 2
        ψs = Angle{F}[]
    else
        ψs = Quaternion{F}[]
    end
    return Polyform{D,T,F}(DirectedDenseNautyGraph(0), PolyEncoder{T}(), zeros(T, 0), zeros(T, 0), SVector{2,F}[], ψs, 0)
end
anatomy(p::Polyform) = p.anatomy
faces(::Type{T}, p::Polyform) where {T} = convert(Vector{T}, p.anatomy.labels)
faces(p::Polyform) = p.anatomy.labels
species(::Type{T}, p::Polyform) where {T} = convert(Vector{T}, p.species)
species(p::Polyform) = p.species
symmetry_number(p::Polyform) = p.σ
positions(p::Polyform) = p.xs
orientations(p::Polyform) = p.ψs

function Base.size(p::Polyform)
    return length(p.xs)
end
nvertices(p::Polyform) = nv(p.anatomy)
function Base.show(io::Core.IO, p::Polyform{D,T,F}) where {D,T,F}
    print(io, "Polyform{$D,$T,$F}[n=$(siez(p))]")
end
Base.show(io::Core.IO, ::Type{Polyform{D,T,F}}) where {D,T,F} = print(io, "Polyform{$D,$T,$F}")

function canonize!(p::Polyform)
    canon_perm, n = NautyGraphs.canonize!(p.anatomy)
    permute!(p.encoder, canon_perm)
    p.bond_partners .= p.bond_partners[canon_perm]

    inv_perm = invperm(canon_perm)
    for (i, v) in enumerate(p.bond_partners)
        if v == 0
            continue
        end
        p.bond_partners[i] = inv_perm[v]
    end
    p.σ = n
    return
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

function Base.copy(p::Polyform)
    return Polyform(copy(p.anatomy), copy(p.encoder), copy(p.bond_partners), copy(p.species), copy(p.xs), copy(p.ψs), p.σ)
end
function Base.copy!(dest::Polyform{D,T,F}, src::Polyform{D,T,F}) where {D,T,F}
    copy!(dest.anatomy, src.anatomy)
    copy!(dest.encoder, src.encoder)
    copy!(dest.bond_partners, src.bond_partners)
    copy!(dest.species, src.species)
    copy!(dest.xs, src.xs)
    copy!(dest.ψs, src.ψs)
    dest.σ = src.σ
    return dest
end

function create_monomer(geometry::AbstractGeometry{T,F},
                        species_idx::Integer,
                        face_labels::AbstractVector{<:Integer}) where {T,F}

    if dimension(geometry) == 2
        xs = [Point{2,F}(0., 0.)]
        ψs = [Angle{F}(0.)]
    elseif dimension(geometry) == 3
        xs = [Point{3,F}(0., 0., 0.)]
        ψs = [Quaternion{F}(1., 0., 0., 0.)]
    else
        error("Polyforms are only implemented for dimension 2 and 3.")
    end

    return Polyform(DirectedDenseNautyGraph(geometry.anatomy, convert(Vector{T}, face_labels)),
                    copy(geometry.encoder),
                    zeros(T, geometry.encoder.n_vertices),
                    [T(species_idx)],
                    xs, ψs)
end