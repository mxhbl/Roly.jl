using LinearAlgebra
using StaticArrays
using Graphs
using NautyGraphs

mutable struct Polyform{D,T<:Integer,F<:AbstractFloat,R<:RotationOperator{F}}
    anatomy::DirectedDenseNautyGraph{Cint}
    translator::EncTranslator{Int}
    species::Vector{T}
    xs::Vector{Point{D,F}}
    ψs::Vector{R}
    σ::T
    function Polyform{D,T,F}(anatomy, translator, species, xs, ψs, σ) where {D,T,F}
        if D == 2
            return new{D,T,F,Angle{F}}(anatomy, translator, species, xs, ψs, σ)
        elseif D == 3
            return new{D,T,F,Quaternion{F}}(anatomy, translator, species, xs, ψs, σ)
        else
            error("Polyforms are only defined in dimension 2 and 3.")
        end
    end
end
Base.show(io::Core.IO, ::Type{Polyform{D,T,F}}) where {D,T,F} = println(io, "Polyform{$D,$T,$F}")
function Polyform{D,T,F}(anatomy, translator, species, xs, ψs) where {D,T,F}
    p = Polyform{D,T,F}(anatomy, translator, species, xs, ψs, 0)
    canonize!(p)
    return p
end 
Polyform(anatomy, translator, species::AbstractVector{T}, xs::AbstractVector{<:Point{D,F}}, ψs, σ) where {D,T,F} = 
Polyform{D,T,F}(anatomy, translator, species, xs, ψs, σ)
Polyform(anatomy, translator, species::AbstractVector{T}, xs::AbstractVector{<:Point{D,F}}, ψs) where {D,T,F} = 
Polyform{D,T,F}(anatomy, translator, species, xs, ψs)

function Polyform{D,T,F}() where {D,T,F}
    if D == 2
        ψs = Angle{F}[]
    else
        ψs = Quaternion{F}[]
    end
    return Polyform{D,T,F}(DirectedDenseNautyGraph(0), EncTranslator{Int}(), zeros(T, 0), SVector{2,F}[], ψs, 0)
end

function canonize!(p::Polyform)
    canon_perm, n = NautyGraphs.canonize!(p.anatomy)
    permute!(p.translator, canon_perm)
    p.σ = n
    return
end

faces(::Type{T}, p::Polyform) where {T} = convert(Vector{T}, p.anatomy.labels)
faces(p::Polyform) = p.anatomy.labels
species(p::Polyform) = p.species

function Base.size(p::Polyform)
    e = ne(p.anatomy)
    v = nv(p.anatomy)
    return length(p.xs), (e - v) ÷ 2, v
end #  TODO: nbonds only valid in 2D
function Base.show(io::IO, p::Polyform{D,T,F}) where {D,T,F}
    let (n, m, _) = size(p)
        println(io, "Polyform{$D,$T,$F}[n=$n, k=$m]")
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

function Base.copy(p::Polyform)
    return Polyform(copy(p.anatomy), copy(p.translator), copy(p.species), copy(p.xs), copy(p.ψs), p.σ)
end
function Base.copy!(dest::Polyform{D,T,F}, src::Polyform{D,T,F}) where {D,T,F}
    copy!(dest.anatomy, src.anatomy)
    copy!(dest.translator, src.translator)
    copy!(dest.species, src.species)
    copy!(dest.xs, src.xs)
    copy!(dest.ψs, src.ψs)
    dest.σ = src.σ
    return dest
end

function create_monomer(geometry::AbstractGeometry{F},
                        species_idx::T,
                        face_labels::AbstractVector{<:Integer}) where {T<:Integer,F<:Real}

    if geometry_dimension(geometry) == 2
        xs = [Point{2,F}(0., 0.)]
        ψs = [Angle{F}(0.)]
    elseif geometry_dimension(geometry) == 3
        xs = [Point{3,F}(0., 0., 0.)]
        ψs = [Quaternion{F}(1., 0., 0., 0.)]
    else
        error("Polyforms are only implemented for dimension 2 and 3.")
    end

    return Polyform(DirectedDenseNautyGraph(geometry.anatomy, face_labels),
                    EncTranslator([collect(1:length(geometry))]),
                    [species_idx],
                    xs, ψs)
end