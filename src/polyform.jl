using LinearAlgebra
using StaticArrays
using Graphs
using NautyGraphs

mutable struct Polyform{D,T<:Integer,F<:AbstractFloat,R<:RotationOperator{F}}
    anatomy::NautyDiGraph
    bond_partners::Vector{T}
    canonical_order::Vector{T}
    species::Vector{T}
    xs::Vector{Point{D,F}}
    ψs::Vector{R}
    σ::T
    function Polyform{D,T,F}(anatomy, bond_status, canonical_order, species, xs, ψs, σ) where {D,T,F}
        if D == 2
            return new{D,T,F,Angle{F}}(anatomy, bond_status, canonical_order, species, xs, ψs, σ)
        elseif D == 3
            return new{D,T,F,Quaternion{F}}(anatomy, bond_status, canonical_order, species, xs, ψs, σ)
        else
            error("Polyforms are only defined in dimension 2 and 3.")
        end
    end
end
function Polyform{D,T,F}(anatomy, bond_status, canonical_order, species, xs, ψs) where {D,T,F}
    p = Polyform{D,T,F}(anatomy, bond_status, canonical_order, species, xs, ψs, 0)
    canonize!(p)
    return p
end 
Polyform(anatomy, bond_status, canonical_order, species::AbstractVector{T}, xs::AbstractVector{<:Point{D,F}}, ψs, σ) where {D,T,F} = 
Polyform{D,T,F}(anatomy, bond_status, canonical_order, species, xs, ψs, σ)
Polyform(anatomy, bond_status, canonical_order, species::AbstractVector{T}, xs::AbstractVector{<:Point{D,F}}, ψs) where {D,T,F} = 
Polyform{D,T,F}(anatomy, bond_status, canonical_order, species, xs, ψs)
function Polyform{D,T,F}() where {D,T,F}
    if D == 2
        ψs = Angle{F}[]
    else
        ψs = Quaternion{F}[]
    end
    return Polyform{D,T,F}(NautyDiGraph(0), zeros(T, 0), zeros(T, 0), zeros(T, 0), SVector{2,F}[], ψs, 0)
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
nparticles(p::Polyform) = size(p)
nvertices(p::Polyform) = nv(p.anatomy)
function Base.show(io::Core.IO, p::Polyform{D,T,F}) where {D,T,F}
    print(io, "Polyform{$D,$T,$F}[n=$(size(p))]")
end
Base.show(io::Core.IO, ::Type{Polyform{D,T,F}}) where {D,T,F} = print(io, "Polyform{$D,$T,$F}")

function canonize!(p::Polyform)
    n, _, canon_perm, _ = NautyGraphs._nautyhash(p.anatomy)
    p.canonical_order .= canon_perm
    p.σ = n
    return
end

function rhash(p::Polyform)
    return ghash(p.anatomy)
end
is_isomorphic(p::Polyform, h::Polyform) = p ≃ h
≃(p::Polyform, h::Polyform) = NautyGraphs.is_isomorphic(p.anatomy, h.anatomy)


function Base.copy(p::Polyform)
    return Polyform(copy(p.anatomy), copy(p.bond_partners), copy(p.canonical_order), copy(p.species), copy(p.xs), copy(p.ψs), p.σ)
end
function Base.copy!(dest::Polyform{D,T,F}, src::Polyform{D,T,F}) where {D,T,F}
    copy!(dest.anatomy, src.anatomy)
    copy!(dest.bond_partners, src.bond_partners)
    copy!(dest.canonical_order, src.canonical_order)
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

    return Polyform(NautyDiGraph(geometry.anatomy, convert(Vector{T}, face_labels)),
                                 zeros(T, nvertices(geometry)),
                                 convert(Vector{T}, 1:nvertices(geometry)),
                                 [T(species_idx)],
                                 xs, ψs)
end

function particle2vertex(p::Polyform, assembly_system, particle::Integer, site::Integer)
    @assert particle <= nparticles(p)
    @assert site <= nsites(assembly_system.geometries[p.species[particle]])
    vps = (vertices_per_site(assembly_system.geometries[i]) for i in @view p.species[1:particle])
    return 1 + sum(sum(nvs) for nvs in Iterators.take(vps, particle-1); init=0) + sum(@view last(vps)[1:site-1])
end
function particle2vertex(p::Polyform, assembly_system, particle::Integer)
    @assert particle <= nparticles(p)
    vps = (vertices_per_site(assembly_system.geometries[i]) for i in @view p.species[1:particle])
    vertices = sum(sum(nvs) for nvs in Iterators.take(vps, particle-1); init=0) .+ cumsum((1, last(vps)[1:end-1]...))
    return vertices
end

function particle2multivertex(p::Polyform, assembly_system, particle::Integer, site::Integer)
    @assert particle <= nparticles(p)
    @assert site <= nsites(assembly_system.geometries[p.species[particle]])
    vps = (vertices_per_site(assembly_system.geometries[i]) for i in @view p.species[1:particle])
    start_vertex = 1 + sum(sum(nvs) for nvs in Iterators.take(vps, particle-1); init=0) + sum(@view last(vps)[1:site-1])
    end_vertex = start_vertex + last(vps)[site] - 1
    return start_vertex:end_vertex
end
function particle2multivertex(p::Polyform, assembly_system, particle::Integer)
    @assert particle <= nparticles(p)
    vps = (vertices_per_site(assembly_system.geometries[i]) for i in @view p.species[1:particle])
    start_vertex = 1 + sum(sum(nvs) for nvs in Iterators.take(vps, particle-1); init=0)
    end_vertex = start_vertex + sum(last(vps)) - 1
    return start_vertex:end_vertex
end
function vertex2particle(p::Polyform, assembly_system, vertex::Integer)
    @assert vertex <= nvertices(p)
    vps = (vertices_per_site(assembly_system.geometries[i]) for i in p.species)
    particle, site = 0, 0
    for vp in vps
        particle += 1
        for v in vp
            site += 1
            vertex -= v
            if vertex <= 0
                return particle, site
            end
        end
        site = 0
    end
    return
end
function vertex2label(p::Polyform, assembly_system, vertex::Integer)
    particle, site = vertex2particle(p, assembly_system, vertex)
    spcs = species(p)[particle]
    return spcssite2label(spcs, site, assembly_system)
end