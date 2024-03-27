using LinearAlgebra
using StaticArrays
using Graphs
using NautyGraphs

mutable struct EncTranslator{T}
    fwd::Vector{Vector{T}}
    bwd::Matrix{T}
    n::T
    m::T

    function EncTranslator(fwd::AbstractVector{<:AbstractVector{T}}) where {T}
        n = length(fwd)
        m = sum(length(f) for f in fwd; init=0)

        bwd = zeros(2, m)
        for i in 1:m
            for (j, f) in enumerate(fwd)
                if i ∈ f
                    bwd[:, i] .= [T(j), T(findfirst(x->x==i, f))]
                    break
                end
            end
        end

        return new{T}(fwd, bwd, n, m)
    end
    function EncTranslator(fwd::AbstractVector{<:AbstractVector{T}}, bwd::AbstractMatrix{T}) where {T}
        return new{T}(fwd, bwd, length(fwd), size(bwd, 2))
    end
end
EncTranslator{T}() where {T} = EncTranslator(Vector{T}[], zeros(T, (0, 0)))

function Base.permute!(t::EncTranslator, p::AbstractVector{<:Integer})
    @assert length(p) == t.m

    t.bwd .= @view t.bwd[:, p]

    inv_p = invperm(p)
    t.fwd .= [inv_p[f] for f in t.fwd]
    return
end


function Base.vcat(t::EncTranslator{T}, h::EncTranslator) where {T}
    fwd = vcat(t.fwd, [hf .+ t.m for hf in h.fwd]) #TODO: reduce allocation
    bwd = hcat(t.bwd, h.bwd .+ t.n*[T(1), T(0)])
    return EncTranslator(fwd, bwd)
end
function Base.deleteat!(t::EncTranslator, i::Integer)
    del_ai = t.fwd[i]

    ai_shifter(ai) = ai - sum(x->x < ai, del_ai) #TODO:check gt or gte
    k_shifter(k) = k - (k > i)

    deleteat!(t.fwd, i)
    for j in eachindex(t.fwd)
        t.fwd[j] .= ai_shifter.(t.fwd[j])
    end
    t.bwd = t.bwd[:, setdiff(1:end, del_ai)]
    t.bwd[1, :] .= k_shifter.(t.bwd[1, :])

    t.n -= 1
    t.m -= length(del_ai)
    return
end
Base.copy(t::EncTranslator) = EncTranslator([copy(f) for f in t.fwd], copy(t.bwd))
function Base.copy!(dest::EncTranslator{T}, src::EncTranslator{T}) where {T}
    # TODO: make this non-allocating
    # copy!(dest.fwd, src.fwd)
    resize!(dest.fwd, src.n)
    dest.fwd .= copy.(src.fwd)
    dest.bwd = copy(src.bwd)
    dest.n = src.n
    dest.m = src.m
    return dest
end

mutable struct Structure{T<:Integer,F<:AbstractFloat}
    anatomy::DirectedDenseNautyGraph{Cint}
    translator::EncTranslator{Int}
    species::Vector{T}
    positions::Coordinates2D{F}
    σ::T

    function Structure{T,F}(anatomy, translator, species, positions) where {T,F}
        s = new{T,F}(anatomy, translator, species, positions, 0)
        canonize!(s)
        return s
    end

    function Structure{T,F}(anatomy, translator, species, positions, σ) where {T,F}
        return new{T,F}(anatomy, translator, species, positions, σ)
    end
end
Structure(args...) = Structure{DefInt,DefFloat}(args...)
Structure{T,F}() where {T,F} = Structure{T,F}(DirectedDenseNautyGraph(0), EncTranslator{Int}(), zeros(T, 0), Coordinates2D{F}([], []), 0)

function canonize!(s::Structure)
    canon_perm, n = NautyGraphs.canonize!(s.anatomy)
    permute!(s.translator, canon_perm)
    s.σ = n
    return
end

faces(::Type{T}, s::Structure) where {T} = convert(Vector{T}, s.anatomy.labels)
faces(s::Structure) = faces(Int, s)
species(s::Structure) = s.species

function Base.size(s::Structure)
    e = ne(s.anatomy)
    v = nv(s.anatomy)
    return length(s.positions), (e - v) ÷ 2, v
end #  TODO: nbonds only valid in 2D
function Base.show(io::IO, s::Structure{T,F}) where {T,F}
    let (n, m, _) = size(s)
        println(io, "Structure{$T, $F}[n=$n, k=$m]")
    end
end

function Base.hash(s::Structure)
    return hash(s.anatomy)
end
function Base.hash(s::Structure, h::UInt64)
    return hash(s.anatomy, h)
end
function Base.:(==)(si::Structure, sj::Structure)
    return hash(si.anatomy) == hash(sj.anatomy)
end

function Base.copy(s::Structure{T,F}) where {T,F}
    return Structure{T,F}(copy(s.anatomy), copy(s.translator), copy(s.species), deepcopy(s.positions), s.σ)
end
function Base.copy!(dest::Structure{T,F}, src::Structure{T,F}) where {T,F}
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

    return Structure{T,F}(DirectedDenseNautyGraph(geometry.anatomy, face_labels),
                          EncTranslator([collect(1:length(geometry))]),
                          [species_idx],
                        #   Coordinates2D{F}([SA[0.0, 0.0, 0.0]], [SA[1.0, 0.0, 0.0, 0.0]]))
                          Coordinates2D{F}([SA[0.0, 0.0]], [0.0]))

end