using NautyGraphs
using StaticArrays
using ElasticArrays

HashType = UInt
DefInt, DefFloat = Int16, Float32
SideLoc{T} = Tuple{T,T}  # Location of a binding site within a structure, in the form (particle, side)
CVec{N,F} = Union{SVector{N,F},MVector{N,F}}


function cart2pol(x::F, y::Real) where {F}
    y = convert(F, y)
    return [sqrt(x^2 + y^2), atan(y, x) / π]
end
function pol2cart(r::F, ψ::Real) where {F}
    ψ = convert(F, ψ)
    return [r * cos(π * ψ),  r * sin(π * ψ)]
end

function irg_flatten(a::Integer, b::Integer, intervals::AbstractVector{<:Integer})
    if a == 1
        return b
    else
        return intervals[a-1] + b
    end
end
function irg_unflatten(i::Integer, intervals::AbstractVector{<:Integer})
    if i <= intervals[1]
        return (1, i)
    else
        a = searchsortedlast(intervals, i, lt=≤) + 1
        b = i - intervals[a - 1]
        return (a, b)
    end
end

mutable struct PolyEncoder{T}
    fwd::Vector{Vector{Vector{T}}}
    bwd::Matrix{T}
    n_particles::T
    n_vertices::T

    function PolyEncoder(fwd::AbstractVector{<:AbstractVector{<:AbstractVector{T}}}, bwd::AbstractMatrix{T}) where {T}
        return new{T}(fwd, bwd, length(fwd), size(bwd, 2))
    end
end
PolyEncoder{T}() where {T} = PolyEncoder(Vector{Vector{T}}[], zeros(T, (3, 0)))
function PolyEncoder(fwd::AbstractVector{<:AbstractVector{<:AbstractVector{T}}}) where {T}
    n_vertices = sum(length(s) for f in fwd for s in f; init=0)

    bwd = zeros(T, 3, n_vertices)
    for i in 1:n_vertices
        for (j, f) in enumerate(fwd)
            for (k, s) in enumerate(f)
                if i ∈ s
                    bwd[:, i] .= [j, k, findfirst(x->x==i, s)]
                    break
                end
            end
        end
    end

    return PolyEncoder(fwd, bwd)
end

function Base.permute!(enc::PolyEncoder, p::AbstractVector{<:Integer})
    @assert length(p) == enc.n_vertices

    enc.bwd .= @view enc.bwd[:, p]

    inv_p = invperm(p)
    for part in enc.fwd
        for i in eachindex(part)
            for j in eachindex(part[i])
                part[i][j] = inv_p[part[i][j]]
            end
        end
    end
    return
end

# TODO: the following two methods should be implemented with better memory management,
# With as few allocations as possible
function concatenate(enc::PolyEncoder{T}, h::PolyEncoder) where {T}
    ne = enc.n_particles
    me = enc.n_vertices
    nh = h.n_particles

    fwd = similar(enc.fwd, ne+nh)
    for i in eachindex(enc.fwd)
        fwd[i] = [copy(s) for s in enc.fwd[i]]
    end
    for i in eachindex(h.fwd)
        fwd[i+ne] = [copy(s) .+ me for s in h.fwd[i]]
    end

    bwd = hcat(enc.bwd, h.bwd .+ ne * [one(T), zero(T), zero(T)])
    return PolyEncoder(fwd, bwd)
end
function Base.deleteat!(enc::PolyEncoder, i::Integer)
    #Flatten all to-be-deleted indices into one vector
    del_vs = [v for verts in enc.fwd[i] for v in verts]

    vertex_shift(v) = sum(x -> x < v, del_vs)
    particle_shift(k) = (k > i)

    deleteat!(enc.fwd, i)
    for part in enc.fwd
        for j in eachindex(part)
            @views part[j] .-= vertex_shift.(part[j])
        end
    end
    enc.bwd = enc.bwd[:, setdiff(1:end, del_vs)] #TODO: optimize this
    @views enc.bwd[1, :] .-= particle_shift.(enc.bwd[1, :])

    enc.n_particles -= 1
    enc.n_vertices -= length(del_vs)
    return
end
# Base.copy(enc::PolyEncoder) = PolyEncoder(deepcopy(enc.fwd), copy(enc.bwd)) # TODO: watch the deepcopy
function Base.copy!(dest::PolyEncoder{T}, src::PolyEncoder{T}) where {T}
    resize!(dest.fwd, src.n_particles)
    for i in eachindex(dest.fwd)
        if isassigned(dest.fwd, i)
            for j in eachindex(dest.fwd[i])
                dest.fwd[i][j] .= src.fwd[i][j]
            end
        else
            dest.fwd[i] = copy.(src.fwd[i])
        end
    end
    # dest.fwd .= [copy.(part) for part in src.fwd]
    
    # resize!(dest.bwd, size(src.bwd))
    dest.bwd = copy(src.bwd)
    dest.n_particles = src.n_particles
    dest.n_vertices = src.n_vertices
    return dest
end


function are_bridge(g::AbstractNautyGraph, vs::AbstractVector{<:Integer})
    neighs = zeros(Bool, nv(g))
    buffer = zeros(Int, nv(g))

    for v in vs
        k = NautyGraphs.outneighbors!(buffer, g, v)
        for neigh in @view buffer[1:k]
            neighs[neigh] = true
        end
    end
    neighs[vs] .= false

    if length(neighs) < 2
        return false
    end

    forbidden = zeros(Bool, nv(g))
    forbidden[vs] .= true

    return !vertices_connected(g, findfirst(neighs), findall(neighs), forbidden)
end
function vertices_connected(g::DirectedDenseNautyGraph, v0::Integer, 
    targets::AbstractVector{<:Integer}, forbidden::AbstractVector{<:Integer})

    n = nv(g)
    neighs = zeros(Int, n)
    explored = zeros(Bool, n)
    
    queue = zeros(Cint, n)
    queue[1] = v0
    q_start = 1
    q_end = 2

    explored[v0] = true
    while q_end > q_start
        v = queue[q_start]
        q_start += 1
        n_neighs = NautyGraphs.outneighbors!(neighs, g, v)
        for neigh in @view neighs[1:n_neighs]
            if forbidden[neigh]
                continue
            end

            if !explored[neigh]
                explored[neigh] = true

                if all(@view explored[targets])
                    return true
                end

                queue[q_end] = neigh
                q_end += 1                
            end
        end
    end

    return false
end