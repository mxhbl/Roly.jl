using NautyGraphs
using StaticArrays

HashType = UInt
DefInt, DefFloat = Int16, Float32
SideLoc{T} = Tuple{T,T}  # Location of a binding site within a structure, in the form (particle, side)
CVec{N,F} = Union{SVector{N,F},MVector{N,F}}

function quaternion_multiply(ψi::AbstractVector, ψj::AbstractVector)
    a, b, c, d = ψi
    x, y, z, w = ψj

    o = a*x - b*y - c*z - d*w
    i = a*y + b*x + c*w - d*z
    j = a*z - b*w + c*x + d*y
    k = a*w + b*z - c*y + d*x
    return SVector(o, i, j, k)
end
quaternion_inv(ψ::AbstractVector) = SVector(ψ[1], -ψ[2], -ψ[3], -ψ[4])

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