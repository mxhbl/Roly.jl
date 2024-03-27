using StaticArrays
using Graphs
using LinearAlgebra

abstract type AbstractCoordinates{F<:Real} end

struct Coordinates2D{F<:AbstractFloat} <: AbstractCoordinates{F}
    xs::Vector{MVector{2,F}}
    ψs::Vector{F}
end
function Coordinates2D(xs::AbstractVector{SVector{2,F}}, ψs::Vector) where {F}
    return Coordinates2D{F}(xs, convert.(F, ψs))
end

struct Coordinates3D{F<:AbstractFloat} <: AbstractCoordinates{F}
    xs::Vector{MVector{3,F}}
    ψs::Vector{MVector{4,F}}
end
function Coordinates3D(xs::AbstractVector{SVector{3,F}}, ψs::AbstractVector) where {F}
    return Coordinates2D{F}(xs, convert.(F, ψs))
end

Base.length(cp::AbstractCoordinates) = length(cp.xs)
Base.getindex(cp::AbstractCoordinates, i::Integer) = (cp.xs[i], cp.ψs[i])
function Base.iterate(cp::AbstractCoordinates, state=1)
    return state > length(cp) ? nothing : (cp[state], state + 1)
end
# Base.copy(qs::Coordinates2D{F}) where {F} = Coordinates2D{F}([copy(x) for x in qs.xs], copy(qs.ψs))
function Base.copy!(dest::AbstractCoordinates{F}, src::AbstractCoordinates{F}) where {F}
    copy!(dest.xs, src.xs)
    copy!(dest.ψs, src.ψs)
    return dest
end
Base.deepcopy(qs::C) where {C<:AbstractCoordinates} = C([copy(x) for x in qs.xs], [copy(ψ) for ψ in qs.ψs])
Base.deepcopy(qs::Coordinates2D{F}) where {F} = Coordinates2D{F}([copy(x) for x in qs.xs], copy(qs.ψs))

#TODO: make a non-deep version of vcat to save some memory
#TODO: watch out for type stability here
function Base.vcat(qi::C, qj::D) where {C<:AbstractCoordinates, D<:AbstractCoordinates}
    return C(vcat(qi.xs, qj.xs), vcat(qi.ψs, qj.ψs))
end

function rotate(v::CVec{2,F}, ψ::F) where {F}
    c, s = cos(π * ψ), sin(π * ψ)
    return @SVector [c*v[1] - s*v[2], s*v[1] + c*v[2]]
end
function rotate(v::CVec{3,F}, ψ::CVec{4,F}) where {F}
    vq = SVector(0, v[1], v[2], v[3])
    return quaternion_multiply(quaternion_multiply(ψ, vq), quaternion_inv(ψ))[2:end]
end

function rotate!(qs::Coordinates2D{F}, ψ::F) where {F}
    c, s = cos(π * ψ), sin(π * ψ)
    for i in eachindex(qs.xs)
        x = qs.xs[i][1]
        y = qs.xs[i][2]
        qs.xs[i][1] = c * x - s * y
        qs.xs[i][2] = s * x + c * y

        qs.ψs[i] = mod.(qs.ψs[i] .+ ψ, F(2.0))
    end
end
function rotate!(qs::Coordinates3D{F}, ψ::CVec{4, F}) where {F}
    for i in eachindex(qs.xs)
        qs.xs[i] .= rotate(qs.xs[i], ψ)
        qs.ψs[i] = quaternion_multiply(ψ, qs.ψs[i])
    end
end

function shift!(qs::AbstractCoordinates, Δx::AbstractVector{F}) where {F}
    for i in eachindex(qs.xs)
        @views qs.xs[i] .+= Δx
    end
    return
end

function grab_at!(qs::AbstractCoordinates, i::Integer) where {F}
    # Transforms qs such that particle i it at the origin with reference orientation
    shift!(qs, deepcopy(-qs.xs[i]))
    rotate!(qs, deepcopy(-qs.ψs[i]))
    return
end

abstract type AbstractGeometry{T<:Integer,F<:AbstractFloat} end
Base.length(geom::AbstractGeometry) = length(geom.xs)

struct PolygonGeometry{T<:Integer,F<:AbstractFloat} <: AbstractGeometry{T,F}
    xs::Vector{SVector{2,F}}                # Displacement vectors from the center to the binding sites
    rs::Vector{F}                           # Distances of the side midpoints from center
    θs::Vector{Rational{T}}                 # Angles of the side midpoints from reference axis (in units of π)
    ϕs::Vector{F}                           # Angle between side and side-midpoint line (in units of π)
    corners::Vector{SVector{2,F}}
    anatomy::DirectedDenseNautyGraph{Cint}  # Oriented graph of the dual polyhedron associated with the bb symmetry group
    R_max::F

    function PolygonGeometry(n::T, a::F) where {T<:Integer,F<:AbstractFloat}
        r_in = convert(F, 0.5a * cot(π / n))
        r_out = convert(F, 0.5a * csc(π / n))

        θs = [mod(-T(2) * i // n - T(1) // T(2), T(2)) for i in T(0):(n - T(1))]
        ϕs = [T(1) // T(2) for _ in T(1):n]

        rs = r_in * ones(T, n)
        corners = [convert.(F, r_out * [cos(π * (θ + T(1) // n)), sin(π * (θ + T(1) // n))])
                   for θ in θs]

        R_max = maximum(norm.(corners))

        anatomy = DirectedDenseNautyGraph(cycle_digraph(n))

        xs = [SVector{2,F}(pol2cart(r, θ)) for (r, θ) in zip(rs, θs)]
        return new{T,F}(xs, rs, θs, ϕs, corners, anatomy, R_max)
    end
end

const UnitTriangleGeometry = PolygonGeometry(DefInt(3), DefFloat(1.0))
const UnitSquareGeometry = PolygonGeometry(DefInt(4), DefFloat(1.0))
const UnitPentagonGeometry = PolygonGeometry(DefInt(5), DefFloat(1.0))
const UnitHexagonGeometry = PolygonGeometry(DefInt(6), DefFloat(1.0))


# TODO:Remove T from all geometries
struct PolyhedronGeometry{T<:Integer,F<:AbstractFloat} <: AbstractGeometry{T,F}
    xs::Vector{SVector{3,F}}                # Displacement vectors from the center to the binding sites
    rs::Vector{F}                           # Distances of the side midpoints from center
    θs_ref::Vector{SVector{4,F}}             
    Δθs::Vector{SVector{4,F}}              
    corners::Vector{SVector{3,F}}
    anatomy::DirectedDenseNautyGraph{Cint}  # Oriented graph of the dual polyhedron associated with the bb symmetry group
    R_max::F

    function PolyhedronGeometry(shape::Symbol, a::F) where {F<:AbstractFloat}
        if shape == :cube
            xs = a * [SVector{3,F}(0, 0, -1),
                      SVector{3,F}(0, 0, 1),
                      SVector{3,F}(-1, 0, 0),
                      SVector{3,F}(1, 0, 0),
                      SVector{3,F}(0, -1, 0),
                      SVector{3,F}(0, 1, 0)]

            rs = a * ones(F, 6)
            corners = vec([a * SVector(i, j, k) for i in F[-1, 1], j in F[-1, 1], k in F[-1, 1]])
            θs_ref = [SVector{4, F}(1, 0, 0, 0),
                      SVector{4, F}(0, 1, 0, 0),
                      SVector{4, F}(1, 0, 1, 0) / √2,
                      SVector{4, F}(1, 0, -1, 0) / √2,
                      SVector{4, F}(1, 1, 0, 0) / √2,
                      SVector{4, F}(1, -1, 0, 0) / √2]
            Δθs = [SVector{4, F}(0, 1, 0, 0),
                   SVector{4, F}(1, 0, 0, 0),
                   SVector{4, F}(1, 0, -1, 0) / √2,
                   SVector{4, F}(1, 0, 1, 0) / √2,
                   SVector{4, F}(1, -1, 0, 0) / √2,
                   SVector{4, F}(1, 1, 0, 0) / √2]

            adjmx = [0 0 1 1 0 0;
                     0 0 1 1 0 0;
                     0 0 0 0 1 1;
                     0 0 0 0 1 1;
                     1 1 0 0 0 0;
                     1 1 0 0 0 0]
            anatomy = DirectedDenseNautyGraph(adjmx)
            R_max = maximum(rs)
        else
            error("Shape $shape not supported.")
        end
        
        return new{DefInt,F}(xs, rs, θs_ref, Δθs, corners, anatomy, R_max)
    end
end

const UnitCubeGeometry = PolyhedronGeometry(:cube, DefFloat(1.))
#TODO Implement basic 3d shapes: Platonic solids and polygon extrusions

# struct MeshGeometry{T<:Int, R<:Real}<:AbstractGeometry
#     rs::Vector{R}               # Distances of the side midpoints from center
#     θs::Vector{R}               # Angles of the side midpoints from reference axis (in units of π)
#     ϕs::Vector{R}               # Angle between side and side-midpoint line (in units of π)
#     mesh::Vector{SVector{2,R}}  # Mesh of the building block (used for overlap checks)
#     anatomy::DiGraph{T}           # Oriented graph of the dual polyhedron associated with the bb symmetry group
# end




function attachment_offset(face_i::Integer, face_j::Integer,
                           geom_i::PolygonGeometry{T,F},
                           geom_j::PolygonGeometry) where {T,F}
    x_i, θ_i, ϕ_i = geom_i.xs[face_i], geom_i.θs[face_i], geom_i.ϕs[face_i]
    x_j, θ_j, ϕ_j = geom_j.xs[face_j], geom_j.θs[face_j], geom_j.ϕs[face_j]

    Δx = x_i + rotate(x_j, θ_i - θ_j + ϕ_i - ϕ_j)
    Δψ = mod(θ_i - θ_j + ϕ_i - ϕ_j + T(1), T(2))
    return Δx, Δψ
end
function attachment_offset(face_i::Integer, face_j::Integer,
                           geom_i::PolyhedronGeometry{T,F},
                           geom_j::PolyhedronGeometry) where {T,F}
    x_i, Δθ = geom_i.xs[face_i], geom_i.Δθs[face_i]
    x_j, θ_ref = geom_j.xs[face_j], geom_j.θs_ref[face_j]

    Δx = x_i - rotate(x_j, quaternion_multiply(Δθ, θ_ref))
    Δψ = quaternion_multiply(Δθ, θ_ref)
    return Δx, Δψ
end

# TODO: this function is susceptible to roundoff errors, we very careful with thresh
function face_pairs(Δx::AbstractVector{F}, ψ_i, ψ_j, geom_i::AbstractGeometry{T,F}, geom_j::AbstractGeometry,
                    convex=false, thresh=1e-2) where {T,F}
    pairs = Tuple{T,T}[]

    for (i, z_i) in enumerate(geom_i.xs), (j, z_j) in enumerate(geom_j.xs)
        z_i = rotate(z_i, ψ_i)
        z_j = rotate(z_j, ψ_j)

        matching = maximum(abs, z_i - z_j - Δx) < F(thresh)
        if matching
            push!(pairs, (T(i), T(j)))
            if convex
                break
            end
        end
    end
    return pairs
end

function overlap(Δx, ψ_i, ψ_j, geom_i, geom_j; buffer=0.1)
    # Assumes lattice!!
    return norm(Δx) < 2geom_i.rs[1] - buffer
    # return __check_overlap(geom_i.mesh, transform(geom_j.mesh, Δx, Δψ))
end

function contact_status(Δx::AbstractVector{F}, ψ_i, ψ_j, geom_i::AbstractGeometry{T,F}, geom_j::AbstractGeometry{T,F}) where {T,F}
    if norm(Δx) > geom_i.R_max + geom_j.R_max
        return true, Tuple{T,T}[]
    end

    if overlap(Δx, ψ_i, ψ_j, geom_i, geom_j)
        return false, Tuple{T,T}[]
    end

    pairs = face_pairs(Δx, ψ_i, ψ_j, geom_i, geom_j)
    return true, pairs
end