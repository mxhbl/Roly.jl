using StaticArrays
using Graphs
using LinearAlgebra

Point{D,F} = SVector{D,F}
abstract type RotationOperator{F<:Real} end

struct Angle{F} <: RotationOperator{F}
    θ::F
    Angle{F}(θ) where {F} = new{F}(mod(θ, 2))
end
Angle(θ::F) where {F<:Real} = Angle{F}(θ)
Base.:(==)(a::Angle, b::Angle) = a.θ == b.θ
Base.isapprox(a::Angle, b::Angle) = isapprox(a.θ, b.θ)
Base.:*(ai::Angle{F}, aj::Angle) where {F} = Angle{F}(ai.θ + aj.θ)
Base.inv(a::Angle{F}) where {F} = Angle{F}(-a.θ)
Base.convert(::Type{Angle{F}}, a::Real) where {F} = Angle{F}(a)
Base.convert(::Type{Angle{F}}, a::Angle) where {F} = Angle{F}(a.θ)

struct Quaternion{F} <: RotationOperator{F}
    w::F
    x::F
    y::F
    z::F
end
Quaternion(w::F, x::F, y::F, z::F) where {F<:Real} = Quaternion{F}(w, x, y, z)
Quaternion(w::Real, x::Real, y::Real, z::Real) = Quaternion(promote(w, x, y, z)...)
Base.convert(::Type{Quaternion{F}}, v::AbstractVector) where {F} = Quaternion{F}(v[1:4]...)
Base.convert(::Type{Quaternion{F}}, q::Quaternion) where {F} = Quaternion{F}(q.w, q.x, q.y, q.z)
# Base.:(==)(a::Quaternion, b::Quaternion) = ...
# Base.isapprox(a::Quaternion, b::Quaternion) = ...
Base.inv(q::Quaternion) = Quaternion(q.w, -q.x, -q.y, -q.z)
LinearAlgebra.norm(q::Quaternion) = sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z)

function Base.:*(qi::Quaternion, qj::Quaternion)
    a, b, c, d = qi.w, qi.x, qi.y, qi.z
    w, x, y, z = qj.w, qj.x, qj.y, qj.z

    o = a*w - b*x - c*y - d*z
    i = a*x + b*w + c*z - d*y
    j = a*y - b*z + c*w + d*x
    k = a*z + b*y - c*x + d*w
    return Quaternion(o, i, j, k)
end
Base.:*(α::Real, q::Quaternion) = Quaternion(α * q.w, α * q.x, α * q.y, α * q.z)
Base.:/(q::Quaternion, α::Real) = inv(α) * q

normalized(q::Quaternion) = inv(norm(q)) * q

function rotate(x::Point{2,F}, ϕ::Angle{F}) where {F}
    c, s = cospi(ϕ.θ), sinpi(ϕ.θ) #TODO: maybe refactor this into a method
    return Point{2,F}(c * x[1] - s * x[2], s * x[1] + c * x[2])
end
function rotate(x::Point{3,F}, ϕ::Quaternion{F}) where {F}
    xq = Quaternion(0, x[1], x[2], x[3])
    xq = ϕ * xq * inv(ϕ)
    return Point{3,F}(xq.x, xq.y, xq.z)
end

function rotate!(xs::AbstractVector{<:Point{2,F}}, ϕ::Angle{F}) where {F}
    c, s = cospi(ϕ.θ), sinpi(ϕ.θ)
    for i in eachindex(xs)
        x, y = xs[i]
        xs[i] = Point{2,F}(c * x - s * y, s * x + c * y)
    end
end
function rotate!(xs::AbstractVector{<:Point{3,F}}, ϕ::Quaternion{F}) where {F}
    for i in eachindex(xs)
        xs[i] = rotate(xs[i], ϕ)
    end
end

function shift!(xs::AbstractVector{<:Point}, Δx::Point)
    for i in eachindex(xs)
        @views xs[i] .+= Δx
    end
    return
end

function grab_at!(xs::AbstractVector{<:Point}, ψs::AbstractVector{<:RotationOperator}, i::Integer)
    # Shifts and rotates a list of xs and ψs such that particle i is at the origin in reference orientation
    Δx = -xs[i]
    ϕ = inv(ψs[i])

    shift!(xs, Δx)
    rotate!(xs, ϕ)
    ψs .= ϕ .* ψs
    return
end

abstract type AbstractGeometry{F<:AbstractFloat} end
Base.length(geom::AbstractGeometry) = length(geom.xs)

struct PolygonGeometry{F<:AbstractFloat} <: AbstractGeometry{F}
    xs::Vector{Point{2,F}}                # Displacement vectors from the center to the binding sites
    θs_ref::Vector{Angle{F}}              # Rotation necessary to rotate side i into reference position
    Δθs::Vector{Angle{F}}                 # Rotation to be applied to the added particle (assumed to be in reference orientation) once attached to side i 
    corners::Vector{Point{2,F}}
    anatomy::DirectedDenseNautyGraph{Cint}  # Oriented graph of the dual polyhedron associated with the bb symmetry group
    R_min::F
    R_max::F

    function PolygonGeometry(n::Integer, a::F) where {F}
        r_in = convert(F, 0.5a * cot(π / n))
        r_out = convert(F, 0.5a * csc(π / n))

        # θs = [mod(-T(2) * i / n - T(1) / T(2), T(2)) for i in T(0):(n - T(1))]
        # ϕs = [T(1) / T(2) for _ in T(1):n]
        θs_ref = [Angle{F}(2/n * i) for i in 0:(n-1)]
        Δθs = [Angle{F}(1 - 2/n*i) for i in 0:(n-1)]

        rs = r_in * ones(F, n)
        corners = [convert.(F, r_out * [cos(π * (θ.θ + 1/n)), sin(π * (θ.θ + 1/n))]) for θ in θs_ref] #TODO: clean up the θ.θ

        xs = [SVector{2,F}(pol2cart(r, -θ.θ - 1/2)) for (r, θ) in zip(rs, θs_ref)]
        anatomy = DirectedDenseNautyGraph(cycle_digraph(n))

        R_min = minimum(norm.(xs))
        R_max = maximum(norm.(corners))

        # θs = [Angle(θ) for θ in θs]
        # ϕs = [Angle(ϕ) for ϕ in ϕs]

        return new{F}(xs, θs_ref, Δθs, corners, anatomy, R_min, R_max)
    end
end

const UnitTriangleGeometry = PolygonGeometry(3, DefFloat(1.0))
const UnitSquareGeometry = PolygonGeometry(4, DefFloat(1.0))
const UnitPentagonGeometry = PolygonGeometry(5, DefFloat(1.0))
const UnitHexagonGeometry = PolygonGeometry(6, DefFloat(1.0))


struct PolyhedronGeometry{F<:AbstractFloat} <: AbstractGeometry{F}
    xs::Vector{Point{3,F}}                # Displacement vectors from the center to the binding sites
    θs_ref::Vector{Quaternion{F}}         # Rotation necessary to rotate side i into reference position   
    Δθs::Vector{Quaternion{F}}            # Rotation to be applied to the added particle (assumed to be in reference orientation) once attached to side i 
    corners::Vector{Point{3,F}}
    anatomy::DirectedDenseNautyGraph{Cint}  # Oriented graph of the dual polyhedron associated with the bb symmetry group
    R_min::F
    R_max::F

    function PolyhedronGeometry(shape::Symbol, a::F) where {F<:AbstractFloat}
        if shape == :cube
            xs = a * [Point{3,F}(0, 0, -1),
                      Point{3,F}(0, 0, 1),
                      Point{3,F}(-1, 0, 0),
                      Point{3,F}(1, 0, 0),
                      Point{3,F}(0, -1, 0),
                      Point{3,F}(0, 1, 0)]

            corners = vec([a * SVector(i, j, k) for i in F[-1, 1], j in F[-1, 1], k in F[-1, 1]])
            θs_ref = [Quaternion{F}(1, 0, 0, 0),
                      Quaternion{F}(0, 1, 0, 0),
                      inv(√2) * Quaternion{F}(1, 0, 1, 0),
                      inv(√2) * Quaternion{F}(1, 0, -1, 0),
                      inv(√2) * Quaternion{F}(1, 1, 0, 0),
                      inv(√2) * Quaternion{F}(1, -1, 0, 0)]
            Δθs = [Quaternion{F}(0, 1, 0, 0),
                   Quaternion{F}(1, 0, 0, 0),
                   inv(√2) * Quaternion{F}(1, 0, -1, 0),
                   inv(√2) * Quaternion{F}(1, 0, 1, 0),
                   inv(√2) * Quaternion{F}(1, -1, 0, 0),
                   inv(√2) * Quaternion{F}(1, 1, 0, 0)]

            adjmx = [0 0 1 1 0 0;
                     0 0 1 1 0 0;
                     0 0 0 0 1 1;
                     0 0 0 0 1 1;
                     1 1 0 0 0 0;
                     1 1 0 0 0 0]
            anatomy = DirectedDenseNautyGraph(adjmx)

            R_min = minimum(norm.(xs))
            R_max = maximum(norm.(corners))
        else
            error("Shape $shape not supported.")
        end
        
        return new{F}(xs, θs_ref, Δθs, corners, anatomy, R_min, R_max)
    end
end

const UnitCubeGeometry = PolyhedronGeometry(:cube, DefFloat(1.))
#TODO Implement basic 3d shapes: Platonic solids and polygon extrusions

geometry_dimension(::PolygonGeometry) = 2
geometry_dimension(::PolyhedronGeometry) = 3


function attachment_offset(face_i::Integer, face_j::Integer,
                           geom_i::AbstractGeometry,
                           geom_j::AbstractGeometry)
    x_i, Δθ = geom_i.xs[face_i], geom_i.Δθs[face_i]
    x_j, θ_ref = geom_j.xs[face_j], geom_j.θs_ref[face_j]

    Δψ = Δθ * θ_ref
    Δx = x_i - rotate(x_j, Δψ)
    return Δx, Δψ
end

# TODO: this function is susceptible to roundoff errors, be very careful with thresh
function face_pairs(Δx::Point{D,F}, ψi::RotationOperator{F}, ψj::RotationOperator{F}, geom_i::AbstractGeometry{F}, geom_j::AbstractGeometry, convex=false, thresh=1e-2) where {D,F}
    pairs = Tuple{Int,Int}[]

    for (i, z_i) in enumerate(geom_i.xs), (j, z_j) in enumerate(geom_j.xs)
        z_i = rotate(z_i, ψi)
        z_j = rotate(z_j, ψj)

        matching = maximum(abs, z_i - z_j - Δx) < F(thresh)
        if matching
            push!(pairs, (i, j))
            if convex
                break
            end
        end
    end
    return pairs
end

function overlap(Δx, ψ_i, ψ_j, geom_i, geom_j; buffer=0.1)
    # Assumes lattice!!
    return norm(Δx) < 2geom_i.R_min - buffer
    # return __check_overlap(geom_i.mesh, transform(geom_j.mesh, Δx, Δψ))
end

function contact_status(Δx::Point, ψi::RotationOperator, ψj::RotationOperator, geom_i::AbstractGeometry, geom_j::AbstractGeometry)
    if norm(Δx) > geom_i.R_max + geom_j.R_max
        return true, Tuple{Int,Int}[]
    end

    if overlap(Δx, ψi, ψj, geom_i, geom_j)
        return false, Tuple{Int,Int}[]
    end

    pairs = face_pairs(Δx, ψi, ψj, geom_i, geom_j)
    return true, pairs
end