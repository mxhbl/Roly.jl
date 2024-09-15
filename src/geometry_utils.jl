using StaticArrays
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
# TODO: take symmetry into account
Base.:(==)(a::Quaternion, b::Quaternion) = a.w == b.w && a.x == b.x && a.y == b.y && a.z == b.z
Base.isapprox(a::Quaternion, b::Quaternion; kwargs...) = isapprox(a.w, b.w; kwargs...) && isapprox(a.x, b.x; kwargs...) && isapprox(a.y, b.y; kwargs...) && isapprox(a.z, b.z; kwargs...)
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
function rotate!(ψs::AbstractVector{<:RotationOperator{F}}, ϕ::RotationOperator{F}) where {F}
    ψs .= Ref(ϕ) .* ψs
    return
end

function shift!(xs::AbstractVector{<:Point}, Δx::Point)
    xs .+= Ref(Δx)
    return
end

function grab!(xs::AbstractVector{<:Point}, ψs::AbstractVector{<:RotationOperator}, i::Integer, Δx=nothing, Δψ=nothing)
    # Shifts and rotates a list of xs and ψs such that particle i is at the origin in reference orientation
    Δψ = isnothing(Δψ) ? inv(ψs[i]) : Δψ * inv(ψs[i])

    shift!(xs, -xs[i])
    rotate!(xs, Δψ)
    rotate!(ψs, Δψ)

    !isnothing(Δx) && shift!(xs, Δx)
    return
end


function cart2pol(x::F, y::Real) where {F}
    y = convert(F, y)
    return [sqrt(x^2 + y^2), atan(y, x) / π]
end
function pol2cart(r::F, ψ::Real) where {F}
    ψ = convert(F, ψ)
    return [r * cos(π * ψ),  r * sin(π * ψ)]
end