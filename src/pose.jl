"""
    Pose{D,F,R<:Rotation{D,F}}

A `Pose` describes the configuration of a `D`-dimensional rigid body.
It consists of translational coordinates `pose.x` and orientational
coordinates `pose.psi`.
"""
struct Pose{D,F,R<:Rotation{D,F}}
    x::SVector{D,F}
    psi::R
end

"""
    Pose(x, psi)

Create a pose from a position vector `x` and orientation `psi`.
"""
Pose(x::AbstractVector, psi::R) where {R<:Rotation} = Pose{length(x),eltype(x),R}(SVector{length(x),eltype(x)}(x), psi)
Pose{D,F}(x::AbstractVector, psi::R) where {D,F,R<:Rotation{D,F}} = Pose{D,F,R}(SVector{D,F}(x), psi)

# Rotation types are not closed under multiplication (`RotXYZ * RotXYZ` is a `RotMatrix3`)
Base.convert(::Type{Pose{D,F,R}}, p::Pose{D,F}) where {D,F,R<:Rotation{D,F}} = Pose{D,F,R}(p.x, convert(R, p.psi))

"""
    one(::Type{Pose{D,F,R}})
    one(p::Pose)

Return the identity pose (origin, identity rotation) of the given type.
"""
Base.one(::Type{Pose{D,F,R}}) where {D,F,R} = Pose{D,F,R}(zero(SVector{D,F}), one(R))
Base.one(p::Pose) = one(typeof(p))

"""
    Pose{2,F}(; orientationtype=Angle2d)
    Pose{3,F}(; orientationtype=RotMatrix3)
    Pose{D}(; ...)

Return the identity pose of eltype `F` in `D` dimensions, using the default rotation type
(`Angle2d` in 2D, `RotMatrix3` in 3D). The rotation type can be overridden with
`orientationtype`, but see the note on `normal_pose` before picking an Euler one.
"""
Pose{2,F}(; orientationtype=Angle2d) where {F} = one(Pose{2,F,orientationtype{F}})
Pose{3,F}(; orientationtype=RotMatrix3) where {F} = one(Pose{3,F,orientationtype{F}})
Pose{D,F}(; kwargs...) where {D,F} = throw(ArgumentError("no default pose for D=$D, use D=2 or D=3"))
Pose{D}(; kwargs...) where {D} = Pose{D,Float64}(; kwargs...)

"""
    Pose{3,F}(p::Pose{2,F})

Lift a 2D pose into 3D by appending a zero z-coordinate and wrapping the 2D rotation angle
into the z-axis rotation of the 3D rotation type (default `RotXYZ`).
"""
function Pose{3,F}(p::Pose{2,F}; orientationtype=RotXYZ) where {F}
    return Pose(SVector{3,F}(p.x..., F(0)), orientationtype(F(0), F(0), rotation_angle(p.psi)))
end
function Pose{3}(p::Pose{2,F}; orientationtype=RotXYZ) where {F}
    return Pose{3,F}(p; orientationtype)
end

Base.:(==)(p1::Pose, p2::Pose) = p1.x == p2.x && p1.psi == p2.psi
Base.isapprox(p1::Pose, p2::Pose; kwargs...) = isapprox(p1.x, p2.x; kwargs...) && isapprox(p1.psi, p2.psi; kwargs...)

Base.:*(p1::Pose, p2::Pose) = Pose(p1.psi * p2.x + p1.x, p1.psi * p2.psi)
Base.:*(p::Pose, v::AbstractVector) = p.psi * v + p.x
Base.:*(R::Rotation, p::Pose) = typeof(p)(R * p.x, R * p.psi)
Base.:*(p::Pose, R::Rotation) = typeof(p)(p.x, p.psi * R)

"""
    translate(x, v)

Return `x` shifted by the vector `v`, its orientation unchanged.

Defined for a [`Pose`](@ref), a `Particle`, a [`BindingSite`](@ref) and an
[`AbstractPolyform`](@ref). A general rigid motion is `g * x` instead, where `g` is a `Pose`: `*`
is the group action, and its two sides differ, `g * x` moving `x` in the world frame and `x * g`
in its own.
"""
translate(p::Pose, v) = typeof(p)(p.x + v, p.psi)

Base.inv(p::Pose) =
    let ψinv = inv(p.psi)
        typeof(p)(-ψinv * p.x, ψinv)
    end
Base.:/(p1::Pose, p2) = p1 * inv(p2)
Base.:\(p1, p2::Pose) = inv(p1) * p2

Base.Matrix(p::Pose{D,F}) where {D,F} = [p.psi p.x; zeros(F, D)' F(1)]

Base.eltype(::Type{<:Pose{D,F}}) where {D,F} = F

dimension(::Pose{D}) where {D} = D
dimension(::Type{<:Pose{D}}) where {D} = D

function Base.show(io::Core.IO, p::Pose{D,F,R}) where {D,F,R}
    if D == 2
        print(io, "[x: $(p.x); ψ: $(rotation_angle(p.psi) / π)π]")
    elseif D == 3
        print(io, "[x: $(p.x); ψ: $(rotation_angle(p.psi) / π)π, $(rotation_axis(p.psi))]")
    end
end
function Base.show(io::Core.IO, ::MIME"text/plain", p::Pose{D,F,R}) where {D,F,R}
    print(io, "$D-dimensional Pose{$D,$F,$(string(nameof(R)))}:\n")
    print(io, " - x: $(p.x)\n")
    if D == 2
        print(io, " - psi: $(rotation_angle(p.psi) / π)π")
    elseif D == 3
        print(io, " - psi: $(rotation_angle(p.psi) / π)π, $(rotation_axis(p.psi))")
    end
end

function cart2pol(x::F, y::Real) where {F}
    y = convert(F, y)
    return SVector(sqrt(x^2 + y^2), atan(y, x))
end
function pol2cart(r::F, psi::Real) where {F}
    psi = convert(F, psi)
    return SVector(r * cos(psi), r * sin(psi))
end

"""
    normal_pose(x, twist=0; orientationtype=RotMatrix3)

Create a pose at point `x`, whose local x axis also points outward along `x`. `twist` rotates
the y,z axes about x, which is how a patch's bond twist is chosen: the twist is the site's
reference direction, and turning it turns the partner it holds.
"""
# The 3D default is `RotMatrix3`, not an Euler parameterisation, because rotation types are not
# closed under multiplication: `RotXYZ * RotXYZ` is a `RotMatrix3`. The pose type is threaded
# through `BindingSite`, `Particle`, `Polyform` and `BindingRules`, so a species whose sites are
# pinned to `RotXYZ` would need every one of those wrappers to convert on the way past. Using
# the type that composes to itself -- which `PatchySphere` and `PolyhedronParticleSpecies`
# already build directly -- means nothing has to, and all 3D species agree on one pose type and
# can share a `BindingRules`.
function normal_pose(x, twist=0; orientationtype=length(x) == 2 ? Angle2d : RotMatrix3)
    n = normalize(x)
    φ = atan(n[2], n[1])
    length(x) == 2 && return Pose(x, orientationtype(φ))

    θ = acos(clamp(n[3], -1, 1))
    R = RotZ(φ) * RotY(θ - π/2) * RotX(twist)
    return Pose(x, orientationtype(R))
end