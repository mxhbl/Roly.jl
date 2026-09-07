"""
    BindingSite{P<:Pose}

A `BindingSite` describes an anchor point at which particles are attached together.
The binding site color, together with an interaction matrix, determines whether two binding sites may bind.

`sitesym` is the order of the site's own rotational symmetry about its outward normal, disregarding the rest
of the particle.
`stab` size of the stabilizer of this site in the rotational symmetry group of the containing particle.
For example, consider a square side of a 3-prism (i.e. an extruded triangle). The side has 4-fold symmetry
(`sitesym == 4`), but only two of these rotations are symmetries of the entire particle (`stab == 2`).

A binding site's pose carries three things. Its position and its outward normal (the pose's local x axis)
are physical. The remaining freedom in the frame, the *twist*, is generally not: 
a site with a `sitesym`-fold symmetry is unchanged by turns of `2π/sitesym` about its normal,
so `psi`, `psi·Rx(2π/sitesym)`, ... all describe the same site.

`locking` determines if the binding site fixes the twist between it and its binding partners. A locking
site holds its partner fixed in the relative orientation defined by its frame (see [`standard_twist`](@ref)), 
up to symmetry transformations of the particle, so the number of the possible twists is equal to `stab`. 
If `locking==false`, the site allows all twists that are compatible with the site symmetry alone, disregarding
the symmetry of the particle, so the number of twists is equal to `sitesym`.
The two choices of locking coincide whenever the stabilizer of the site is equal to its individual symmetry group.
"""
struct BindingSite{P<:Pose,F<:Real}
    pose::P
    color::Int
    vertices::UnitRange{Int}
    touching_tolerance::F
    alignment_tolerance::F
    sitesym::Int
    stab::Int
    locking::Bool
end
function BindingSite(
    pose::P,
    color::Integer,
    vertices::UnitRange{<:Integer},
    touching_tolerance::Real,
    alignment_tolerance::Real,
    sitesym::Integer=1,
    stab::Integer=1,
    locking::Bool=true,
) where {P<:Pose}
    F = eltype(P)
    return BindingSite{P,F}(
        pose, color, vertices, convert(F, touching_tolerance), convert(F, alignment_tolerance), sitesym, stab, locking
    )
end

@inline BindingSite(
    site::BindingSite;
    pose=site.pose,
    color=site.color,
    vertices=site.vertices,
    sitesym=site.sitesym,
    stab=site.stab,
    locking=site.locking,
) = typeof(site)(pose, color, vertices, site.touching_tolerance, site.alignment_tolerance, sitesym, stab, locking)

Base.:(==)(b1::BindingSite, b2::BindingSite) = b1.vertices == b2.vertices && b1.color == b2.color
Base.hash(b::BindingSite, h::UInt) = hash(b.color, hash(b.vertices, h))
Base.isless(b1::BindingSite, b2::BindingSite) = isless(b1.vertices, b2.vertices)
Base.:*(p, site::BindingSite) = BindingSite(site; pose=p * site.pose)
Base.:*(site::BindingSite, p) = BindingSite(site; pose=site.pose * p)
translate(site::BindingSite, v) = BindingSite(site; pose=translate(site.pose, v))
@inline shift_vertices(site::BindingSite, v::Integer) = BindingSite(site; vertices=site.vertices .+ v)
@inline shift_color(site::BindingSite, c::Integer) = BindingSite(site; color=site.color + c)
@inline setcolor(site::BindingSite, c::Integer) = BindingSite(site; color=c)
@inline setstab(site::BindingSite, s::Integer) = BindingSite(site, stab=s)
posetype(::BindingSite{P}) where {P} = P
posetype(::Type{<:BindingSite{P}}) where {P} = P

@inline standard_rotation(::Type{F}, ::Val{2}, args...) where {F} = Angle2d{F}(π)
@inline standard_rotation(::Type{F}, ::Val{3}, t::Integer=0, ntwists::Integer=1) where {F} =
    RotXYZ{F}(2F(π) * t / ntwists, 0, π)

"""
    standard_twist(b::BindingSite, t=0, ntwists=1)
    standard_twist(p::Pose, t=0, ntwists=1)

Return the pose a partner of binding site `b` takes when the bond is in twist `t` of `ntwists`.

We employ the convention that attached binding sites are "facing each other": their poses are
related via a 180 degree rotation around their (shared) z-axes.

A bond in 3D additionally admits a *twist* about the bond (x-)axis.
This twist is generally fixed by the two sites' frames, except that the frame of a symmetric particle
can only be determined up to symmetry transformations.
If there are `ntwists` admissable twists, then any twist angle that is a multiple of `2π/ntwists`
is valid, and `t` indexes them from 0.
"""
@inline standard_twist(b::BindingSite{<:Pose{D,F}}, t::Integer=0, ntwists::Integer=1) where {D,F} =
    b.pose * standard_rotation(F, Val(D), t, ntwists)
@inline standard_twist(p::Pose{D,F}, t::Integer=0, ntwists::Integer=1) where {D,F} =
    p * standard_rotation(F, Val(D), t, ntwists)

"""
    twistfreedom(b::BindingSite)

Return how many possible twists binding site `b` allows. Equal to `stab` if the site is locking,
and equal to `sitesym` otherwise. See [`BindingSite`](@ref).
"""
@inline twistfreedom(b::BindingSite) = b.locking ? b.stab : b.sitesym

"""
    twistfreedom(b1::BindingSite, b2::BindingSite)

Return `ntwists`, the number of twists about the bond axis that a bond between `b1` and `b2`
admits. Equal to the least common multiple of the two sites' own twist freedoms, since the
admissible set has to be closed under both. See [`standard_twist`](@ref).
"""
@inline twistfreedom(b1::BindingSite, b2::BindingSite) = lcm(twistfreedom(b1), twistfreedom(b2))

"""
    _ndistincttwists(bbody::BindingSite, battach::BindingSite)

Return how many of the `twistfreedom(bbody, battach)` twists give *distinct* structures when
binding site `battach` is attached to a binding site `bbody` on a polyform.

This count excludes twists that differ by a symmetry of the particle being attached, of which
there are exactly `battach.stab`.
Note the asymmetry: only the *attached* binding site's stabilizer is quotiented out. The symmetries of the
polyform carrying `bbody` can merge structures too, but these are not knowable from one site,
and need to be caught through canonization.
"""
@inline _ndistincttwists(bbody::BindingSite, battach::BindingSite) = twistfreedom(bbody, battach) ÷ battach.stab

"""
    Contact(vs1, vs2, twist, ntwists)

Two binding sites found in contact, and how they meet.

`vs1` and `vs2` are the [`graphrep`](@ref) vertex ranges the two sites occupy. `twist` is which of the `ntwists`
twists the bond admits they are in. See also [`contact_pairing`](@ref),
"""
struct Contact
    vs1::UnitRange{Int}
    vs2::UnitRange{Int}
    twist::Int
    ntwists::Int
end

Base.show(io::Core.IO, c::Contact) = print(io, "Contact[", c.vs1, "-", c.vs2, ", twist ", c.twist, "/", c.ntwists, "]")

"""
    contact_pairing(vs1::UnitRange, vs2::UnitRange, t=0, ntwists=1)

Return the pairs of graph vertices that should be joined when two bonding binding sites occupying the vertex
ranges `vs1` and `vs2` meet in twist `t` of `ntwists` (see [`standard_twist`](@ref)).

Imagine the vertices corresponds to corners of a regular n-gon attached to the site, number the vertices of each site
starting from zero and denote `kᵢ = length(vsᵢ)`. Vertex `a` of site 1 then sits at angle `2πa/k₁` in site 1's frame,
and vertex `b` of site 2  (whose frame is turned by `2πt/ntwists` and then flipped, so its cyclic order runs the other
way around) sits at `2πt/ntwists - 2πb/k₂`. Writing `K = lcm(k₁, k₂)`, `s₁ = K/k₁`, `s₂ = K/k₂`, the coincidences are
the solutions of

    a*s₁ + b*s₂ = m   (mod K),      m = t*K/ntwists

`m` is always an integer: a site's `sitesym` `gᵢ` divides its vertex count, so
`ntwists = lcm(g₁, g₂)` divides `K`. Since `gcd(s₁, s₂) = 1`, the congruence has exactly
`gcd(k₁, k₂)` solutions for every `m`, and distinct twists give disjoint solution sets, so the
pairing count never depends on the twist, and the graph records which twist a bond is in.
"""
@inline contact_pairing(c::Contact) = contact_pairing(c.vs1, c.vs2, c.twist, c.ntwists)

@inline function contact_pairing(vs1::UnitRange{<:Integer}, vs2::UnitRange{<:Integer}, t::Integer=0, ntwists::Integer=1)
    k1, k2 = length(vs1), length(vs2)
    G = gcd(k1, k2)
    K = k1 * k2 ÷ G
    a0, b0 = _twist_shift(k1, k2, K, t, ntwists)
    return (vs1[1 + mod(a0 + j * (k1 ÷ G), k1)] => vs2[1 + mod(b0 - j * (k2 ÷ G), k2)] for j in 0:(G - 1))
end

# One solution of a*s1 + b*s2 = t*K/ntwists (mod K); the rest follow by stepping the two ranges
# in opposite directions. Zero-based, and `(0, 0)` for the untwisted bond.
@inline function _twist_shift(k1::Int, k2::Int, K::Int, t::Integer, ntwists::Integer)
    t == 0 && return 0, 0
    K % ntwists == 0 || throw(
        ArgumentError(
            "a bond in $ntwists twists cannot be encoded by sites of $k1 and $k2 graph vertices: " *
            "$ntwists must divide lcm($k1, $k2) = $K. Give the sites more vertices, i.e. a graph built " *
            "by `dartencoding` rather than `cycleencoding`, or declare a smaller sitesym.",
        ),
    )
    s1, s2 = K ÷ k1, K ÷ k2
    m = t * K ÷ ntwists
    # u*s1 = 1 (mod s2), so a0 = m*u solves the congruence mod s2 and hence mod K once b0
    # absorbs the remainder. gcd(s1, s2) = 1 always, so the inverse exists.
    a0 = s2 == 1 ? 0 : mod(m * invmod(s1, s2), k1)
    return a0, mod((m - a0 * s1) ÷ s2, k2)
end

"""
    color(b::BindingSite)

Return the binding site's color.

If the binding site comes from an `BindingRules`, its interactions with other
binding sites are determined by `interactionmatrix(rules)[:, color(b)]`.
"""
color(b::BindingSite) = b.color

"""
    istouching(b1::BindingSite, b2::BindingSite)

Check whether the translational components of the binding sites' poses
are approximately equal.
"""
function istouching(b1::BindingSite, b2::BindingSite)
    return isapprox(b1.pose.x, b2.pose.x; atol=b1.touching_tolerance + b2.touching_tolerance, rtol=0)
end

"""
    twist(b1::BindingSite, b2::BindingSite)

Return which of the bond's twists the two sites are in, or `nothing` if their
orientations are not related by any of them.

This function is symmetric in its arguments.
"""
function twist(b1::BindingSite, b2::BindingSite)
    atol = b1.alignment_tolerance + b2.alignment_tolerance
    ntwists = twistfreedom(b1, b2)
    for t in 0:(ntwists - 1)
        isapprox(b1.pose.psi, standard_twist(b2, t, ntwists).psi; atol, rtol=0) && return t
    end
    return nothing
end

"""
    isaligned(b1::BindingSite, b2::BindingSite)

Check whether the orientation components of the binding sites' poses
differ by a [`standard_twist`](@ref), in any twist the bond admits.
"""
isaligned(b1::BindingSite, b2::BindingSite) = !isnothing(twist(b1, b2))

"""
    isincontact(b1::BindingSite, b2::BindingSite)

Check whether the binding sites are touching and aligned.
"""
function isincontact(b1::BindingSite, b2::BindingSite)
    return istouching(b1, b2) && isaligned(b1, b2)
end

function Base.show(io::Core.IO, b::BindingSite)
    print(io, "BindingSite[c=$(b.color), vs=$(b.vertices)]")
end
function Base.show(io::Core.IO, ::MIME"text/plain", b::BindingSite)
    println(io, "$(dimension(b.pose))-dimensional BindingSite:")
    print(io, " - color: \t")
    println(io, b.color)
    print(io, " - vertices:\t$(b.vertices)\n")
    print(io, " - pose: \t$(b.pose)")
end