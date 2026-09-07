"""
    PolygonParticleSpecies{F}

A 2D convex regular polygon with `n` edges and one binding site per edge.
"""
struct PolygonParticleSpecies{F,B<:BindingSite} <: ParticleSpecies{2,B}
    g::NautyDiGraph
    sites::Vector{B}
    corners::Vector{SVector{2,F}}
    rmin::F
    rmax::F
    skin::F
end

"""
    PolygonParticleSpecies(n, a=1.0; colors=1:n)

Construct a regular `n`-gon with edge length `a`. Each edge carries one binding site.

`colors` assigns interaction colors to the binding sites, which is what the interaction matrix
uses to decide which sites bond.
"""
function PolygonParticleSpecies(n::Integer, a::F=1.0; colors=1:n) where {F<:Real}
    r_in = convert(F, 0.5a * cot(π / n))
    r_out = convert(F, 0.5a * csc(π / n))

    tol = sqrt(eps(F)) * r_out
    poses = map(1:n) do i
        ψ = Angle2d{F}(-F(π) * (1 / 2 + 2 / n * (i - 1)))
        Pose(SVector{2,F}(pol2cart(r_in, rotation_angle(ψ))), ψ)
    end
    # A 2D site has no turn about its in-plane normal, so its sitesym is 1 throughout, and with
    # it the stabilizer: the only rotation fixing a site is the identity.
    sitesyms = ones(Int, n)
    labels = siteorbits(poses, sitesyms, collect(colors))
    stabs = stabilizerorders(poses, sitesyms, labels)

    g, ranges = cycleencoding(n; labels)
    sites = BindingSite{Pose{2,F,Angle2d{F}},F}[]
    for (i, c) in enumerate(colors)
        push!(sites, BindingSite(poses[i], c, ranges[i], tol, tol / r_in, 1, stabs[i]))
    end

    corners = [
        SVector{2,F}(r_out * cos(-π / 2 - (2k - 1) * π / n), r_out * sin(-π / 2 - (2k - 1) * π / n)) for k in 1:n
    ]
    return check_encoding(PolygonParticleSpecies{F,eltype(sites)}(g, sites, corners, r_in, r_out, tol))
end

function Base.show(io::Core.IO, ps::PolygonParticleSpecies)
    return print(io, "$(dimension(ps))d PolygonParticleSpecies with $(nsites(ps)) sites")
end

function Base.copy(pps::PolygonParticleSpecies)
    return PolygonParticleSpecies(copy(pps.g), copy(pps.sites), copy(pps.corners), pps.rmin, pps.rmax, pps.skin)
end

graphrep(p::PolygonParticleSpecies) = p.g
nsites(p::PolygonParticleSpecies) = length(p.sites)
bindingsite(p::PolygonParticleSpecies, i::Integer) = p.sites[i]
isconvex(::PolygonParticleSpecies) = true

"""
    _isregular(corners)

Whether `corners` is a regular polygon, i.e. whether its own `2π/n` turn carries it onto itself
corner for corner. Accepts either winding.
"""
function _isregular(corners::AbstractVector{SVector{2,F}}) where {F}
    n = length(corners)
    atol = sqrt(eps(F)) * maximum(norm, corners)
    return any((1, -1)) do sgn
        R = Angle2d{F}(sgn * 2F(π) / n)
        all(i -> isapprox(R * corners[i], corners[mod1(i + 1, n)]; atol), 1:n)
    end
end

# The three regular polygons that tile the plane. `rmin` is the inradius r, and a regular n-gon
# has edge length 2r*tan(pi/n).
function _tilingcell(ps::PolygonParticleSpecies)
    nsites(ps) in (3, 4, 6) && _isregular(ps.corners) ? (nsites(ps), 2 * ps.rmin * tan(π / nsites(ps))) : nothing
end

function overlap(p1::SpeciesAndPose{<:PolygonParticleSpecies}, p2::SpeciesAndPose{<:PolygonParticleSpecies}; kwargs...)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    skin = spcs1.skin + spcs2.skin
    d = norm(pose1.x - pose2.x)

    # check in and out radii first
    d >= spcs1.rmax + spcs2.rmax && return false
    d < (spcs1.rmin + spcs2.rmin) - skin && return true

    return sat_overlap(
        Iterators.flatten((edgenormals(spcs1.corners, pose1), edgenormals(spcs2.corners, pose2))),
        spcs1.corners,
        pose1,
        spcs2.corners,
        pose2,
        skin,
    )
end

bounding_radius(ps::PolygonParticleSpecies) = ps.rmax

"""
    UnitNgon(n, a=1.0; kwargs...)

A regular `n`-gon with edge length `a`, with one binding site per edge. Rotation group `C_n`.

Spelled to match [`UnitPrism`](@ref) and the other body constructors; it is
[`PolygonParticleSpecies`](@ref) under a name that reads the same way.
"""
UnitNgon(n::Integer, a::Real=1.0; kwargs...) = PolygonParticleSpecies(n, a; kwargs...)

"""
    UnitTriangle

An equilateral triangle with unit-length edges and one binding site per edge.
"""
const UnitTriangle = PolygonParticleSpecies(3)

"""
    UnitSquare

A square with unit-length edges and one binding site per edge.
"""
const UnitSquare = PolygonParticleSpecies(4)

"""
    UnitHexagon

A regular hexagon with unit-length edges and one binding site per edge.
"""
const UnitHexagon = PolygonParticleSpecies(6)


"""
    SymmetricUnitTriangle

An equilateral triangle with unit-length edges, every edge the same color.

Where [`UnitTriangle`](@ref) gives each edge a color of its own, this one leaves the body its full
rotation group, so `symmetrynumber` is 3 rather than 1 and edges bond interchangeably.
"""
const SymmetricUnitTriangle = PolygonParticleSpecies(3; colors=fill(1, 3))

"""
    SymmetricUnitSquare

A square with unit-length edges, every edge the same color.

Where [`UnitSquare`](@ref) gives each edge a color of its own, this one leaves the body its full
rotation group, so `symmetrynumber` is 4 rather than 1 and edges bond interchangeably.
"""
const SymmetricUnitSquare = PolygonParticleSpecies(4; colors=fill(1, 4))

"""
    SymmetricUnitHexagon

A regular hexagon with unit-length edges, every edge the same color.

Where [`UnitHexagon`](@ref) gives each edge a color of its own, this one leaves the body its full
rotation group, so `symmetrynumber` is 6 rather than 1 and edges bond interchangeably.
"""
const SymmetricUnitHexagon = PolygonParticleSpecies(6; colors=fill(1, 6))
