"""
    PolyhedronParticleSpecies{F}

A 3D convex polyhedron with one binding site per face.
"""
struct PolyhedronParticleSpecies{F,B<:BindingSite} <: ParticleSpecies{3,B}
    g::NautyDiGraph
    sites::Vector{B}
    polyhedron::Polyhedron{F}
    normals::Vector{SVector{3,F}}
    edgedirections::Vector{SVector{3,F}}
    rmin::F
    rmax::F
    skin::F
end

"""
    PolyhedronParticleSpecies(p::Polyhedron; colors=1:nfaces(p), locking=true, twists=0)

Build a particle species from the polyhedron `p`, with one binding site at each face centroid.

`colors` assigns interaction colors to the binding sites, which is what the interaction matrix
uses to decide which sites bond.

`locking` says whether a bond locks orientation, or is free to bond in all geometrically permitted
relative orientations.

`twists` turns a face's binding site about its own normal, by an angle in radians.
"""
function PolyhedronParticleSpecies(p::Polyhedron; colors=1:nfaces(p), locking=true, twists=0)
    _polyhedronspecies(p, colors, locking, twists, nothing)
end

# `usecycle` forces an encoding; see `_facesites`. Internal, and the reason the public
# constructor above is a one-liner.
function _polyhedronspecies(p::Polyhedron{F}, colors, locking, twists, usecycle::Union{Nothing,Bool}) where {F}
    n = nfaces(p)
    length(colors) == n || throw(ArgumentError("expected $n colors, one per face, got $(length(colors))"))

    rmin = inradius(p)
    rmax = bounding_radius(p)
    tol = sqrt(eps(F)) * rmax

    g, sites = _facesites(
        p,
        fs -> _faceposes(p, fs),
        colors,
        _expandperface(locking, n, "locking flags"),
        _expandperface(twists, n, "twists"),
        usecycle,
        tol,
        tol / rmin,
    )

    return check_encoding(
        PolyhedronParticleSpecies{F,eltype(sites)}(g, sites, p, facenormals(p), _edgedirections(p), rmin, rmax, tol)
    )
end

"""
    _faceposes(p::Polyhedron, fs)

The binding site frames of `p`'s faces, given face corner lists `fs`: local x along the
outward normal, local z pointing at the midpoint of each face's first edge.
"""
function _faceposes(p::Polyhedron{F}, fs::Vector{Vector{Int}}) where {F}
    P = Pose{3,F,RotMatrix3{F}}
    return map(eachindex(fs)) do i
        x = facecentroid(p, i)
        ex = facenormal(p, i)
        ez = normalize((corners(p)[fs[i][1]] + corners(p)[fs[i][2]]) / 2 - x)
        P(x, RotMatrix3{F}(hcat(ex, cross(ez, ex), ez)))
    end
end

# Separating axis candidates contributed by the edges. Only the direction matters, and only
# up to sign, so parallel edges are collapsed: a cube contributes only 3 instead of 12.
function _edgedirections(p::Polyhedron{F}) where {F}
    dirs = SVector{3,F}[]
    for f in faces(p), k in eachindex(f)
        d = normalize(corners(p)[f[mod1(k + 1, length(f))]] - corners(p)[f[k]])
        any(e -> isapprox(abs(dot(e, d)), 1; atol=sqrt(eps(F))), dirs) && continue
        push!(dirs, d)
    end
    return dirs
end

function Base.show(io::Core.IO, ps::PolyhedronParticleSpecies)
    return print(io, "$(dimension(ps))d PolyhedronParticleSpecies with $(nsites(ps)) sites")
end

function Base.copy(ps::PolyhedronParticleSpecies)
    return typeof(ps)(
        copy(ps.g),
        copy(ps.sites),
        copy(ps.polyhedron),
        copy(ps.normals),
        copy(ps.edgedirections),
        ps.rmin,
        ps.rmax,
        ps.skin,
    )
end

graphrep(p::PolyhedronParticleSpecies) = p.g
nsites(p::PolyhedronParticleSpecies) = length(p.sites)
bindingsite(p::PolyhedronParticleSpecies, i::Integer) = p.sites[i]
isconvex(::PolyhedronParticleSpecies) = true
bounding_radius(ps::PolyhedronParticleSpecies) = ps.rmax

"""
    polyhedron(ps::PolyhedronParticleSpecies)

Return the [`Polyhedron`](@ref) the species was built from.
"""
polyhedron(ps::PolyhedronParticleSpecies) = ps.polyhedron
corners(ps::PolyhedronParticleSpecies) = corners(ps.polyhedron)

function overlap(
    p1::SpeciesAndPose{<:PolyhedronParticleSpecies}, p2::SpeciesAndPose{<:PolyhedronParticleSpecies}; kwargs...
)
    spcs1, pose1 = p1
    spcs2, pose2 = p2
    skin = spcs1.skin + spcs2.skin
    t = pose2.x - pose1.x

    # check in and out radii first
    d = norm(t)
    d >= spcs1.rmax + spcs2.rmax && return false
    d < (spcs1.rmin + spcs2.rmin) - skin && return true

    # separating axes
    axes = Iterators.flatten((
        (pose1.psi * nrm for nrm in spcs1.normals),
        (pose2.psi * nrm for nrm in spcs2.normals),
        (cross(pose1.psi * e1, pose2.psi * e2) for e1 in spcs1.edgedirections, e2 in spcs2.edgedirections),
    ))
    return sat_overlap(axes, corners(spcs1), pose1, corners(spcs2), pose2, skin)
end

"""
    UnitTetrahedron

A regular tetrahedron with unit-length edges and one binding site per face.
"""
const UnitTetrahedron = PolyhedronParticleSpecies(Tetrahedron())

"""
    UnitCube

A cube with unit-length edges and one binding site per face.
"""
const UnitCube = PolyhedronParticleSpecies(Cube())

"""
    UnitOctahedron

A regular octahedron with unit-length edges and one binding site per face.
"""
const UnitOctahedron = PolyhedronParticleSpecies(Octahedron())

"""
    UnitDodecahedron

A regular dodecahedron with unit-length edges and one binding site per face.
"""
const UnitDodecahedron = PolyhedronParticleSpecies(Dodecahedron())

"""
    UnitIcosahedron

A regular icosahedron with unit-length edges and one binding site per face.
"""
const UnitIcosahedron = PolyhedronParticleSpecies(Icosahedron())

"""
    UnitPyramid(n, a=1.0; h=a, kwargs...)

A pyramid over a regular `n`-gon with edge length `a`, with one binding site per face.
Rotation group `C_n`.
"""
UnitPyramid(n::Integer, a::Real=1.0; h::Real=a, kwargs...) = PolyhedronParticleSpecies(Pyramid(n, a; h); kwargs...)

"""
    UnitPrism(n, a=1.0; h=a, kwargs...)

A prism over a regular `n`-gon with edge length `a`, with one binding site per face.
Rotation group `D_n`, except for `UnitPrism(4)` whose default height makes it a cube.
"""
UnitPrism(n::Integer, a::Real=1.0; h::Real=a, kwargs...) = PolyhedronParticleSpecies(Prism(n, a; h); kwargs...)

"""
    UnitAntiprism(n, a=1.0; kwargs...)

A uniform antiprism over a regular `n`-gon with edge length `a`, with one binding site per
face. Rotation group `D_n`, except for `UnitAntiprism(3)` which is a regular octahedron.
"""
UnitAntiprism(n::Integer, a::Real=1.0; kwargs...) = PolyhedronParticleSpecies(Antiprism(n, a); kwargs...)


"""
    SymmetricUnitTetrahedron

A regular tetrahedron with unit-length edges, every face the same color.

Where [`UnitTetrahedron`](@ref) gives each face a color of its own, this one leaves the body its full
rotation group, so `symmetrynumber` reports that group's order rather than 1 and faces bond
interchangeably.
"""
const SymmetricUnitTetrahedron = PolyhedronParticleSpecies(Tetrahedron(); colors=fill(1, nfaces(Tetrahedron())))

"""
    SymmetricUnitCube

A cube with unit-length edges, every face the same color.

Where [`UnitCube`](@ref) gives each face a color of its own, this one leaves the body its full
rotation group, so `symmetrynumber` reports that group's order rather than 1 and faces bond
interchangeably.
"""
const SymmetricUnitCube = PolyhedronParticleSpecies(Cube(); colors=fill(1, nfaces(Cube())))

"""
    SymmetricUnitOctahedron

A regular octahedron with unit-length edges, every face the same color.

Where [`UnitOctahedron`](@ref) gives each face a color of its own, this one leaves the body its full
rotation group, so `symmetrynumber` reports that group's order rather than 1 and faces bond
interchangeably.
"""
const SymmetricUnitOctahedron = PolyhedronParticleSpecies(Octahedron(); colors=fill(1, nfaces(Octahedron())))

"""
    SymmetricUnitDodecahedron

A regular dodecahedron with unit-length edges, every face the same color.

Where [`UnitDodecahedron`](@ref) gives each face a color of its own, this one leaves the body its full
rotation group, so `symmetrynumber` reports that group's order rather than 1 and faces bond
interchangeably.
"""
const SymmetricUnitDodecahedron = PolyhedronParticleSpecies(Dodecahedron(); colors=fill(1, nfaces(Dodecahedron())))

"""
    SymmetricUnitIcosahedron

A regular icosahedron with unit-length edges, every face the same color.

Where [`UnitIcosahedron`](@ref) gives each face a color of its own, this one leaves the body its full
rotation group, so `symmetrynumber` reports that group's order rather than 1 and faces bond
interchangeably.
"""
const SymmetricUnitIcosahedron = PolyhedronParticleSpecies(Icosahedron(); colors=fill(1, nfaces(Icosahedron())))
