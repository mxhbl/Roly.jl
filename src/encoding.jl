"""
    RotationGroup

A proper (or *chiral*) point group: the rotations carrying a rigid body onto itself, without reflections.

The five families are [`Cyclic`](@ref), [`Dihedral`](@ref), [`Tetrahedral`](@ref),
[`Octahedral`](@ref) and [`Icosahedral`](@ref). Use [`rotationgroup`](@ref) to get the rotational
symmetries of a body.
"""
abstract type RotationGroup end

"""
    Cyclic(n)

The rotation group `C_n`: `n` turns about a single axis. Order `n`, realized by [`Pyramid`](@ref).
"""
struct Cyclic{N} <: RotationGroup end
Cyclic(n::Integer) = Cyclic{Int(n)}()

"""
    Dihedral(n)

The rotation group `D_n`: the `n` turns about a principal axis together with the `n` half turns
about axes perpendicular to it. Order `2n`, realized by [`Prism`](@ref) and [`Antiprism`](@ref).
"""
struct Dihedral{N} <: RotationGroup end
Dihedral(n::Integer) = Dihedral{Int(n)}()

"""
    Tetrahedral()

The rotation group `T` of a regular tetrahedron, of order 12.
"""
struct Tetrahedral <: RotationGroup end

"""
    Octahedral()

The rotation group `O` shared by the cube and the regular octahedron, of order 24.
"""
struct Octahedral <: RotationGroup end

"""
    Icosahedral()

The rotation group `I` shared by the regular dodecahedron and icosahedron, of order 60.
"""
struct Icosahedral <: RotationGroup end

Base.show(io::Core.IO, ::Cyclic{N}) where {N} = print(io, "Cyclic($N)")
Base.show(io::Core.IO, ::Dihedral{N}) where {N} = print(io, "Dihedral($N)")

"""
    grouporder(group::RotationGroup)

Return the number of elements in `group`.
"""
grouporder(::Cyclic{N}) where {N} = N
grouporder(::Dihedral{N}) where {N} = 2N
grouporder(::Tetrahedral) = 12
grouporder(::Octahedral) = 24
grouporder(::Icosahedral) = 60

"""
    Polyhedron{F}

A convex polyhedron, stored as a list of `corners` (vertices) and a list of `faces`.

Each face is a list of indices into `corners`, wound counter-clockwise as seen from
*outside* the body. This winding is what fixes the orientation of the graph encoding.
On construction, each face is rotated to start at a canonical corner.
"""
struct Polyhedron{F<:AbstractFloat}
    corners::Vector{SVector{3,F}}
    faces::Vector{Vector{Int}}

    function Polyhedron(corners::Vector{SVector{3,F}}, faces::Vector{Vector{Int}}) where {F}
        _check_winding(length(corners), faces)
        corners = _recenter(corners)
        _check_convex(corners, faces)
        return new{F}(corners, _canonical_faces(corners, faces))
    end
end

Base.copy(p::Polyhedron) = Polyhedron(copy(p.corners), [copy(f) for f in p.faces])

"""
    Polyhedron(corners, faces)
    Polyhedron(corners)

Construct a convex polyhedron from its `corners` and `faces`. 

Faces must be wound counter-clockwise as seen from outside the body.
If `faces` is omitted, they are derived from the corners by finding all supporting planes of
the convex hull.
"""
function Polyhedron(corners::AbstractVector{<:AbstractVector}, faces)
    F = float(eltype(first(corners)))
    return Polyhedron(SVector{3,F}[SVector{3,F}(c) for c in corners], Vector{Int}[collect(Int, f) for f in faces])
end
function Polyhedron(corners::AbstractVector{<:AbstractVector})
    F = float(eltype(first(corners)))
    cs = _recenter(SVector{3,F}[SVector{3,F}(c) for c in corners])
    return Polyhedron(cs, _derive_faces(cs))
end

"""
    _check_winding(ncorners, faces)

Verify that `faces` describes a closed, consistently oriented surface: every directed edge
occurs exactly once, its reverse occurs in exactly one other face, and all corners are part
of some face.
"""
function _check_winding(ncorners::Integer, faces::AbstractVector{<:AbstractVector{Int}})
    length(faces) < 3 && throw(ArgumentError("a polyhedron needs at least 3 faces"))
    used = falses(ncorners)
    for f in faces, v in f
        checkindex(Bool, eachindex(used), v) && (used[v] = true)
    end
    all(used) || throw(
        ArgumentError(
            "corner$(count(!, used) > 1 ? "s" : "") $(join(findall(!, used), ", ")) " *
            "$(count(!, used) > 1 ? "are" : "is") used by no face",
        ),
    )

    seen = Dict{Tuple{Int,Int},Int}()
    for (i, f) in enumerate(faces)
        length(f) < 3 && throw(ArgumentError("face $i has fewer than 3 corners"))
        allunique(f) || throw(ArgumentError("face $i repeats a corner"))
        all(in(1:ncorners), f) || throw(ArgumentError("face $i indexes a nonexistent corner"))
        for k in eachindex(f)
            e = (f[k], f[mod1(k + 1, length(f))]) # edge
            haskey(seen, e) && throw(
                ArgumentError(
                    "directed edge $e occurs in faces $(seen[e]) and $i; faces must all be wound counter-clockwise seen from outside",
                ),
            )
            seen[e] = i
        end
    end
    for (e, i) in seen
        haskey(seen, reverse(e)) ||
            throw(ArgumentError("edge $e of face $i is not shared with a second face; the surface is not closed"))
    end
    return nothing
end

"""
    _check_convex(corners, faces)

Verify that every corner lies on the inner side of every face's plane, so that the body is
the intersection of its faces' half-spaces.
"""
function _check_convex(corners::Vector{SVector{3,F}}, faces::AbstractVector{<:AbstractVector{Int}}) where {F}
    atol = sqrt(eps(F)) * maximum(norm, corners)
    for (i, f) in enumerate(faces)
        centroid = sum(corners[v] for v in f) / length(f)
        nrm = zero(SVector{3,F})
        for k in eachindex(f)
            nrm += cross(corners[f[k]] - centroid, corners[f[mod1(k + 1, length(f))]] - centroid)
        end
        nrm = normalize(nrm)
        d = dot(nrm, centroid)
        for (j, x) in enumerate(corners)
            dot(nrm, x) <= d + atol ||
                throw(ArgumentError("corner $j lies outside the plane of face $i, so the body is not convex."))
        end
    end
    return nothing
end

"""
    _canonical_faces(corners, faces)

Return `faces` with each corner list cyclically permuted so that the face begins at an
intrinsically chosen corner, fixing the twist reference of that face's binding site.

Convention: a site's local z axis points to the midpoint of its face's first edge.
Faces that are related by a rotation of the body must pick compatible orderings.
For example, a rectangle has degree 4 but only 2-fold symmetry, so corners fall into two classes a 90deg turn apart.
Picking as a reference edge both long edges on some faces and a short edges on others hides the symmetry that relates them.

The choice of twist reference is made from the geometry of the faces: rotate to the lexicographically
least cyclic word of `(edge length, interior angle)`, which ensures compatible references throughout.
We also impose that opposing faces should have aligned orientations.
"""
function _canonical_faces(corners::Vector{SVector{3,F}}, faces::Vector{Vector{Int}}) where {F}
    atol = sqrt(eps(F)) * maximum(norm, corners)
    return map(faces) do f
        k = length(f)
        edge(m) = norm(corners[f[mod1(m + 1, k)]] - corners[f[mod1(m, k)]])
        function interior_angle(m)
            a = corners[f[mod1(m - 1, k)]] - corners[f[mod1(m, k)]]
            b = corners[f[mod1(m + 1, k)]] - corners[f[mod1(m, k)]]
            return acos(clamp(dot(a, b) / (norm(a) * norm(b)), -one(F), one(F)))
        end
        # Does the word starting at s come before the one starting at t?
        function isbefore(s, t)
            for j in 0:(k - 1), (x, y) in ((edge(s + j), edge(t + j)), (interior_angle(s + j), interior_angle(t + j)))
                isapprox(x, y; atol) || return x < y
            end
            return false
        end
        best = 1
        for s in 2:k
            isbefore(s, best) && (best = s)
        end
        circshift(f, 1 - best)
    end |> fs -> _mateopposing(corners, fs)
end

"""
    _mateopposing(corners, faces)

Return `faces` with the corner list of each face that has an opposite turned so that the two meet
squarely rather than askew.

The second half of the convention [`_canonical_faces`](@ref) sets. The first half picks each
face's starting corner from that face alone, which leaves two faces that face each other starting
wherever their own geometry says, and a bond between them turning the neighbour by whatever their
two references happen to differ by.

This is a choice among bonds that are all already sound, not a repair. A bond admits only the
twists the face's own symmetry allows, and those carry the polygon onto itself, so two bonded
faces line up whichever one it takes. What the starting corner settles is *which* of them the
bond is, and that decides whether the neighbour ends up a translate.

Call the face's own symmetry angle, `2π/sitesym`, one step -- and note that turning the corner
list is the only way to move, so only whole numbers of steps are available. That is also what
makes this free: a whole step carries a frame onto an equivalent one, so it marks nothing and
costs the body no symmetry.

What decides the case is how the two faces sit in the body, not how symmetric they are:

  - they are turned alike, so the turn that makes bonded copies *translates* of each other is a
    whole number of steps: take it, and rules over those faces can tile by translation. Every
    prism, whatever its parity, and the cube
  - they are staggered by half a step, so that turn is a half-integer number of them. Half a turn
    is then itself half-integer when the symmetry is odd (`sitesym/2`), so adding it lands on a
    whole number: take that, and copies meet turned by 180 degrees. Octahedra, dodecahedra,
    icosahedra, odd antiprisms
  - staggered, with even symmetry, so neither is whole: leave the face alone. Nothing is given up
    by doing so -- the stagger is in the body itself, so no rotation about the bond makes two of
    them translates, and the bond they already have is the one their geometry asks for. Only the
    caps of an even antiprism

A face with no symmetry at all has no whole step but zero, so it is never turned, and opposing
faces are under no obligation to be congruent.
"""
function _mateopposing(cs::Vector{SVector{3,F}}, faces::Vector{Vector{Int}}) where {F}
    tol = sqrt(eps(F)) * maximum(norm, cs)
    centroid(f) = sum(cs[v] for v in f) / length(f)
    function frame(f)
        x = centroid(f)
        ex = normalize(cross(cs[f[2]] - cs[f[1]], cs[f[3]] - cs[f[2]]))
        ez = normalize((cs[f[1]] + cs[f[2]]) / 2 - x)
        return RotMatrix3{F}(hcat(ex, cross(ez, ex), ez))
    end
    function symmetry(f)
        k = length(f)
        c = centroid(f)
        rel = [cs[v] - c for v in f]
        nrm = frame(f)[:, 1]
        return count(0:(k - 1)) do s
            R = AngleAxis(2F(π) * s / k, nrm[1], nrm[2], nrm[3])
            all(j -> isapprox(R * rel[j], rel[mod1(j + s, k)]; atol=tol), 1:k)
        end
    end

    out = [copy(f) for f in faces]
    for i in eachindex(faces)
        j = findfirst(m -> isapprox(dot(frame(faces[i])[:, 1], frame(faces[m])[:, 1]), -1; atol=tol),
                      eachindex(faces))
        (isnothing(j) || j <= i) && continue
        f = faces[j]
        # the frame a partner of face `i` must wear, against the one face `j` wears already
        turn = transpose(frame(f)) * (frame(faces[i]) * standard_rotation(F, Val(3)))
        mate = atan(turn[3, 2], turn[2, 2])

        # a whole step of the face's own symmetry, counted in corners
        perstep = length(f) ÷ symmetry(f)
        function wholesteps(θ)
            steps = θ * length(f) / (2F(π))
            m = round(Int, steps)
            return abs(steps - m) < sqrt(eps(F)) && iszero(mod(m, perstep)) ? m : nothing
        end
        m = something(wholesteps(mate), wholesteps(mate + F(π)), 0)
        out[j] = circshift(f, -m)
    end
    return out
end

"""
    _propagate_faces(corners, faces, labels)

Return `faces` with each corner list cyclically permuted so that faces related by a
label-preserving rotation of the body have the compatible ordering of corners.

For example, consider a triangular prism with square sides. All square sides need to be
oriented in such a way that rotations around the major axis map them to each other.
Otherwise, spurrous 90deg rotations might appear. The issue comes from a mismatch of the
faces' sitesym (== 4) and stabilizer (== 2). `_propagate_faces` fixes this by propagating the
orientation of a reference face across its symmetry orbit.

It only makes an orbit internally consistent; which corner the orbit's own reference face starts
at is left as it arrived, from [`_canonical_faces`](@ref). Both passes are needed, and in that
order: the intrinsic one has no group to consult, and this one has no way to pick a reference.
"""
function _propagate_faces(cs::Vector{SVector{3,F}}, faces::Vector{Vector{Int}}, labels) where {F}
    atol = sqrt(eps(F)) * maximum(norm, cs)
    centroid(f) = sum(cs[v] for v in f) / length(f)
    midpoint(f, k) = (cs[f[k]] + cs[f[mod1(k + 1, length(f))]]) / 2

    # Build the rotation subgroup that carries every face onto one of the same label, i.e. the 
    # symmetry group of the labeled polyhedron
    faceof(x) = findfirst(i -> isapprox(centroid(faces[i]), x; atol), eachindex(faces))
    group = filter(_rotationgroup(cs, faces)) do Q
        all(eachindex(faces)) do i
            j = faceof(Q * centroid(faces[i]))
            !isnothing(j) && labels[j] == labels[i]
        end
    end

    # non-symmetric
    length(group) == 1 && return faces

    faces = [copy(f) for f in faces]
    for j in 2:length(faces)
        cj = centroid(faces[j])
        for i in 1:(j - 1)
            labels[i] == labels[j] || continue
            # A symmetry carrying face i's centroid onto face j's carries the face itself, and
            # so carries its first edge midpoint onto one of face j's.

            # find a group element that maps face face i's midpoint to face j
            k = findfirst(Q -> isapprox(Q * centroid(faces[i]), cj; atol), group)
            isnothing(k) && continue

            # find the edge of face j that the first edge of face i is mapped to
            target = group[k] * midpoint(faces[i], 1)
            d = findfirst(m -> isapprox(midpoint(faces[j], m), target; atol), eachindex(faces[j]))
            isnothing(d) && continue
            # (there can only be a single group element that maps an edge of i to an edge of j)

            # permute face j to make it congruent to face i under the action of the group
            faces[j] = circshift(faces[j], 1 - d)
            break
        end
    end
    return faces
end

"""
    _derive_faces(corners)

Find the faces of the convex hull of `corners` by testing every corner triple for a
supporting plane, then winding each face counter-clockwise about its outward normal.

`corners` must be centered around the origin, so that the sign of a plane's offset orients its
normal outward.

Every corner on a face's plane joins that face, including one sitting mid-edge, which allows a
face containing three collinear corners.
"""
function _derive_faces(corners::Vector{SVector{3,F}}) where {F}
    n = length(corners)
    atol = sqrt(eps(F)) * maximum(norm, corners)

    faces = Vector{Int}[]
    planes = Tuple{SVector{3,F},F}[]
    for i in 1:(n - 2), j in (i + 1):(n - 1), k in (j + 1):n
        nrm = cross(corners[j] - corners[i], corners[k] - corners[i])
        # collinear triple: no plane through it, and any plane it lies on is found by another.
        norm(nrm) < atol && continue
        nrm = normalize(nrm)
        d = dot(nrm, corners[i])
        # point the normal away from the centroid, which is the origin
        d < 0 && ((nrm, d) = (-nrm, -d))
        # supporting plane of the hull?
        all(c -> dot(nrm, c) <= d + atol, corners) || continue
        # already found this plane?
        any(pl -> isapprox(pl[1], nrm; atol) && isapprox(pl[2], d; atol), planes) && continue

        push!(planes, (nrm, d))
        push!(faces, _wind_counterclockwise!(corners, findall(c -> abs(dot(nrm, c) - d) <= atol, corners), nrm))
    end
    return faces
end

"""
    _wind_counterclockwise!(corners, idxs, nrm)

In-place sort `idxs`, the indices of `corners` lying on a supporting plane, into counter-clockwise order seen
from outside, i.e. looking down `nrm`.
"""
function _wind_counterclockwise!(corners::Vector{SVector{3,F}}, idxs::Vector{Int}, nrm::SVector{3,F}) where {F}
    c = sum(corners[i] for i in idxs) / length(idxs)
    u = normalize(corners[first(idxs)] - c)
    v = cross(nrm, u)
    return sort!(idxs; by=i -> atan(dot(corners[i] - c, v), dot(corners[i] - c, u)))
end

"""
    corners(p::Polyhedron)

Return the corner positions of `p`.
"""
corners(p::Polyhedron) = p.corners

"""
    faces(p::Polyhedron)

Return the faces of `p` as lists of corner indices, each wound counter-clockwise seen from
outside.
"""
faces(p::Polyhedron) = p.faces

"""
    facevertices(p::Polyhedron, i)

Return the corner indices of the `i`th face of `p`.
"""
facevertices(p::Polyhedron, i::Integer) = p.faces[i]

ncorners(p::Polyhedron) = length(p.corners)

"""
    nfaces(p::Polyhedron)

Return the number of faces of `p`.
"""
nfaces(p::Polyhedron) = length(p.faces)

"""
    nedges(p::Polyhedron)

Return the number of edges of `p`.
"""
nedges(p::Polyhedron) = sum(length, p.faces) ÷ 2

facedegree(p::Polyhedron, i::Integer) = length(p.faces[i])

Base.eltype(::Type{<:Polyhedron{F}}) where {F} = F
Base.eltype(p::Polyhedron) = eltype(typeof(p))

"""
    facecentroid(p::Polyhedron, i)

Return the centroid of the `i`th face of `p`.
"""
facecentroid(p::Polyhedron, i::Integer) = sum(p.corners[v] for v in p.faces[i]) / facedegree(p, i)

"""
    facecentroids(p::Polyhedron)

Return the centroids of every face of `p`.
"""
facecentroids(p::Polyhedron) = [facecentroid(p, i) for i in 1:nfaces(p)]

"""
    facenormal(p::Polyhedron, i)

Return the outward unit normal of the `i`th face of `p`.
"""
function facenormal(p::Polyhedron{F}, i::Integer) where {F}
    # Summed over the whole face to handle degenerate edges
    f = p.faces[i]
    c = facecentroid(p, i)
    nrm = zero(SVector{3,F})
    for k in eachindex(f)
        nrm += cross(p.corners[f[k]] - c, p.corners[f[mod1(k + 1, length(f))]] - c)
    end
    return normalize(nrm)
end

"""
    facenormals(p::Polyhedron)

Return the outward unit normals of every face of `p`.
"""
facenormals(p::Polyhedron) = [facenormal(p, i) for i in 1:nfaces(p)]

"""
    edgemidpoint(p::Polyhedron, i, k)

Return the midpoint of edge `k` of face `i` of `p`.
"""
function edgemidpoint(p::Polyhedron, i::Integer, k::Integer)
    f = p.faces[i]
    return (p.corners[f[k]] + p.corners[f[mod1(k + 1, length(f))]]) / 2
end

"""
    bounding_radius(p::Polyhedron)

Return the distance from the origin to the farthest corner of `p`.
"""
bounding_radius(p::Polyhedron) = maximum(norm, p.corners)

"""
    inradius(p::Polyhedron)

Return the distance from the origin to the closest face centroid of `p`.
"""
inradius(p::Polyhedron) = minimum(norm(facecentroid(p, i)) for i in 1:nfaces(p))

"""
    minedgelength(p::Polyhedron)

Return the length of the shortest edge of `p`.
"""
function minedgelength(p::Polyhedron)
    return minimum(norm(p.corners[f[k]] - p.corners[f[mod1(k + 1, length(f))]]) for f in p.faces for k in eachindex(f))
end

"""
    rotationgroup(p::Polyhedron)

Return the rotations that map `p` onto itself, as `RotMatrix3`es about the corner
centroid.
"""
rotationgroup(p::Polyhedron) = _rotationgroup(corners(p), faces(p))

function _rotationgroup(cs::Vector{SVector{3,F}}, faces::Vector{Vector{Int}}) where {F}
    atol = sqrt(eps(F)) * maximum(norm, cs)

    # An orthonormal frame attached to dart k of face i, i.e. to that face's kth edge: along
    # the edge, along the outward normal, and their cross product.
    function dartframe(i, k)
        f = faces[i]
        e1 = normalize(cs[f[mod1(k + 1, length(f))]] - cs[f[k]])
        e2 = _facenormal(cs, f)
        return hcat(e1, e2, cross(e1, e2))
    end

    Mref = dartframe(1, 1)
    group = RotMatrix3{F}[]

    # Every symmetry maps the first dart to some dart, and that correspondence determines the
    # rotation, so it suffices to test the `2 * nedges(p)` candidates this generates.
    for i in eachindex(faces), k in eachindex(faces[i])
        # every potential symmetry rotation maps the ref frame to some other dartframe, and is therefore
        # given by the matrix
        R = dartframe(i, k) * Mref'
        # R is a symmetry if every corners is mapped to some other corner
        all(c -> any(c2 -> isapprox(R * c, c2; atol), cs), cs) || continue
        push!(group, RotMatrix3{F}(R))
    end
    return group
end

function _facenormal(cs::Vector{SVector{3,F}}, f::Vector{Int}) where {F}
    c = sum(cs[v] for v in f) / length(f)
    nrm = zero(SVector{3,F})
    for k in eachindex(f)
        nrm += cross(cs[f[k]] - c, cs[f[mod1(k + 1, length(f))]] - c)
    end
    return normalize(nrm)
end

"""
    faceorbits(p::Polyhedron)

Group the faces into the orbits of [`rotationgroup`](@ref), and return one orbit index per
face.

Passing these to [`dartencoding`](@ref) yields the body's true rotation group, whereas
`labels=fill(1, nfaces(p))` yields the *combinatorial* symmetry of the face lattice, which
can be larger. For example, `Pyramid(3)` is combinatorially a tetrahedron and would report 12 instead of
its actual 3.
"""
function faceorbits(p::Polyhedron)
    cents = facecentroids(p)
    atol = sqrt(eps(eltype(p))) * maximum(norm, cents)

    labels = zeros(Int, nfaces(p))
    group = rotationgroup(p)
    next = 0
    for i in eachindex(labels)
        labels[i] == 0 || continue
        next += 1
        for R in group
            j = findfirst(c -> isapprox(R * cents[i], c; atol), cents)
            isnothing(j) || (labels[j] = next)
        end
    end
    return labels
end

"""
    facesym(p::Polyhedron, i)
    facesym(p::Polyhedron)

Return the order of face `i`'s own rotational symmetry about its outward normal.

This is the `sitesym` of a binding site placed on that face: a face invariant under turns of
`2π/sitesym` has that many equally good twist references, so nothing about the particle in
isolation may depend on which was picked. It divides the face degree and equals it only for a
regular face; a rectangular face has degree 4 but only 2-fold symmetry.

It is a property of the face alone, not of the body it belongs to. A triangular prism is only
2-fold about a square side face, but the face is still a square, so `facesym` is 4 there.
"""
function facesym(p::Polyhedron{F}, i::Integer) where {F}
    f = facevertices(p, i)
    k = length(f)
    c = facecentroid(p, i)
    nrm = facenormal(p, i)
    rel = [corners(p)[v] - c for v in f]
    atol = sqrt(eps(F)) * maximum(norm, rel)
    # A rotational symmetry of a k-gon shifts its corner ring cyclically, and a shift by s is
    # realized by the turn 2π/k * s about the normal
    return count(0:(k - 1)) do s
        R = AngleAxis(2F(π) * s / k, nrm[1], nrm[2], nrm[3])
        all(j -> isapprox(R * rel[j], rel[mod1(j + s, k)]; atol), 1:k)
    end
end
facesym(p::Polyhedron) = [facesym(p, i) for i in 1:nfaces(p)]

function Base.show(io::Core.IO, p::Polyhedron)
    return print(io, "Polyhedron[v=$(ncorners(p)), e=$(nedges(p)), f=$(nfaces(p))]")
end
function Base.show(io::Core.IO, ::MIME"text/plain", p::Polyhedron{F}) where {F}
    println(io, "Polyhedron{$F}:")
    println(io, " - corners: \t$(ncorners(p))")
    println(io, " - edges: \t$(nedges(p))")
    print(io, " - faces: \t$(nfaces(p)) $(_facesummary(p))")
end
function _facesummary(p::Polyhedron)
    degs = sort!(unique(facedegree(p, i) for i in 1:nfaces(p)))
    return "(" * join(("$(count(i -> facedegree(p, i) == d, 1:nfaces(p)))×$d-gon" for d in degs), ", ") * ")"
end

"""
    dartencoding(p::Polyhedron; labels=1:nfaces(p))
    dartencoding(faces::Vector{Vector{Int}}; labels=1:length(faces))

Return `(g, ranges)`, the dart encoding of polyhedron `p` and the graph vertices belonging to each face.

A *dart* is a half-edge: one of the two directed traversals of a polyhedron edge, the one
belonging to a given face. Each edge therefore has two darts, one per adjoining face, a face of
degree `k` owns `k` of them, and the dart encoding graph has `2 * nedges(p)` vertices in total. 
The encoding consists of

1. a directed `k`-cycle through each face's own darts, following the face's counter-clockwise
   winding, and
2. a bidirectional edge joining the two darts that sit on each shared polyhedron edge.

The automorphism group of the resulting `g` is the polyhedron's rotational symmetry group.
`labels` assigns one symmetry label per face, inherited by that face's darts: all labels
distinct gives a symmetry number of 1, all labels equal gives the full rotation group, and
merging some faces gives the subgroup preserving that labeling. The dart encoding is closely
related to the Cayley graph of the body's symmetry group.
"""
dartencoding(p::Polyhedron; labels=1:nfaces(p)) = dartencoding(faces(p); labels)

function dartencoding(fs::Vector{Vector{Int}}; labels=1:length(fs))
    length(labels) == length(fs) ||
        throw(ArgumentError("expected $(length(fs)) labels, one per face, got $(length(labels))"))

    # graph vertices for each face of the polyhedron
    ranges = Vector{UnitRange{Int}}(undef, length(fs))
    ndarts = 0
    for (i, f) in enumerate(fs)
        ranges[i] = (ndarts + 1):(ndarts + length(f))
        ndarts += length(f)
    end

    vertex_labels = Cint[labels[i] for (i, f) in enumerate(fs) for _ in f]
    g = NautyDiGraph(ndarts; vertex_labels)

    # Dart k of face i sits on the directed edge f[k] -> f[k+1]. Its partner is the dart of
    # the adjacent face carrying the reversed edge
    dart_of_edge = Dict{Tuple{Int,Int},Int}()
    for (i, f) in enumerate(fs), k in eachindex(f)
        dart_of_edge[(f[k], f[mod1(k + 1, length(f))])] = first(ranges[i]) + k - 1
    end

    for (i, f) in enumerate(fs)
        v0 = first(ranges[i])
        for k in eachindex(f)
            v = v0 + k - 1
            add_edge!(g, v, v0 + mod(k, length(f)))
            partner = dart_of_edge[(f[mod1(k + 1, length(f))], f[k])]
            if v < partner
                add_edge!(g, v, partner)
                add_edge!(g, partner, v)
            end
        end
    end
    return g, ranges
end

"""
    cycleencoding(nsites; labels=1:nsites)

Return `(g, ranges)`, the cycle encoding of a particle with `nsites` binding sites: a single
directed cycle carrying one vertex per site.

This is the 2D polygon encoding. It is also valid in 3D when every site's twist freedom is 1
and all labels are distinct, and not otherwise.
"""
function cycleencoding(nsites::Integer; labels=1:nsites)
    length(labels) == nsites || throw(ArgumentError("expected $nsites labels, one per site, got $(length(labels))"))
    nsites < 1 && throw(ArgumentError("a particle needs at least one binding site"))

    g = NautyDiGraph(cycle_digraph(nsites); vertex_labels=collect(Cint, labels))
    return g, [i:i for i in 1:nsites]
end

"""
    _expandperface(x, n, what)

Expand a per-face keyword to a length-`n` vector, a scalar meaning the same for every face.
`what` names the argument in the length-mismatch error.
"""
function _expandperface(x, n::Integer, what::AbstractString)
    x isa Union{Bool,Real} && return fill(x, n)
    length(x) == n || throw(ArgumentError("expected $n $what, one per face, got $(length(x))"))
    return collect(x)
end

"""
    _facesites(p, poseof, colors, locking, twists, usecycle, touching_tol, alignment_tol)

Return `(g, sites)`, the graph encoding and the binding sites of a species carrying one site
per face of the polyhedron `p`.

  - `poseof(fs)` gives the site poses implied by the face corner lists `fs`. It is a function
    rather than a vector of poses because `_facesites` re-winds the faces itself, and a site's
    frame is read off its face's *first* corner, so the poses have to be taken again afterwards.
  - `colors`, `locking` and `twists` hold one entry per face. `twists` is an angle in radians,
    turning that face's site about its own normal.
  - `usecycle` forces an encoding; `nothing` picks [`cycleencoding`](@ref) whenever it is
    equivalent to [`dartencoding`](@ref), see [`_cycle_suffices`](@ref).
  - `touching_tol` and `alignment_tol` become the sites' tolerances.

## Why the order of the steps is forced

A face's first corner is its site's twist reference, and settling it takes three steps:

 1. `p` arrives from the [`Polyhedron`](@ref) constructor with each face already started at an
    intrinsically chosen corner, which pins the references up to each face's `sitesym`.
 2. That is enough to derive the labeling, since [`siteorbits`](@ref) compares frames only up
    to `sitesym`. Knowing the labeling gives the symmetry group.
 3. [`_propagate_faces`](@ref) re-winds along that group, pinning the references up to `stab`.

`twists` is applied last, once the references are fixed. It is folded into the key
[`siteorbits`](@ref) groups by, alongside the color, so that giving two faces of one orbit
different twists splits the orbit.
"""
# What a twist marks a face with, for the labels to tell faces apart by.
#
# A turn by the site's own symmetry angle carries its frame onto an equivalent one, so it marks
# nothing and has to reduce away: without this a face turned by a whole `2π` -- the identity --
# would read as distinct from its neighbours and cost the body symmetry it has. What survives the
# reduction is the part the face's own symmetry cannot absorb, which does mark it.
function _twistmarks(twists, sitesyms, ::Type{F}) where {F}
    return map(zip(twists, sitesyms)) do (t, sitesym)
        step = 2F(π) / sitesym
        r = mod(F(t), step)
        return min(r, abs(r - step)) < sqrt(eps(F)) ? zero(F) : r
    end
end

function _facesites(
    p::Polyhedron{F},
    poseof,
    colors,
    locking,
    twists,
    usecycle::Union{Nothing,Bool},
    touching_tol::Real,
    alignment_tol::Real,
) where {F}
    n = nfaces(p)
    sitesyms = facesym(p)
    labels = siteorbits(poseof(faces(p)), sitesyms, collect(zip(colors, _twistmarks(twists, sitesyms, F))))

    fs = _propagate_faces(corners(p), faces(p), labels)
    # A twist is an angle about the site's own normal. Whole dart steps are taken by rotating
    # the face's corner list, which moves the frame and the vertex numbering together and keeps
    # `contact_pairing`'s anchor on a pair of coincident darts; the remainder is applied to the
    # frame alone. That remainder is fine: the anchor has to be consistent and to tell
    # twists apart, not to mark a physical coincidence, and both survive.
    steps = [round(Int, t * length(f) / (2F(π))) for (f, t) in zip(fs, twists)]
    fs = [circshift(f, -mod(m, length(f))) for (f, m) in zip(fs, steps)]
    poses = [pose * RotX(F(t) - 2F(π) * m / length(f)) for (pose, t, m, f) in zip(poseof(fs), twists, steps, fs)]
    stabs = stabilizerorders(poses, sitesyms, labels)

    freedoms = [l ? s : g for (g, s, l) in zip(sitesyms, stabs, locking)]
    cyclic = something(usecycle, _cycle_suffices(freedoms, labels))
    g, ranges = cyclic ? cycleencoding(n; labels) : dartencoding(fs; labels)

    sites = [
        BindingSite(poses[i], colors[i], ranges[i], touching_tol, alignment_tol, sitesyms[i], stabs[i], locking[i]) for
        i in 1:n
    ]
    return g, sites
end

"""
    _cycle_suffices(twistfreedoms, labels)

Return `true` if one graph vertex per site suffices to encode the particle, so that
[`cycleencoding`](@ref) can stand in for [`dartencoding`](@ref).

Two things have to fit. The symmetry number must come out right, which needs all `labels`
distinct: no rotation can then preserve the labeling, so the answer is 1 whatever structure
the graph has internally. The bonds must be distinguishable, which means every site's twist 
freedom must be 1.
"""
_cycle_suffices(twistfreedoms, labels) = allunique(labels) && all(isone, twistfreedoms)

############### Polyhedra
_recenter(cs) = (c0=sum(cs) / length(cs); [c - c0 for c in cs])

# Scale a solid to edge length `a`, after construction so that the faces exist and the divisor
# is the true minimum edge length.
_scaleto(p::Polyhedron, a) = Polyhedron(corners(p) * (a / minedgelength(p)), faces(p))

"""
    Tetrahedron(a=1.0)

A regular tetrahedron with edge length `a`. Proper rotation group `T`, of order 12.
"""
function Tetrahedron(a::Real=1.0)
    F = float(typeof(a))
    cs = SVector{3,F}[(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
    return _scaleto(Polyhedron(cs), a)
end

"""
    Cube(a=1.0)

A cube with edge length `a`. Proper rotation group `O`, of order 24.
"""
function Cube(a::Real=1.0)
    F = float(typeof(a))
    cs = SVector{3,F}[(x, y, z) for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
    return _scaleto(Polyhedron(cs), a)
end

"""
    Octahedron(a=1.0)

A regular octahedron with edge length `a`. Proper rotation group `O`, of order 24.
"""
function Octahedron(a::Real=1.0)
    F = float(typeof(a))
    cs = SVector{3,F}[(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    return _scaleto(Polyhedron(cs), a)
end

"""
    Icosahedron(a=1.0)

A regular icosahedron with edge length `a`. Proper rotation group `I`, of order 60.
"""
function Icosahedron(a::Real=1.0)
    F = float(typeof(a))
    φ = (1 + sqrt(F(5))) / 2
    cs = SVector{3,F}[]
    for s1 in (-1, 1), s2 in (-1, 1)
        push!(cs, SVector{3,F}(0, s1, s2 * φ), SVector{3,F}(s1, s2 * φ, 0), SVector{3,F}(s2 * φ, 0, s1))
    end
    return _scaleto(Polyhedron(cs), a)
end

"""
    Dodecahedron(a=1.0)

A regular dodecahedron with edge length `a`. Proper rotation group `I`, of order 60.
"""
function Dodecahedron(a::Real=1.0)
    F = float(typeof(a))
    φ = (1 + sqrt(F(5))) / 2
    cs = SVector{3,F}[SVector{3,F}(x, y, z) for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
    for s1 in (-1, 1), s2 in (-1, 1)
        push!(cs, SVector{3,F}(0, s1 / φ, s2 * φ), SVector{3,F}(s1 / φ, s2 * φ, 0), SVector{3,F}(s2 * φ, 0, s1 / φ))
    end
    return _scaleto(Polyhedron(cs), a)
end

"""
    Pyramid(n, a=1.0; h=a)

A pyramid over a regular `n`-gon with edge length `a` and apex height `h`. Proper rotation
group `C_n`, of order `n`.
"""
function Pyramid(n::Integer, a::Real=1.0; h::Real=a)
    n >= 3 || throw(ArgumentError("a pyramid needs n >= 3, got n=$n"))
    F = float(promote_type(typeof(a), typeof(h)))
    r = F(a) / (2 * sin(F(π) / n))
    cs = SVector{3,F}[SVector{3,F}(r * cos(2F(π) * k / n), r * sin(2F(π) * k / n), 0) for k in 0:(n - 1)]
    push!(cs, SVector{3,F}(0, 0, F(h)))
    return Polyhedron(cs)
end

"""
    Prism(n, a=1.0; h=a)

A prism over a regular `n`-gon with edge length `a` and height `h`. Proper rotation group
`D_n`, of order `2n`. `Prism(4)` is a cube.
"""
function Prism(n::Integer, a::Real=1.0; h::Real=a)
    n >= 3 || throw(ArgumentError("a prism needs n >= 3, got n=$n"))
    F = float(promote_type(typeof(a), typeof(h)))
    r = F(a) / (2 * sin(F(π) / n))
    cs = SVector{3,F}[]
    for s in (-1, 1), k in 0:(n - 1)
        push!(cs, SVector{3,F}(r * cos(2F(π) * k / n), r * sin(2F(π) * k / n), s * F(h) / 2))
    end
    return Polyhedron(cs)
end

"""
    Antiprism(n, a=1.0)

A uniform antiprism over a regular `n`-gon with edge length `a`. Proper rotation group `D_n`,
of order `2n`. `Antiprism(3)` is an octahedron.
"""
function Antiprism(n::Integer, a::Real=1.0)
    n >= 3 || throw(ArgumentError("an antiprism needs n >= 3, got n=$n"))
    F = float(typeof(a))
    r = F(a) / (2 * sin(F(π) / n))
    h = F(a) * sqrt(1 - 1 / (4 * cos(F(π) / 2n)^2))
    cs = SVector{3,F}[]
    for k in 0:(n - 1)
        push!(cs, SVector{3,F}(r * cos(2F(π) * k / n), r * sin(2F(π) * k / n), -h / 2))
    end
    for k in 0:(n - 1)
        θ = 2F(π) * (k + F(1) / 2) / n
        push!(cs, SVector{3,F}(r * cos(θ), r * sin(θ), h / 2))
    end
    return Polyhedron(cs)
end

"""
    _sitetwists(psi, sitesym)

The orientations of a site that are equivalent by symmetry of the face: turns by `2π/sitesym` 
about the sites own outward normal, which `normal_pose` puts on the site's local x axis.
"""
_sitetwists(psi::Rotation{3,F}, sitesym::Integer) where {F} = (psi * RotX(F(2π) * m / sitesym) for m in 0:(sitesym - 1))
_sitetwists(psi::Rotation{2}, ::Integer) = (psi,)

"""
    sitelabel(ps::ParticleSpecies, i::Integer)

Return the symmetry label of site `i` of `ps`, the graph label all of that site's vertices carry.
"""
sitelabel(ps::ParticleSpecies, i::Integer) = Int(label(graphrep(ps), first(bindingsite(ps, i).vertices)))

"""
    permutationgroup(ps::ParticleSpecies)

Return the rotational symmetry group of `ps` as site permutations: one permutation per rotation
about the particle origin that maps every binding site onto a site with the same symmetry label,
matching orientation as well as position.

The counterpart of [`rotationgroup`](@ref), which returns the same group as rotations. Compare
[`symmetrynumber`](@ref), which reads the order of this group off the graph encoding instead;
[`check_encoding`](@ref) is the check that the two agree.
"""
function permutationgroup(ps::ParticleSpecies)
    perms = Vector{Int}[]
    _eachsitesymmetry((_, perm) -> push!(perms, perm), _sitegeometry(ps)...;
                      group=_rotationcandidates(ps))
    return perms
end

"""
    rotationgroup(ps::ParticleSpecies)

Return the rotational symmetry group of `ps` as rotations about the particle origin: those that
map every binding site onto a site with the same symmetry label, matching orientation as well as
position.

The counterpart of [`permutationgroup`](@ref), which returns the same group as site
permutations.
"""
function rotationgroup(ps::ParticleSpecies)
    rots = _rotationtype(posetype(ps))[]
    _eachsitesymmetry((Q, _) -> push!(rots, Q), _sitegeometry(ps)...; group=_rotationcandidates(ps))
    return rots
end

"""
    _rotationcandidates(ps::ParticleSpecies)

The rotations to test when asking for `ps`'s symmetry group, or `nothing` to derive them from its
sites; the `group` argument of [`_eachsitesymmetry`](@ref).

Only ever a *superset* of the answer, since every consumer drops the ones that move a site. A
species whose sites do not determine it overrides this, see the warning on
[`_eachsitesymmetry`](@ref).
"""
_rotationcandidates(::ParticleSpecies) = nothing

_rotationtype(::Type{<:Pose{D,F,R}}) where {D,F,R} = R

# Poses, site symmetries and the label each site is matched by. The group functions match on the
# graph's labels, since they ask what the graph claims; `siteorbits` matches on colors, since it
# asks what the arrangement is.
function _sitegeometry(ps::ParticleSpecies)
    n = nsites(ps)
    sites = [bindingsite(ps, i) for i in 1:n]
    return ([s.pose for s in sites], [s.sitesym for s in sites], [sitelabel(ps, i) for i in 1:n])
end

"""
    _eachsitesymmetry(f, poses, sitesyms, sitelabels; group=nothing)

Call `f(Q, perm)` on every rotation `Q` about the particle (= site poses, syms, and labels) origin that carries each site onto
one with the same `sitelabels` entry, matching position and orientation, together with the site
permutation `perm` it induces. Return how many there were.

`sitelabels` is whatever a symmetry has to preserve: graph labels when asking what the graph
claims, colors when asking what the arrangement is. Pass `f = Returns(nothing)` to only count.

`group` is where the candidate rotations come from. Given one, each of its elements is tried once
and those moving a site off the set are dropped, which is how a caller holding the particle's
symmetry group by other means reaches the same permutations.

Given `nothing`, the candidates are derived from the sites, by orbit and stabilizer: a rotation
is fixed by where it sends site 1's frame, it must send site 1 to a site of the same label, and
it may leave the frame in any of that site's `sitesym` equivalent turns. Every group element
appears exactly once, since the outer loop runs over site 1's orbit and the inner one over its
stabilizer.

!!! warning "When the group has to be supplied"
    Deriving it rests on two things. The sites must have *distinct positions*, which makes a
    rotation determined by where it sends site 1 and the site map injective for free. And the
    sites must determine the body, which an assembled [`Polyform`](@ref) does not: two open sites
    can sit in symmetric poses with different structures behind them. A
    [`MetaParticleSpecies`](@ref) fails both and supplies its `Polyform`'s own group instead.
"""
function _eachsitesymmetry(f, poses, sitesyms, sitelabels; group=nothing)
    n = length(poses)
    tol = sqrt(eps(eltype(typeof(first(poses)))))
    atol = tol * maximum(norm(p.x) for p in poses)

    function permutation(Q)
        perm = zeros(Int, n)
        for i in 1:n
            j = findfirst(1:n) do j
                sitelabels[j] == sitelabels[i] &&
                    isapprox(Q * poses[i].x, poses[j].x; atol) &&
                    any(psi -> isapprox(Q * poses[i].psi, psi; atol=tol), _sitetwists(poses[j].psi, sitesyms[j]))
            end
            isnothing(j) && return nothing
            perm[i] = j
        end
        return perm
    end

    count = 0
    if group !== nothing
        for Q in group
            perm = permutation(Q)
            isnothing(perm) && continue
            f(Q, perm)
            count += 1
        end
        return count
    end

    for a in 1:n
        sitelabels[a] == sitelabels[1] || continue
        for psi in _sitetwists(poses[a].psi, sitesyms[a])
            Q = psi * inv(poses[1].psi)
            perm = permutation(Q)
            isnothing(perm) && continue
            # Keep only the `a` that `Q` really sends site 1 to, so an element cannot be reached
            # twice. Two sites of one particle never share a frame, so this is a no-op here, but
            # it is what makes the orbit-stabilizer count exact rather than incidental.
            perm[1] == a || continue
            f(Q, perm)
            count += 1
        end
    end
    return count
end

"""
    _sitesymmetries(poses, sitesyms, sitelabels; group=nothing)

Return the site permutations of [`_eachsitesymmetry`](@ref), for callers that hold the site
geometry rather than a [`ParticleSpecies`](@ref).
"""
function _sitesymmetries(poses, sitesyms, sitelabels; group=nothing)
    perms = Vector{Int}[]
    _eachsitesymmetry((_, perm) -> push!(perms, perm), poses, sitesyms, sitelabels; group)
    return perms
end

"""
    _canonicalpartition(xs)

Renumber the classes of `xs` to `1:k` in the order their first member appears, so that two
partitions of the same sites compare equal with `==` exactly when they are the same partition.
[`siteorbits`](@ref) already returns this form.
"""
function _canonicalpartition(xs)
    seen = Dict{eltype(xs),Int}()
    return [get!(seen, x, length(seen) + 1) for x in xs]
end

"""
    _graphsiteorbits(ps::ParticleSpecies, vertexorbits)

Return one orbit index per site of `ps`, collapsed from the graph vertex orbits `vertexorbits`
that nauty reports for `graphrep(ps)`.

Sites are joined whenever they share a vertex orbit, which is what an automorphism carrying one
site onto another forces. A site's own vertices need *not* share an orbit: a square side face of
a triangular prism owns four darts, and the prism's six rotations split them into two orbits, so
reading the orbit of a site's first vertex would be wrong.
"""
function _graphsiteorbits(ps::ParticleSpecies, vertexorbits)
    n = nsites(ps)
    orbit = collect(1:n)
    firstsite = Dict{eltype(vertexorbits),Int}()
    for i in 1:n, v in bindingsite(ps, i).vertices
        j = get!(firstsite, vertexorbits[v], i)
        lo, hi = minmax(orbit[i], orbit[j])
        hi == lo || replace!(orbit, hi => lo)
    end
    return _canonicalpartition(orbit)
end

"""
    siteorbits(poses, sitesyms, colors; group=nothing)

Group the sites into the orbits of the rotations that preserve the colored arrangement, and
return one orbit index per site.
"""
function siteorbits(poses, sitesyms, colors; group=nothing)
    n = length(poses)
    orbit = collect(1:n)
    # we need sitesyms to quotient out the site symmetry groups
    for perm in _sitesymmetries(poses, sitesyms, colors; group)
        for i in 1:n
            lo, hi = minmax(orbit[i], orbit[perm[i]])
            hi == lo && continue
            replace!(orbit, hi => lo)
        end
    end
    # compact to 1:k in place, so the result can be used as graph labels directly.
    ids = sort!(unique(orbit))
    for i in eachindex(orbit)
        orbit[i] = searchsortedfirst(ids, orbit[i])
    end
    return orbit
end

"""
    stabilizerorders(ps::ParticleSpecies)

Return the order of each site's stabilizer: how many of the particle's own
symmetries leave the site where it is.
"""
stabilizerorders(ps::ParticleSpecies) = [bindingsite(ps, i).stab for i in 1:nsites(ps)]

"""
    stabilizerorders(poses, sitesyms, sitelabels; group=nothing)

Return the order of each site's stabilizer: how many of the particle's own
symmetries leave the site where it is.
"""
function stabilizerorders(poses, sitesyms, sitelabels; group=nothing)
    perms = _sitesymmetries(poses, sitesyms, sitelabels; group)
    return [count(perm -> perm[i] == i, perms) for i in eachindex(poses)]
end

"""
    _recolor!(ps::ParticleSpecies, sites, colors)

Give `sites` the interaction colors `colors`, and update `ps`' labeling and the stabilizers.

Throws if the recoloring asks for a symmetry the graph encoding cannot express.
"""
function _recolor!(ps::ParticleSpecies, sites::AbstractVector{<:BindingSite}, colors)
    # store old labeling and restore on error
    oldsites, oldlabels = copy(sites), copy(labels(graphrep(ps)))
    _recolor!(sites, graphrep(ps), colors)
    try
        check_encoding(ps)
    catch err
        err isa ArgumentError || rethrow()
        copy!(sites, oldsites)
        setlabels!(graphrep(ps), oldlabels)
        throw(
            ArgumentError(
                "this recoloring changes the particle's symmetry by more than its graph encoding " *
                "can express, so it cannot be applied. Build the species with its dartencoding instead ($(err.msg))",
            ),
        )
    end
    return nothing
end

"""
    _recolor!(sites, g::NautyDiGraph, colors)

Give `sites` the interaction colors `colors`, and update `g`'s labeling and the stabilizers.
"""
function _recolor!(sites::AbstractVector{<:BindingSite}, g::NautyDiGraph, colors)
    length(colors) == length(sites) || throw(ArgumentError("incorrect number of colors"))
    poses = [s.pose for s in sites]
    sitesyms = [s.sitesym for s in sites]
    orbits = siteorbits(poses, sitesyms, collect(colors))
    stabs = stabilizerorders(poses, sitesyms, orbits)

    for i in eachindex(sites)
        sites[i] = setstab(setcolor(sites[i], colors[i]), stabs[i])
        for v in sites[i].vertices
            setlabel!(g, v, orbits[i])
        end
    end
    return nothing
end

"""
    _check_labeling(ps::ParticleSpecies)

Throw unless `ps`'s labeling is at least as fine as its coloring, and is exactly its symmetry
orbits: two sites carry the same label if and only if a rotation of the particle carries one
onto the other.
"""
function _check_labeling(ps::ParticleSpecies)
    seen = Dict{Int,Int}()
    for i in 1:nsites(ps)
        b = bindingsite(ps, i)
        l = sitelabel(ps, i)
        c = get!(seen, l, color(b))
        c == color(b) || throw(
            ArgumentError(
                "Sites sharing symmetry label $l have colors $c and $(color(b)). A labeling must " *
                "be at least as fine as the coloring.",
            ),
        )
    end

    labeling = _canonicalpartition([sitelabel(ps, i) for i in 1:nsites(ps)])
    orbits = siteorbits(_sitegeometry(ps)...; group=_rotationcandidates(ps))
    labeling == orbits || throw(
        ArgumentError(
            "This labeling groups the sites as $labeling, but the rotations preserving it group " *
            "them as $orbits. A labeling has to be exactly the symmetry orbits: " *
            (
                if length(unique(labeling)) < length(unique(orbits))
                    "here it declares sites equivalent that no rotation maps onto each other. "
                else
                    "here it splits sites that a rotation does map onto each other. "
                end
            ) *
            "Derive it with `siteorbits` rather than writing it out.",
        ),
    )
    return ps
end

"""
    check_encoding(ps::ParticleSpecies)

Throw if `ps`'s graph claims a different symmetry from its geometry and colors, or if its
labeling is not exactly its symmetry orbits. Return `ps`.

Three things have to agree:

 1. the labeling is at least as fine as the coloring, and is exactly the symmetry orbits
 2. the graph's automorphism group and the geometry's rotation group have the same order
 3. they also have the same orbits, since two groups of equal order can act differently

Every built-in species runs this in its constructor. Call it in your own if you write graph
labels by hand instead of deriving them with [`siteorbits`](@ref).
"""
function check_encoding(ps::ParticleSpecies)
    _check_labeling(ps)
    geom, group = _sitegeometry(ps), _rotationcandidates(ps)
    geometric = _eachsitesymmetry(Returns(nothing), geom...; group)
    # Not `canonize`: a species' graph must stay in construction order, see `symmetrynumber`.
    _, autg = nauty(graphrep(ps))
    graph = convert(Int, autg.n)
    graph == geometric || throw(
        ArgumentError(
            "Graph encoding claims a symmetry number of $graph, but the geometry of this " *
            "$(dimension(ps))d species have a rotational symmetry of $geometric. " *
            (
                if graph > geometric
                    "The labels declare sites equivalent that no rotation maps onto each other."
                else
                    "The graph distinguishes sites that a rotation does map onto each other; the " *
                    "encoding does not describe this site arrangement."
                end
            ),
        ),
    )

    orbits = siteorbits(geom...; group)
    fromgraph = _graphsiteorbits(ps, autg.orbits)
    orbits == fromgraph || throw(
        ArgumentError(
            "Graph encoding and geometry agree that this $(dimension(ps))d species has a " *
            "rotational symmetry of $graph, but not on which sites it relates: the geometry " *
            "groups them as $orbits and the graph as $fromgraph. Two groups of the same order " *
            "can act differently, and the enumeration reads the graph's action.",
        ),
    )
    return ps
end
