"""
    AbstractPolyform{D}

An aggregate of particles in `D` dimensions, connected at binding sites and represented by a directed graph.
"""
abstract type AbstractPolyform{D} end

"""
    Polyform

A `Polyform` is an aggregate of particles connected at binding
sites, represented by a directed graph.
"""
mutable struct Polyform{D,P<:Particle,S<:BindingRules,G<:AbstractNautyGraph} <: AbstractPolyform{D}
    graphrep::G
    sigma::Int
    canon2orig::Vector{Int}
    orig2canon::Vector{Int}
    particles::Vector{P}
    bindingrules::S
end

"""
    Polyform(rules::BindingRules{D}) where {D}

Create an empty polyform containing no particles.
"""
function Polyform(rules::BindingRules{D}) where {D}
    P = posetype(rules)
    g = NautyDiGraph(0)
    return Polyform{D,Particle{P},typeof(rules),typeof(g)}(g, 1, Int[], Int[], Particle{P}[], rules)
end

"""
    Polyform(rules::BindingRules{D}, i::Integer) where {D}

Create a single-particle polyform, consisting of species `i` of `rules`.
"""
function Polyform(rules::BindingRules{D}, i::Integer) where {D}
    P = posetype(rules)
    ps = species(rules, i)
    g = copy(graphrep(ps))
    part = Particle(rules, i; leadingvertex=1)
    perm = first(nauty(g; canonize=true))
    cvs = convert(Vector{Int}, perm)
    return Polyform{D,Particle{P},typeof(rules),typeof(g)}(g, symmetrynumber(ps), cvs, invperm(cvs), [part], rules)
end

"""
    translate!(p::AbstractPolyform, v)

Shift every particle of `p` by the vector `v`, leaving its graph and orientations untouched.
"""
function translate!(p::AbstractPolyform, v)
    map!(q -> translate(q, v), p.particles, p.particles)
    return p
end

"""
    translate(p::AbstractPolyform, v)

Return a copy of `p` shifted by the vector `v`.
"""
translate(p::AbstractPolyform, v) = translate!(copy(p), v)

function Base.copy(p::Polyform)
    return typeof(p)(
        copy(p.graphrep), p.sigma, copy(p.canon2orig), copy(p.orig2canon), copy(p.particles), p.bindingrules
    )
end
function Base.copy!(dst::Polyform, src::Polyform)
    copy!(dst.graphrep, src.graphrep)
    dst.sigma = src.sigma
    copy!(dst.canon2orig, src.canon2orig)
    copy!(dst.orig2canon, src.orig2canon)
    copy!(dst.particles, src.particles)
    dst.bindingrules = src.bindingrules
    return dst
end

Base.show(io::Core.IO, p::Polyform{D}) where {D} = print(io, "Polyform{$D}[n=$(nparticles(p))]")
Base.show(io::Core.IO, ::Type{Polyform{D}}) where {D} = print(io, "Polyform{$D}")
Base.:(==)(p1::Polyform, p2::Polyform) = bindingrules(p1) === bindingrules(p2) && graphrep(p1) == graphrep(p2)

"""
    ParticleSite(particle, site)

Which binding site of which particle inside a [`Polyform`](@ref): where a site actually is in an
assembled structure.

Distinct from [`SpeciesSite`](@ref), which names a site of a *species* and is what a set of
[`BindingRules`](@ref) speaks in. Both were `NTuple{2,Int}` once.

Iterates and indexes like the pair it replaces, so `(p, k) = ps` still works.
"""
struct ParticleSite
    particle::Int
    site::Int
end

Base.iterate(l::ParticleSite, i::Int=1) = i > 2 ? nothing : (getfield(l, i), i + 1)
Base.length(::ParticleSite) = 2
Base.getindex(l::ParticleSite, i::Integer) = getfield(l, Int(i))
Base.show(io::Core.IO, l::ParticleSite) = print(io, "ParticleSite(", l.particle, ", ", l.site, ")")
Base.isless(a::ParticleSite, b::ParticleSite) = (a.particle, a.site) < (b.particle, b.site)

"""
    nparticles(p::AbstractPolyform)

Return the number of particles in `p`.
"""
@inline nparticles(p::AbstractPolyform) = length(p.particles)

"""
    nsites(p::AbstractPolyform)

Return the total number of binding sites across all particles in `p`, including bound sites.
"""
@inline nsites(p::AbstractPolyform) = sum(prt -> nsites(prt, p.bindingrules), p.particles; init=0)

"""
    bindingrules(p::AbstractPolyform)

Return the `BindingRules` that `p` belongs to.
"""
@inline bindingrules(p::AbstractPolyform) = p.bindingrules

"""
    symmetrynumber(p::AbstractPolyform)

Return the symmetry number of `p`, i.e. the size of its automorphism group.
"""
@inline symmetrynumber(p::AbstractPolyform) = p.sigma

@inline graphrep(p::AbstractPolyform) = p.graphrep
@inline dimension(::AbstractPolyform{D}) where {D} = D
@inline posetype(p::AbstractPolyform) = posetype(bindingrules(p))
@inline posetype(::Type{<:Polyform{D,<:Particle{<:P}}}) where {D,P} = P
@inline numtype(p::AbstractPolyform) = eltype(posetype(p))

"""
    particletype(p::AbstractPolyform)

The concrete `Particle` type of every particle in `p`, for sizing a container that holds them.
"""
@inline particletype(p::AbstractPolyform) = particletype(bindingrules(p))

"""
    sitetype(p::AbstractPolyform)

The concrete [`BindingSite`](@ref) type of every site of `p`, for sizing a container that holds
them.
"""
@inline sitetype(p::AbstractPolyform) = sitetype(bindingrules(p))
@inline numtype(::Type{<:Polyform{D,<:Particle{<:P}}}) where {D,P} = eltype(P)

"""
    tocanon(p::AbstractPolyform, v::Integer)

Convert an original vertex index `v` (as stored in `Particle.leadingvertex` or
`BindingSite.vertices`) to the corresponding canonical graph vertex index in `graphrep(p)`.
"""
@inline tocanon(p::AbstractPolyform, v::Integer) = p.orig2canon[v]

"""
    toorig(p::AbstractPolyform, v::Integer)

Convert a canonical graph vertex index `v` (as returned by iterating `graphrep(p)`)
to the corresponding stable original vertex index.
"""
@inline toorig(p::AbstractPolyform, v::Integer) = p.canon2orig[v]

function _apply_perm!(poly::Polyform, perm)
    # `orig2canon` is rebuilt from the result three lines below, 
    # so we borrow it as scratch buffer
    resize!(poly.orig2canon, length(poly.canon2orig))
    @inbounds for i in eachindex(poly.orig2canon, perm)
        poly.orig2canon[i] = poly.canon2orig[perm[i]]
    end
    @inbounds for i in eachindex(poly.canon2orig, poly.orig2canon)
        poly.canon2orig[i] = poly.orig2canon[i]
    end
    @inbounds for i in eachindex(poly.canon2orig)
        poly.orig2canon[poly.canon2orig[i]] = i
    end
    return nothing
end

@inline is_leadingvertex(p::AbstractPolyform, v::Integer) = any(pt -> pt.leadingvertex == v, p.particles)

@inline particles(p::AbstractPolyform, i::Integer) = p.particles[i]

# Look up the particle whose leadingvertex equals the original vertex v (O(n) scan).
@inline function particle_from_leadingvertex(p::AbstractPolyform, v::Integer)
    i = findfirst(pt -> pt.leadingvertex == v, p.particles)
    return isnothing(i) ? nothing : p.particles[i]
end

"""
    interior_edges(p::AbstractPolyform)

Return a lazy iterator over the internal particle edges of `graphrep(p)`.
"""
@inline interior_edges(p::AbstractPolyform) = _filter_edges(p, Val(false))

"""
    exterior_edges(p::AbstractPolyform)

Return a lazy iterator over the external edges of `graphrep(p)`, those joining two particles.

A bond contributes *several* of these, one per vertex pair `contact_pairing` makes, so use
[`bonds`](@ref) to iterate bonds.
"""
@inline exterior_edges(p::AbstractPolyform) = (e for e in _filter_edges(p, Val(true)) if e.src < e.dst)

"""
    _same_particle(p::AbstractPolyform, u::Integer, v::Integer; canonidxs)

Return `true` if the graph vertices `u` and `v` belong to the same particle.

`canonidxs` says which numbering `u` and `v` are in, and has no default: a caller that does not
say is a caller that has not thought about it, which is how a bond scan once read canonical
indices as original ones and silently lost a class of environment.
"""
@inline function _same_particle(p::AbstractPolyform, u::Integer, v::Integer; canonidxs::Bool)
    # Each particle owns a contiguous block of original vertices starting at its leading vertex, so
    # `u` and `v` are split apart exactly when some leading vertex falls between them.
    canonidxs && ((u, v) = (toorig(p, u), toorig(p, v)))
    lo, hi = minmax(u, v)
    return !any(pt -> lo < leadingvertex(pt) <= hi, p.particles)
end

# An edge is a bond exactly when its endpoints belong to different particles
function _filter_edges(p::AbstractPolyform, ::Val{exterior}) where {exterior}
    return Iterators.filter(edges(graphrep(p))) do (; src, dst)
        same = _same_particle(p, src, dst; canonidxs=true)
        return exterior ? !same : same
    end
end

"""
    _isbound_vertex(p::AbstractPolyform, part::Particle, v::Integer; canonidxs)

Return `true` if the graph vertex `v` of particle `part` is bonded to another particle, i.e. if
it has a neighbor outside `part`'s own block of vertices.

`canonidxs` says which numbering `v` is in, the same way it does for [`bondindex`](@ref) and
[`_same_particle`](@ref).
"""
function _isbound_vertex(p::Polyform, part::Particle, v::Integer; canonidxs::Bool)
    own = graphvertices(part, bindingrules(p))
    neighs = NautyGraphs.adjrow(graphrep(p), canonidxs ? v : tocanon(p, v))
    for w in eachindex(neighs)
        neighs[w] || continue
        toorig(p, w) in own || return true
    end
    return false
end

"""
    bondindex(poly::AbstractPolyform, src::Integer, dst::Integer; canonidxs=true)

Return the index into `bonded_colors(bindingrules(poly))` for the bond between the graph
vertices `src` and `dst`, or `nothing` if they don't form a valid bond type.

`canonidxs` says which numbering `src` and `dst` are in: `true`, the default, for the canonical
one that `graphrep(poly)` and its `edges` are in, and `false` for the stable original one that
`BindingSite.vertices` and `Particle.leadingvertex` are in.
"""
function bondindex(poly::AbstractPolyform, src::Integer, dst::Integer; canonidxs::Bool=true)
    a = _vertex_to_particle_site(poly, src; canonidxs)
    b = _vertex_to_particle_site(poly, dst; canonidxs)
    return bondindex(poly, a, b)
end

"""
    bondindex(poly::AbstractPolyform, a::ParticleSite, b::ParticleSite)

Return the index into `bonded_colors(bindingrules(poly))` for a bond between the sites `a` and
`b` of `poly`, or `nothing` if their colors form no valid bond type.
"""
function bondindex(poly::AbstractPolyform, a::ParticleSite, b::ParticleSite)
    c1, c2 = color(bindingsite(poly, a)), color(bindingsite(poly, b))
    return findfirst(==(minmax(c1, c2)), bonded_colors(bindingrules(poly)))
end

# Map a graph vertex back to (particleindex, siteindex).
function _vertex_to_particle_site(p::AbstractPolyform, v::Integer; canonidxs::Bool)
    rules = bindingrules(p)
    orig_v = canonidxs ? toorig(p, v) : v
    for (i, part) in enumerate(p.particles)
        orig_v in graphvertices(part, rules) || continue
        for j in 1:nsites(part, rules)
            orig_v in bindingsite(part, rules, j).vertices && return ParticleSite(i, j)
        end
    end
    return nothing
end

"""
    bindingsite(p::AbstractPolyform, loc::ParticleSite)

The binding site `loc` names: site `loc.site` of particle `loc.particle` of `p`.
"""
@inline function bindingsite(p::AbstractPolyform, loc::ParticleSite)
    return bindingsite(particles(p, loc.particle), bindingrules(p), loc.site)
end

"""
    bindingsite(p::AbstractPolyform, i::Integer)

Return the `i`-th binding site of `p`, counting through `p`'s particles in the order they are
stored and through each particle's own sites, exactly as on a [`ParticleSpecies`](@ref).

The ordering depends on how `p` was assembled. Use [`canonbindingsite`](@ref) for iterating through
binding sites in canonical order.
"""
function bindingsite(p::AbstractPolyform, i::Integer)
    rules = bindingrules(p)
    k = 0
    for prtcl in p.particles
        for j in 1:nsites(prtcl, rules)
            k += 1
            k == i && return bindingsite(prtcl, rules, j)
        end
    end
    return nothing
end

"""
    bindingsites(p::AbstractPolyform)

Return a lazy iterator over all binding sites of `p`, counting through `p`'s particles in the order they are
stored and through each particle's own sites, exactly as on a [`ParticleSpecies`](@ref).

The ordering depends on how `p` was assembled. Use [`canonbindingsite`](@ref) for iterating through
binding sites in canonical order.
"""
bindingsites(p::AbstractPolyform) = (bindingsite(p, i) for i in 1:nsites(p))

"""
    siteindex(p::AbstractPolyform, siteloc::ParticleSite)

Return the index `i` of site `siteloc`, inverting `bindingsite(p, i)`.
"""
function siteindex(p::AbstractPolyform, siteloc::ParticleSite)
    rules = bindingrules(p)
    return sum(nsites(p.particles[q], rules) for q in 1:(siteloc.particle - 1); init=0) + siteloc.site
end

"""
    canonbindingsite(p::AbstractPolyform, i::Integer)

Return the `i`-th binding site of `p` in canonical order, which follows the canonical graph
labeling and is therefore the same for any two isomorphic polyforms.
"""
function canonbindingsite(p::AbstractPolyform, i::Integer)
    rules = bindingrules(p)
    k = 0
    for v in p.canon2orig
        prtcl = particle_from_leadingvertex(p, v)
        isnothing(prtcl) && continue
        for j in 1:nsites(prtcl, rules)
            k += 1
            k == i && return bindingsite(prtcl, rules, j)
        end
    end
    return nothing
end

"""
    canonbindingsites(p::AbstractPolyform)

Return a lazy iterator over all binding sites of `p` in canonical order.
"""
canonbindingsites(p::AbstractPolyform) = (canonbindingsite(p, i) for i in 1:nsites(p))

"""
    rotationcenter(p::Polyform)

The point every rotational symmetry of `p` fixes: the centroid of its particle positions.

A symmetry permutes the particles, so it fixes their centroid, and [`rotationgroup`](@ref
rotationgroup(::Polyform)) turns about this point rather than about `p`'s pose origin.
"""
rotationcenter(p::Polyform) = sum(pt.pose.x for pt in p.particles) / nparticles(p)

"""
    permutationgroup(p::Polyform)

Return the rotational symmetry group of `p` as particle permutations: one permutation per
rotation about the centroid of `p`'s particle positions that carries every particle onto a
particle of the same species, matching position and orientation.

The counterpart of [`rotationgroup(::Polyform)`](@ref), which returns the same group as
rotations, and the assembled analogue of [`permutationgroup(::ParticleSpecies)`](@ref), which
permutes a single particle's sites.

`length(permutationgroup(p))` should equal [`symmetrynumber`](@ref)`(p)`, which reads the same
order off the graph encoding. Use it to check an encoding rather than to compute the symmetry
number, which nauty gives far more cheaply.

The list holds one entry per rotation, so it may repeat: the action on particles is not
faithful. A straight chain of cubes turned about the chain axis leaves every particle where it
was and still moves the graph. Call `unique` on the result for the quotient group.
"""
function permutationgroup(p::Polyform)
    perms = Vector{Int}[]
    _eachpolyformsymmetry((_, perm) -> push!(perms, perm), p)
    return perms
end

"""
    rotationgroup(p::Polyform)

Return the rotational symmetry group of `p` as rotations about the centroid of `p`'s particle
positions: those carrying every particle onto a particle of the same species, matching position
and orientation.

The counterpart of [`permutationgroup(::Polyform)`](@ref), which returns the same group as
particle permutations.
"""
function rotationgroup(p::Polyform)
    rots = _rotationtype(posetype(p))[]
    _eachpolyformsymmetry((Q, _) -> push!(rots, Q), p)
    return rots
end

"""
    _eachpolyformsymmetry(f, p::Polyform)

Call `f(Q, perm)` on every rotation `Q` about the centroid of `p`'s particle positions that
carries each particle onto one of the same species, matching position and orientation, together
with the particle permutation `perm` it induces. Return how many there were.

A particle's frame need only match up to its species' own
[`rotationgroup`](@ref rotationgroup(::ParticleSpecies)), so the candidates are `Q = Rₐ S R₁⁻¹`
for each particle `a` of particle 1's species and each `S` in that group. Every element of the
group appears exactly once.
"""
function _eachpolyformsymmetry(f, p::Polyform)
    F = numtype(p)
    n = nparticles(p)
    # An empty polyform constrains nothing; report the identity, so the count still matches the
    # symmetry number its graph carries.
    n == 0 && (f(one(_rotationtype(posetype(p))), Int[]); return 1)

    rules = bindingrules(p)
    poses = [pt.pose for pt in p.particles]
    spcs = [speciesindex(pt) for pt in p.particles]
    own = Dict(s => rotationgroup(species(rules, s)) for s in unique(spcs))

    xs = [q.x - rotationcenter(p) for q in poses]
    tol = sqrt(eps(F))
    atol = tol * max(one(F), maximum(norm, xs))

    function permutation(Q)
        perm = zeros(Int, n)
        for i in 1:n
            j = findfirst(1:n) do j
                spcs[j] == spcs[i] &&
                    isapprox(Q * xs[i], xs[j]; atol) &&
                    any(S -> isapprox(Q * poses[i].psi, poses[j].psi * S; atol=tol), own[spcs[j]])
            end
            isnothing(j) && return nothing
            perm[i] = j
        end
        return perm
    end

    nfound = 0
    for a in 1:n
        spcs[a] == spcs[1] || continue
        for S in own[spcs[1]]
            Q = poses[a].psi * S * inv(poses[1].psi)
            perm = permutation(Q)
            isnothing(perm) && continue
            # A candidate is built from frames alone, and two particles of one species often
            # carry the same frame, so the same `Q` is reached once per such particle. Keeping
            # only the `a` that `Q` really sends particle 1 to picks each element out once.
            perm[1] == a || continue
            f(Q, perm)
            nfound += 1
        end
    end
    return nfound
end

"""
    _bondedges(p::Polyform)

Return one exterior edge per bond of `p`.

A bond joins two binding sites, and [`contact_pairing`](@ref) pairs `gcd(k₁, k₂)` of their
vertices, so a bond between two dart-encoded faces reaches `graphrep(p)` as several edges. Two
sites are joined by at most one bond, so the pair of sites an edge lands on names the bond.
"""
function _bondedges(p::Polyform)
    seen = Set{NTuple{2,ParticleSite}}()
    return filter(collect(exterior_edges(p))) do (; src, dst)
        key = minmax(_vertex_to_particle_site(p, src; canonidxs=true),
                     _vertex_to_particle_site(p, dst; canonidxs=true))
        key in seen && return false
        push!(seen, key)
        return true
    end
end

"""
    bonds(p::Polyform)

Return a lazy iterator of bonds in `p`, each one once, as [`ParticleSite`](@ref) pairs.
"""
function bonds(p::Polyform)
    return (
        let (lhs, rhs) = minmax(_vertex_to_particle_site(p, e.src; canonidxs=true),
                                _vertex_to_particle_site(p, e.dst; canonidxs=true))
            lhs => rhs
        end for e in _bondedges(p)
    )
end

"""
    nbonds(p::AbstractPolyform)

Return the number of bonds in `p`. Note that `nbonds(::BindingRules)` instead counts how many
*kinds* of bond a set of rules allows.
"""
nbonds(p::AbstractPolyform) = count(Returns(true), bonds(p))

"""
    bondtypes(p::AbstractPolyform)

Return the bond type of every bond of `p`, each one once, indexing `bonded_colors(bindingrules(p))`.

For a [`Tiling`](@ref) these are the bonds of one cell of the infinite structure, the ones a
translate closes counted alongside the ones inside the cell.
"""
function bondtypes(p::AbstractPolyform)
    return map(bonds(p)) do (a, b)
        i = bondindex(p, a, b)
        isnothing(i) && error("Internal error: a bond has no bond type. Please file an issue.")
        return i
    end
end

"""
    composition(p::AbstractPolyform)

Return the composition vector of `p`: counts of each particle species (indices
`1:nspecies(rules)`) followed by counts of each bond type (indices `nspecies+1:end`).
Bond types are ordered as in `bonded_colors(bindingrules(p))`.
"""
function composition(p::AbstractPolyform)
    rules = bindingrules(p)
    ns = nspecies(rules)
    nb = nbonds(rules)
    comp = zeros(Int, ns + nb)

    for part in p.particles
        comp[part.speciesindex] += 1
    end

    for (a, b) in bonds(p)
        i = bondindex(p, a, b)
        isnothing(i) || (comp[ns + i] += 1)
    end

    return comp
end

"""
    raise!(poly::Polyform, site::BindingSite, loc::SpeciesSite, t=0)

Attach a new particle to `poly` at the open binding site `site`, with the species and site index given by `loc`,
in twist `t` of the bond (see [`standard_twist`](@ref)).

Returns `poly` on success, or `missing` if the attachment is geometrically forbidden (overlap or misaligned contact).
"""
function raise!(poly::Polyform, site::BindingSite, loc::SpeciesSite, t::Integer=0; kwargs...)
    rules = bindingrules(poly)
    particle_species = species(rules, loc.species)
    leadingvertex = nv(graphrep(poly)) + 1
    mate = bindingsite(particle_species, loc.site)
    particle_pose = standard_twist(site, t, twistfreedom(site, mate)) * inv(mate.pose)
    attached_particle = Particle(particle_pose, leadingvertex, loc.species)

    overlap, contacts = overlap_and_contacts(poly, attached_particle; kwargs...)
    overlap && return missing

    push!(poly.particles, attached_particle)

    g_attach = graphrep(particle_species)
    add_vertices!(graphrep(poly); vertex_labels=labels(g_attach))
    for (; src, dst) in edges(g_attach)
        add_edge!(graphrep(poly), src + leadingvertex - 1, dst + leadingvertex - 1)
    end

    for contact in contacts
        for (v1, v2) in contact_pairing(contact)
            # the vertices corresponding to the particles on poly need to be transformed to canonical order
            add_edge!(graphrep(poly), tocanon(poly, v1), v2)
            add_edge!(graphrep(poly), v2, tocanon(poly, v1))
        end
    end

    append!(poly.canon2orig, graphvertices(attached_particle, rules))
    perm, autg = nauty(graphrep(poly); canonize=true)
    _apply_perm!(poly, perm)
    poly.sigma = convert(Int, autg.n)
    return poly
end

"""
    lower!(poly::Polyform)

Remove the last particle from `poly`, following canonical ordering.

Returns `poly` on success, `nothing` if `poly` is already empty.
"""
function lower!(poly::Polyform)
    n = nparticles(poly)
    n == 0 && return nothing
    if n == 1
        pop!(poly.particles)
        resize!(poly.canon2orig, 0)
        resize!(poly.orig2canon, 0)
        rem_vertices!(graphrep(poly), vertices(graphrep(poly)))
        poly.sigma = 1
        return poly
    end

    rules = bindingrules(poly)
    nv_g = nv(graphrep(poly))
    target = zeros(Bool, nv_g)
    dist = zeros(Int, nv_g)
    queue = zeros(Int, nv_g)

    # A particle owns >=1 vertex per site, so the scan will revisit vertices many times.
    # -> store for performance.
    tested = Int[]
    part = nothing
    for v in Iterators.reverse(poly.canon2orig)
        # Walk backward from v to find the leading vertex of its particle.
        for k in v:-1:1
            is_leadingvertex(poly, k) || continue
            v = k
            break
        end
        v in tested && continue
        push!(tested, v)

        part = particle_from_leadingvertex(poly, v)
        is_cutset(graphrep(poly), @view(poly.orig2canon[graphvertices(part, rules)]); target, dist, queue) || break
        part = nothing
    end
    isnothing(part) && error("Internal error: no removable particle found in connected polyform. Please file an issue.")
    return _remove_particle!(poly, part)
end

"""
    _remove_particle!(poly::Polyform, part::Particle)

Remove particle `part` from `poly`. Does not check for connectedness before.
"""
function _remove_particle!(poly::Polyform, part::Particle)
    rules = bindingrules(poly)
    vs0 = graphvertices(part, rules)
    vs = sort!(poly.orig2canon[vs0])

    lv = leadingvertex(part)
    idx = findfirst(pt -> pt.leadingvertex == lv, poly.particles)
    last_idx = lastindex(poly.particles)
    idx < last_idx && (poly.particles[idx] = poly.particles[last_idx])
    pop!(poly.particles)

    rem_vertices!(graphrep(poly), vs)
    deleteat!(poly.canon2orig, vs)

    # shift remaining indices and leading vertices
    for i in eachindex(poly.canon2orig)
        poly.canon2orig[i] > last(vs0) && (poly.canon2orig[i] -= length(vs0))
    end
    for i in eachindex(poly.particles)
        leadingvertex(poly.particles[i]) > last(vs0) || continue
        poly.particles[i] = shift_leadingvertex(poly.particles[i], -length(vs0))
    end

    perm, autg = nauty(graphrep(poly); canonize=true)
    _apply_perm!(poly, perm)
    poly.sigma = convert(Int, autg.n)
    return poly
end

"""
    subpolyform(poly::Polyform, particleids)

Return the sub-polyform induced by the particles `particleids`, renumbered in the given order.

The particle subset is expected to be connected; bonds to removed particles are dropped and their
sites become unbound.
"""
function subpolyform(poly::Polyform, particleids)
    rules = bindingrules(poly)
    newparticles = particletype(poly)[]
    verts = Int[]
    lv = 1
    for p in particleids
        part = particles(poly, p)
        for ov in graphvertices(part, rules)
            push!(verts, tocanon(poly, ov))
        end
        push!(newparticles, typeof(part)(part.pose, lv, part.speciesindex))
        lv += nv(graphrep(species(rules, part.speciesindex)))
    end

    g = graphrep(poly)[verts]
    # `g` starts out in new-original vertex order, so the canon maps are the permutation itself.
    perm, autg = nauty(g; canonize=true)
    return typeof(poly)(g, convert(Int, autg.n), collect(Int, perm), invperm(perm), newparticles, rules)
end

# Species `i` of `from` and species `i` of `rules` should be a "compatible" species. What that means exactly can be
# refined later. Right now, we check the vertex count and site poses.
# TODO: this is overly restrictive
function _checkspecies(from::BindingRules, rules::BindingRules, i::Integer)
    ok = i <= nspecies(rules)
    if ok
        a, b = species(from, i), species(rules, i)
        ok =
            nsites(a) == nsites(b) &&
            nv(graphrep(a)) == nv(graphrep(b)) &&
            all(k -> _siteoverlap(bindingsite(a, k), bindingsite(b, k)), 1:nsites(a))
    end
    ok || throw(
        ArgumentError(
            "species $i of `rules` is not the species that would stand in its place. `rules` has " *
            "to list the corresponding species in the same order.",
        ),
    )
    return nothing
end

# Check every species of `poly`, against the ones `rules` has in those places.
function _checkspecies(poly::Polyform, rules::BindingRules)
    foreach(q -> _checkspecies(bindingrules(poly), rules, speciesindex(q)), poly.particles)
    return poly
end

# check if the sites are at the same location and pose
function _siteoverlap(sa::BindingSite, sb::BindingSite)
    return isapprox(sa.pose.x, sb.pose.x; atol=sa.touching_tolerance + sb.touching_tolerance, rtol=0) &&
           isapprox(sa.pose.psi, sb.pose.psi; atol=sa.alignment_tolerance + sb.alignment_tolerance, rtol=0)
end

# Return the polyforms each species of `src` is replaced by, one per species, as a polyform of `rules`.
# Return `given` if present, return the underlying polyform of a meta species, or simply return the corresponding species.
function _resolvesubstitutions(src::BindingRules, rules::BindingRules, given)
    return map(collect(enumerate(species(src)))) do (i, ps)
        sub = get(given, i, nothing)
        if !isnothing(sub)
            bindingrules(sub) === rules || throw(ArgumentError("A substitution has to be a polyform of `rules`."))
            return sub
        end
        return _defaultsubstitution(ps, src, rules, i)
    end
end

# The substitution for species `i` of `src` when none was given; returns species `i` of `rules` and checks if valid
function _defaultsubstitution(::ParticleSpecies, src::BindingRules, rules::BindingRules, i::Integer)
    _checkspecies(src, rules, i)
    return Polyform(rules, i)
end

# Substitute every particle of `poly` of species `i` with the particles of polyform `subs[i]`.
function _substitute_particles(poly::Polyform, subs)
    P = particletype(first(subs))
    parts = P[]
    off = 0
    for part in poly.particles
        sub = subs[speciesindex(part)]
        subrules = bindingrules(sub)
        for q in sub.particles
            si = speciesindex(q)
            push!(parts, P(part.pose * q.pose, off + 1, si))
            off += nv(graphrep(species(subrules, si)))
        end
    end
    return parts
end

"""
    recast(poly::Polyform, rules::BindingRules; substitutions=Dict())

Recast `poly` as a [`Polyform`](@ref) of `rules`, substituting every particle with one or multiple particles
from `rules`.

`substitutions` maps a species index of `poly`'s own rules to the `Polyform` that species is
replaced by. It has to be a polyform of `rules`. A [`MetaParticleSpecies`](@ref) is substituted by the polyform it
wraps and needs no entry. Any other species type not in `substitutions`is relplaced by the species of `rules`
with the same species index.
"""
function recast(poly::Polyform{D}, rules::BindingRules; substitutions=Dict()) where {D}
    src = bindingrules(poly)
    subs = _resolvesubstitutions(src, rules, substitutions)

    parts = _substitute_particles(poly, subs)
    contacts = Contact[]
    g = NautyDiGraph(0)
    for (i, sp) in enumerate(parts)
        placed = view(parts, 1:(i - 1))
        overlap, cts = _overlap_and_contacts(placed, sp, rules)
        overlap && _recastfailed(placed, sp, rules)
        append!(contacts, cts)
        blockdiag!(g, graphrep(species(rules, speciesindex(sp))))
    end

    for contact in contacts
        for (v1, v2) in contact_pairing(contact)
            add_edge!(g, v1, v2)
            add_edge!(g, v2, v1)
        end
    end

    # `g` is built in original vertex order, so the canonical permutation is `canon2orig` itself.
    perm, autg = nauty(g; canonize=true)
    cvs = collect(Int, perm)
    return Polyform{D,particletype(rules),typeof(rules),typeof(g)}(
        g, convert(Int, autg.n), cvs, invperm(cvs), parts, rules
    )
end

# `_overlap_and_contacts` refuses for three different reasons and reports all of them the same
# way, so ask it again to find out which. Only ever reached on the way to an error.
function _recastfailed(parts, part, rules)
    refuses(; kwargs...) = first(_overlap_and_contacts(parts, part, rules; kwargs...))

    refuses(; allow_noninteracting=true, allow_misaligned=true) &&
        throw(ArgumentError("two of the particles overlap, the recast polyform is invalid."))
    why = if refuses(; allow_noninteracting=true)
        "at a twist `rules` does not allow"
    else
        "at a pair of sites `rules` leaves inert"
    end
    return throw(
        ArgumentError("Two of the particles touch $why, so recasting does not result in a polyform valid under `rules`")
    )
end

function _overlap_and_contacts(
    polyparticles::AbstractVector{<:Particle},
    part::Particle,
    rules::BindingRules;
    allow_noninteracting=false,
    allow_misaligned=false,
    kwargs...,
)
    intmat = interactionmatrix(rules)
    contacts = Contact[]

    for polypart in polyparticles
        could_contact(polypart, part, rules) || continue
        overlap(polypart, part, rules) && return true, nothing

        for b1 in bindingsites(polypart, rules), b2 in bindingsites(part, rules)
            istouching(b1, b2) || continue
            interacting = intmat[color(b1), color(b2)]
            !allow_noninteracting && !interacting && return true, nothing
            twst = twist(b1, b2)
            !allow_misaligned && isnothing(twst) && return true, nothing
            push!(contacts, Contact(b1.vertices, b2.vertices, something(twst, 0), twistfreedom(b1, b2)))
        end
    end
    return false, contacts
end

function overlap_and_contacts(poly::Polyform, part::Particle; kwargs...)
    return _overlap_and_contacts(poly.particles, part, bindingrules(poly); kwargs...)
end

# Walk the unbound binding sites of `poly` in canonical order, yielding each one's
# `(particle, site)` location together with the site itself. The four accessors below project it.
function _exposed(poly::AbstractPolyform)
    rules = bindingrules(poly)
    index = Dict(leadingvertex(p) => i for (i, p) in enumerate(poly.particles))
    out = Tuple{ParticleSite,sitetype(rules)}[]
    for orig_v in poly.canon2orig
        part = particle_from_leadingvertex(poly, orig_v)
        isnothing(part) && continue
        for k in 1:nsites(part, rules)
            s = bindingsite(part, rules, k)
            _isbound_vertex(poly, part, first(s.vertices); canonidxs=false) && continue
            push!(out, (ParticleSite(index[leadingvertex(part)], k), s))
        end
    end
    return out
end

"""
    exposedsites(poly::AbstractPolyform)

The [`ParticleSite`](@ref) of every *unbound* binding site of `poly`, in canonical order.

Bound sites are consumed by the bonds holding `poly` together and are never listed. Sites whose
color takes part in no rule are listed, even though nothing can attach through them as `poly`
stands; [`opensites`](@ref) is this list without them.

Addresses rather than the sites themselves, since `bindingsite(poly, loc)` gets the site from
the address but nothing gets the address back from a site.

These are the sites a [`MetaParticleSpecies`](@ref) may expose. It exposes the open ones by
default, and an inert one becomes usable simply by being named and given a live color.
"""
exposedsites(poly::AbstractPolyform) = [l for (l, _) in _exposed(poly)]

"""
    opensites(poly::AbstractPolyform)

The [`ParticleSite`](@ref) of every binding site of `poly` a partner can still attach through:
the unbound ones whose color some rule uses, in canonical order.

See [`exposedsites`](@ref), which lists the inert ones too.
"""
function opensites(poly::AbstractPolyform)
    rules = bindingrules(poly)
    return [l for (l, s) in _exposed(poly) if !isinert(rules, color(s))]
end

"""
    _deletable_species(poly; target, dist, queue)

Return `(top, runnerup, top_lv)`: the two largest species indices among the particles `lower!`
could delete from `poly`, and the top leading vertex of the particle achieving `top`.

`lower!` deletes from the highest label class holding a removable particle. Canonical position
runs with vertex label (nauty orders classes by color, `vertexlabels2labptn` by label), and
`_adjust_labels_and_colors` puts species `s`'s labels above species `s-1`'s, so species index
order is label order. Removability is read off `poly` rather than the child: attaching a
particle cannot disconnect what removing a different one leaves behind.

The runner-up lets the anchor particle exclude itself in O(1); see
[`collect_attachments!`](@ref).
"""
function _deletable_species(poly::Polyform; target, dist, queue)
    rules = bindingrules(poly)
    n = nparticles(poly)
    top, runnerup, top_lv = 0, 0, 0
    for part in poly.particles
        n > 1 &&
            is_cutset(graphrep(poly), @view(poly.orig2canon[graphvertices(part, rules)]); target, dist, queue) &&
            continue
        s = speciesindex(part)
        if s > top
            top, runnerup, top_lv = s, top, leadingvertex(part)
        elseif s > runnerup
            runnerup = s
        end
    end
    return top, runnerup, top_lv
end

"""
    collect_attachments!(attachments, poly::Polyform)

Fill `attachments` with every way of growing `poly` by one particle, as triples
`(site, loc, t)`: an open binding site of `poly`, the species and site index of the
particle to attach, and which of the bond's distinct twists to attach it in.

Two filters keep the list down to what reverse search can use. Only one partner site per
symmetry orbit is kept, since the others give the same child. And a candidate is dropped when
`lower!` would not undo it, i.e. when the child holds a removable particle from a species
index above the one being attached: the child's parent would then be a different polyform and
reverse search would reject the pair anyway.

The second filter is conservative about the anchor particle, the one carrying `site`. The
anchor is excluded from the removable set, because a new particle bonded to it alone leaves it
a cut vertex of the child. When the attachment also closes a ring the anchor does stay
removable, and excluding it then only lets a few extra candidates through, which reverse search
rejects.
"""
function collect_attachments!(attachments, poly::Polyform)
    rules = bindingrules(poly)
    empty!(attachments)
    nv_g = nv(graphrep(poly))
    # `dist` carries distances and the -1/-2 markers `is_cutset` needs, so it cannot be a Bool
    target, dist, queue = zeros(Bool, nv_g), zeros(Int, nv_g), zeros(Int, nv_g)

    top, runnerup, top_lv = _deletable_species(poly; target, dist, queue)
    for orig_v in poly.canon2orig
        anchor = particle_from_leadingvertex(poly, orig_v)
        isnothing(anchor) && continue

        # The anchor is excluded, so it takes the runner-up when it is itself the top scorer.
        deletable = leadingvertex(anchor) == top_lv ? runnerup : top
        for k in 1:nsites(anchor, rules)
            site = bindingsite(anchor, rules, k)
            _isbound_vertex(poly, anchor, first(site.vertices); canonidxs=false) && continue
            isinert(rules, color(site)) && continue

            for loc in distinct_attachments(rules, color(site))
                # The new particle is always removable, so `lower!` stops at its label class
                # unless a higher one survives.
                loc.species < deletable && continue
                mate = bindingsite(rules, loc)
                for t in 0:(_ndistincttwists(site, mate) - 1)
                    push!(attachments, (site, loc, t))
                end
            end
        end
    end
    return attachments
end

"""
    collect_attachments(poly::Polyform)

Return every way of growing `poly` by one particle, allocating the vector.
See [`collect_attachments!`](@ref).
"""
function collect_attachments(poly::Polyform)
    rules = bindingrules(poly)
    BS = sitetype(rules)
    return collect_attachments!(Tuple{BS,SpeciesSite,Int}[], poly)
end
