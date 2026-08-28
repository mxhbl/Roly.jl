"""
    MetaParticleSpecies{D,F,B,G,PF}

A `Polyform` wrapped as a `ParticleSpecies`, so that meta-assemblies can be built out of whole
polyforms.

Selected exposed sites of the polyform become the species' sites (by default are all non-inert sites are forwarded).
"""
struct MetaParticleSpecies{D,F,B<:BindingSite,G<:AbstractNautyGraph,PF} <: ParticleSpecies{D,B}
    g::G
    sites::Vector{B}
    poly::PF
    rmax::F
    skin::F
end

"""
    MetaParticleSpecies(poly::Polyform; colors=nothing)
    MetaParticleSpecies(poly::Polyform, sites; colors=nothing)

Wrap `poly` as a particle species whose sites are the ones named by `sites`, given as `(particle, site)`
pairs and exposed in that order. Without them every open non-inert site is exposed, in
[`exposedsites`](@ref) order.

  - `colors`: one interaction color per exposed site; by default each keeps the color it has
    inside `poly`
  - `exposeinert`: whether to also expose previously inert sites to the new binding rules

Each site's `sitesym` and `locking` carry over from the polyform, while its `stab` is recomputed
against the polyform by [`siteorbits`](@ref) and [`stabilizerorders`](@ref).
"""
function MetaParticleSpecies(poly::Polyform; colors=nothing, exposeinert::Bool=false)
    return MetaParticleSpecies(poly, exposeinert ? exposedsites(poly) : opensites(poly); colors)
end

function MetaParticleSpecies(poly::Polyform{D}, sites; colors=nothing) where {D}
    rules = bindingrules(poly)
    picks = [(Int(p), Int(k)) for (p, k) in sites]  # materialize and narrow to Int
    n = length(picks)
    n > 0 || throw(ArgumentError("a meta-species needs at least one exposed site"))
    allunique(picks) || throw(ArgumentError("`sites` names the same site twice"))

    exposable = Set(exposedsites(poly))
    open = map(picks) do (p, k)
        1 <= p <= nparticles(poly) ||
            throw(ArgumentError("`poly` has $(nparticles(poly)) particles, so ($p, $k) is out of range"))
        part = poly.particles[p]
        1 <= k <= nsites(part, rules) ||
            throw(ArgumentError("particle $p has $(nsites(part, rules)) sites, so ($p, $k) is out of range"))
        ParticleSite(p, k) in exposable || throw(
            ArgumentError(
                "site ($p, $k) is bound inside the polyform; a bond already consumes it, so it cannot be exposed"
            ),
        )
        return bindingsite(part, rules, k)
    end

    cols = colors === nothing ? [color(s) for s in open] : collect(Int, colors)
    length(cols) == n || throw(DimensionMismatch("$n sites were exposed but $(length(cols)) colors were given"))

    # The polyform's graph in its original vertex order, where every particle owns a contiguous
    # block of vertices and each site keeps the vertex range it has inside the polyform.
    g = graphrep(poly)[poly.orig2canon]

    # `sitesym` and `locking` belong to the site and carry over. The orbits and the stabilizers
    # are the ordinary ones, measured against the polyform's own symmetry group rather than one
    # derived from just the site poses, which would be unaware of shape of the polyform.
    group = rotationgroup(poly)
    c = rotationcenter(poly)
    poses, sitesyms = [s.pose + (-c) for s in open], [s.sitesym for s in open]
    orbits = siteorbits(poses, sitesyms, cols; group)
    stabs = stabilizerorders(poses, sitesyms, cols; group)
    metasites = map(1:n) do i
        return BindingSite(
            open[i].pose,
            cols[i],
            open[i].vertices,
            open[i].touching_tolerance,
            open[i].alignment_tolerance,
            open[i].sitesym,
            stabs[i],
            open[i].locking,
        )
    end
    _labelsites!(g, metasites, orbits)

    F = numtype(poly)
    rmax = maximum(poly.particles; init=zero(F)) do part
        norm(part.pose.x) + bounding_radius(species(rules, speciesindex(part)))
    end

    return check_encoding(
        MetaParticleSpecies{D,F,eltype(metasites),typeof(g),typeof(poly)}(
            g, metasites, copy(poly), convert(F, rmax), sqrt(eps(F)) * convert(F, rmax)
        ),
    )
end

# candidate rotation symmetries, which is the symmetry of the underlying polyform.
# Exposing only some of the polyform's open sites or recoloring some sites will lower the symmetry
_rotationcandidates(ps::MetaParticleSpecies) = rotationgroup(polyform(ps))

function _sitegeometry(ps::MetaParticleSpecies)
    c = rotationcenter(polyform(ps))
    return ([s.pose + (-c) for s in ps.sites], [s.sitesym for s in ps.sites], [sitelabel(ps, i) for i in 1:nsites(ps)])
end

# Label every open-site vertex by its symmetry orbit, placed above every interior label.
function _labelsites!(g, sites, orbits)
    sitevertices = Set(v for s in sites for v in s.vertices)
    base = maximum((label(g, v) for v in 1:nv(g) if v ∉ sitevertices); init=0)
    for (i, s) in enumerate(sites)
        for v in s.vertices
            setlabel!(g, v, base + orbits[i])
        end
    end
    return g
end

function Base.show(io::Core.IO, ps::MetaParticleSpecies)
    print(io, "$(dimension(ps))d MetaParticleSpecies[", nparticles(ps.poly), " particles, ", nsites(ps), " sites]")
end

Base.copy(ps::MetaParticleSpecies) = typeof(ps)(copy(ps.g), copy(ps.sites), copy(ps.poly), ps.rmax, ps.skin)

graphrep(ps::MetaParticleSpecies) = ps.g
nsites(ps::MetaParticleSpecies) = length(ps.sites)
bindingsite(ps::MetaParticleSpecies, i::Integer) = ps.sites[i]
bounding_radius(ps::MetaParticleSpecies) = ps.rmax

"""
    polyform(ps::MetaParticleSpecies)

The `Polyform` that `ps` wraps.
"""
polyform(ps::MetaParticleSpecies) = ps.poly

# The generic `setcolors!` re-derives labels and stabilizers from the site poses, which do not
# describe a polyform. Recoloring here keeps the polyform's interior labeling and rewrites only the
# open-site labels and stabilizers, both against the polyform's own symmetry group.
function setcolors!(ps::MetaParticleSpecies, colors::AbstractVector{<:Integer})
    n = nsites(ps)
    length(colors) == n || throw(ArgumentError("expected $n colors, got $(length(colors))"))
    poses, sitesyms, _ = _sitegeometry(ps)
    group = rotationgroup(ps)
    orbits = siteorbits(poses, sitesyms, colors; group)
    stabs = stabilizerorders(poses, sitesyms, colors; group)
    for i in eachindex(ps.sites)
        ps.sites[i] = setstab(setcolor(ps.sites[i], colors[i]), stabs[i])
    end
    _labelsites!(ps.g, ps.sites, orbits)
    return nothing
end

function could_contact(p1::SpeciesAndPose{<:MetaParticleSpecies}, p2::SpeciesAndPose{<:MetaParticleSpecies}; kwargs...)
    (s1, pose1), (s2, pose2) = p1, p2
    return norm(pose1.x - pose2.x) < bounding_radius(s1) + bounding_radius(s2) + s1.skin + s2.skin
end

"""
    overlap(p1::SpeciesAndPose{<:MetaParticleSpecies}, p2::SpeciesAndPose{<:MetaParticleSpecies})

Whether two posed meta-particles overlap, i.e. whether any of their constituent particles do.
"""
function overlap(p1::SpeciesAndPose{<:MetaParticleSpecies}, p2::SpeciesAndPose{<:MetaParticleSpecies}; kwargs...)
    (s1, pose1), (s2, pose2) = p1, p2
    rules1, rules2 = bindingrules(s1.poly), bindingrules(s2.poly)
    for a in s1.poly.particles
        pa = species(rules1, speciesindex(a)) => pose1 * a.pose
        for b in s2.poly.particles
            pb = species(rules2, speciesindex(b)) => pose2 * b.pose
            could_contact(pa, pb) || continue
            overlap(pa, pb; kwargs...) && return true
        end
    end
    return false
end

"""
    BindingRules(ps::MetaParticleSpecies)

Lift the interactions of the wrapped polyform to the meta-species itself: the rules under which
two meta-particles bond exactly where the sites they expose would have bonded as ordinary sites.

Only meaningful while the meta-species still carries the colors it inherited from its polyform,
which is the default.
"""
BindingRules(ps::MetaParticleSpecies) = BindingRules([ps])

"""
    BindingRules(pss::AbstractVector{<:MetaParticleSpecies})

Lift the interactions of the wrapped polyforms to the meta-species themselves: the rules under
which two meta-particles bond exactly where the sites they expose would have bonded as ordinary
sites, across the species as well as within each.

Only works if the meta-species still carry the colors they inherited from their
polyforms, which is the default. All of them must wrap polyforms built under the same rules.
"""
function BindingRules(pss::AbstractVector{<:MetaParticleSpecies})
    isempty(pss) && throw(ArgumentError("The list of particle species is empty."))
    base = bindingrules(polyform(first(pss)))
    all(ps -> bindingrules(polyform(ps)) === base, pss) ||
        throw(ArgumentError("The meta particle species wrap polyforms built under different binding rules."))

    intmat = interactionmatrix(base)
    cols = [[color(bindingsite(ps, i)) for i in 1:nsites(ps)] for ps in pss]
    all(c -> 1 <= c <= ncolors(base), Iterators.flatten(cols)) ||
        throw(ArgumentError("These meta particle species carry colors the original polyforms' rules do not have."))

    # Each unordered pair of sites once: within a species from `i` on, across species all of them.
    rows = [
        [a, i, b, j] for a in eachindex(pss) for i in eachindex(cols[a]) for b in a:length(pss) for
        j in (a == b ? i : 1):length(cols[b]) if intmat[cols[a][i], cols[b][j]]
    ]
    bonds = isempty(rows) ? zeros(Int, 0, 4) : permutedims(reduce(hcat, rows))
    return BindingRules(bonds, pss)
end

"""
    MetaBindingRules{D}
    MetaPolyform{D,P}

Rules whose species are all [`MetaParticleSpecies`](@ref), and a `Polyform` built under such
rules -- the two things the meta level is made of, named so that dispatch can say so. Rules
mixing meta-species with plain ones are neither, since their element type is then the abstract
supertype, which is what makes a meta-particle bonded to a plain one a `MethodError`.
"""
const MetaBindingRules{D} = BindingRules{D,<:MetaParticleSpecies}
const MetaPolyform{D,P} = Polyform{D,P,<:MetaBindingRules{D}}

"""
    originalrules(ps::MetaParticleSpecies)
    originalrules(rules::MetaBindingRules)
    originalrules(poly::MetaPolyform)

Return the rules the meta level was lifted from: the ones the wrapped polyforms are built under.

Distinct from [`inducedrules`](@ref), which projects the meta-rules down and may add bonds the
original rules did not allow. These are the original rules themselves, unchanged.
"""
originalrules(ps::MetaParticleSpecies) = bindingrules(polyform(ps))

"""
    originalrules(rules::MetaBindingRules)

Return the rules every species of `rules` wraps a polyform of, or throw if they disagree.
"""
function originalrules(rules::MetaBindingRules)
    spcs = species(rules)
    out = originalrules(first(spcs))
    all(ps -> originalrules(ps) === out, spcs) ||
        throw(ArgumentError("The meta-species wrap polyforms built under different binding rules."))
    return out
end

"""
    originalrules(poly::MetaPolyform)

Return the rules the meta-particles of `poly` were lifted from.
"""
originalrules(poly::MetaPolyform) = originalrules(bindingrules(poly))

# The color each site of `ps` carries inside the polyform it was taken from. A site keeps the
# vertex range it has there whatever it is recolored to, so any one of its vertices names it.
function _underlyingcolors(ps::MetaParticleSpecies)
    poly = polyform(ps)
    return map(1:nsites(ps)) do i
        # look up the original polyform binding site through the graph vertex of the species' site
        v = first(bindingsite(ps, i).vertices)
        return color(bindingsite(poly, _vertex_to_particle_site(poly, v; canonidxs=false)))
    end
end

"""
    inducedrules(rules::BindingRules)

"Project" the meta-rules `rules` to down to the original binding rules of the underlying particles.
Returns the binding rules of the original particles with additional bonds that make all contacts in the meta-polyform
valid.
"""
function inducedrules(rules::MetaBindingRules)
    spcs = species(rules)
    origrules = originalrules(rules)

    intmat = Matrix(interactionmatrix(origrules))
    cols = map(_underlyingcolors, spcs)
    # add bonds to the original binding rules
    for (sites1, sites2) in bonded_sites(rules)
        for l1 in sites1, l2 in sites2
            c1, c2 = cols[l1.species][l1.site], cols[l2.species][l2.site]
            intmat[c1, c2] = intmat[c2, c1] = true
        end
    end
    return BindingRules(intmat, species(origrules))
end

# Expand the meta-species; return all underyling particles of the underlying binding rules
function _underlying_particles(meta::MetaPolyform)
    subs = [polyform(ps) for ps in species(bindingrules(meta))]
    return _substitute_particles(meta, subs)
end

# A meta-species substitutes the polyform it wraps
function _defaultsubstitution(ps::MetaParticleSpecies, ::BindingRules, rules::BindingRules, ::Integer)
    return _checkspecies(polyform(ps), rules)
end
