"""
    Tiling{D}

A periodic structure: particles bonded to each other and to their own translates under a lattice
of rank at most `D`. Returned by [`tilings`](@ref).

An [`AbstractPolyform`](@ref) like [`Polyform`](@ref), and read the same way -- [`bonds`](@ref),
[`bondtypes`](@ref), [`composition`](@ref) and [`opensites`](@ref) all speak about one cell of the
infinite structure. The two differ in what an edge of the graph means: in a polyform it joins two
particles that sit side by side, in a tiling it may instead join a particle to a translate of one,
so a polyform is the tiling whose lattice is trivial.

The particles are one cell's worth, and which cell that is carries no meaning: cutting the same
structure along different lattice translations gives the same graph, and two tilings are equal
exactly when their graphs are. [`latticevectors`](@ref) is a basis of the lattice the graph
implies, kept because a basis is what its readers want and the search already has one, and left
out of the identity for the same reason the cell is.
"""
struct Tiling{D,P<:Particle,S<:BindingRules,G<:AbstractNautyGraph,V} <: AbstractPolyform{D}
    graphrep::G
    sigma::Int
    canon2orig::Vector{Int}
    orig2canon::Vector{Int}
    particles::Vector{P}
    bindingrules::S
    vectors::Vector{V}
end

"""
    Tiling(cell::Polyform, periodic)

Assemble the tiling whose cell is `cell` and whose translates close the bonds `periodic`, a
[`Contact`](@ref) per bond in `cell`'s own vertex numbering.

Neither the cell nor the lattice is taken as given. Both are read back off the graph once it is
canonically labeled: where the structure is cut, which frame it sits in, and which translations
generate it all follow from that labeling and not from how the search arrived.
"""
function Tiling(cell::Polyform{D}, periodic) where {D}
    rules = bindingrules(cell)
    g = NautyDiGraph(0)
    for part in cell.particles
        blockdiag!(g, graphrep(species(rules, speciesindex(part))))
    end

    # the cell's own bonds, read back into original vertex order, and then the periodic ones,
    # which arrive in that order already. Nothing marks which is which: that a bond crosses the
    # cut is a fact about the cut, not about the structure
    pairs = [(toorig(cell, e.src), toorig(cell, e.dst)) for e in exterior_edges(cell)]
    for contact in periodic, (v1, v2) in contact_pairing(contact)
        push!(pairs, (v1, v2))
    end
    for (u, v) in pairs
        _addmarker!(g, rules, u, v)
    end

    perm, autg = nauty(g; canonize=true)
    cvs = collect(Int, perm)
    P, S, G, V = particletype(rules), typeof(rules), typeof(g), SVector{D,numtype(rules)}
    t = Tiling{D,P,S,G,V}(
        g, convert(Int, autg.n), cvs, invperm(cvs), copy(cell.particles), rules, V[]
    )
    _canonicalcut!(t)
    _canonicalframe!(t)
    append!(t.vectors, _latticebasis(t))
    return t
end

# Put the cell in the frame of its first particle: that particle at the origin, unrotated. The
# canonical labeling names the same particle in any two equal tilings, so this leaves both in one
# frame rather than in frames related by a rigid motion, and their poses, cells and lattice
# vectors come out equal outright. A pose carries a full frame on its own, so one particle is
# enough in any dimension, with no case to make for cells too small to offer a second.
function _canonicalframe!(t::Tiling)
    parts = t.particles
    q = inv(parts[_rootparticle(t)].pose)
    for i in eachindex(parts)
        parts[i] = q * parts[i]
    end
    return t
end

# Re-choose which translate of each particle the cell holds, so that where the structure is cut
# follows from the canonical graph rather than from the cell the search happened to grow.
#
# Walking the canonical spanning tree and pulling each particle onto the bond that reaches it
# leaves every tree bond meeting inside the cell. Which bonds stay periodic is then a function of
# the graph alone, and so is the cell `unitcell` reads back off the poses. Two equal tilings that
# arrived cut differently come out cut the same, up to the lattice translation and rigid motion
# that carry one onto the other and that no cell can distinguish anyway.
function _canonicalcut!(t::Tiling)
    parts = t.particles
    placed = falses(length(parts))
    placed[_rootparticle(t)] = true
    while true
        step = _reachbond(t, placed)
        isnothing(step) && break
        j, delta = step
        parts[j] = parts[j] + delta
        placed[j] = true
    end
    return t
end

# The particle holding the first vertex in canonical order, which is where the walk starts.
function _rootparticle(t::Tiling)
    for v in vertices(graphrep(t))
        loc = _vertex_to_particle_site(t, v; canonidxs=true)
        isnothing(loc) || return loc.particle
    end
    return 1
end

# The first bond in canonical order that reaches an unplaced particle from a placed one, as that
# particle and the translation laying it against its partner. `nothing` once all are placed.
function _reachbond(t::Tiling, placed)
    for (_, (u, v)) in _markers(t)
        a = _vertex_to_particle_site(t, u; canonidxs=false)
        b = _vertex_to_particle_site(t, v; canonidxs=false)
        (isnothing(a) || isnothing(b)) && continue
        placed[a.particle] == placed[b.particle] && continue
        sa, sb = bindingsite(t, a), bindingsite(t, b)
        placed[a.particle] && return b.particle, sa.pose.x - sb.pose.x
        return a.particle, sb.pose.x - sa.pose.x
    end
    return nothing
end

"""
    _markerlabel(rules::BindingRules)

A graph vertex label no species of `rules` uses, worn by the vertices that stand in for bonds.
"""
_markerlabel(rules::BindingRules) =
    maximum(maximum(labels(graphrep(species(rules, i)))) for i in 1:nspecies(rules)) + 1

# Record the bond between graph vertices `u` and `v` as a vertex of its own joined to both, rather
# than as the edge `u -- v`.
#
# A particle of a tiling bonds to its own translate, so a bond can join two sites of one particle
# -- and the species may already carry an edge between exactly those two vertices, as a cycle
# encoding does for the cube's opposite faces. Writing the bond as an edge would then write
# nothing, and the tiling would be indistinguishable from one without the bond. A fresh vertex
# collides with nothing. Every bond is marked, the cell's own alongside the periodic ones, so
# what the marks record is that a bond exists and not where the structure was cut.
function _addmarker!(g, rules::BindingRules, u::Integer, v::Integer)
    add_vertices!(g, 1; vertex_labels=[_markerlabel(rules)])
    m = nv(g)
    add_edge!(g, u, m)
    add_edge!(g, m, u)
    add_edge!(g, v, m)
    add_edge!(g, m, v)
    return g
end

# The marker vertices of `t`, in canonical order, each with the two site vertices it joins, in
# original numbering. Markers are added after every particle's block, so an original index beyond
# the last block is a marker and no particle owns it.
function _markers(t::Tiling)
    g = graphrep(t)
    ml = _markerlabel(bindingrules(t))
    ls = labels(g)
    return ((m, _markerends(t, m)) for m in vertices(g) if ls[m] == ml)
end

function _markerends(t::Tiling, m::Integer)
    neighs = NautyGraphs.adjrow(graphrep(t), m)
    ends = (0, 0)
    for w in eachindex(neighs)
        neighs[w] || continue
        ends = first(ends) == 0 ? (toorig(t, w), 0) : (first(ends), toorig(t, w))
    end
    return ends
end

"""
    bonds(t::Tiling)

Return a lazy iterator of the bonds of one cell of `t`, each one once, as [`ParticleSite`](@ref)
pairs. The bonds a translate closes are listed alongside the ones inside the cell.
"""
function bonds(t::Tiling)
    seen = Set{NTuple{2,ParticleSite}}()
    out = Pair{ParticleSite,ParticleSite}[]
    for (_, (u, v)) in _markers(t)
        key = minmax(_vertex_to_particle_site(t, u; canonidxs=false),
                     _vertex_to_particle_site(t, v; canonidxs=false))
        key in seen && continue
        push!(seen, key)
        push!(out, first(key) => last(key))
    end
    return out
end

# A site of a tiling is bound when a marker hangs off it, wherever the bond leads: to another
# particle, or to a translate of its own.
function _isbound_vertex(t::Tiling, ::Particle, v::Integer; canonidxs::Bool)
    g = graphrep(t)
    ml = _markerlabel(bindingrules(t))
    ls = labels(g)
    neighs = NautyGraphs.adjrow(g, canonidxs ? v : tocanon(t, v))
    for w in eachindex(neighs)
        neighs[w] && ls[w] == ml && return true
    end
    return false
end

function Base.show(io::Core.IO, t::Tiling{D}) where {D}
    return print(io, "Tiling{$D}[n=$(nparticles(t)), d=$(length(latticevectors(t)))]")
end

Base.:(==)(a::Tiling, b::Tiling) = bindingrules(a) === bindingrules(b) && graphrep(a) == graphrep(b)
Base.hash(t::Tiling, h::UInt) = hash(graphrep(t), h)

"""
    latticevectors(t::Tiling)

The shortest basis of the lattice `t` repeats under, at most one vector per dimension. Fewer means
a partial closure: a column in the plane, a column or a sheet in space.

Read off `t`'s own periodic bonds rather than off the search that found it, in the frame `t` is
canonically placed in, so two tilings that are equal carry equal vectors.
"""
latticevectors(t::Tiling) = t.vectors

"""
    iscomplete(t::Tiling)

Whether `t` leaves no open site: every site of every particle is closed, by a neighbor in the
cell or by a translate. A tiling that is not complete closes along fewer directions than the
space has, i.e. a 1d line in the 2d plane.
"""
iscomplete(t::Tiling) = isempty(opensites(t))

"""
    tilingorder(t::Tiling)

How many particles one cell of `t` holds.
"""
tilingorder(t::Tiling) = nparticles(t)

"""
    unitcell(t::Tiling)

One cell of `t` as a [`Polyform`](@ref): its particles, keeping only the bonds that do not cross
into a translate.

Which bonds those are depends on where the structure is cut, and `t` is cut canonically, so equal
tilings give equal cells. It is still one cell among the several that describe `t` equally well,
chosen rather than distinguished.
"""
function unitcell(t::Tiling{D}) where {D}
    rules = bindingrules(t)
    g = NautyDiGraph(0)
    for part in t.particles
        blockdiag!(g, graphrep(species(rules, speciesindex(part))))
    end
    for (_, (u, v)) in _markers(t)
        _isperiodic(t, u, v; canonidxs=false) && continue
        add_edge!(g, u, v)
        add_edge!(g, v, u)
    end

    perm, autg = nauty(g; canonize=true)
    cvs = collect(Int, perm)
    P, S, G = particletype(rules), typeof(rules), typeof(g)
    return Polyform{D,P,S,G}(g, convert(Int, autg.n), cvs, invperm(cvs), copy(t.particles), rules)
end

# Whether the bond between two graph vertices crosses into a translate, which is to say the two
# sites it joins do not meet inside the cell. `_translation` is then the one that carries the
# second onto the first, and the lattice is what those generate.
#
# Position alone decides it: two sites that mate are antiparallel rather than equally posed, so
# their poses never agree, and `_siteoverlap` would call every bond periodic.
function _isperiodic(t::Tiling, src::Integer, dst::Integer; canonidxs::Bool=true)
    a = bindingsite(t, _vertex_to_particle_site(t, src; canonidxs))
    b = bindingsite(t, _vertex_to_particle_site(t, dst; canonidxs))
    return !isapprox(a.pose.x, b.pose.x; atol=a.touching_tolerance + b.touching_tolerance, rtol=0)
end

function _translation(t::Tiling, src::Integer, dst::Integer; canonidxs::Bool=true)
    a = bindingsite(t, _vertex_to_particle_site(t, src; canonidxs))
    b = bindingsite(t, _vertex_to_particle_site(t, dst; canonidxs))
    return a.pose.x - b.pose.x
end

# The translations the periodic bonds name, in canonical order, each once. A bond names the same
# translation read either way round, so they are oriented; two bonds can still name one and the
# same, which is what a cell touching one neighbor through two pairs of sites does.
function _translations(t::Tiling)
    out = eltype(t.vectors)[]
    for (_, (u, v)) in _markers(t)
        _isperiodic(t, u, v; canonidxs=false) || continue
        w = _orient(_translation(t, u, v; canonidxs=false))
        any(x -> x ≈ w, out) || push!(out, w)
    end
    return out
end

"""
    _latticebasis(t::Tiling)

The shortest basis of the lattice `t` repeats under.

The translations its periodic bonds name generate that lattice: cutting along a spanning tree
leaves the remaining bonds a free basis of the graph's cycles, which the translations carry onto
the whole group. There can be more of them than the lattice has rank, though, since a cell may
touch at `v`, at `w` and again at `v + w`, so a subset is taken -- and it has to be one that
generates the rest, not merely an independent one, since independent vectors can span a coarser
sublattice. Every subset is tried, so a basis is found whenever one exists.

That basis is then shortened, since the shortest vectors of the lattice need not carry bonds at
all. Both steps are settled by canonical order where lengths tie, so the result is a function of
the graph.
"""
function _latticebasis(t::Tiling)
    vs = _translations(t)
    isempty(vs) && return vs
    r = rank(reduce(hcat, vs))
    for idx in _indexsubsets(length(vs), r)
        _generates(vs[idx], vs) && return _shortestbasis(vs[idx], vs)
    end
    return error(
        "Internal error: a tiling's periodic bonds name no basis of its lattice. Please file an issue."
    )
end

# The shortest basis of the lattice `basis` spans. Every lattice vector no longer than the longest
# of `basis` is enumerated -- the rows of the pseudoinverse bound how far the coefficients can
# reach -- and the shortest independent ones are taken, which is a basis for rank at most three.
# Ties in length go to a bonded translation, in canonical order, before anything else.
function _shortestbasis(basis, preferred)
    r = length(basis)
    B = reduce(hcat, basis)
    reach = maximum(norm, basis)
    P = pinv(B)
    bound = [max(1, ceil(Int, reach * norm(view(P, i, :)))) for i in 1:r]

    cands = eltype(basis)[]
    for c in Iterators.product(((-bound[i]):bound[i] for i in 1:r)...)
        all(iszero, c) && continue
        v = _orient(sum(c[i] * basis[i] for i in 1:r))
        norm(v) < _tol(v) && continue
        any(x -> x ≈ v, cands) || push!(cands, v)
    end
    sort!(cands; by=v -> _shortestfirst(v, preferred))

    out = eltype(basis)[]
    for v in cands
        push!(out, v)
        rank(reduce(hcat, out)) == length(out) || pop!(out)
        length(out) == r && break
    end
    _generates(out, basis) || error(
        "Internal error: the shortest vectors of a tiling's lattice are not a basis of it. Please file an issue."
    )
    return out
end

# Shortest first; among equally short ones a bonded translation before an unbonded lattice vector,
# and earlier bonds before later. Coordinates settle what is left, which only arises between
# vectors the structure gives no reason to tell apart.
function _shortestfirst(v, preferred)
    i = findfirst(w -> w ≈ v, preferred)
    return (round(norm(v); digits=7), isnothing(i) ? length(preferred) + 1 : i,
            Tuple(round.(v; digits=7)))
end

# Whether every vector of `vs` is an integer combination of the independent vectors `basis`.
function _generates(basis, vs)
    B = reduce(hcat, basis)
    rank(B) == length(basis) || return false
    return all(vs) do v
        c = B \ v
        norm(B * c - v) < _tol(v) || return false
        return all(x -> abs(x - round(x)) < _tol(v), c)
    end
end

# Every `r`-element subset of `1:n`, as index vectors, in increasing order.
function _indexsubsets(n::Integer, r::Integer)
    out = Vector{Int}[]
    r == 0 && return push!(out, Int[])
    grow(prefix, from) =
        for i in from:n
            length(prefix) + 1 == r ? push!(out, [prefix; i]) : grow([prefix; i], i + 1)
        end
    grow(Int[], 1)
    return out
end

"""
    tilingenum(f, poly::Polyform; nreps=2, maxorder=1)

Enumerate the periodic closures of `poly`, streaming each one to `f` as a [`Tiling`](@ref).

`f(t)` returns one of three signals:

  - `ACCEPT` (or `true`): enumeration continues as normal.
  - `REJECT` (or `false`): no further vector is added to `t`'s lattice.
  - `BREAK`: the enumeration terminates immediately.

Each tiling is visited exactly once: a cell can close the same way along several choices of
vectors, and the repeats never reach `f`. Keyword arguments are as in [`tilings`](@ref).
"""
function tilingenum(f::F, poly::Polyform; nreps::Integer=2, maxorder::Integer=1) where {F}
    isempty(opensites(poly)) && return
    # `exposeinert=true` to catch overlaps of inert sites
    metarules = BindingRules(MetaParticleSpecies(poly; exposeinert=true))

    # keep track of duplicates and only call f on newly found tilings
    # still need to descent to offspring of seen tilings
    seen = Set{_tilingtype(poly)}()
    once(t) = t in seen ? ACCEPT : (push!(seen, t); f(t))

    polyenum(metarules; maxsize=maxorder) do cell, _
        return _celltilings(once, cell, nreps)
    end
    return
end

"""
    tilings(poly::Polyform; nreps=2, maxorder=1)

Return the periodic closures of `poly`: choices of up to `dimension` translation vectors under
which copies of a cell form valid bonds with each other and never overlap, each vector
contributing at least one bond.

  - `nreps`: neighbor shells placed and checked per vector.
  - `maxorder`: let a cell hold up to `maxorder` copies of `poly`, grown as a meta-polyform and
    then checked for translation tilings. Above 1 this reaches tilings in which `poly` appears in
    rotated configuraitons.
  - returns a vector of [`Tiling`](@ref)

Candidate vectors connect interacting, aligned pairs of open sites, so a closure that needs a
turn is found only through the cell it makes: growing cells with [`MetaParticleSpecies`](@ref)
puts the turn inside the cell, leaving a lattice that is a pure translation again.

Use [`tilingenum`](@ref) to process the closures as they are found, and to stop early.
"""
function tilings(poly::Polyform; kwargs...)
    out = _tilingtype(poly)[]
    tilingenum(t -> (push!(out, t); ACCEPT), poly; kwargs...)
    return out
end

function _tilingtype(poly::Polyform{D}) where {D}
    rules = bindingrules(poly)
    V = SVector{D,numtype(rules)}
    return Tiling{D,particletype(rules),typeof(rules),typeof(graphrep(poly)),V}
end

# The first tiling of `poly` that `pred` accepts, or `nothing`, leaving the rest unenumerated.
function _findtiling(pred::F, poly::Polyform; kwargs...) where {F}
    hit = Ref{Union{Nothing,_tilingtype(poly)}}(nothing)
    tilingenum(poly; kwargs...) do t
        pred(t) || return ACCEPT
        hit[] = t
        return BREAK
    end
    return hit[]
end

# What every placement is measured against: the cell being translated, the candidate translations,
# and what the cell brings to the search before a single copy is laid down.
struct _ShellSearch{PF,V}
    cell::PF                 # the candidate cell, as the meta-polyform it was grown as
    vectors::Vector{V}
    ofvertex::Vector{Int}    # first vertex of a site -> its index among the cell's, else 0
    nreps::Int
    spent::BitVector         # sites no translate can use: bound inside the cell, or inert
end

# The rules the cell's meta-particles were lifted from: each wraps a polyform that carries them.
_baserules(cell::Polyform) = bindingrules(polyform(first(species(bindingrules(cell)))))

# Every tiling one candidate cell admits, streamed to `f`.
function _celltilings(f::F, cell::Polyform, nreps::Integer) where {F}
    sites = collect(bindingsites(cell))
    # a translate can use a site that is unbound inside the cell and not inert; everything else is
    # spent before the search starts, and a `complete` closure is one that uses up what is left
    spent = trues(length(sites))
    for l in opensites(cell)
        spent[siteindex(cell, l)] = false
    end
    free = [s for (i, s) in enumerate(sites) if !spent[i]]
    isempty(free) && return ACCEPT
    metarules = bindingrules(cell)
    vectors = _candidatelatticevectors(free, metarules)
    isempty(vectors) && return ACCEPT

    # only a cell that got this far is worth reading back as a polyform of the underlying rules.
    # The meta cell and its recast share a vertex numbering, so the contacts the search records
    # against the one are edges of the other
    cellpoly = recast(cell, _baserules(cell))
    search = _ShellSearch(cell, vectors, _siteofvertex(sites), Int(nreps), spent)
    emit(_, contacts) = f(Tiling(cellpoly, contacts))
    return _tilings!(emit, search, dimension(metarules), Int[], 1)
end

# Translations that lay an open site onto a compatible, already-aligned partner site. A site pair
# is seen both ways round, so every direction arrives twice, once with each sign; orienting each
# one keeps a single representative, the shells covering the other sign.
function _candidatelatticevectors(sites, rules)
    intmat = interactionmatrix(rules)
    vecs = typeof(first(sites).pose.x)[]
    for s1 in sites, s2 in sites
        intmat[color(s1), color(s2)] || continue
        isaligned(s1, s2) || continue
        v = s1.pose.x - s2.pose.x
        norm(v) < _tol(v) && continue
        v = _orient(v)
        any(u -> u ≈ v, vecs) || push!(vecs, v)
    end
    return vecs
end

# Depth-first search over strictly growing vector index sets, streaming every set that closes.
# A set that fails cannot be rescued by adding a vector to it -- the added vector only lays down
# more copies, and every objection is to a copy -- so a failure prunes the whole subtree, as does
# an `f` that rejects.
function _tilings!(f::F, s::_ShellSearch, maxvecs::Integer, chosen::Vector{Int},
    from::Integer) where {F}
    length(chosen) == maxvecs && return ACCEPT
    for idx in from:length(s.vectors)
        push!(chosen, idx)
        closure = _closure(s, chosen)
        if closure !== nothing
            signal = f(chosen, closure)
            signal == BREAK && return BREAK
            signal != REJECT && _tilings!(f, s, maxvecs, chosen, idx + 1) == BREAK && return BREAK
        end
        pop!(chosen)
    end
    return ACCEPT
end

# Lay a copy of the cell at every lattice point `chosen` reaches within `nreps` shells, and return
# the bonds the copies form with it, one [`Contact`](@ref) per bond. `nothing` if the placement is
# no tiling at all -- copies overlap, a contact is not a valid bond or does not join two sites of
# the cell, a site is claimed twice, or one of the chosen vectors buys no bond, which would only
# stack disconnected copies.
function _closure(s::_ShellSearch, chosen::Vector{Int})
    parts = s.cell.particles
    rules = bindingrules(s.cell)    # the *meta* rules: the cell is translated as meta-particles
    # what each of the cell's sites is bonded to: 0 while it is free, -1 for one no translate can
    # use. A site carries at most one bond, so a second claim on one is a contradiction
    partner = [s.spent[i] ? -1 : 0 for i in eachindex(s.spent)]
    contacts = Contact[]
    bought = zeros(Int, length(chosen))
    placed = eltype(parts)[]

    # the shell runs both ways along every vector, so one choice of vectors stands for the whole
    # lattice it generates rather than for one cone of it
    for m in Iterators.product(ntuple(_ -> (-s.nreps):(s.nreps), length(chosen))...)
        all(iszero, m) && continue
        t = sum(m[i] * s.vectors[chosen[i]] for i in eachindex(chosen))
        for part in parts
            sp = part + t
            first(_overlap_and_contacts(placed, sp, rules)) === true && return nothing
            ov, cts = _overlap_and_contacts(parts, sp, rules)
            ov && return nothing
            for contact in cts
                # both endpoints must be sites of the cell
                i1, i2 = _lookup(s, contact.vs1), _lookup(s, contact.vs2)
                (i1 == 0 || i2 == 0) && return nothing
                # the shell runs both ways, so the copies on either side both report the bond
                # between them; the second sighting is that same bond, not a second claim
                partner[i1] == i2 && partner[i2] == i1 && continue
                (partner[i1] == 0 && partner[i2] == 0) || return nothing
                partner[i1], partner[i2] = i2, i1
                push!(contacts, contact)
                # the copy is a translate along the last vector it moves on, so credit that one
                bought[findlast(!iszero, m)] += 1
            end
            push!(placed, sp)
        end
    end
    all(>(0), bought) || return nothing
    return contacts
end

# A contact comes back as the vertex range its site occupies, and a translated copy keeps the
# ranges of the cell it was translated from, so the first vertex identifies the site. An array
# rather than a `Dict`: the vertices are a contiguous range, so indexing is the whole lookup.
function _siteofvertex(sites)
    ofvertex = zeros(Int, maximum(s -> last(s.vertices), sites; init=0))
    for (i, s) in enumerate(sites)
        ofvertex[first(s.vertices)] = i
    end
    return ofvertex
end
_lookup(s::_ShellSearch, vs) = first(vs) <= length(s.ofvertex) ? s.ofvertex[first(vs)] : 0

# spatial units are O(1), so use sqrt(eps) for absolute tolerance
_tol(::AbstractVector{F}) where {F} = sqrt(eps(F))

# A vector and its negative generate the same translations, so pick the one whose first
# significant component is positive.
function _orient(v::AbstractVector)
    i = findfirst(x -> abs(x) > _tol(v), v)
    return isnothing(i) || v[i] > 0 ? v : -v
end

"""
    isunitcell(poly::Polyform; kwargs...)

Whether `poly` tiles space by translations with every open site closed — see [`tilings`](@ref).
"""
isunitcell(poly::Polyform; kwargs...) = _findtiling(iscomplete, poly; kwargs...) !== nothing

"""
    tilelatticevectors(poly::Polyform; kwargs...)

The lattice vectors of the first complete tiling, or `nothing` — see [`tilings`](@ref).
"""
function tilelatticevectors(poly::Polyform; kwargs...)
    t = _findtiling(iscomplete, poly; kwargs...)
    return isnothing(t) ? nothing : latticevectors(t)
end

"""
    cantile(rules::BindingRules; maxtilesize, kwargs...)

Search the structures of `rules` up to `maxtilesize` particles for a unit cell; returns the
first one found, or `nothing`.
"""
function cantile(rules::BindingRules; maxtilesize, kwargs...)
    hit = Ref{Any}(nothing)
    f(s, _) = isunitcell(s; kwargs...) ? (hit[]=copy(s); BREAK) : ACCEPT
    polyenum(f, rules; maxsize=maxtilesize)
    return hit[]
end

"""
    canchain(rules::BindingRules; maxlength=chainstatebound(rules) + 1, kwargs...)

Grow linear chains of `rules`, attaching only at their two ends, and return the first one that
closes periodically — or `nothing` if none does within `maxlength` particles.

  - `maxlength`: longest chain to try. The default is long enough that no periodic chain can hide
    below it, see [`chainstatebound`](@ref); lower it only to trade certainty for speed
  - other keyword arguments go to [`tilings`](@ref)

A chain that closes proves the rules admit arbitrarily large structures: truncating its periodic
continuation gives a valid polyform of every length. This is far cheaper than [`cantile`](@ref),
which sweeps every structure, because chains branch only at their ends, and a system whose
structures all close runs out of chains on its own well before the bound.

Restricting the search to chains costs nothing in reach. Any connected part of a valid structure
is valid here, and an infinite structure's bond graph is connected and locally finite, so it
contains an infinite path: branching growth always implies chain growth.

**In 2D** `nothing` at the default `maxlength` therefore means no infinite structure exists. A
chain long enough to repeat one of its finitely many states repeats the rigid motion between the
two occurrences forever, and in the plane that motion is a translation or a rotation — a rotation
keeps every particle on a circle, so its iterates must eventually collide. Unbounded growth
leaves only the translation, which is a periodic chain.

**In 3D** it does not: the motion can be a screw, whose translation along the axis never
collides. A screw by `2πp/q` closes into a translation over `q` copies and is found with
`maxorder ≥ q`, but an irrational one never closes at all, and no periodicity test can see it.
See [`isunbounded`](@ref).
"""
function canchain(rules::BindingRules;
    maxlength::Integer=chainstatebound(rules) + 1, kwargs...)
    frontier = [Polyform(rules, i) for i in 1:nspecies(rules)]
    seen = Set(hash(graphrep(p)) for p in frontier)
    while !isempty(frontier)
        poly = popfirst!(frontier)
        _findtiling(Returns(true), poly; kwargs...) === nothing || return poly
        nparticles(poly) < maxlength || continue
        for child in _extendends(poly)
            hash(graphrep(child)) in seen && continue
            push!(seen, hash(graphrep(child)))
            push!(frontier, child)
        end
    end
    return nothing
end

"""
    isunbounded(rules::BindingRules; maxlength=chainstatebound(rules) + 1)

Whether `rules` admits arbitrarily large structures. See [`growthwitness`](@ref), which returns
the repeating motion behind a `true`.
"""
isunbounded(rules::BindingRules; kwargs...) = growthwitness(rules; kwargs...) !== nothing

"""
    growthwitness(rules::BindingRules; maxlength=chainstatebound(rules) + 1)

Find a rigid motion `g` and a polyform `C` such that `C, g(C), g²(C), …` is an infinite valid
structure, or `nothing` if none exists.

  - returns `(; generator, cell, period)`: the motion, the particles it repeats, and how many
    particles that is

Chains are grown one particle at a time, tracking each particle's state — its species, the site
its incoming bond uses, and that bond's phase. There are [`chainstatebound`](@ref) such states, so
a longer chain repeats one, and the motion `g` between the two occurrences repeats forever.
Restricting to chains loses nothing: any connected part of a valid structure is valid, and an
infinite structure's bond graph contains an infinite path.

By Chasles' theorem `g` is a screw — a rotation by `θ` about an axis together with a translation
`h` along it, degenerating to a pure rotation (`h = 0`) or a pure translation (`θ = 0`). The two
cases settle the question:

  - `h = 0`: every copy stays the same distance from the axis, so all of them lie in one bounded
    region. Infinitely many particles of positive volume cannot, so growth stops. In 2D this is
    every rotation, which is why the plane has no unbounded non-periodic growth
  - `h ≠ 0`: copies `m` periods apart sit `m·h` apart along the axis, so once `|m·h|` exceeds the
    cell's own axial extent they cannot touch. Only `m ≤ ⌈L/|h|⌉` placements can collide, and
    checking those settles the whole infinite chain

Both bounds are exact rather than cutoffs, so `false` is a proof and not a budget running out.
The shapes never enter: a particle's reach along the axis is bounded by its own radius, so
non-convex pieces can interlock only within the same window. What the guarantee does inherit is
the exactness of pairwise `overlap` for the species in play.
"""
function growthwitness(rules::BindingRules; maxlength::Integer=chainstatebound(rules) + 1)
    for i in 1:nspecies(rules)
        w = _walkchain(Polyform(rules, i), [(i, 0, 0)], maxlength)
        w === nothing || return w
    end
    return nothing
end

"""
    chainstatebound(rules::BindingRules)

How long a chain of `rules` can get before it must repeat itself.

What a chain can do next depends only on the species at its end, which of that species' sites
carries the incoming bond, and with which phase — finitely many states. A longer chain visits one
twice, and the stretch between the two visits is a cell that repeats forever, so searching past
this length can find no periodic chain that a shorter one would have missed.
"""
function chainstatebound(rules::BindingRules)
    total = 0
    for i in 1:nspecies(rules)
        ps = species(rules, i)
        for k in 1:nsites(ps)
            site = bindingsite(ps, k)
            phases = 1
            for loc in possible_attachments(rules, color(site))
                mate = bindingsite(rules, loc)
                phases = max(phases, _ndistincttwists(mate, site))
            end
            total += phases
        end
    end
    return total
end

# Depth-first over directed chains, extending only at the particle last added, so the states seen
# so far are exactly the chain read from its start.
function _walkchain(poly::Polyform, states, maxlength)
    rules = bindingrules(poly)
    n = nparticles(poly)
    part = poly.particles[end]
    for k in 1:nsites(part, rules)
        site = bindingsite(part, rules, k)
        _isbound_vertex(poly, part, first(site.vertices); canonidxs=false) && continue
        isinert(rules, color(site)) && continue
        for loc in possible_attachments(rules, color(site))
            mate = bindingsite(rules, loc)
            for r in 0:(_ndistincttwists(site, mate)-1)
                child = copy(poly)
                ismissing(raise!(child, site, loc, r)) && continue
                # `raise!` forms every geometric contact, so growth may cross-link back onto the
                # chain. That is a valid structure and often the interesting one -- the walk only
                # insists that one particle was added, and always extends the newest
                nparticles(child) == n + 1 || continue

                state = (loc.species, loc.site, Int(r))
                j = findlast(==(state), states)
                if j !== nothing
                    w = _screwwitness(child, j, n + 1)
                    w === nothing || return w
                end
                if n + 1 < maxlength
                    w = _walkchain(child, push!(copy(states), state), maxlength)
                    w === nothing || return w
                end
            end
        end
    end
    return nothing
end

# The stretch between two occurrences of a state repeats under `g`. Whether that goes on forever
# is decided by the screw's pitch, and if it can, by finitely many placements.
function _screwwitness(poly::Polyform, i::Integer, j::Integer)
    rules = bindingrules(poly)
    cell = poly.particles[i:(j-1)]
    g = poly.particles[j].pose * inv(poly.particles[i].pose)

    axis, h = _screwpitch(g)
    h === nothing && return nothing               # a pure rotation cannot escape its own annulus

    us = (dot(p.pose.x, axis) for p in cell)
    rmax = maximum(bounding_radius(species(rules, speciesindex(p))) for p in cell)
    extent = maximum(us) - minimum(us) + 2 * rmax
    gm = g
    for _ in 1:ceil(Int, extent/abs(h))
        for p in cell
            moved = gm * p
            first(_overlap_and_contacts(cell, moved, rules)) === true && return nothing
        end
        gm = gm * g
    end
    return (; generator=g, cell=cell, period=length(cell))
end

# `(axis, pitch)` of a rigid motion, or `(nothing, nothing)` when it has no translation along its
# axis: a pure rotation in 3D, any rotation in 2D, or the identity.
function _screwpitch(g::Pose{D,F}) where {D,F}
    tol = sqrt(eps(F))
    # a composed motion accumulates its angle, so a full turn comes back as 2π rather than 0;
    # wrap to (-π, π] before asking whether it turns at all
    turning = abs(rem(rotation_angle(g.psi), 2 * F(π), RoundNearest)) > tol
    if !turning
        n = norm(g.x)
        n > tol || return nothing, nothing        # the identity: no repeat at all
        return g.x / n, n                         # a pure translation is a screw of zero angle
    end
    D == 3 || return nothing, nothing             # in the plane a turn has no axis to climb
    axis = rotation_axis(g.psi)
    pitch = dot(g.x, axis)
    return axis, abs(pitch) > tol ? pitch : nothing
end

# How many of a particle's sites are bonded; a chain's ends are the particles with at most one.
function _bonddegree(poly::Polyform, part)
    rules = bindingrules(poly)
    return count(1:nsites(part, rules)) do k
        _isbound_vertex(poly, part, first(bindingsite(part, rules, k).vertices); canonidxs=false)
    end
end

# Every chain one particle longer, grown at an end. `raise!` forms all geometric contacts, so an
# attachment can branch the chain or close it into a ring; those are dropped, leaving paths.
function _extendends(poly::Polyform)
    rules = bindingrules(poly)
    out = typeof(poly)[]
    for part in poly.particles
        nparticles(poly) == 1 || _bonddegree(poly, part) <= 1 || continue
        for k in 1:nsites(part, rules)
            site = bindingsite(part, rules, k)
            _isbound_vertex(poly, part, first(site.vertices); canonidxs=false) && continue
            isinert(rules, color(site)) && continue
            for loc in possible_attachments(rules, color(site))
                mate = bindingsite(rules, loc)
                for r in 0:(_ndistincttwists(site, mate)-1)
                    child = copy(poly)
                    ismissing(raise!(child, site, loc, r)) && continue
                    _ischain(child) && push!(out, child)
                end
            end
        end
    end
    return out
end

# A path: every particle bonded to at most two others, and no cycle.
function _ischain(poly::Polyform)
    nbonds = sum(part -> _bonddegree(poly, part), poly.particles; init=0) ÷ 2
    return nbonds == nparticles(poly) - 1 &&
           all(part -> _bonddegree(poly, part) <= 2, poly.particles)
end
