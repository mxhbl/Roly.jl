"""
    PolyformEnvironment{N}

The local environment of `N` root particles inside a [`Polyform`](@ref): all particles within
graph distance `depth` of a root, as an isomorphism class. See [`ParticleEnvironment`](@ref) and
[`BondEnvironment`](@ref).

  - `graph`: canonized graph representation with the root particles distinguished, in order
  - `rootvertices`: one canonical site vertex per root particle
  - `depth`: crop radius around each root
  - `rules`: the `BindingRules` the labels refer to
"""
struct PolyformEnvironment{N,G<:AbstractNautyGraph,S<:BindingRules}
    graph::G
    rootvertices::NTuple{N,Int}
    depth::Int
    rules::S
end

"""
    ParticleEnvironment

[`PolyformEnvironment`](@ref) of a single particle: the ball of all particles within graph distance
`depth` of the root.
"""
const ParticleEnvironment{G,S} = PolyformEnvironment{1,G,S}

"""
    BondEnvironment

[`PolyformEnvironment`](@ref) of a bond: all particles within graph distance `depth` of either endpoint,
with the endpoints as roots, in order.
"""
const BondEnvironment{G,S} = PolyformEnvironment{2,G,S}

Base.hash(e::PolyformEnvironment, h::UInt) = hash(e.graph, hash(e.depth, h))
function Base.:(==)(a::PolyformEnvironment{N}, b::PolyformEnvironment{N}) where {N}
    return a.depth == b.depth && a.rules === b.rules && a.graph == b.graph
end
Base.:(==)(::PolyformEnvironment, ::PolyformEnvironment) = false

_envname(N) = N == 1 ? "ParticleEnvironment" : N == 2 ? "BondEnvironment" : "PolyformEnvironment{$N}"
function Base.show(io::Core.IO, e::PolyformEnvironment{N}) where {N}
    return print(io, _envname(N), "[k=$(e.depth), nv=$(nv(e.graph))]")
end

# Root marks are added in multiples of this, so marked labels never collide with site labels.
function _markoffset(rules::BindingRules)
    m = 0
    # find the largest label of any site in `rules`
    for i in 1:nspecies(rules)
        m = max(m, maximum(labels(graphrep(species(rules, i)))))
    end
    return m
end

"""
    _orig2canon!(h)

Canonize `h` in place and return its `orig2canon` map, the direction [`Polyform`](@ref) stores
under that name: `orig2canon[v]` is where the vertex that was `v` ended up.

NautyGraphs runs the other way. `canonize!` hands back `canon2orig`, so this inverts once, here,
rather than at each call site -- every caller in this file wants to follow a vertex it already
holds into the canonical graph.
"""
_orig2canon!(h::AbstractNautyGraph) = invperm(collect(Int, canonize!(h)))

# Mark and canonize: bump the labels of the i-th group of vertices by i*offset. Marking every vertex
# of a root particle (rather than a single one) preserves the particle's internal symmetry.
# Returns the marked graph and its `orig2canon`, see `_orig2canon!`.
function _canonmarked(g::AbstractNautyGraph, groups, offset)
    h = copy(g)
    for (i, vs) in enumerate(groups), v in vs
        setlabel!(h, v, labels(h)[v] + i * offset)
    end
    return h, _orig2canon!(h)
end

# Buffers for the particle-hop search, reusable across calls and polyforms.
struct EnvironmentBuffers
    dist::Vector{Int}               # graph-distance of the particles
    queue::Vector{Int}              # interal buffer
    vertex2particle::Vector{Int}    # original vertex -> owning particle
    neighbors::Vector{Int}          # neighbor scratch for `bfs!`
end
EnvironmentBuffers() = EnvironmentBuffers(Int[], Int[], Int[], Int[])
function Base.copy(bufs::EnvironmentBuffers)
    return EnvironmentBuffers(copy(bufs.dist), copy(bufs.queue), copy(bufs.vertex2particle),
        copy(bufs.neighbors))
end

# Graph distances from the `roots` via `bfs!`: vertices of the same particle are free, bonds
# cost one hop. Capped at `maxdepth`. Fills and returns `bufs.dist`.
#
# `vertex2particle` rather than polyform.jl's `_same_particle`: the search needs the neighbor's
# particle, not just whether it is this one, and both that and `_same_particle` scan the particle
# list per edge where the table is a lookup.
function _particledists!(bufs::EnvironmentBuffers, poly::Polyform, roots; maxdepth)
    rules = bindingrules(poly)
    g = graphrep(poly)

    resize!(bufs.dist, nparticles(poly))
    resize!(bufs.queue, nparticles(poly))
    resize!(bufs.vertex2particle, nv(g))
    for (i, part) in enumerate(poly.particles), ov in graphvertices(part, rules)
        bufs.vertex2particle[ov] = i
    end

    fill!(bufs.dist, -1)
    bfs!(bufs.dist, bufs.queue, roots; maxdepth) do p
        # Any edge leaving the particle is a bond; interior edges stay within it.
        empty!(bufs.neighbors)
        for ov in graphvertices(particles(poly, p), rules)
            outneighs = NautyGraphs.adjrow(g, tocanon(poly, ov))
            for w in eachindex(outneighs)
                outneighs[w] || continue
                q = bufs.vertex2particle[toorig(poly, w)]
                q == p || push!(bufs.neighbors, q)
            end
        end
        bufs.neighbors
    end
    return bufs.dist
end

# Induced marked subgraph on the particles with dist in [0, depth], with the `roots` particles
# marked in order. `rootvertices` gives one original site vertex per root; the returned tuple holds
# their canonical indices in the new graph.
function _envgraph(poly::Polyform, dist::AbstractVector{<:Integer}, depth::Integer, roots, rootvertices)
    rules = bindingrules(poly)
    verts = Int[]
    offsets = zeros(Int, nparticles(poly))   # particle -> position of its first vertex in `verts`
    for p in 1:nparticles(poly)
        0 <= dist[p] <= depth || continue
        offsets[p] = length(verts) + 1
        for ov in graphvertices(particles(poly, p), rules)
            push!(verts, tocanon(poly, ov))
        end
    end

    h = graphrep(poly)[verts]
    groups = map(p -> offsets[p]:(offsets[p]+nsites(particles(poly, p), rules)-1), roots)
    hm, orig2canon = _canonmarked(h, groups, _markoffset(rules))

    # Follow each root vertex through the two renumberings the crop imposed. `ov` is an original
    # vertex of particle `p`, so `ov - leadingvertex` is its offset inside that particle's block,
    # `offsets[p]` is where the block starts in `h`, and `orig2canon` carries the result into `hm`.
    canonrootvertices = map(roots, rootvertices) do p, ov
        orig2canon[offsets[p]+ov-leadingvertex(particles(poly, p))]
    end
    return hm, canonrootvertices
end

"""
    ParticleEnvironment(poly::Polyform, root::Integer; depth::Integer)

Extract the environment of particle `root` in `poly`: the ball of all particles within graph
distance `depth`.
"""
function ParticleEnvironment(poly::Polyform, root::Integer; depth::Integer, bufs=EnvironmentBuffers())
    dist = _particledists!(bufs, poly, (root,); maxdepth=depth)
    rootvertex = first(graphvertices(particles(poly, root), bindingrules(poly)))
    hm, rootvertices = _envgraph(poly, dist, depth, (root,), (rootvertex,))
    return PolyformEnvironment(hm, rootvertices, Int(depth), bindingrules(poly))
end

"""
    BondEnvironment(poly::Polyform, bond; depth::Integer)

Extract the environment of `bond` (an element of `bonds(poly)`): all particles within graph
distance `depth` of either endpoint. The endpoints become the environment's roots, in the order
given.
"""
function BondEnvironment(poly::Polyform, bond::Pair; depth::Integer, bufs=EnvironmentBuffers())
    rules = bindingrules(poly)
    (p1, s1), (p2, s2) = bond.first, bond.second
    p1 == p2 && throw(ArgumentError("bond endpoints must be distinct particles"))

    dist = _particledists!(bufs, poly, (p1, p2); maxdepth=depth)
    sitevertex(p, s) = first(bindingsite(poly, ParticleSiteLoc(p, s)).vertices)
    hm, rootvertices = _envgraph(poly, dist, depth, (p1, p2),
        (sitevertex(p1, s1), sitevertex(p2, s2)))
    return PolyformEnvironment(hm, rootvertices, Int(depth), rules)
end

"""
    particleenvironments(poly::Polyform; depth::Integer)

Return the environment of every particle of `poly`, in particle order.
"""
function particleenvironments(poly::Polyform; depth::Integer)
    bufs = EnvironmentBuffers()
    return [ParticleEnvironment(poly, p; depth, bufs) for p in 1:nparticles(poly)]
end

"""
    reverse(env::BondEnvironment)

Return the environment with the two roots swapped.
"""
function Base.reverse(env::BondEnvironment)
    offset = _markoffset(env.rules)
    h = copy(env.graph)
    for v in vertices(h)
        l = labels(h)[v]
        if l > 2offset
            setlabel!(h, v, l - offset)
        elseif l > offset
            setlabel!(h, v, l + offset)
        end
    end
    orig2canon = _orig2canon!(h)
    return PolyformEnvironment(h, (orig2canon[env.rootvertices[2]], orig2canon[env.rootvertices[1]]), env.depth,
        env.rules)
end

### Enumeration of all particle environments of a rule set, by reverse search over rooted balls.

mutable struct EnvironmentState{P<:Polyform,G<:AbstractNautyGraph}
    poly::P               # working polyform; the root is always particle 1
    key::G                # root-marked canonical graph, the reverse-search identity
    rootvertex::Int       # canonical index of the root's first site vertex in `key`
end

function EnvironmentState(rules::BindingRules)
    poly = Polyform(rules)
    return EnvironmentState(poly, copy(graphrep(poly)), 0)
end
Base.copy(s::EnvironmentState) = EnvironmentState(copy(s.poly), copy(s.key), s.rootvertex)
function Base.copy!(dst::EnvironmentState, src::EnvironmentState)
    copy!(dst.poly, src.poly)
    copy!(dst.key, src.key)
    dst.rootvertex = src.rootvertex
    return dst
end

_rootgroup(poly::Polyform) =
    [tocanon(poly, ov) for ov in graphvertices(particles(poly, 1), bindingrules(poly))]

# Recompute the root-marked canonical key of a state; the whole polyform is its own ball.
function _rekey!(s::EnvironmentState)
    poly = s.poly
    hm, orig2canon = _canonmarked(graphrep(poly), (_rootgroup(poly),), _markoffset(bindingrules(poly)))
    s.key = hm
    s.rootvertex = orig2canon[tocanon(poly, first(graphvertices(particles(poly, 1), bindingrules(poly))))]
    return s
end

mutable struct EnvironmentEnumAux{BS<:BindingSite,G<:AbstractNautyGraph}
    seen::Set{G}                              # offspring keys of the current parent
    pairs::Vector{Tuple{BS,SpeciesSiteLoc}}   # attachments that stay within the ball
    bufs::EnvironmentBuffers
    depth::Int
end
function Base.copy(aux::EnvironmentEnumAux)
    return EnvironmentEnumAux(copy(aux.seen), copy(aux.pairs), copy(aux.bufs), aux.depth)
end

# Attachment options whose owner particle is at distance <= depth-1 from the root, so the new
# particle lands inside the ball.
function _collectballpairs!(aux::EnvironmentEnumAux, poly::Polyform)
    rules = bindingrules(poly)
    empty!(aux.pairs)
    dist = _particledists!(aux.bufs, poly, (1,); maxdepth=aux.depth)
    for (p, part) in enumerate(poly.particles)
        0 <= dist[p] <= aux.depth - 1 || continue
        for k in 1:nsites(part, rules)
            site = bindingsite(part, rules, k)
            _isbound_vertex(poly, part, first(site.vertices); canonidxs=false) && continue
            isinert(rules, color(site)) && continue
            for loc in possible_attachments(rules, color(site))
                push!(aux.pairs, (site, loc))
            end
        end
    end
    return aux.pairs
end

function _adjenv!(u::EnvironmentState, v::EnvironmentState, j::Integer, aux::EnvironmentEnumAux)
    rules = bindingrules(v.poly)
    if nparticles(v.poly) == 0
        j > nspecies(rules) && return nothing
        copy!(u.poly, Polyform(rules, j))
        return _rekey!(u)
    end

    j == 1 && _collectballpairs!(aux, v.poly)
    j > length(aux.pairs) && return nothing

    site, loc = aux.pairs[j]
    copy!(u.poly, v.poly)
    ismissing(raise!(u.poly, site, loc)) && return missing
    _rekey!(u)
    u.key in aux.seen && return missing
    push!(aux.seen, copy(u.key))
    return u
end

# The parent drops the canonically-last particle among those at maximal distance from the root.
# Removing a maximal-distance particle keeps a ball a ball with every other distance unchanged: a
# shortest root-path passes only through vertices of strictly smaller distance, so it never routes
# through a maximal-distance one.
function _lsenv!(w::EnvironmentState, v::EnvironmentState, bufs::EnvironmentBuffers)
    n = nparticles(v.poly)
    rules = bindingrules(v.poly)
    if n <= 1
        copy!(w.poly, Polyform(rules))
        copy!(w.key, graphrep(w.poly))
        w.rootvertex = 0
        return w
    end

    dist = _particledists!(bufs, v.poly, (1,); maxdepth=typemax(Int))
    D = maximum(dist)
    _, orig2canon = _canonmarked(graphrep(v.poly), (_rootgroup(v.poly),), _markoffset(rules))
    canonpos(p) = minimum(orig2canon[tocanon(v.poly, ov)]
                          for ov in graphvertices(particles(v.poly, p), rules))
    drop = argmax(canonpos, (p for p in 1:n if dist[p] == D))

    w.poly = subpolyform(v.poly, [p for p in 1:n if p != drop])
    return _rekey!(w)
end

"""
    particleenvironments(f, rules::BindingRules; depth, maxsize=Inf, maxstrs=Inf)

Enumerate every particle environment of `rules` at radius `depth`, streaming each to
`f(env, nparticles)`. Each class is visited exactly once; `f` may return `REJECT` to prune all
extensions of an environment, or `BREAK` to stop.
"""
function particleenvironments(f, rules::BindingRules; depth::Integer, maxsize=Inf, maxstrs=Inf,
    kwargs...)
    v₀ = EnvironmentState(rules)
    BS = sitetype(rules)
    G = typeof(graphrep(v₀.poly))
    aux = EnvironmentEnumAux(Set{G}(), Tuple{BS,SpeciesSiteLoc}[], EnvironmentBuffers(), Int(depth))
    lsbufs = EnvironmentBuffers()
    rsys = RSSystem((w, v) -> _lsenv!(w, v, lsbufs), _adjenv!, v₀;
        compare=(a, b) -> a.key == b.key, aux)

    frs = (s, _) -> f(PolyformEnvironment(copy(s.key), (s.rootvertex,), Int(depth), rules), nparticles(s.poly))
    return reversesearch(frs, rsys; maxdepth=maxsize, maxverts=maxstrs + 1, kwargs...)
end

"""
    particleenvironments(rules::BindingRules; depth, kwargs...)

Return all particle environments of `rules` at radius `depth`.
"""
function particleenvironments(rules::BindingRules; depth::Integer, kwargs...)
    G = typeof(graphrep(Polyform(rules)))
    envs = ParticleEnvironment{G,typeof(rules)}[]
    particleenvironments((e, _) -> (push!(envs, e); ACCEPT), rules; depth, kwargs...)
    return envs
end

# ------------------------------------------------------------------------------------------------
# Bond environments of a particle environment, by pure graph surgery: particles are the components
# of the interior edges, root-incident bonds are the bidirectional edges leaving the root
# component, and the crop is a component-level breadth-first search. No poses are involved.

# Partition the vertices of an environment graph into particles: bond edges are bidirectional,
# interior edges are not. (Same assumption as `exterior_edges`: no double edges inside one
# particle's graph.)
function _components(g::AbstractNautyGraph)
    nvg = nv(g)
    comp = zeros(Int, nvg)
    ncomp = 0
    stack = Int[]
    for v₀ in 1:nvg
        comp[v₀] == 0 || continue
        ncomp += 1
        comp[v₀] = ncomp
        push!(stack, v₀)
        while !isempty(stack)
            v = pop!(stack)
            for w in outneighbors(g, v)
                (has_edge(g, w, v) || comp[w] != 0) && continue
                comp[w] = ncomp
                push!(stack, w)
            end
            for w in inneighbors(g, v)
                (has_edge(g, v, w) || comp[w] != 0) && continue
                comp[w] = ncomp
                push!(stack, w)
            end
        end
    end
    compverts = [Int[] for _ in 1:ncomp]
    for v in 1:nvg
        push!(compverts[comp[v]], v)
    end
    return comp, compverts
end

"""
    bondenvironments(env::ParticleEnvironment)

Return the environments of every root-incident bond of `env`, cropped at `depth - 1`, with the
root as the first endpoint.

The crop inside the ball equals the crop inside any larger structure containing it: every path of
length at most `depth - 1` from a root neighbor stays within the ball.
"""
function bondenvironments(env::ParticleEnvironment)
    env.depth >= 1 || throw(ArgumentError("bond environments need `depth >= 1`"))
    g = env.graph
    offset = _markoffset(env.rules)
    comp, compverts = _components(g)
    ncomp = length(compverts)
    rootcomp = comp[env.rootvertices[1]]

    # particle-level adjacency; duplicate entries are harmless to `bfs!`
    adj = [Int[] for _ in 1:ncomp]
    for v in vertices(g), w in outneighbors(g, v)
        (has_edge(g, w, v) && comp[v] != comp[w]) || continue
        push!(adj[comp[v]], comp[w])
    end

    dist = zeros(Int, ncomp)
    queue = zeros(Int, ncomp)
    pos = zeros(Int, nv(g))
    out = BondEnvironment{typeof(env.graph),typeof(env.rules)}[]
    for u in compverts[rootcomp], w in outneighbors(g, u)
        (has_edge(g, w, u) && comp[w] != rootcomp) || continue
        partner = comp[w]

        # crop: every particle within depth-1 of either endpoint
        fill!(dist, -1)
        bfs!(c -> adj[c], dist, queue, (rootcomp, partner); maxdepth=env.depth - 1)

        verts = Int[]
        for c in 1:ncomp
            dist[c] >= 0 || continue
            for v in compverts[c]
                push!(verts, v)
                pos[v] = length(verts)
            end
        end

        h = g[verts]
        for (k, v) in enumerate(verts)
            l = mod1(labels(h)[k], offset)   # strip the old root mark
            l += comp[v] == rootcomp ? offset : comp[v] == partner ? 2offset : 0
            setlabel!(h, k, l)
        end
        orig2canon = _orig2canon!(h)
        push!(out, PolyformEnvironment(h, (orig2canon[pos[u]], orig2canon[pos[w]]), env.depth - 1, env.rules))
    end
    return out
end

"""
    crop(env::ParticleEnvironment, depth)

The environment of the same root at a smaller `depth`: the radius-`depth` sub-ball.

Distances inside a ball equal the ambient distances to the root, so the crop of the ball is the
ball the ambient structure would have produced at the smaller radius.
"""
function crop(env::ParticleEnvironment, depth::Integer)
    0 <= depth <= env.depth || throw(ArgumentError("`depth` must lie in [0, $(env.depth)]"))
    depth == env.depth && return env
    g = env.graph
    comp, compverts = _components(g)
    ncomp = length(compverts)
    rootcomp = comp[env.rootvertices[1]]

    adj = [Int[] for _ in 1:ncomp]
    for v in vertices(g), w in outneighbors(g, v)
        (has_edge(g, w, v) && comp[v] != comp[w]) || continue
        push!(adj[comp[v]], comp[w])
    end

    dist = fill(-1, ncomp)
    queue = zeros(Int, ncomp)
    bfs!(c -> adj[c], dist, queue, (rootcomp,); maxdepth=depth)

    verts = Int[]
    pos = zeros(Int, nv(g))
    for c in 1:ncomp
        dist[c] >= 0 || continue
        for v in compverts[c]
            push!(verts, v)
            pos[v] = length(verts)
        end
    end

    # the root's mark survives the crop unchanged, so no relabeling is needed
    h = g[verts]
    orig2canon = _orig2canon!(h)
    return PolyformEnvironment(h, (orig2canon[pos[env.rootvertices[1]]],), Int(depth), env.rules)
end

"""
    crop(env::BondEnvironment, depth)

The environment of the same bond at a smaller `depth`: particles within `depth` of either root.
"""
function crop(env::BondEnvironment, depth::Integer)
    0 <= depth <= env.depth || throw(ArgumentError("`depth` must lie in [0, $(env.depth)]"))
    depth == env.depth && return env
    g = env.graph
    comp, compverts = _components(g)
    adj = _componentadjacency(g, comp, length(compverts))

    dist = fill(-1, length(compverts))
    queue = zeros(Int, length(compverts))
    bfs!(c -> adj[c], dist, queue,
        (comp[env.rootvertices[1]], comp[env.rootvertices[2]]); maxdepth=depth)

    verts, pos = _keptvertices(compverts, dist, nv(g))
    # both roots keep their marks, so the labels carry over unchanged
    h = g[verts]
    orig2canon = _orig2canon!(h)
    return PolyformEnvironment(h, (orig2canon[pos[env.rootvertices[1]]], orig2canon[pos[env.rootvertices[2]]]),
        Int(depth), env.rules)
end

"""
    rootenvironment(env::PolyformEnvironment, i, depth)

The [`ParticleEnvironment`](@ref) of root `i` of `env` at radius `depth`, extracted from `env`.

Valid whenever the radius-`depth` ball of that root lies inside `env` (for a bond environment,
any `depth <= env.depth`).
"""
function rootenvironment(env::PolyformEnvironment{N}, i::Integer, depth::Integer) where {N}
    1 <= i <= N || throw(ArgumentError("`env` has $N roots"))
    0 <= depth <= env.depth || throw(ArgumentError("`depth` must lie in [0, $(env.depth)]"))
    g = env.graph
    offset = _markoffset(env.rules)
    comp, compverts = _components(g)
    adj = _componentadjacency(g, comp, length(compverts))
    rootcomp = comp[env.rootvertices[i]]

    dist = fill(-1, length(compverts))
    queue = zeros(Int, length(compverts))
    bfs!(c -> adj[c], dist, queue, (rootcomp,); maxdepth=depth)

    verts, pos = _keptvertices(compverts, dist, nv(g))
    h = g[verts]
    for (k, v) in enumerate(verts)
        # strip every root mark, then re-mark the requested root as the sole root
        l = mod1(labels(h)[k], offset)
        setlabel!(h, k, l + (comp[v] == rootcomp ? offset : 0))
    end
    orig2canon = _orig2canon!(h)
    return PolyformEnvironment(h, (orig2canon[pos[env.rootvertices[i]]],), Int(depth), env.rules)
end

# Particle-level adjacency of an environment graph (bond edges are bidirectional).
function _componentadjacency(g, comp, ncomp)
    adj = [Int[] for _ in 1:ncomp]
    for v in vertices(g), w in outneighbors(g, v)
        (has_edge(g, w, v) && comp[v] != comp[w]) || continue
        push!(adj[comp[v]], comp[w])
    end
    return adj
end

# Vertices of the components reached by a crop search, with their positions in the kept order.
function _keptvertices(compverts, dist, nvg)
    verts = Int[]
    pos = zeros(Int, nvg)
    for c in eachindex(compverts)
        dist[c] >= 0 || continue
        for v in compverts[c]
            push!(verts, v)
            pos[v] = length(verts)
        end
    end
    return verts, pos
end
