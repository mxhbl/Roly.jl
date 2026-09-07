"""
    Particle{P<:Pose}

A `Particle` describes an indivisible subunit that makes up a `Polyform`.

Every particle is a member of a `ParticleSpecies` (identified by `speciesindex` within
an `BindingRules`) and is located in space at a specific `Pose`.
"""
struct Particle{P<:Pose}
    pose::P
    leadingvertex::Int
    speciesindex::Int
end

"""
    Particle(rules::BindingRules, speciesindex::Integer, pose=nothing; leadingvertex::Integer)

Create a particle of species `speciesindex` within `rules` at a given `pose`.
"""
function Particle(rules::BindingRules, speciesindex::Integer, pose=nothing; leadingvertex::Integer)
    if isnothing(pose)
        pose = one(posetype(rules))
    end
    return Particle(pose, leadingvertex, speciesindex)
end

Base.:*(p, part::Particle) = typeof(part)(p * part.pose, part.leadingvertex, part.speciesindex)
Base.:*(part::Particle, p) = typeof(part)(part.pose * p, part.leadingvertex, part.speciesindex)
translate(part::Particle, v) =
    typeof(part)(translate(part.pose, v), part.leadingvertex, part.speciesindex)

@inline graphvertices(p::Particle, rules::BindingRules) =
    (1:nv(graphrep(species(rules, p.speciesindex)))) .+ (p.leadingvertex - 1)

"""
    leadingvertex(p::Particle)

Return the leading vertex of particle `p`.
"""
leadingvertex(p::Particle) = p.leadingvertex

"""
    shift_leadingvertex(p::Particle, v::Integer)

Return a copy of `p` whose block of graph vertices starts `v` further along.

Used when removing a particle compacts the original vertex numbering; compare
[`shift_vertices`](@ref) for binding sites.
"""
@inline shift_leadingvertex(p::Particle, v::Integer) = typeof(p)(p.pose, p.leadingvertex + v, p.speciesindex)

"""
    speciesindex(p::Particle)

Return the species index of particle `p` within its `BindingRules`.
"""
speciesindex(p::Particle) = p.speciesindex

"""
    nsites(p::Particle, rules::BindingRules)

Return the number of binding sites of particle `p`.
"""
nsites(p::Particle, rules::BindingRules) = nsites(species(rules, p.speciesindex))

"""
    bindingsite(p::Particle, rules::BindingRules, i::Integer)

Return the `i`th binding site of particle `p`.
"""
function bindingsite(p::Particle, rules::BindingRules, i::Integer)
    return p.pose * shift_vertices(bindingsite(species(rules, p.speciesindex), i), leadingvertex(p) - 1)
end

"""
    bindingsites(p::Particle, rules::BindingRules)

Return an iterator over the binding sites of particle `p`.
"""
function bindingsites(p::Particle, rules::BindingRules)
    return (bindingsite(p, rules, i) for i in 1:nsites(p, rules))
end

"""
    could_contact(p1::Particle, p2::Particle, rules::BindingRules)

Return `true` if the particles could potentially be in contact.
"""
function could_contact(p1::Particle, p2::Particle, rules::BindingRules)
    return could_contact(species(rules, p1.speciesindex) => p1.pose, species(rules, p2.speciesindex) => p2.pose)
end

"""
    overlap(p1::Particle, p2::Particle, rules::BindingRules)

Return `true` if the particles are overlapping.
"""
function overlap(p1::Particle, p2::Particle, rules::BindingRules)
    spcs1, spcs2 = species(rules, p1.speciesindex), species(rules, p2.speciesindex)

    # Binding rules can override the overlap check
    if rules._onlattice
        atol = sqrt(eps(numtype(spcs1))) * (bounding_radius(spcs1) + bounding_radius(spcs2))
        return norm(p1.pose.x - p2.pose.x) < atol
    end
    return overlap(spcs1 => p1.pose, spcs2 => p2.pose)
end

"""
    isconvex(p::Particle, rules::BindingRules)

Return true if the particle has a convex shape, which enables minor optimizations when
checking for overlaps.
"""
isconvex(p::Particle, rules::BindingRules) = isconvex(species(rules, p.speciesindex))

function Base.show(io::Core.IO, p::Particle)
    print(io, "Particle[si=$(p.speciesindex), lv=$(p.leadingvertex)]")
end
function Base.show(io::Core.IO, ::MIME"text/plain", p::Particle)
    print(io, "Particle:\n")
    print(io, " - speciesindex: $(p.speciesindex)\n")
    print(io, " - leadingvertex: $(p.leadingvertex)\n")
    print(io, " - pose: $(p.pose)")
end