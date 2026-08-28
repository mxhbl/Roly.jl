# A meta-particle has no shape of its own: it is the polyform it wraps. So it is drawn by drawing
# that polyform's particles, with the colors of the meta-species painted back onto the sites it
# exposes. Everything else -- sites a bond inside the cluster already consumes, and open sites the
# species chose not to expose -- is drawn in `INERT_COLOR`, the same neutral the rest of the
# package uses for a site nothing can bond to.
#
# Neutral rather than a tint of the species color: a species palette ramps from pale to dark, and
# a tinted species color lands squarely on its pale end. On a block of squares the first exposed
# site came out within 0.1 of the interior and the two were indistinguishable. Grey has no hue to
# collide with.

"""
    _metasites(spcs::MetaParticleSpecies, speciesindex, rules, sitecolor)

Return `(speciesindex, exposed)`, where `exposed` maps the [`ParticleSiteLoc`](@ref) of each site
`spcs` exposes to the color it should be drawn in.

The colors come from [`_resolve_colors`](@ref) on the meta-species itself, so they follow the
meta `rules`: a site the new rules leave inert is drawn inert, whatever it was inside the
polyform.
"""
function _metasites(spcs::MetaParticleSpecies, speciesindex, rules, sitecolor)
    # Untinted, in 3D as well as 2D. A polyhedron tints its own faces because its whole surface is
    # bulk, but here the bulk is the neutral interior and the exposed sites are the few markers on
    # it, so they carry the weight. Tinting them towards white would walk them into the interior.
    si, _, colors, _ = _resolve_colors(spcs, speciesindex, rules, sitecolor)

    # Match the species' sites back to the polyform's by vertex range, not by `==`: a recolored
    # site is no longer equal to the one it was taken from.
    poly = polyform(spcs)
    byvertex = Dict(first(bindingsite(poly, l).vertices) => l for l in exposedsitelocs(poly))
    exposed = Dict(byvertex[first(bindingsite(spcs, i).vertices)] => colors[i] for i in 1:nsites(spcs))
    return si, exposed
end

# The constituent particles of `spcs` at `pose`, each with the species to draw it as, its world
# pose, and a `sitecolor` callback that knows which particle it belongs to. The callback is what
# `_resolve_colors` down the line asks for every site, so the inner rules never enter: colors are
# settled here.
function _metaparts(spcs::MetaParticleSpecies, pose, speciesindex, rules, sitecolor, interiorcolor)
    _, exposed = _metasites(spcs, speciesindex, rules, sitecolor)
    poly = polyform(spcs)
    inner = bindingrules(poly)
    wash = isnothing(interiorcolor) ? RGBf(INERT_COLOR) : RGBf(interiorcolor)
    return [
        (
            species(inner, part.speciesindex),
            pose * part.pose,
            part.speciesindex,
            (_, k) -> get(exposed, ParticleSiteLoc(p, k), wash),
        ) for (p, part) in enumerate(poly.particles)
    ]
end

"""
    particlemesh(spcs::MetaParticleSpecies, pose; kwargs...)

Merge the meshes of the particles `spcs` wraps into one, so a meta-assembly depth-sorts as a
whole. `nothing` if any constituent has no mesh form, as in 2D.
"""
function particlemesh(
    spcs::MetaParticleSpecies{D,F},
    pose::Pose=Pose{D,F}();
    sitecolor=nothing,
    speciesindex=nothing,
    rules=nothing,
    interiorcolor=nothing,
    kwargs...,
) where {D,F}
    pts, tris, cols = Point3f[], NTuple{3,Int}[], RGBAf[]
    for (ps, ppose, si, sc) in _metaparts(spcs, pose, speciesindex, rules, sitecolor, interiorcolor)
        geom = particlemesh(ps, ppose; sitecolor=sc, speciesindex=si, rules=nothing, kwargs...)
        isnothing(geom) && return nothing
        p, t, c = geom
        offset = length(pts)
        append!(pts, p)
        append!(cols, c)
        for (a, b, d) in t
            push!(tris, (a + offset, b + offset, d + offset))
        end
    end
    return pts, tris, cols
end

"""
    plot_particlespecies!(ax, spcs::MetaParticleSpecies, pose; kwargs...)

Draw a meta-particle as the polyform it wraps, recolored to the sites it exposes.

`interiorcolor` overrides the color the rest of the polyform is drawn in; by default it is
`INERT_COLOR`, since nothing can attach through those sites at this level.
"""
function plot_particlespecies!(
    ax,
    spcs::MetaParticleSpecies{D,F},
    pose::Pose=Pose{D,F}();
    sitecolor=nothing,
    speciesindex=nothing,
    rules=nothing,
    interiorcolor=nothing,
    kwargs...,
) where {D,F}
    parts = _metaparts(spcs, pose, speciesindex, rules, sitecolor, interiorcolor)
    geom = particlemesh(spcs, pose; sitecolor, speciesindex, rules, interiorcolor, kwargs...)
    if !isnothing(geom)
        pts, tris, cols = geom
        return mesh!(ax, pts, _facematrix(tris); color=cols, shading=NoShading)
    end
    out = nothing
    for (ps, ppose, si, sc) in parts
        out = plot_particlespecies!(ax, ps, ppose; sitecolor=sc, speciesindex=si, rules=nothing, kwargs...)
    end
    return out
end
