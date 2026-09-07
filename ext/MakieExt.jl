module MakieExt

using Roly
using Makie
using LinearAlgebra: dot, norm, normalize
import Roly: species, bindingrules, polyformplot, polyformplot!, render
import Roly: Tiling, unitcell, latticevectors, Pose
import Roly: corners, faces, facevertices, nfaces, facecentroid, facenormal, polyhedron, bounding_radius

@recipe PolyformPlot (poly, ) begin
    bindingrules = nothing
    pose = nothing
    "How many cells to draw along each lattice vector, when the thing drawn is a `Tiling`."
    repeats = 2
    "How solid the repeated cells are drawn against the one cell, for a `Tiling`. `1` draws them alike."
    ghost = 0.4
    "Whether to draw a `Tiling`'s lattice vectors, from the cell they repeat."
    showlattice = true
    Makie.mixin_generic_plot_attributes()...
end

function Makie.plot!(p::PolyformPlot{<:Tuple{<:Polyform}})
    plot_polyform!(p, p.poly[], p.pose[])
    return p
end
function Makie.plot!(p::PolyformPlot{<:Tuple{<:ParticleSpecies}})
    # With no pose given, let each species method supply its own identity-pose default. The
    # rules are optional too, and decide which sites are drawn as bonding.
    pose = p.pose[]
    args = isnothing(pose) ? (p.poly[],) : (p.poly[], pose)
    plot_particlespecies!(p, args...; rules=p.bindingrules[])
    return p
end

# A tiling is drawn as a patch of itself: its cell, repeated along each lattice vector. One cell
# with `repeats = 1`, and enough to read the pattern by default.
function Makie.plot!(p::PolyformPlot{<:Tuple{<:Tiling}})
    t = p.poly[]
    cell = unitcell(t)
    vs = latticevectors(t)
    base = p.pose[]
    n = max(1, p.repeats[])
    origin = sum(q.pose.x for q in cell.particles) / nparticles(cell)
    rot = one(typeof(first(cell.particles).pose.psi))
    place(offset) = (shift = Pose(offset, rot); isnothing(base) ? shift : base * shift)

    # the cell at full strength and its translates behind it, so that the tile itself reads out
    # of the pattern rather than dissolving into it
    for m in Iterators.product(ntuple(_ -> 0:(n - 1), length(vs))...)
        offset = isempty(vs) ? zero(origin) : sum(m[i] * vs[i] for i in eachindex(vs))
        alpha = all(iszero, m) ? 1 : p.ghost[]
        plot_polyform!(p, cell, place(offset); rules=bindingrules(cell), alpha)
    end

    # and the translations themselves, drawn from the cell they carry
    if p.showlattice[] && !isempty(vs)
        anchor = isnothing(base) ? origin : (base * Pose(origin, rot)).x
        D = length(origin)
        segs = Point{D,Float32}[]
        for v in vs
            append!(segs, _openarrow(anchor, anchor + v))
        end
        linesegments!(p, segs; color=:black, linewidth=3)
    end
    return p
end

# The segments of an open arrow from `tail` to `head`: the shaft, and two barbs swept back from
# the point. Drawn rather than asked for, so that the head is a bare `>` and not a filled wedge,
# which reads as another particle among particles.
function _openarrow(tail, head; sweep=π / 7, fraction=0.18)
    D = length(head)
    u = normalize(head - tail)
    # any direction across the shaft will do to sweep the barbs into: take the axis least
    # aligned with it, so that what is left after projecting the shaft out is never small
    k = argmin(abs.(u))
    e = typeof(u)(ntuple(i -> i == k ? oneunit(eltype(u)) : zero(eltype(u)), D))
    w = normalize(e - dot(e, u) * u)
    len = fraction * norm(head - tail)
    barb(s) = head - len * (cos(sweep) * u + s * sin(sweep) * w)
    P = Point{D,Float32}
    return [P(tail), P(head), P(head), P(barb(1)), P(head), P(barb(-1))]
end

# Documented on the stub in `Roly`, so that Documenter finds it without loading this extension.
function render(p; hidedecorations=true, kwargs...)
    f = Figure()
    ax = dimension(p) == 3 ? Axis3(f[1, 1], aspect=:data) : Axis(f[1, 1], aspect=DataAspect())
    hidedecorations && hidedecorations!(ax)
    polyformplot!(ax, p; kwargs...)
    return f
end


"""
    plot_particlespecies!(ax, spcs::ParticleSpecies, pose::Pose; sitecolor, kwargs...)

Draw the particle species `spcs` at `pose` onto `ax`.

`sitecolor` is an optional callback `(speciesindex, siteindex) -> color` that controls
per-site coloring. When omitted, the species' palette is used, with sites that no bond can
use greyed out.
"""
function plot_particlespecies!(ax, spcs::ParticleSpecies, args...; kwargs...)
    error("$(typeof(spcs)) does not implement `plot_particlespecies`!")
end

"""
    particlemesh(spcs::ParticleSpecies, pose; kwargs...)

Return `(points, triangles, colors)` for `spcs` drawn at `pose`, or `nothing` for a species
with no triangle-mesh form.

This exists so that a whole polyform can go into a *single* `mesh!` call. A backend that
composites plot by plot rather than depth-testing — CairoMakie — draws whole plots in
insertion order, so one plot per particle means particles occlude each other by construction
order rather than by depth. Within one mesh CairoMakie does sort the triangles, so merging
fixes it. GLMakie has a depth buffer and is unaffected either way.
"""
particlemesh(::ParticleSpecies, args...; kwargs...) = nothing

# Triangle triples as the `ntriangles × 3` index matrix `mesh!` wants.
_facematrix(tris::Vector{NTuple{3,Int}}) = [tris[i][j] for i in eachindex(tris), j in 1:3]

"""
    plot_polyform!(ax, poly::Polyform, pose=nothing; kwargs...)

Draw all particles of `poly` onto `ax`, coloring inert sites with `INERT_COLOR`.

Species that provide a [`particlemesh`](@ref) are merged into one mesh so they depth-sort
against each other; the rest are drawn one plot per particle.

`alpha` scales the transparency of everything drawn, for showing one polyform behind another.
"""
function plot_polyform!(ax, poly::Polyform, pose=nothing; alpha=1, kwargs...)
    rules = bindingrules(poly)
    pts, tris, cols = Point3f[], NTuple{3,Int}[], RGBAf[]

    for part in poly.particles
        ps = species(rules, part.speciesindex)
        part_pose = isnothing(pose) ? part.pose : pose * part.pose
        geom = particlemesh(ps, part_pose; rules, kwargs...)
        if isnothing(geom)
            # a species with no mesh form draws itself, so the transparency has to go with it
            faded = alpha == 1 ? kwargs : (; kwargs..., alpha)
            plot_particlespecies!(ax, ps, part_pose; rules, faded...)
            continue
        end
        p, t, c = geom
        offset = length(pts)
        append!(pts, p)
        append!(cols, alpha == 1 ? c : [RGBAf(x.r, x.g, x.b, x.alpha * alpha) for x in c])
        for (a, b, d) in t
            push!(tris, (a + offset, b + offset, d + offset))
        end
    end

    isempty(tris) || mesh!(ax, pts, _facematrix(tris); color=cols, shading=NoShading)
    return
end

include("palette.jl")
include("plot_polygonparticle.jl")
include("plot_polyhedronparticle.jl")
include("plot_patchyparticle.jl")
include("plot_metaparticle.jl")

end # module
