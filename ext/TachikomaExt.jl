module TachikomaExt

using Tachikoma
using LinearAlgebra: cross, dot, normalize
using StaticArrays: SVector
using Roly: Roly, BindingRules, BindingSite, ParticleSpecies, Pose
using Roly: PolygonParticleSpecies, PolyhedronParticleSpecies
using Roly:
    bindingrules,
    bindingsite,
    bindingsites,
    bonded_colors,
    bonded_sites,
    bounding_radius,
    color,
    composition,
    corners,
    dimension,
    facenormal,
    facevertices,
    interactionmatrix,
    isaligned,
    isinert,
    istouching,
    nbonds,
    ncolors,
    nfaces,
    nparticles,
    nsites,
    nspecies,
    overlap,
    polygen,
    polyhedron,
    posetype,
    species,
    standard_twist,
    twistfreedom

import Tachikoma: view, update!, should_quit
import Roly: ruleeditor

### Palette

# Species hues, the same ones the Makie extension identifies a species by (`species_basecolor`
# in `ext/palette.jl`), so a species looks like itself in the terminal and in a figure.
const SPECIES_HEX = (
    0x1A78C6, 0xE75451, 0xFFB12F, 0xB571C4, 0x43CFE2, 0x19B684, 0xFF8835, 0x7C85FF, 0xB9D63E, 0xFF5C83, 0x21AB53
)

# Ends of the ramp a species' sites are shaded along.
const PALE = ColorRGB(0xFF, 0xFF, 0xFF)
const DEEP = ColorRGB(0x1C, 0x1C, 0x1C)

hexrgb(h::UInt32) = ColorRGB(UInt8((h >> 16) & 0xFF), UInt8((h >> 8) & 0xFF), UInt8(h & 0xFF))
rgb256(c::ColorRGB) = hex_to_color256(UInt32(c.r) << 16 | UInt32(c.g) << 8 | UInt32(c.b))

speciesrgb(i::Integer) = hexrgb(SPECIES_HEX[mod1(i, length(SPECIES_HEX))])
speciescolor(i::Integer) = rgb256(speciesrgb(i))
speciesstyle(i::Integer; kwargs...) = Style(; fg=speciescolor(i), kwargs...)

"""
    siterbg(spidx, k, n)

Return the color of site `k` of `n` on species `spidx`: a shade of that species' own hue.

A site's color has to say which particle it belongs to. Indexing one flat palette by the site's
color number instead, as this used to, gave species 2's first site the hue of species 5, which
told you nothing.
"""
function siterbg(spidx::Integer, k::Integer, n::Integer)
    base = speciesrgb(spidx)
    t = n <= 1 ? 0.5 : (k - 1) / (n - 1)
    return color_lerp(color_lerp(base, PALE, 0.45), color_lerp(base, DEEP, 0.4), t)
end

"""
    bondrgb(a, b)

Return the color of a bond between two sites: their colors mixed. A bond between two species
then reads as a blend of the two hues, and one within a species keeps that hue.
"""
bondrgb(a::ColorRGB, b::ColorRGB) = color_lerp(a, b, 0.5)

# Color indices are drawn as single characters so that a site label fits in one cell.
colorlabel(c::Integer) = c <= 9 ? Char('0' + c) : Char('a' + mod(c - 10, 26))

### Camera

# The isometric viewpoint, looking down on the origin from the direction (1, 1, 1), and how far
# the elevation may ever be driven from the horizon.
const ISO_AZIMUTH = π / 4
const ISO_ELEVATION = atan(1 / sqrt(2))
const MAX_ELEVATION = 85π / 180
const WORLD_UP = SVector(0.0, 0.0, 1.0)

"""
    Camera(azimuth, elevation)

Return the orthographic camera that looks at the origin from the direction given in spherical
coordinates: `azimuth` measured in the world's xy plane from the x axis, `elevation` above that
plane. `view` is the unit vector from the scene toward the camera, and `right` and `up` are the
axes of the projection plane, pointing right and up on screen.

An orthographic projection preserves lengths perpendicular to the view, so a particle's bounding
radius still bounds it on screen and the camera can be fitted with the same arithmetic in 2D and
3D. Parallel projection is also what makes depth sorting enough: nothing changes size with
distance, so a nearer convex particle simply covers what is behind it.
"""
struct Camera
    view::SVector{3,Float64}
    right::SVector{3,Float64}
    up::SVector{3,Float64}
end

function Camera(azimuth::Real, elevation::Real)
    v = SVector(cos(elevation) * cos(azimuth), cos(elevation) * sin(azimuth), sin(elevation))
    r = normalize(cross(v, WORLD_UP))
    return Camera(v, r, cross(r, v))
end

const ISOCAM = Camera(ISO_AZIMUTH, ISO_ELEVATION)

"""
    plane(cam, p)

Return the coordinates of a point in the camera's projection plane. A 2D point is already in it,
so the camera is ignored and the drawing of a 2D species is unchanged.
"""
plane(::Camera, p::SVector{2}) = (p[1], p[2])
plane(cam::Camera, p::SVector{3}) = (dot(p, cam.right), dot(p, cam.up))

"""
    facing(cam, site)

Whether a binding site's outward normal, which is the local x axis of its pose, points toward
the camera. Always true in 2D, where nothing is hidden.
"""
facing(::Camera, ::BindingSite{<:Pose{2}}) = true
facing(cam::Camera, s::BindingSite{<:Pose{3}}) = dot(s.pose.psi * SVector(1.0, 0.0, 0.0), cam.view) > 0

"""
    interiorbonds(spcs)

Whether a bond between two placed particles is worth labelling on the drawing.

In 2D it is, the two sites meeting along an edge of the outline. In 3D the two faces meet inside
the solid, where nothing of them shows, so a label there would sit on whichever particle happens
to stand in front of it and say nothing about either.
"""
interiorbonds(spcs::ParticleSpecies) = dimension(spcs) == 2

### Placements

# A placement is a species index paired with an absolute pose. The lattice the previous editor
# used is gone: poses come from seating one binding site against another, which is what
# `Roly.raise!` does when it grows a polyform.
const Placement{P} = Tuple{Int,P}

function absolutesites(placements, species)
    return [[pose * s for s in bindingsites(species[i])] for (i, pose) in placements]
end

"""
    seatpose(anchor, mate, t)

Return the pose that puts binding site `mate`, given in its species' body frame, against the
already placed site `anchor`, in twist `t` of the bond.
"""
function seatpose(anchor::BindingSite, mate::BindingSite, t::Integer)
    return standard_twist(anchor, t, twistfreedom(anchor, mate)) * inv(mate.pose)
end

"""
    freesites(placements, species, cam)

Return `(placement index, site index)` for every site that no other placement is touching,
i.e. every site still available as an attachment anchor, ordered by [`sortfree`](@ref).
"""
function freesites(placements, species, cam::Camera)
    sites = absolutesites(placements, species)
    free = Tuple{Int,Int}[]
    for i in eachindex(placements), k in eachindex(sites[i])
        taken = false
        for j in eachindex(placements)
            j == i && continue
            if any(s2 -> istouching(sites[i][k], s2), sites[j])
                taken = true
                break
            end
        end
        taken || push!(free, (i, k))
    end
    return sortfree(free, placements, species, cam)
end

"""
    sortfree(free, placements, species, cam)

Return the free-site list ordered clockwise around the structure's center as the camera sees it,
so that stepping through it walks the perimeter instead of jumping between particles in the
order they happened to be placed. Sites at the same bearing are ordered outermost first, which
keeps the walk moving outward before it turns.
"""
function sortfree(free, placements, species, cam::Camera)
    isempty(free) && return free
    sites = absolutesites(placements, species)
    cu, cv = plane(cam, sum(pose.x for (_, pose) in placements) / length(placements))
    # Descending bearing puts increasing indices clockwise on screen, since the projection puts
    # the plane's second axis up in the drawing.
    return sort!(free; by=((i, k),) -> begin
        u, v = plane(cam, sites[i][k].pose.x)
        du, dv = u - cu, v - cv
        (-atan(dv, du), -hypot(du, dv))
    end)
end

"""
    contacts(placements, species)

Return `(position, species1, site1, species2, site2)` for every pair of touching sites between
two placements, so that a bond can be labelled on the drawing with the two colors it joins.
"""
function contacts(placements, species)
    out = Tuple{typeof(first(placements)[2].x),Int,Int,Int,Int}[]
    length(placements) > 1 || return out
    sites = absolutesites(placements, species)
    for i in eachindex(placements), j in eachindex(placements)
        j > i || continue
        for (k1, s1) in enumerate(sites[i]), (k2, s2) in enumerate(sites[j])
            istouching(s1, s2) && push!(out, (s1.pose.x, placements[i][1], k1, placements[j][1], k2))
        end
    end
    return out
end

### Bond inference
# Attachment decides where particles go; the rules are read off every pair of placements
# afterwards. Recording a bond at attachment time instead would miss the contacts a placement
# makes with particles it was not attached to, which are exactly the ones worth discovering.

function inferred_bonds(placements, species)
    sites = absolutesites(placements, species)
    bondpairs = Set{Tuple{Int,Int}}()
    for i in eachindex(placements), j in eachindex(placements)
        j > i || continue
        for s1 in sites[i], s2 in sites[j]
            if istouching(s1, s2) && isaligned(s1, s2)
                push!(bondpairs, minmax(color(s1), color(s2)))
            end
        end
    end
    return sort!(collect(bondpairs))
end

function bonds_matrix(placements, species)
    sites = absolutesites(placements, species)
    seen = Set{NTuple{4,Int}}()
    rows = NTuple{4,Int}[]
    for i in eachindex(placements), j in eachindex(placements)
        j > i || continue
        sp_i = placements[i][1]
        sp_j = placements[j][1]
        for (k1, s1) in enumerate(sites[i]), (k2, s2) in enumerate(sites[j])
            if istouching(s1, s2) && isaligned(s1, s2)
                a = (sp_i, k1, sp_j, k2)
                b = (sp_j, k2, sp_i, k1)
                if a ∉ seen && b ∉ seen
                    push!(rows, a)
                    push!(seen, a)
                    push!(seen, b)
                end
            end
        end
    end
    return isempty(rows) ? zeros(Int, 0, 4) : reduce(vcat, [collect(r)' for r in rows])
end

"""
    usedspecies(placements)

Return the species indices that carry at least one placement, sorted.
"""
usedspecies(placements) = sort!(collect(Set(i for (i, _) in placements)))

"""
    remappedbonds(placements, species)

Return the `nx4` bonds matrix the placements imply, with species renumbered to drop any that
were never placed.
"""
function remappedbonds(placements, species)
    bonds = bonds_matrix(placements, species)
    remap = Dict(sp => i for (i, sp) in enumerate(usedspecies(placements)))
    for row in axes(bonds, 1)
        bonds[row, 1] = remap[bonds[row, 1]]
        bonds[row, 3] = remap[bonds[row, 3]]
    end
    return bonds
end

"""
    buildrules(placements, species)

Return the `BindingRules` the placements imply on their own, or `nothing` if no bond has been
formed. This is the geometry alone; the editor's own rules come from [`buildrules(::EditorModel)`](@ref),
which lays the hand-edited bonds over this.
"""
function buildrules(placements, species)
    isempty(placements) && return nothing
    bonds = remappedbonds(placements, species)
    isempty(bonds) && return nothing
    kept = [species[sp] for sp in usedspecies(placements)]
    return length(kept) == 1 ? BindingRules(bonds, kept[1]) : BindingRules(bonds, kept)
end

"""
    bondsfromrules(rules)

Return the `nx4` bonds matrix that reproduces `rules`, one representative site pair per bonding
color pair. Used for `output=:bonds`, which has to reflect bonds added by hand as well as the
ones read off the geometry.
"""
function bondsfromrules(rules::BindingRules)
    rows = NTuple{4,Int}[]
    for (locs1, locs2) in bonded_sites(rules)
        (sp1, k1) = first(locs1)
        (sp2, k2) = first(locs2)
        push!(rows, (sp1, k1, sp2, k2))
    end
    return isempty(rows) ? zeros(Int, 0, 4) : reduce(vcat, [collect(r)' for r in rows])
end

### Geometry for drawing

"""
    outline(spcs, pose)

Return the world-space corners of the species' boundary, closed, at the given pose. Species
that expose no corners are drawn as their bounding circle.
"""
outline(spcs::PolygonParticleSpecies, pose::Pose) = [pose * c for c in spcs.corners]
function outline(spcs::ParticleSpecies{2}, pose::Pose)
    r = bounding_radius(spcs)
    return [pose * SVector(r * cos(2π * k / 24), r * sin(2π * k / 24)) for k in 0:23]
end

"""
    World(scale, cx, cy, w, h, cam=ISOCAM)

Map from world coordinates to braille dot coordinates on a canvas `w` cells wide and `h` cells
tall. `cam` projects the point onto its plane first, which does nothing in 2D. `scale` is dots
per world unit and `(cx, cy)` is the point of the projection plane placed at the center. Braille
dots are 2 per cell across and 4 down, which on a terminal cell roughly twice as tall as it is
wide makes them square, so one scale serves both axes.
"""
struct World
    scale::Float64
    cx::Float64
    cy::Float64
    w::Int
    h::Int
    cam::Camera
end

World(scale::Real, cx::Real, cy::Real, w::Int, h::Int) = World(scale, cx, cy, w, h, ISOCAM)

function todotsf(v::World, p)
    u, w = plane(v.cam, p)
    return ((u - v.cx) * v.scale + v.w, (v.cy - w) * v.scale + 2v.h)
end

todots(v::World, p) = round.(Int, todotsf(v, p))

"""
    tocells(v, rect, p)

Return the position of a point in fractional cell coordinates within `rect`, where cell `(i, j)`
covers `[i, i+1) x [j, j+1)`. Used by the filled 3D drawing, which paints whole cells rather
than braille dots.
"""
function tocells(v::World, rect::Rect, p)
    dx, dy = todotsf(v, p)
    return (rect.x + dx / 2, rect.y + dy / 4)
end

inrect(r::Rect, col::Int, row::Int) = r.x <= col <= right(r) && r.y <= row <= bottom(r)

"""
    worldbox(placements, species, cam)

Return `(xmin, xmax, ymin, ymax)`, the extent the placements cover in the camera's projection
plane, including each particle's bounding radius.
"""
function worldbox(placements, species, cam::Camera)
    xs = Float64[]
    ys = Float64[]
    for (i, pose) in placements
        r = bounding_radius(species[i])
        u, w = plane(cam, pose.x)
        push!(xs, u - r, u + r)
        push!(ys, w - r, w + r)
    end
    return (minimum(xs), maximum(xs), minimum(ys), maximum(ys))
end

"""
    boxscale(box, w, h)

Return the largest scale at which `box` still fits a `w` by `h` cell canvas, with a cell of
margin on each side.
"""
function boxscale((xmin, xmax, ymin, ymax), w::Int, h::Int)
    spanx = max(xmax - xmin, eps())
    spany = max(ymax - ymin, eps())
    return min((2w - 4) / spanx, (4h - 4) / spany)
end

"""
    fitworld(placements, species, w, h; cam=ISOCAM)

Return the `World` that fits every placement into a `w` by `h` cell canvas. Used for the rules
gallery and the enumeration preview, where each structure is drawn on its own and should fill
the space it is given, from the fixed isometric viewpoint in 3D.
"""
function fitworld(placements, species, w::Int, h::Int; cam::Camera=ISOCAM)
    (isempty(placements) || w < 1 || h < 1) && return World(4.0, 0.0, 0.0, max(w, 1), max(h, 1), cam)
    box = worldbox(placements, species, cam)
    cx = (box[1] + box[2]) / 2
    cy = (box[3] + box[4]) / 2
    return World(boxscale(box, w, h), cx, cy, w, h, cam)
end

function drawoutline!(canvas::Canvas, v::World, pts)
    n = length(pts)
    for k in 1:n
        x0, y0 = todots(v, pts[k])
        x1, y1 = todots(v, pts[mod1(k + 1, n)])
        line!(canvas, x0, y0, x1, y1)
    end
    return canvas
end

# A particle carries a filled triangle at its center pointing away from site 1, which is what
# tells two otherwise identical particles apart when they are turned differently. Below this many
# dots it would be a smudge rather than a direction, so it is left off.
const ARROW_MIN_DOTS = 4

"""
    drawreject!(canvas, v, pts)

Draw a cross through the center of a shape, marking a placement the editor will not accept.

A cross rather than the heading triangle, since a refused placement has no orientation worth
reading, and drawn by its caller in a neutral gray rather than red: red is one of the species
colors, so a refused ghost drawn in it reads as a different species instead of as a refusal.
"""
function drawreject!(canvas::Canvas, v::World, pts)
    dots = [todots(v, p) for p in pts]
    xs, ys = first.(dots), last.(dots)
    cx = round(Int, (minimum(xs) + maximum(xs)) / 2)
    cy = round(Int, (minimum(ys) + maximum(ys)) / 2)
    arm = round(Int, clamp(0.4 * min(maximum(xs) - minimum(xs), maximum(ys) - minimum(ys)) / 2, 2.0, 8.0))
    line!(canvas, cx - arm, cy - arm, cx + arm, cy + arm)
    line!(canvas, cx - arm, cy + arm, cx + arm, cy - arm)
    return canvas
end

"""
    drawheading!(canvas, v, pose, spcs, pts)

Draw a small filled triangle at a particle's center, pointing away from its first binding site.

A triangle rather than an arrow: a shaft with a head needs more dots than a particle at this
scale has to give, whereas a triangle is a direction and a marker at once, and being symmetric
about its axis it sits centered without any of the half-cell correction a stroke needs.

`pts` is the outline already computed, so the triangle centers on the shape as drawn rather than
on the pose, which rounding can leave half a cell away from it.
"""
function drawheading!(canvas::Canvas, v::World, pose::Pose, spcs, pts)
    nsites(spcs) < 1 && return canvas
    dots = [todots(v, p) for p in pts]
    cx = (minimum(first, dots) + maximum(first, dots)) / 2
    cy = (minimum(last, dots) + maximum(last, dots)) / 2
    sx, sy = todots(v, (pose * bindingsite(spcs, 1)).pose.x)
    dx, dy = cx - sx, cy - sy
    len = hypot(dx, dy)
    len < ARROW_MIN_DOTS && return canvas

    ux, uy = dx / len, dy / len
    height = clamp(len * 0.8, 3.0, 12.0)
    halfbase = height * 0.45
    apex = (cx + ux * height * 0.55, cy + uy * height * 0.55)
    back = (cx - ux * height * 0.45, cy - uy * height * 0.45)
    # Fill by fanning lines from the apex across the base, which is cheap and leaves no gaps at
    # the sizes involved.
    steps = max(2, ceil(Int, 2 * halfbase))
    for t in 0:steps
        f = 2 * (t / steps) - 1
        bx = round(Int, back[1] - uy * halfbase * f)
        by = round(Int, back[2] + ux * halfbase * f)
        line!(canvas, round(Int, apex[1]), round(Int, apex[2]), bx, by)
    end
    return canvas
end

### Solid drawing

# The gray a refused placement is filled in, and how far a pending one is faded toward the
# background so that it reads as provisional without changing hue.
const REJECT_GRAY = ColorRGB(0x78, 0x78, 0x80)
const GHOST_FADE = 0.45

# Ramps are cached because finding the nearest entry of the 256-color cube searches all of it,
# and every face of every particle asks for one on every frame.
const RAMPS = Dict{ColorRGB,NTuple{3,Color256}}()

"""
    faceramp(base)

Return three shades of `base`, dark to pale, that stay distinct after quantization to the
256-color cube.

The cube has six levels per channel, so a plausible-looking lighting model produces shades that
collapse onto the same entry and a cube comes out in two colors instead of three. The spread is
therefore widened until the three codes actually differ.
"""
function faceramp(base::ColorRGB)
    return get!(RAMPS, base) do
        for spread in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
            shades = (rgb256(color_lerp(base, DEEP, spread)), rgb256(base), rgb256(color_lerp(base, PALE, spread)))
            length(unique(c.code for c in shades)) == 3 && return shades
        end
        return (rgb256(color_lerp(base, DEEP, 0.9)), rgb256(base), rgb256(color_lerp(base, PALE, 0.9)))
    end
end

"""
    faceshade(cam, nrm)

Return which of the three shades a face with outward normal `nrm` takes: the palest for a face
turned upward, and otherwise one of two, split by which side of the screen it looks toward.

Three buckets rather than a continuous lighting term, which is both what the drawing style asks
for and what survives quantization, a smooth ramp having put a cube's top and one of its sides
on the same entry of the color cube.
"""
function faceshade(cam::Camera, nrm)
    dot(nrm, WORLD_UP) > 0.5 && return 3
    return dot(nrm, cam.right) >= 0 ? 2 : 1
end

"""
    visiblefaces(spcs, pose, cam)

Return the polygons of a particle's surface that face the camera, each paired with the index of
the shade it is drawn in. Faces turned away are dropped, which for a convex body leaves exactly
the silhouette and leaves nothing to sort among the faces that remain.

A 3D species with no polyhedron behind it is returned as the silhouette of its bounding sphere,
in the same way the 2D drawing falls back to a circle.
"""
function visiblefaces(spcs::PolyhedronParticleSpecies, pose::Pose, cam::Camera)
    body = polyhedron(spcs)
    cs = [pose * c for c in corners(body)]
    out = Tuple{Vector{eltype(cs)},Int}[]
    for i in 1:nfaces(body)
        nrm = pose.psi * facenormal(body, i)
        dot(nrm, cam.view) > 1e-9 || continue
        push!(out, ([cs[k] for k in facevertices(body, i)], faceshade(cam, nrm)))
    end
    return out
end

function visiblefaces(spcs::ParticleSpecies{3}, pose::Pose, cam::Camera)
    r = bounding_radius(spcs)
    disk = [pose.x + r * (cos(2π * k / 24) * cam.right + sin(2π * k / 24) * cam.up) for k in 0:23]
    return [(disk, 2)]
end

"""
    fillpolygon!(paint, pts)

Call `paint(column, row)` for every cell whose center falls inside the polygon `pts`, given in
fractional cell coordinates.

A scanline fill: a horizontal line through the polygon crosses its edges an even number of
times, so the spans between successive crossings are its interior.
"""
function fillpolygon!(paint, pts)
    n = length(pts)
    n > 2 || return nothing
    ys = last.(pts)
    xs = Float64[]
    for row in floor(Int, minimum(ys)):floor(Int, maximum(ys))
        yc = row + 0.5
        empty!(xs)
        for i in 1:n
            x1, y1 = pts[i]
            x2, y2 = pts[mod1(i + 1, n)]
            (y1 <= yc < y2 || y2 <= yc < y1) || continue
            push!(xs, x1 + (yc - y1) / (y2 - y1) * (x2 - x1))
        end
        sort!(xs)
        for k in 1:2:(length(xs) - 1), col in ceil(Int, xs[k] - 0.5):floor(Int, xs[k + 1] - 0.5)
            paint(col, row)
        end
    end
    return nothing
end

"""
    crossout!(buf, rect, cells)

Draw a cross over the extent a particle covers, marking a placement the editor will not accept.

The 3D counterpart of [`drawreject!`](@ref). It writes characters over the fill rather than
drawing into a canvas, a filled particle leaving no free braille dots to draw into.
"""
function crossout!(buf::Buffer, rect::Rect, cells)
    isempty(cells) && return buf
    x0, x1 = extrema(floor(Int, p[1]) for p in cells)
    y0, y1 = extrema(floor(Int, p[2]) for p in cells)
    steps = max(x1 - x0, y1 - y0, 1)
    style = Style(; fg=rgb256(PALE), bold=true)
    for t in 0:steps
        f = t / steps
        col = round(Int, x0 + f * (x1 - x0))
        for row in (round(Int, y0 + f * (y1 - y0)), round(Int, y1 - f * (y1 - y0)))
            inrect(rect, col, row) && set_char!(buf, col, row, '╳', style)
        end
    end
    return buf
end

"""
    drawparticles!(buf, rect, v, placements, species; ghost=nothing, blocked=false, wire=false)

Draw every placement into `rect`, followed by the pending particle when `ghost` is given as a
`(species index, pose)` pair. `blocked` says that the pending particle is one the editor will
refuse.

2D particles are drawn as outlines on a braille canvas, one canvas per species so that each
keeps its own color. 3D particles are drawn as filled faces shaded by their normal, far to near,
so that a nearer particle covers what stands behind it.

`wire` asks for the 3D drawing to be an outline instead, the edges of the faces that turn toward
the camera, on the same braille canvas the 2D drawing uses. A filled face is a whole cell, four
times the height and twice the width of a braille dot, so it needs a large pane to read; the
small drawings, the rules gallery and the enumeration thumbnails, take the wireframe.
"""
function drawparticles!(
    buf::Buffer, rect::Rect, v::World, placements, spcs; ghost=nothing, blocked::Bool=false, wire::Bool=false
)
    (isempty(placements) && ghost === nothing) && return buf
    P = ghost === nothing ? typeof(placements[1][2]) : typeof(ghost[2])
    return drawparticles!(Val(dimension(P)), buf, rect, v, placements, spcs; ghost, blocked, wire)
end

function drawparticles!(
    ::Val{2}, buf::Buffer, rect::Rect, v::World, placements, spcs; ghost=nothing, blocked::Bool=false, wire::Bool=false
)
    canvases = Dict{Int,Canvas}()
    for (i, pose) in placements
        cv = get!(() -> Canvas(rect.width, rect.height; style=speciesstyle(i)), canvases, i)
        pts = outline(spcs[i], pose)
        drawoutline!(cv, v, pts)
        drawheading!(cv, v, pose, spcs[i], pts)
    end
    for (_, cv) in sort!(collect(canvases); by=first)
        render(cv, rect, buf)
    end
    if ghost !== nothing
        i, pose = ghost
        cv = Canvas(rect.width, rect.height; style=blocked ? tstyle(:text_dim) : speciesstyle(i; dim=true))
        pts = outline(spcs[i], pose)
        drawoutline!(cv, v, pts)
        blocked ? drawreject!(cv, v, pts) : drawheading!(cv, v, pose, spcs[i], pts)
        render(cv, rect, buf)
    end
    return buf
end

function drawparticles!(
    ::Val{3}, buf::Buffer, rect::Rect, v::World, placements, spcs; ghost=nothing, blocked::Bool=false, wire::Bool=false
)
    wire && return drawwires!(buf, rect, v, placements, spcs; ghost, blocked)
    idxs = [i for (i, _) in placements]
    poses = [pose for (_, pose) in placements]
    pending = falses(length(poses))
    if ghost !== nothing
        push!(idxs, ghost[1])
        push!(poses, ghost[2])
        push!(pending, true)
    end
    # Painter's algorithm: the farthest particle first. Within a convex particle the faces turned
    # away are dropped outright, so those that survive never overlap and need no ordering.
    for n in sortperm([dot(pose.x, v.cam.view) for pose in poses])
        base = speciesrgb(idxs[n])
        ramp = if !pending[n]
            faceramp(base)
        else
            faceramp(blocked ? REJECT_GRAY : color_lerp(base, DEEP, GHOST_FADE))
        end
        seen = Tuple{Float64,Float64}[]
        for (poly, shade) in visiblefaces(spcs[idxs[n]], poses[n], v.cam)
            cells = [tocells(v, rect, p) for p in poly]
            append!(seen, cells)
            style = Style(; bg=ramp[shade])
            fillpolygon!(cells) do col, row
                inrect(rect, col, row) && set_char!(buf, col, row, ' ', style)
                return nothing
            end
        end
        pending[n] && blocked && crossout!(buf, rect, seen)
    end
    return buf
end

"""
    drawwires!(buf, rect, v, placements, species; ghost=nothing, blocked=false)

Draw 3D particles as braille outlines: the edges of the faces that turn toward the camera.

Hidden edges go in two steps. Within a particle, dropping the faces turned away removes exactly
the edges on its far side, the body being convex. Between particles, the drawing runs farthest
first and each particle erases the dots its silhouette covers from everything already drawn
behind it, so a nearer particle hides a farther one as it would if both were filled.

One canvas per particle, since the erasing has to happen in depth order and a canvas carries a
single color. They are composited at the end, in the same order.
"""
function drawwires!(buf::Buffer, rect::Rect, v::World, placements, spcs; ghost=nothing, blocked::Bool=false)
    idxs = [i for (i, _) in placements]
    poses = [pose for (_, pose) in placements]
    pending = falses(length(poses))
    if ghost !== nothing
        push!(idxs, ghost[1])
        push!(poses, ghost[2])
        push!(pending, true)
    end

    drawn = Canvas[]
    for n in sortperm([dot(pose.x, v.cam.view) for pose in poses])
        style = if !pending[n]
            speciesstyle(idxs[n])
        elseif blocked
            tstyle(:text_dim)
        else
            speciesstyle(idxs[n]; dim=true)
        end
        cv = Canvas(rect.width, rect.height; style)
        faces = visiblefaces(spcs[idxs[n]], poses[n], v.cam)
        # A dot is at the center of its own cell of the dot grid, so the polygon is offset by
        # half a dot to line the scanline's sampling up with where `todots` rounds to.
        for (poly, _) in faces
            fillpolygon!([todotsf(v, p) .+ 0.5 for p in poly]) do dx, dy
                for behind in drawn
                    unset_point!(behind, dx, dy)
                end
                return nothing
            end
        end
        for (poly, _) in faces
            drawoutline!(cv, v, poly)
        end
        pending[n] && blocked && drawreject!(cv, v, [p for (poly, _) in faces for p in poly])
        push!(drawn, cv)
    end
    for cv in drawn
        render(cv, rect, buf)
    end
    return buf
end

### Model

const SIDEBAR_W = 22
const RULES_W = 40
const PANE_W = 44
const PREVIEW_MIN_W = 24
const GALLERY_ROW = 6
const GALLERY_MIN_ROW = 4
const GALLERY_MIN_W = 11
const PAIR_W = 7

# A bond cell is split between the colors of the two sites it joins, using the upper half block:
# the row's color above, the column's below.
bondglyph() = '▀'

# How many braille dots wide one particle is when the editor opens. Small enough that the first
# several attachments land inside the view already on screen, which keeps the camera still. A 3D
# particle needs more of them: the projection puts three faces' worth of edges inside the
# silhouette that a 2D outline covers with one.
const DOTS_PER_PARTICLE = 14
const DOTS_PER_PARTICLE_3D = 20

# What `-` and `=` do to the scale, and how far either may take it from the opening view.
const ZOOM_STEP = 1.3
const ZOOM_RANGE = 16.0

# What the camera keys turn by, how squarely a site has to face the camera for the view to be
# left where it is, and how squarely the view puts a site it does move to bring forward. The
# second is the larger of the two: turning a little past the threshold means a site that has
# only just come into view does not set the camera swinging again on the next keystroke.
const TURN_STEP = π / 12
const SITE_HIDDEN = 0.2
const SITE_MARGIN = 0.3

# Bounds on the enumeration the preview runs. `maxstrs` is what actually keeps it cheap: the
# search stops after that many structures however permissive the rules are.
const MAXSIZE_RANGE = 2:12
const MAXSTRS_RANGE = (25, 50, 100, 200, 500, 1000, 2000, 5000)
const THUMB_W = 13
const THUMB_H = 6
const DETAIL_W = 24

@kwdef mutable struct EditorModel{S,P} <: Model
    base::S
    species::Vector{S}
    placements::Vector{Placement{P}}
    free::Vector{Tuple{Int,Int}} = Tuple{Int,Int}[]
    rules::Union{Nothing,BindingRules} = nothing
    anchor::Int = 1
    incoming::Int = 1
    twist::Int = 0
    active_species::Int = 1
    # Which pane the arrow keys and enter act on.
    focus::Symbol = :rules
    # The rules pane's cursor: the color pair it is pointing at.
    pair::Tuple{Int,Int} = (1, 1)
    # Bonds turned on or off by hand, laid over the ones the geometry implies. Keyed by ordered
    # color pair.
    overrides::Dict{Tuple{Int,Int},Bool} = Dict{Tuple{Int,Int},Bool}()
    # One color per site, shaded from its species' hue. Cached because finding the nearest
    # terminal color searches the whole 256-color cube.
    sitergbs::Vector{ColorRGB} = ColorRGB[]
    palette::Vector{Color256} = Color256[]
    # The structures the last enumeration found, which one is selected, and what the run was
    # bounded by. `stale` marks the rules as having changed since that run.
    polyforms::Vector{Any} = Any[]
    enuminfo::String = ""
    selected::Int = 1
    maxsize::Int = 6
    maxstrs::Int = 200
    stale::Bool = true
    # How many thumbnails the last frame fitted across, so that up and down move a whole row.
    gridcols::Int = 1
    # The construction pane starts hidden: the rules are the deliverable, and placing blocks is
    # one way of arriving at them rather than the only one.
    showconstruction::Bool = false
    message::String = ""
    messagekind::Symbol = :info
    quit::Bool = false
    # Camera, held across frames so that attaching a particle does not move the structure that
    # is already on screen. See `worldfor!`.
    scale::Float64 = 1.0
    basescale::Float64 = 1.0
    cx::Float64 = 0.0
    cy::Float64 = 0.0
    refit::Bool = false
    manualzoom::Bool = false
    # Where the camera stands in 3D. Ignored in 2D, where the drawing plane is the world.
    azimuth::Float64 = ISO_AZIMUTH
    elevation::Float64 = ISO_ELEVATION
end

should_quit(m::EditorModel) = m.quit

camera(m::EditorModel) = Camera(m.azimuth, m.elevation)

"""
    poseof(m)

Return the pose type the model's placements are stored in.
"""
poseof(::EditorModel{S,P}) where {S,P} = P

function EditorModel(spcs::S) where {S<:ParticleSpecies}
    P = posetype(spcs)
    dots = dimension(spcs) == 3 ? DOTS_PER_PARTICLE_3D : DOTS_PER_PARTICLE
    opening = dots / (2 * bounding_radius(spcs))
    m = EditorModel{S,P}(;
        base=spcs, species=S[copy(spcs)], placements=Placement{P}[(1, one(P))], scale=opening, basescale=opening
    )
    refresh!(m)
    resetcamera!(m)
    return m
end

"""
    zoom!(m, factor)

Scale the view by `factor`, holding the cursor still on screen so that zooming in follows the
site being worked on rather than the middle of the structure.

Zooming by hand switches the automatic fitting off, since a view the user chose should not be
undone by the next attachment. `0` gives it back.
"""
function zoom!(m::EditorModel, factor::Real)
    new = clamp(m.scale * factor, m.basescale / ZOOM_RANGE, m.basescale * ZOOM_RANGE)
    new == m.scale && return m
    p = anchorposition(m)
    u, w = p === nothing ? (m.cx, m.cy) : plane(camera(m), p)
    m.cx = u - (u - m.cx) * m.scale / new
    m.cy = w - (w - m.cy) * m.scale / new
    m.scale = new
    m.manualzoom = true
    return m
end

"""
    worldfor!(m, placements, w, h)

Return the `World` to draw the construction pane with, adjusting the camera only when it has to.

The camera holds its scale and center between frames, so attaching a particle leaves everything
already on screen where it was. It recenters and zooms out when the structure would otherwise
leave the view, and it never zooms back in on its own: a structure that has been trimmed stays
at the scale it reached, until `0` asks for a refit. Once the view has been zoomed by hand it
stops fitting at all, again until `0`.
"""
function worldfor!(m::EditorModel, placements, w::Int, h::Int)
    cam = camera(m)
    (w < 1 || h < 1) && return World(m.scale, m.cx, m.cy, max(w, 1), max(h, 1), cam)
    isempty(placements) && return World(m.scale, m.cx, m.cy, w, h, cam)
    box = worldbox(placements, m.species, cam)
    xmin, xmax, ymin, ymax = box
    if m.refit
        m.cx = (xmin + xmax) / 2
        m.cy = (ymin + ymax) / 2
        m.scale = boxscale(box, w, h)
        m.manualzoom = false
        m.refit = false
    elseif !m.manualzoom
        halfw = w - 2.0
        halfh = 2h - 2.0
        fits =
            (xmax - m.cx) * m.scale <= halfw &&
            (m.cx - xmin) * m.scale <= halfw &&
            (ymax - m.cy) * m.scale <= halfh &&
            (m.cy - ymin) * m.scale <= halfh
        if !fits
            m.cx = (xmin + xmax) / 2
            m.cy = (ymin + ymax) / 2
            m.scale = min(m.scale, boxscale(box, w, h))
        end
    end
    return World(m.scale, m.cx, m.cy, w, h, cam)
end

"""
    recenter!(m)

Put the middle of the structure back at the middle of the view, keeping the scale.

Called after the camera turns, which redefines the projection plane the center is held in:
carrying the old center over unchanged would slide the structure off the screen.
"""
function recenter!(m::EditorModel)
    isempty(m.placements) && return m
    box = worldbox(m.placements, m.species, camera(m))
    m.cx = (box[1] + box[2]) / 2
    m.cy = (box[3] + box[4]) / 2
    return m
end

"""
    perimetersort!(m)

Reorder the free-site list clockwise as the camera now sees it, keeping the cursor on the site
it was already on. The order is what the arrow keys walk, so in 3D it has to follow the camera.
"""
function perimetersort!(m::EditorModel)
    isempty(m.free) && return m
    at = m.free[m.anchor]
    m.free = sortfree(m.free, m.placements, m.species, camera(m))
    m.anchor = something(findfirst(==(at), m.free), 1)
    return m
end

"""
    pivot!(m)

Turn the camera until the site the cursor is on faces it, and no further. Does nothing in 2D,
where no site is ever hidden.

The cursor is what the user drives and the camera follows it, so stepping onto a face on the far
side brings that face round rather than leaving the cursor somewhere invisible. The turn is the
smallest one that clears the face of the silhouette, taken in the plane the view direction and
the face normal span, which keeps the structure recognizable across the move.
"""
function pivot!(m::EditorModel)
    dimension(m.base) == 3 || return m
    s = anchorsite(m)
    s === nothing && return m
    cam = camera(m)
    nrm = SVector{3,Float64}(s.pose.psi * SVector(1.0, 0.0, 0.0))
    dot(nrm, cam.view) >= SITE_HIDDEN && return m

    # A view exactly opposite the normal spans no plane with it, so it is nudged sideways first
    # and the turn goes whichever way that nudge points.
    v = dot(nrm, cam.view) < -1 + 1e-6 ? normalize(cam.view + 1e-3 * cam.right) : cam.view
    w = normalize(v - nrm * dot(nrm, v))
    θ = acos(SITE_MARGIN)
    target = nrm * cos(θ) + w * sin(θ)

    m.elevation = clamp(asin(clamp(target[3], -1.0, 1.0)), -MAX_ELEVATION, MAX_ELEVATION)
    m.azimuth = atan(target[2], target[1])
    recenter!(m)
    perimetersort!(m)
    return m
end

"""
    resetcamera!(m)

Put the camera back on the isometric viewpoint and the cursor on a site that viewpoint already
shows. Does nothing in 2D.

The editor opens this way, so that the first thing on screen is the canonical view rather than
one the cursor has already swung the camera away from.
"""
function resetcamera!(m::EditorModel)
    dimension(m.base) == 3 || return m
    m.azimuth = ISO_AZIMUTH
    m.elevation = ISO_ELEVATION
    perimetersort!(m)
    sites = absolutesites(m.placements, m.species)
    cam = camera(m)
    n = findfirst(((i, k),) -> facing(cam, sites[i][k]), m.free)
    n === nothing || (m.anchor = n)
    return recenter!(m)
end

"""
    turn!(m, dazimuth, delevation)

Turn the camera by hand, by the given increments in azimuth and elevation. Does nothing in 2D.
"""
function turn!(m::EditorModel, dazimuth::Real, delevation::Real)
    dimension(m.base) == 3 || return m
    m.azimuth += dazimuth
    m.elevation = clamp(m.elevation + delevation, -MAX_ELEVATION, MAX_ELEVATION)
    recenter!(m)
    perimetersort!(m)
    return m
end

function ensurespecies!(m::EditorModel, n::Integer)
    while length(m.species) < n
        push!(m.species, copy(m.base))
    end
    return m
end

"""
    refresh!(m)

Recompute the free-site list and the inferred rules, and clamp every selection into range.
Called after any change to the placements, rather than per frame, since inference is quadratic
in the number of placements.
"""
function refresh!(m::EditorModel; near=nothing)
    m.free = freesites(m.placements, m.species, camera(m))
    m.rules = buildrules(m)
    refreshpalette!(m)
    n = totalcolors(m)
    m.pair = (clamp(m.pair[1], 1, max(n, 1)), clamp(m.pair[2], 1, max(n, 1)))
    m.stale = true
    if isempty(m.free)
        m.anchor = 1
    elseif near === nothing
        m.anchor = mod1(m.anchor, length(m.free))
    else
        # Attaching consumes the anchor site and reorders the list, so carrying the index over
        # would put the cursor somewhere unrelated. Carrying the position over leaves it on the
        # nearest site to the one just used.
        sites = absolutesites(m.placements, m.species)
        m.anchor = argmin(n -> sum(abs2, sites[m.free[n][1]][m.free[n][2]].pose.x - near), eachindex(m.free))
    end
    m.incoming = mod1(m.incoming, nsites(m.species[m.active_species]))
    m.twist = mod(m.twist, max(ntwists(m), 1))
    pivot!(m)
    return m
end

"""
    sitecolors(m)

Return a `Dict` from `(species index, site index)` to the color that site carries in the rules
the editor will return.

`Roly.BindingRules` gives each species a contiguous color range of its own, so a species drawn
with its own site colors would disagree with the row it occupies in the interaction matrix. The
active species is included even before it is placed, so the pending particle is labelled with
the colors it will have once committed.
"""
function sitecolors(m::EditorModel)
    out = Dict{Tuple{Int,Int},Int}()
    c = 1
    for spidx in shownspecies(m)
        spcs = m.species[spidx]
        cols = [color(bindingsite(spcs, k)) for k in 1:nsites(spcs)]
        lo, hi = extrema(cols)
        for k in eachindex(cols)
            out[(spidx, k)] = cols[k] - lo + c
        end
        c += hi - lo + 1
    end
    return out
end

"""
    refreshpalette!(m)

Recompute every site's color, shading each species' hue across the sites it has.
"""
function refreshpalette!(m::EditorModel)
    gc = sitecolors(m)
    m.sitergbs = fill(speciesrgb(1), totalcolors(m))
    for spidx in shownspecies(m)
        n = nsites(m.species[spidx])
        for k in 1:n
            m.sitergbs[gc[(spidx, k)]] = siterbg(spidx, k, n)
        end
    end
    m.palette = [rgb256(c) for c in m.sitergbs]
    return m
end

"""
    sitecolor(m, c)
    colorstyle(m, c; kwargs...)
    bondstyle(m, c1, c2; kwargs...)

The color of one site, a style in it, and the style of a bond between two sites, which is the
two mixed.
"""
sitecolor(m::EditorModel, c::Integer) = 1 <= c <= length(m.palette) ? m.palette[c] : speciescolor(c)
colorstyle(m::EditorModel, c::Integer; kwargs...) = Style(; fg=sitecolor(m, c), kwargs...)

function bondstyle(m::EditorModel, c1::Integer, c2::Integer; kwargs...)
    n = length(m.sitergbs)
    (1 <= c1 <= n && 1 <= c2 <= n) || return colorstyle(m, c1; kwargs...)
    return Style(; fg=rgb256(bondrgb(m.sitergbs[c1], m.sitergbs[c2])), kwargs...)
end

"""
    enumerate!(m)

Enumerate the polyforms the current rules allow, for the preview pane.

Run only when asked, since a permissive rules set can take long enough to be felt between
keystrokes, and it is the one part of the editor whose cost is not bounded by the size of the
drawing. `m.maxstrs` is what actually bounds it: the search stops after that many structures
however many the rules admit. A rules set that fails to enumerate leaves the pane empty with a
message rather than taking the editor down.
"""
function enumerate!(m::EditorModel)
    empty!(m.polyforms)
    m.selected = 1
    m.stale = false
    rules = m.rules
    if rules === nothing || nbonds(rules) == 0
        m.enuminfo = "no bonds"
        return m
    end
    try
        polys = polygen(rules; maxsize=m.maxsize, maxstrs=m.maxstrs)
        filter!(p -> nparticles(p) > 0, polys)
        append!(m.polyforms, polys)
        capped = length(polys) >= m.maxstrs
        m.enuminfo = string(length(polys), capped ? "+" : "", " ≤ ", m.maxsize)
    catch err
        err isa InterruptException && rethrow()
        m.enuminfo = "cannot enumerate"
    end
    return m
end

"""
    shownspecies(m)

Return the species the rules cover, which is every instance the digit keys have created.

Covering only the placed ones would mean that selecting species 2 and then 3, without placing
either, dropped 2 again and left the matrix the size it was, so a species could not be given
bonds before it was used.
"""
shownspecies(m::EditorModel) = eachindex(m.species)

"""
    totalcolors(m)

Return how many colors the rules span, the sum of each shown species' color range.
"""
totalcolors(m::EditorModel) = maximum(values(sitecolors(m)); init=0)

"""
    inferredmatrix(m)

Return the interaction matrix the placements imply on their own: every pair of sites that touch
and align, in the colors the rules use.
"""
function inferredmatrix(m::EditorModel)
    n = totalcolors(m)
    mat = zeros(Bool, n, n)
    n == 0 && return mat
    gc = sitecolors(m)
    sites = absolutesites(m.placements, m.species)
    for i in eachindex(m.placements), j in eachindex(m.placements)
        j > i || continue
        for (k1, s1) in enumerate(sites[i]), (k2, s2) in enumerate(sites[j])
            if istouching(s1, s2) && isaligned(s1, s2)
                c1 = gc[(m.placements[i][1], k1)]
                c2 = gc[(m.placements[j][1], k2)]
                mat[c1, c2] = mat[c2, c1] = true
            end
        end
    end
    return mat
end

"""
    effectivematrix(m)

Return the interaction matrix the editor's rules use: what the geometry implies, with the pairs
toggled by hand laid over it.

Editing does not detach the rules from the drawing. A hand-set bond stays set as the structure
grows, and clearing one the geometry keeps producing leaves it cleared, which is what makes the
two halves of the editor usable together.
"""
function effectivematrix(m::EditorModel)
    mat = inferredmatrix(m)
    n = size(mat, 1)
    for ((c1, c2), on) in m.overrides
        (1 <= c1 <= n && 1 <= c2 <= n) || continue
        mat[c1, c2] = mat[c2, c1] = on
    end
    return mat
end

"""
    buildrules(m::EditorModel)

Return the `BindingRules` the editor currently describes, or `nothing` if nothing is placed.

Unlike the geometry-only form these exist as soon as a particle does, with no bonds at all if
none have been made, so that the rules pane has a matrix to edit from the start.
"""
function buildrules(m::EditorModel)
    isempty(m.placements) && return nothing
    kept = [m.species[sp] for sp in shownspecies(m)]
    return BindingRules(effectivematrix(m), kept)
end

"""
    notify!(m, text; kind=:info)

Put a line in the status bar. `kind` is `:warning` for something the editor declined to do.
"""
function notify!(m::EditorModel, text::AbstractString; kind::Symbol=:info)
    m.message = text
    m.messagekind = kind
    return m
end

"""
    addspecies!(m)

Append another instance of the editor's species and make it the active one.
"""
function addspecies!(m::EditorModel)
    push!(m.species, copy(m.base))
    m.active_species = length(m.species)
    refresh!(m)
    return notify!(m, string("added species ", length(m.species)))
end

"""
    removespecies!(m)

Drop the last species, if nothing depends on it.

The last one only, so that the colors of the others do not shift under bonds already set. A
species carrying a bond or a placement is kept, with a note saying which, rather than silently
taking those with it.
"""
function removespecies!(m::EditorModel)
    n = length(m.species)
    n > 1 || return notify!(m, "one species is the minimum"; kind=:warning)
    n in usedspecies(m.placements) && return notify!(m, "species $n is placed"; kind=:warning)

    gc = sitecolors(m)
    cols = [gc[(n, k)] for k in 1:nsites(m.species[n])]
    mat = effectivematrix(m)
    any(mat[c, c2] for c in cols, c2 in axes(mat, 2)) && return notify!(m, "species $n has bonds"; kind=:warning)

    pop!(m.species)
    m.active_species = min(m.active_species, length(m.species))
    # The colors that species held no longer exist, so any override naming them goes too.
    filter!(kv -> kv.first[1] ∉ cols && kv.first[2] ∉ cols, m.overrides)
    refresh!(m)
    return notify!(m, string("removed species ", n))
end

"""
    togglebond!(m, pair)

Turn the bond between a pair of colors on or off, recording it as an override so that it holds
against whatever the geometry says.
"""
function togglebond!(m::EditorModel, (c1, c2)::Tuple{Int,Int})
    n = totalcolors(m)
    (1 <= c1 <= n && 1 <= c2 <= n) || return m
    mat = effectivematrix(m)
    on = !mat[c1, c2]
    m.overrides[(c1, c2)] = on
    m.overrides[(c2, c1)] = on
    refresh!(m)
    return notify!(m, string(colorlabel(c1), on ? " ─ " : " ╌ ", colorlabel(c2), on ? " bonded" : " cleared"))
end

"""
    anchorposition(m)

Return where the cursor currently sits in world coordinates, or `nothing` if there is no free
site. Read before a mutation so that [`refresh!`](@ref) can put the cursor back near it.
"""
function anchorposition(m::EditorModel)
    s = anchorsite(m)
    return s === nothing ? nothing : s.pose.x
end

"""
    anchorsite(m)

Return the absolute binding site the next attachment would seat against, or `nothing` if the
structure has no free site left.
"""
function anchorsite(m::EditorModel)
    isempty(m.free) && return nothing
    i, k = m.free[m.anchor]
    spidx, pose = m.placements[i]
    return pose * bindingsite(m.species[spidx], k)
end

matesite(m::EditorModel) = bindingsite(m.species[m.active_species], m.incoming)

# A 2D bond fixes the partner's orientation completely: `Roly.standard_twist` turns it by π
# regardless of `t`. Only a 3D bond leaves a twist to choose.
function ntwists(m::EditorModel)
    a = anchorsite(m)
    (a === nothing || dimension(m.base) == 2) && return 1
    return twistfreedom(a, matesite(m))
end

"""
    pendingcontacts(m, preview)

Return `(position, placement index, site, pending site)` for every site of the pending particle
that would touch one already placed.

Seating a particle against the chosen site can put it against several at once, which is how a
ring closes. Those extra bonds are what makes an arrangement reveal rules that were not aimed
at, so the drawing names them before the placement is committed rather than after.
"""
function pendingcontacts(m::EditorModel, preview)
    out = Tuple{typeof(preview.x),Int,Int,Int}[]
    spcs = m.species[m.active_species]
    pending = [preview * s for s in bindingsites(spcs)]
    sites = absolutesites(m.placements, m.species)
    for i in eachindex(m.placements), (k1, s1) in enumerate(sites[i]), (k2, s2) in enumerate(pending)
        istouching(s1, s2) && push!(out, (s1.pose.x, i, k1, k2))
    end
    return out
end

"""
    blockedby(m, pose)

Return the index of the first placement the pending particle would overlap, or `nothing` if it
would sit clear.

Read while drawing as well as on attaching, so that a placement the editor is going to refuse
looks refused beforehand rather than only saying so afterwards.
"""
function blockedby(m::EditorModel, pose)
    spcs = m.species[m.active_species]
    for (i, (j, other)) in enumerate(m.placements)
        overlap(spcs => pose, m.species[j] => other) && return i
    end
    return nothing
end

"""
    previewpose(m)

Return the pose the pending attachment would take, or `nothing` if there is nowhere to attach.
"""
function previewpose(m::EditorModel)
    a = anchorsite(m)
    a === nothing && return nothing
    return seatpose(a, matesite(m), m.twist)
end

function attach!(m::EditorModel)
    pose = previewpose(m)
    if pose === nothing
        return notify!(m, "no free site"; kind=:warning)
    end
    blocked = blockedby(m, pose)
    blocked === nothing || return notify!(m, "overlaps particle $blocked"; kind=:warning)
    near = anchorposition(m)
    push!(m.placements, (m.active_species, pose))
    m.message = ""
    refresh!(m; near)
    return m
end

function seed!(m::EditorModel)
    # A new component goes far enough from everything already placed that it cannot touch it,
    # so its species can be given rules of its own.
    offset = 0.0
    for (i, pose) in m.placements
        offset = max(offset, pose.x[1] + 3 * bounding_radius(m.species[i]))
    end
    origin = one(poseof(m))
    pose = typeof(origin)(origin.x + SVector(offset, zero(offset)), origin.psi)
    push!(m.placements, (m.active_species, pose))
    # The cursor follows onto the new component, which is where the next attachment belongs.
    refresh!(m; near=pose.x)
    return m
end

function undo!(m::EditorModel)
    near = anchorposition(m)
    length(m.placements) > 1 && pop!(m.placements)
    refresh!(m; near)
    return m
end

"""
    reset!(m)

Clear the construction back to a single particle, keeping the rules it produced.

The bonds the geometry currently implies are written into the overrides first, so that emptying
the drawing does not empty the table with it: an arrangement can be built to discover a rule and
then cleared away to build the next one.
"""
function reset!(m::EditorModel)
    mat = effectivematrix(m)
    for c1 in axes(mat, 1), c2 in axes(mat, 2)
        m.overrides[(c1, c2)] = mat[c1, c2]
    end
    P = poseof(m)
    m.placements = Placement{P}[(m.active_species, one(P))]
    m.refit = true
    refresh!(m)
    resetcamera!(m)
    return m
end

### Input

# Screen directions in world coordinates. The drawing puts world +y up, so these are what the
# arrow keys mean geometrically.
const ARROWS = (up=SVector(0.0, 1.0), down=SVector(0.0, -1.0), left=SVector(-1.0, 0.0), right=SVector(1.0, 0.0))

"""
    stepanchor!(m, dir)

Move the cursor one free site along the perimeter, in whichever of the two directions better
matches the screen direction `dir`.

The cursor only ever moves to a neighbor in the clockwise order, so it never steps over a site
on the way, and there is always a neighbor to move to, so no arrow key is ever dead. Which
neighbor it takes is decided by the direction alone, comparing the bearing of each against the
key rather than the distance, so a near neighbor is not passed over for a distant one.
"""
function stepanchor!(m::EditorModel, dir::Symbol)
    length(m.free) > 1 || return m
    d = getfield(ARROWS, dir)
    cam = camera(m)
    sites = absolutesites(m.placements, m.species)
    at(n) = plane(cam, sites[m.free[n][1]][m.free[n][2]].pose.x)
    hu, hv = at(m.anchor)
    function alignment(n)
        u, v = at(n)
        du, dv = u - hu, v - hv
        return (du * d[1] + dv * d[2]) / max(hypot(du, dv), eps())
    end
    nxt = mod1(m.anchor + 1, length(m.free))
    prv = mod1(m.anchor - 1, length(m.free))
    m.anchor = alignment(nxt) >= alignment(prv) ? nxt : prv
    return pivot!(m)
end

"""
    steppair!(m, dir)

Move the rules cursor, which points at one cell of the interaction matrix.

Up and down move the row, left and right the column, so the cursor moves the way the matrix is
laid out. The pair list beside it marks the same cell when that pair is bonded, so it follows
along without needing a cursor of its own.
"""
function steppair!(m::EditorModel, dir::Symbol)
    n = totalcolors(m)
    n > 0 || return m
    c1, c2 = m.pair
    dir === :up && (c1 = mod1(c1 - 1, n))
    dir === :down && (c1 = mod1(c1 + 1, n))
    dir === :left && (c2 = mod1(c2 - 1, n))
    dir === :right && (c2 = mod1(c2 + 1, n))
    m.pair = (c1, c2)
    return m
end

"""
    cyclefocus!(m, step)

Move the focus between the panes that take input. The construction pane drops out of the cycle
while it is hidden.
"""
function cyclefocus!(m::EditorModel, step::Int)
    panes = m.showconstruction ? (:rules, :construction, :enumeration) : (:rules, :enumeration)
    i = something(findfirst(==(m.focus), panes), 1)
    m.focus = panes[mod1(i + step, length(panes))]
    return m
end

"""
    setmaxsize!(m, step)
    setmaxstrs!(m, step)

Move the enumeration's bounds one notch. Changing either marks the last run as out of date
rather than re-running, since the point of bounding a run is to decide before paying for it.
"""
function setmaxsize!(m::EditorModel, step::Int)
    new = clamp(m.maxsize + step, first(MAXSIZE_RANGE), last(MAXSIZE_RANGE))
    new == m.maxsize && return m
    m.maxsize = new
    m.stale = true
    return m
end

function setmaxstrs!(m::EditorModel, step::Int)
    i = something(findfirst(>=(m.maxstrs), MAXSTRS_RANGE), 1)
    new = MAXSTRS_RANGE[clamp(i + step, 1, length(MAXSTRS_RANGE))]
    new == m.maxstrs && return m
    m.maxstrs = new
    m.stale = true
    return m
end

"""
    stepselection!(m, dir)

Move the selection through the enumerated structures, along the grid the pane lays them out in.
"""
function stepselection!(m::EditorModel, dir::Symbol, ncols::Int)
    isempty(m.polyforms) && return m
    step = if dir === :left
        -1
    elseif dir === :right
        1
    elseif dir === :up
        -ncols
    else
        ncols
    end
    m.selected = clamp(m.selected + step, 1, length(m.polyforms))
    return m
end

function update!(m::EditorModel, evt::KeyEvent)
    m.message = ""
    k, c = evt.key, evt.char

    # Keys that mean the same thing wherever the focus is.
    if k === :escape || (k === :char && c == 'q')
        m.quit = true
        return nothing
    elseif k === :tab
        cyclefocus!(m, 1)
        return nothing
    elseif k === :backtab
        cyclefocus!(m, -1)
        return nothing
    elseif k === :char && c == 'b'
        m.showconstruction = !m.showconstruction
        if m.showconstruction
            m.focus = :construction
        elseif m.focus === :construction
            m.focus = :rules
        end
        return nothing
    elseif k === :char && c == 'e'
        enumerate!(m)
        return nothing
    end

    if m.focus === :rules
        if k === :up || k === :down || k === :left || k === :right
            steppair!(m, k)
        elseif k === :enter
            togglebond!(m, m.pair)
        elseif k === :char && c == 'a'
            addspecies!(m)
        elseif k === :char && c == 'd'
            removespecies!(m)
        end
        return nothing
    end

    if m.focus === :enumeration
        if k === :up || k === :down || k === :left || k === :right
            stepselection!(m, k, m.gridcols)
        elseif k === :char && c == 's'
            setmaxsize!(m, -1)
        elseif k === :char && c == 'S'
            setmaxsize!(m, 1)
        elseif k === :char && c == 'x'
            setmaxstrs!(m, -1)
        elseif k === :char && c == 'X'
            setmaxstrs!(m, 1)
        end
        return nothing
    end

    if k === :up || k === :down || k === :left || k === :right
        stepanchor!(m, k)
    elseif k === :char && (c == ',' || c == '.')
        if !isempty(m.free)
            m.anchor = mod1(m.anchor + (c == '.' ? 1 : -1), length(m.free))
            pivot!(m)
        end
    elseif k === :enter
        attach!(m)
    elseif k === :backspace
        undo!(m)
    elseif k === :char && c == 'r'
        m.incoming = mod1(m.incoming + 1, nsites(m.species[m.active_species]))
        m.twist = mod(m.twist, ntwists(m))
    elseif k === :char && c == 'R'
        m.incoming = mod1(m.incoming - 1, nsites(m.species[m.active_species]))
        m.twist = mod(m.twist, ntwists(m))
    elseif k === :char && (c == 't' || c == 'T')
        m.twist = mod(m.twist + (c == 't' ? 1 : -1), ntwists(m))
    elseif k === :char && c == 'n'
        seed!(m)
    elseif k === :char && c == 'c'
        reset!(m)
    elseif k === :char && c == '0'
        m.refit = true
    elseif k === :char && c == '-'
        zoom!(m, 1 / ZOOM_STEP)
    elseif k === :char && (c == '=' || c == '+')
        zoom!(m, ZOOM_STEP)
    elseif k === :char && (c == '[' || c == ']')
        turn!(m, c == ']' ? TURN_STEP : -TURN_STEP, 0.0)
    elseif k === :char && (c == '{' || c == '}')
        turn!(m, 0.0, c == '}' ? TURN_STEP : -TURN_STEP)
    elseif k === :char && isdigit(c)
        d = parse(Int, c)
        ensurespecies!(m, d)
        m.active_species = d
        m.incoming = mod1(m.incoming, nsites(m.species[d]))
        m.twist = mod(m.twist, ntwists(m))
        refresh!(m)
    end
    return nothing
end

### Rendering

"""
    view(m, f)

Draw the editor: keys on the left, then the rules, then the construction, then the enumeration.

The rules are always on screen, next to the keys that edit them. The construction pane is
optional, so that a design can be worked on as a bond table alone, and whatever is left over
goes to the enumeration, which wants the most room. Rules and construction are given the same
width so that a species drawn in one is the size it is in the other.
"""
function view(m::EditorModel, f::Frame)
    buf = f.buffer
    # One row at the foot for the status line, which is where the editor says what it just did
    # or why it declined to.
    body, status = split_layout(Layout(Vertical, Constraint[Fill(), Fixed(1)]), f.area)
    draw_status(m, status, buf)
    f = Frame(buf, body, f.gfx_regions, f.pixel_snapshots)

    room = f.area.width - SIDEBAR_W - PANE_W
    # A matrix too big for the rules pane gets a pane of its own rather than being dropped for
    # the pair list: the whole window height holds a far larger one than a third of a pane does.
    inline = matrixinline(m, PANE_W - 2, f.area.height)
    wanted = m.rules === nothing ? 0 : matrixwidth(m.rules) + 2
    withmatrix = !inline && wanted <= room
    room -= withmatrix ? wanted : 0
    withconstruction = m.showconstruction && room >= PANE_W
    room -= withconstruction ? PANE_W : 0
    withpreview = room >= PREVIEW_MIN_W

    cs = Constraint[Fixed(SIDEBAR_W), withpreview ? Fixed(PANE_W) : Fill()]
    withmatrix && push!(cs, Fixed(wanted))
    withconstruction && push!(cs, withpreview ? Fixed(PANE_W) : Fill())
    withpreview && push!(cs, Fill())

    rects = split_layout(Layout(Horizontal, cs), f.area)
    next = 2
    draw_rules(m, rects[next], buf; withmatrix=inline, elsewhere=withmatrix)
    next += 1
    if withmatrix
        draw_matrixpane(m, rects[next], buf)
        next += 1
    end
    withconstruction && draw_construction(m, rects[next], buf)
    withpreview && draw_preview(m, rects[end], buf)
    draw_sidebar(m, rects[1], buf)
    return nothing
end

"""
    draw_matrixpane(m, rect, buf)

Draw the interaction matrix in a pane of its own, which is where it goes once it outgrows the
rules pane.
"""
function draw_matrixpane(m::EditorModel, rect::Rect, buf::Buffer)
    inner = render(focusblock(m, :rules; title="Matrix"), rect, buf)
    (inner.width < 3 || inner.height < 2) && return buf
    draw_matrix(m, inner, buf)
    return buf
end

"""
    draw_status(m, rect, buf)

Draw the status line: what the editor last did on the left, what it is holding on the right.
"""
function draw_status(m::EditorModel, rect::Rect, buf::Buffer)
    style = m.messagekind === :warning ? tstyle(:warning) : tstyle(:text)
    right = string(
        length(m.species),
        length(m.species) == 1 ? " species  " : " species  ",
        m.rules === nothing ? 0 : nbonds(m.rules),
        " bonds  ",
        length(m.placements),
        " placed",
    )
    bar = StatusBar(; left=[Span(m.message, style)], right=[Span(right, tstyle(:text_dim))])
    render(bar, rect, buf)
    return buf
end

"""
    focusblock(m, pane; kwargs...)

Return a `Block` whose border shows whether `pane` currently has the keys.
"""
function focusblock(m::EditorModel, pane::Symbol; kwargs...)
    focused = m.focus === pane
    return Block(;
        border_style=focused ? tstyle(:border_focus) : tstyle(:border),
        title_style=tstyle(focused ? :title : :text_dim; bold=focused),
        kwargs...,
    )
end

function drawlines!(buf::Buffer, r::Rect, lines)
    for (i, (text, style)) in enumerate(lines)
        i > r.height && break
        set_string!(buf, r.x, r.y + i - 1, text, style; max_x=right(r))
    end
    return buf
end

"""
    draw_speciesstrip(m, rect, buf)

Draw the species as a row of swatches, marking the one the next placement would use.
"""
function draw_speciesstrip(m::EditorModel, rect::Rect, buf::Buffer)
    x = rect.x
    for i in eachindex(m.species)
        x + 2 > right(rect) && break
        active = i == m.active_species
        active && set_char!(buf, x, rect.y, '▶', tstyle(:text))
        set_char!(buf, x + 1, rect.y, '■', speciesstyle(i; bold=active))
        set_string!(buf, x + 2, rect.y, string(i), speciesstyle(i; bold=active); max_x=right(rect))
        x += 4
    end
    return buf
end

"""
    draw_sidebar(m, rect, buf)

Draw the species list, then the keys, in a section per pane.

Each pane keeps its own section, in a fixed order, whether or not it has the focus: only the
header lights up. Moving the sections about as the focus changes would mean hunting for a key
that has not gone anywhere.
"""
function draw_sidebar(m::EditorModel, rect::Rect, buf::Buffer)
    inner = render(Block(; title="Editor"), rect, buf)
    header(text, pane) =
        (rpad(string("── ", text, " "), 20, "─"), m.focus === pane ? tstyle(:title; bold=true) : tstyle(:text_dim))
    key(k, what) = (rpad(k, 6) * what, tstyle(:text_dim))

    lines = Tuple{String,Style}[]
    push!(lines, header("Any pane", :none))
    push!(lines, key("tab", "next pane"))
    push!(lines, key("b", m.showconstruction ? "hide build" : "build"))
    push!(lines, key("q", "accept"))

    c1, c2 = m.pair
    on = m.rules !== nothing && c1 <= ncolors(m.rules) && c2 <= ncolors(m.rules) && interactionmatrix(m.rules)[c1, c2]
    rules = [
        header("Rules", :rules),
        (
            string("  ", colorlabel(c1), on ? " ─ " : " ╌ ", colorlabel(c2), on ? "  bonded" : ""),
            on ? tstyle(:success) : tstyle(:text),
        ),
        key("↑↓←→", "cell"),
        key("enter", "bond"),
        key("a / d", "add/drop"),
    ]

    build = Tuple{String,Style}[]
    if m.showconstruction
        anchor = isempty(m.free) ? "none" : string(m.free[m.anchor][1], ".", m.free[m.anchor][2])
        pose = previewpose(m)
        blocked = pose === nothing ? nothing : blockedby(m, pose)
        closing = pose === nothing || blocked !== nothing ? 0 : length(pendingcontacts(m, pose)) - 1
        note = if blocked !== nothing
            (string("  blocked by ", blocked), tstyle(:warning))
        elseif closing > 0
            (string("  closes ", closing, closing == 1 ? " bond" : " bonds"), tstyle(:success))
        else
            ("", RESET)
        end
        append!(
            build,
            [
                header("Build", :construction),
                (string("  at ", anchor, " site ", m.incoming), tstyle(:text)),
                note,
                key("↑↓←→", "site"),
                key(", .", "step"),
                key("r / R", "turn"),
            ],
        )
        ntwists(m) > 1 && push!(build, key("t / T", "twist"))
        if dimension(m.base) == 3
            append!(build, [key("[ ]", "turn"), key("{ }", "tilt")])
        end
        append!(build, [key("enter", "attach"), key("bksp", "undo"), key("1-9", "species"), key("n / c", "part/clear")])
    end

    enum = [
        header("Enumeration", :enumeration),
        (string("  size ", m.maxsize, "  max ", m.maxstrs), tstyle(m.stale ? :warning : :text)),
        key("e", m.stale ? "run (stale)" : "run"),
        key("s / S", "size"),
        key("x / X", "max"),
        key("↑↓←→", "pick"),
    ]

    for section in (rules, build, enum)
        isempty(section) && continue
        push!(lines, ("", RESET))
        append!(lines, section)
    end

    drawlines!(buf, inner, lines)
    return buf
end

function draw_construction(m::EditorModel, rect::Rect, buf::Buffer)
    title = string(length(m.placements), " particles")
    box = render(focusblock(m, :construction; title="Construction", title_right=title), rect, buf)
    (box.width < 2 || box.height < 3) && return buf

    # The species live here rather than in the sidebar: which one is active only decides what
    # the next placement will be, so it belongs to the pane that does the placing.
    draw_speciesstrip(m, Rect(box.x, box.y, box.width, 1), buf)
    inner = Rect(box.x, box.y + 1, box.width, box.height - 1)

    preview = previewpose(m)
    bounds = preview === nothing ? m.placements : vcat(m.placements, [(m.active_species, preview)])
    v = worldfor!(m, bounds, inner.width, inner.height)

    blocked = preview === nothing ? nothing : blockedby(m, preview)
    ghost = preview === nothing ? nothing : (m.active_species, preview)
    drawparticles!(buf, inner, v, m.placements, m.species; ghost, blocked=blocked !== nothing, wire=true)

    # Labels sit at each site's own position, so the colors read straight off the drawing. A
    # bond gets both of its colors, written side by side, since which pair meets is the thing
    # being decided.
    cell(x) = (inner.x + x[1] ÷ 2, inner.y + x[2] ÷ 4)
    label!(x, style, ch) = set_char!(buf, cell(x)..., ch, style)
    function pair!(x, c1, style1, c2, style2)
        col, row = cell(x)
        set_char!(buf, col, row, colorlabel(c1), style1)
        set_char!(buf, col + 1, row, colorlabel(c2), style2)
        return buf
    end

    # Colors are the ones the returned rules will use, not each species' own, so that a label
    # here names the same row the interaction matrix does.
    gc = sitecolors(m)

    if preview !== nothing
        # Naming the ghost's sites is what shows its orientation: `r` turns it, and the labels
        # move round with it. Sites that would land on a bond are left to the pending pairs
        # below, which write them next to the site they would meet.
        spcs = m.species[m.active_species]
        meeting = Set(k for (_, _, _, k) in pendingcontacts(m, preview))
        for (k, s) in enumerate(bindingsites(spcs))
            k in meeting && continue
            site = preview * s
            facing(v.cam, site) || continue
            c = gc[(m.active_species, k)]
            sitestyle = blocked === nothing ? Style(; fg=sitecolor(m, c), dim=true) : tstyle(:text_dim)
            label!(todots(v, site.pose.x), sitestyle, colorlabel(c))
        end
    end

    if interiorbonds(m.base)
        for (x, sp1, k1, sp2, k2) in contacts(m.placements, m.species)
            c1, c2 = gc[(sp1, k1)], gc[(sp2, k2)]
            pair!(todots(v, x), c1, colorstyle(m, c1), c2, colorstyle(m, c2))
        end
    end

    # Free sites carry their color, so the construction and the rules pane label the same thing
    # the same way. The anchor is the one the next attachment seats against.
    sites = absolutesites(m.placements, m.species)
    for (n, (i, k)) in enumerate(m.free)
        n == m.anchor && continue
        facing(v.cam, sites[i][k]) || continue
        c = gc[(m.placements[i][1], k)]
        label!(todots(v, sites[i][k].pose.x), colorstyle(m, c), colorlabel(c))
    end

    # The pending bonds last, so they win their cells. Seating against the chosen site can put
    # the particle against several at once, which is how a ring closes, so every contact the
    # placement would make is labelled and not just the one that was aimed at.
    a = anchorsite(m)
    if a !== nothing && preview === nothing
        i, k = m.free[m.anchor]
        ca = gc[(m.placements[i][1], k)]
        label!(todots(v, a.pose.x), Style(; fg=sitecolor(m, ca), bold=true, underline=true), colorlabel(ca))
    elseif preview !== nothing && blocked !== nothing
        i, k = m.free[m.anchor]
        ca = gc[(m.placements[i][1], k)]
        label!(todots(v, a.pose.x), tstyle(:text_dim; bold=true), colorlabel(ca))
    elseif preview !== nothing
        aimed = a === nothing ? (0, 0) : m.free[m.anchor]
        for (x, i, k1, k2) in pendingcontacts(m, preview)
            c1 = gc[(m.placements[i][1], k1)]
            c2 = gc[(m.active_species, k2)]
            pair!(
                todots(v, x),
                c1,
                Style(; fg=sitecolor(m, c1), bold=(i, k1) == aimed, underline=true),
                c2,
                Style(; fg=sitecolor(m, c2), dim=true, underline=true),
            )
        end
    end
    return buf
end

"""
    draw_preview(m, rect, buf)

Draw the polyforms the current rules enumerate, smallest first, as many as the pane holds.

Each is drawn with the same outlines as the construction pane, so a structure here is recognisably
made of the blocks over there. The header reports how many the run found and whether it was cut
short, which is how a rules set that explodes announces itself.
"""
function draw_preview(m::EditorModel, rect::Rect, buf::Buffer)
    label = m.stale ? string(m.enuminfo, isempty(m.enuminfo) ? "" : "  ", "e to run") : m.enuminfo
    inner = render(focusblock(m, :enumeration; title="Enumeration", title_right=label), rect, buf)
    (inner.width < 4 || inner.height < 2) && return buf

    # A structure drawn at thumbnail size loses edges shorter than a braille dot, so the selected
    # one is redrawn as large as the pane allows, beside the grid.
    detail = !isempty(m.polyforms) && inner.width >= THUMB_W + DETAIL_W
    detailw = detail ? clamp(inner.width ÷ 3, DETAIL_W, 44) : 0
    grid = detail ? Rect(inner.x, inner.y, inner.width - detailw, inner.height) : inner

    if isempty(m.polyforms)
        m.gridcols = 1
        drawlines!(buf, inner, [(m.stale ? "press e to enumerate" : m.enuminfo, tstyle(:text_dim))])
        return buf
    end

    ncols = max(1, grid.width ÷ THUMB_W)
    nrows = max(1, grid.height ÷ THUMB_H)
    m.gridcols = ncols
    # Scroll by whole rows, so the selection is always on screen.
    firstrow = max(0, (m.selected - 1) ÷ ncols - nrows + 1)
    for n in (firstrow * ncols + 1):length(m.polyforms)
        slot = n - firstrow * ncols
        slot > ncols * nrows && break
        col = (slot - 1) % ncols
        row = (slot - 1) ÷ ncols
        r = Rect(grid.x + col * THUMB_W, grid.y + row * THUMB_H, THUMB_W - 1, THUMB_H - 1)
        draw_thumbnail(m, m.polyforms[n], r, buf; wire=true)
        # Each one carries its position in the list, so a structure can be named when talking
        # about it and the detail box can say which it is showing.
        selected = n == m.selected
        set_string!(buf, r.x, r.y, string(n), selected ? tstyle(:text; bold=true) : tstyle(:text_dim); max_x=right(r))
        selected && set_string!(
            buf,
            r.x,
            bottom(r),
            "─"^r.width,
            m.focus === :enumeration ? tstyle(:border_focus) : tstyle(:border);
            max_x=right(r),
        )
    end

    detail && draw_detail(m, Rect(right(inner) - detailw + 2, inner.y, detailw - 2, inner.height), buf)
    return buf
end

"""
    draw_detail(m, rect, buf)

Draw the selected structure on its own, at whatever scale the box allows, which is where a
thumbnail too small to show every edge can be read properly.
"""
function draw_detail(m::EditorModel, rect::Rect, buf::Buffer)
    i = clamp(m.selected, 1, length(m.polyforms))
    poly = m.polyforms[i]
    n = nparticles(poly)
    title = string("#", i)
    inner = render(Block(; title, title_right=string(n, n == 1 ? " particle" : " particles")), rect, buf)
    (inner.width < 3 || inner.height < 2) && return buf
    if inner.height > 4
        draw_thumbnail(m, poly, Rect(inner.x, inner.y, inner.width, inner.height - 2), buf; labels=true)
        draw_composition(m, poly, Rect(inner.x, bottom(inner) - 1, inner.width, 2), buf)
    else
        draw_thumbnail(m, poly, inner, buf; labels=true)
    end
    return buf
end

"""
    draw_composition(m, poly, rect, buf)

Write the structure's composition vector: how many particles of each species it has, then how
many bonds of each bonding color pair, in the order `Roly.composition` returns them.
"""
function draw_composition(m::EditorModel, poly, rect::Rect, buf::Buffer)
    # A polyform's composition is indexed by the rules it was enumerated under, which are not
    # necessarily the ones on screen: the preview keeps its last run while the rules move on.
    rules = bindingrules(poly)
    comp = composition(poly)
    ns = nspecies(rules)
    pairs = bonded_colors(rules)

    set_string!(buf, rect.x, rect.y, "spc", tstyle(:text_dim); max_x=right(rect))
    x = rect.x + 4
    for i in 1:ns
        x + 1 > right(rect) && break
        set_string!(buf, x, rect.y, lpad(comp[i], 2), speciesstyle(i; bold=comp[i] > 0))
        x += 3
    end

    set_string!(buf, rect.x, rect.y + 1, "bnd", tstyle(:text_dim); max_x=right(rect))
    x = rect.x + 4
    for (j, (c1, c2)) in enumerate(pairs)
        x + 1 > right(rect) && break
        count = comp[ns + j]
        set_string!(buf, x, rect.y + 1, lpad(count, 2), bondstyle(m, c1, c2; bold=count > 0))
        x += 3
    end
    return buf
end

function draw_thumbnail(m::EditorModel, poly, rect::Rect, buf::Buffer; labels::Bool=false, wire::Bool=false)
    (rect.width < 3 || rect.height < 2) && return buf
    # Drawn from the polyform's own rules, since its particles index that species list rather
    # than whatever the editor has moved on to.
    spcs = species(bindingrules(poly))
    placements = [(p.speciesindex, p.pose) for p in poly.particles]
    isempty(placements) && return buf
    v = fitworld(placements, spcs, rect.width, rect.height)
    drawparticles!(buf, rect, v, placements, spcs; wire)
    labels && drawsitelabels!(m, poly, spcs, placements, v, rect, buf)
    return buf
end

"""
    drawsitelabels!(m, poly, spcs, placements, v, rect, buf)

Name every binding site of a drawn structure with its color, in that color, and write both
colors side by side wherever two sites meet.

The same labelling the construction pane uses, so a bond in an enumerated structure can be read
off against the matrix without counting sites round a particle.
"""
function drawsitelabels!(m::EditorModel, poly, spcs, placements, v::World, rect::Rect, buf::Buffer)
    cell(x) = (rect.x + x[1] ÷ 2, rect.y + x[2] ÷ 4)
    sites = [[pose * s for s in bindingsites(spcs[i])] for (i, pose) in placements]

    touching = Set{Tuple{Int,Int}}()
    for i in eachindex(placements), j in eachindex(placements)
        j > i || continue
        for (k1, s1) in enumerate(sites[i]), (k2, s2) in enumerate(sites[j])
            istouching(s1, s2) || continue
            push!(touching, (i, k1))
            push!(touching, (j, k2))
            interiorbonds(spcs[placements[i][1]]) || continue
            c1, c2 = color(s1), color(s2)
            col, row = cell(todots(v, s1.pose.x))
            set_char!(buf, col, row, colorlabel(c1), colorstyle(m, c1))
            set_char!(buf, col + 1, row, colorlabel(c2), colorstyle(m, c2))
        end
    end
    for i in eachindex(placements), (k, s) in enumerate(sites[i])
        (i, k) in touching && continue
        facing(v.cam, s) || continue
        c = color(s)
        col, row = cell(todots(v, s.pose.x))
        set_char!(buf, col, row, colorlabel(c), colorstyle(m, c))
    end
    return buf
end

"""
    matrixfits(rules, w)

Whether the interaction matrix is narrow enough to draw in a `w` column pane. Each color costs
two columns, plus two for the row labels.
"""
matrixfits(rules::BindingRules, w::Int) = 2 + 2 * ncolors(rules) <= w

"""
    matrixwidth(rules)
    matrixinline(m, w, h)

How wide the matrix needs to be, and whether the rules pane can hold it beside the pair list in
`w` by `h` cells.
"""
matrixwidth(rules::BindingRules) = 2 + 2 * ncolors(rules)

function matrixinline(m::EditorModel, w::Int, h::Int)
    rules = m.rules
    rules === nothing && return true
    pairsw = min(PAIR_W + 1, w ÷ 2)
    reserved = min(length(shownspecies(m)) + 1, GALLERY_MIN_ROW)
    return matrixfits(rules, w - pairsw) && ncolors(rules) + 1 <= max(h - 1 - reserved, 1)
end

"""
    gallerylayout(n, w, h)

Return `(ncols, nrows, rowheight)` for drawing `n` species in a `w` by `h` region, or `nothing`
if they cannot be drawn large enough to read.

The width is filled first and the height spent afterwards: a species drawn in a tall thin column
of its own wastes most of the pane, so as many fit across as the boxes have room for and the rows
take what is left. Once the rows would be too short to read the gallery gives up and the caller
falls back to a compact list.
"""
function gallerylayout(n::Int, w::Int, h::Int)
    (n < 1 || w < GALLERY_MIN_W || h < GALLERY_MIN_ROW) && return nothing
    ncols = min(n, max(1, w ÷ GALLERY_MIN_W))
    nrows = cld(n, ncols)
    rowheight = h ÷ nrows
    rowheight < GALLERY_MIN_ROW && return nothing
    return (ncols, nrows, min(rowheight, GALLERY_ROW))
end

"""
    pairrows(rules, w)

Return how many rows the pair list needs in a `w` column pane, laying pairs out in columns.
"""
pairrows(rules::BindingRules, w::Int) = cld(max(length(bonded_colors(rules)), 1), max(1, w ÷ PAIR_W))

function draw_rules(m::EditorModel, rect::Rect, buf::Buffer; withmatrix::Bool=true, elsewhere::Bool=false)
    rules = m.rules
    if rules === nothing
        inner = render(focusblock(m, :rules; title="Rules"), rect, buf)
        (inner.width < 2 || inner.height < 2) && return buf
        drawlines!(buf, inner, [("no bonds yet", tstyle(:text_dim))])
        return buf
    end

    # The matrix is either here, in a pane of its own, or nowhere: say which, since an absent
    # matrix with no explanation looks like a fault.
    label = withmatrix || elsewhere ? "" : "matrix needs a wider window"
    inner = render(focusblock(m, :rules; title="Rules", title_right=label), rect, buf)
    (inner.width < 2 || inner.height < 2) && return buf
    used = shownspecies(m)

    # The species summary always keeps a few rows, so the tables can never squeeze it out.
    reserved = min(length(used) + 1, GALLERY_MIN_ROW)
    budget = max(inner.height - 1 - reserved, 1)
    pairsw = min(PAIR_W + 1, inner.width ÷ 2)
    tableheight = if withmatrix
        max(ncolors(rules) + 1, min(pairrows(rules, pairsw), budget))
    else
        min(pairrows(rules, inner.width), budget)
    end
    room = max(inner.height - tableheight - 1, 0)
    layout = gallerylayout(length(used), inner.width, room)
    galleryheight = layout === nothing ? room : layout[2] * layout[3]

    parts = split_layout(Layout(Vertical, Constraint[Fixed(galleryheight), Fixed(1), Fill()]), inner)
    if layout === nothing
        draw_specieslist(m, rules, used, parts[1], buf)
    else
        draw_gallery(m, rules, used, layout, parts[1], buf)
    end
    set_string!(buf, parts[2].x, parts[2].y, "─"^parts[2].width, tstyle(:border); max_x=right(parts[2]))

    # The matrix says what is possible and the list says what is set. When the matrix has been
    # moved to a pane of its own the list takes the whole width here.
    if withmatrix
        tables = split_layout(Layout(Horizontal, Constraint[Fill(), Fixed(pairsw)]), parts[3])
        draw_matrix(m, tables[1], buf)
        draw_pairs(m, tables[2], buf)
    else
        draw_pairs(m, parts[3], buf)
    end
    return buf
end

"""
    draw_specieslist(m, rules, used, rect, buf)

List the species and their site colors on one line each, the fallback for when there are too
many species to draw them. Carries the same information as the gallery apart from the shapes.
"""
function draw_specieslist(m::EditorModel, rules::BindingRules, used, rect::Rect, buf::Buffer)
    (isempty(used) || rect.height < 1) && return buf
    gc = sitecolors(m)
    entryw = 5 + 2 * maximum(nsites(m.species[i]) for i in used)
    ncols = max(1, rect.width ÷ entryw)
    capacity = ncols * rect.height
    shown = length(used) > capacity ? capacity - 1 : length(used)
    for n in 1:shown
        spidx = used[n]
        x = rect.x + entryw * ((n - 1) % ncols)
        y = rect.y + (n - 1) ÷ ncols
        set_string!(buf, x, y, string(rpad(n, 2), "■"), speciesstyle(spidx); max_x=right(rect))
        for k in 1:nsites(m.species[spidx])
            cx = x + 3 + 2 * k
            cx > right(rect) && break
            c = gc[(spidx, k)]
            inert = c <= ncolors(rules) ? isinert(rules, c) : true
            set_char!(buf, cx, y, colorlabel(c), inert ? tstyle(:text_dim) : colorstyle(m, c; bold=true))
        end
    end
    if shown < length(used)
        x = rect.x + entryw * (shown % ncols)
        y = rect.y + shown ÷ ncols
        set_string!(buf, x, y, string("+", length(used) - shown), tstyle(:text_dim); max_x=right(rect))
    end
    return buf
end

"""
    draw_gallery(m, rules, used, rect, buf)

Draw each placed species once, with every site marked by its color index. Sites whose color
bonds nothing render dim, matching how the Makie extension greys inert sites. Species are
labelled by their position in the returned rules, but keep the color the construction pane
draws them in.

`layout` is the `(ncols, nrows, rowheight)` that [`gallerylayout`](@ref) chose for the region, so
that a few species are drawn large and many are drawn small.
"""
function draw_gallery(m::EditorModel, rules::BindingRules, used, layout, rect::Rect, buf::Buffer)
    (isempty(used) || rect.height < 2) && return buf
    gc = sitecolors(m)
    ncols, nrows, rowheight = layout
    # Columns are measured out rather than split with `Fill`, which shares space unevenly when
    # several of equal weight compete: five over 64 cells come out 12, 10, 8, 6, 28.
    colwidth = rect.width ÷ ncols
    for (n, spidx) in enumerate(used)
        rowidx, colidx = fldmod1(n, ncols)
        rowidx > nrows && break
        r = Rect(rect.x + (colidx - 1) * colwidth, rect.y + (rowidx - 1) * rowheight, colwidth, rowheight)
        r.height < 2 && continue
        set_string!(buf, r.x, r.y, string("species ", n), speciesstyle(spidx); max_x=right(r))
        cr = Rect(r.x, r.y + 1, r.width, r.height - 1)
        spcs = m.species[spidx]
        pose = one(posetype(spcs))
        v = fitworld([(spidx, pose)], m.species, cr.width, cr.height)
        drawparticles!(buf, cr, v, [(spidx, pose)], m.species; wire=true)
        # The two colors the rules cursor names are marked here as well, so that a cell of the
        # matrix can be read as the sites it stands for. Every site is named, those on the far
        # side of a 3D particle included: this is the reference drawing of the species, so a
        # color that appears in the matrix has to appear here too, only dimmed to say that the
        # face carrying it is turned away.
        for (k, s) in enumerate(bindingsites(spcs))
            c = gc[(spidx, k)]
            dx, dy = todots(v, s.pose.x)
            inert = c <= ncolors(rules) ? isinert(rules, c) : true
            style = if c == m.pair[1] || c == m.pair[2]
                Style(; fg=sitecolor(m, c), bold=true, underline=true)
            elseif inert || !facing(v.cam, s)
                tstyle(:text_dim)
            else
                colorstyle(m, c; bold=true)
            end
            set_char!(buf, cr.x + dx ÷ 2, cr.y + dy ÷ 4, colorlabel(c), style)
        end
    end
    return buf
end

"""
    draw_matrix(m, rect, buf)

Draw the color interaction matrix as a grid, `■` where two colors bond and `·` where they do
not. This is the `output=:matrix` return value, shown as it is being built.

The cell the rules cursor points at is picked out, along with its row and column labels, since
which pair a cell stands for is otherwise hard to trace across a wide matrix. A bond is drawn in
its two site colors mixed, so it reads as belonging to both.
"""
function draw_matrix(m::EditorModel, rect::Rect, buf::Buffer)
    rules = m.rules
    (rules === nothing || rect.width < 3 || rect.height < 2) && return buf
    n = ncolors(rules)
    imat = interactionmatrix(rules)
    active = m.focus === :rules
    r1, c1 = m.pair
    for c in 1:n
        x = rect.x + 2 + 2 * (c - 1)
        x > right(rect) && break
        set_char!(buf, x, rect.y, colorlabel(c), c == c1 ? colorstyle(m, c; bold=true) : tstyle(:text_dim))
    end
    for c in 1:n
        y = rect.y + c
        y > bottom(rect) && break
        set_char!(buf, rect.x, y, colorlabel(c), colorstyle(m, c; bold=c == r1))
        for c2 in 1:n
            x = rect.x + 2 + 2 * (c2 - 1)
            x > right(rect) && break
            on = imat[c, c2]
            if !on
                set_char!(buf, x, y, '·', (c, c2) == m.pair ? tstyle(:text; underline=active) : tstyle(:text_dim))
            else
                # A bond cell is split between its two site colors rather than painted in one,
                # the row's as foreground and the column's as background, so a cell names both
                # ends of the bond it stands for. The matrix being symmetric, the mirrored cell
                # shows the same pair the other way round.
                set_char!(
                    buf,
                    x,
                    y,
                    bondglyph(),
                    Style(;
                        fg=sitecolor(m, c),
                        bg=sitecolor(m, c2),
                        bold=(c, c2) == m.pair,
                        underline=(c, c2) == m.pair && active,
                    ),
                )
            end
        end
    end
    return buf
end

"""
    draw_pairs(m, rect, buf)

List the bonding color pairs beside the matrix.

The matrix says what is possible and the list says what is set, which is what you want when the
relation is sparse and the grid is mostly dots. Pairs fill as many columns as the space allows,
and any that still do not fit are counted.
"""
function draw_pairs(m::EditorModel, rect::Rect, buf::Buffer)
    rules = m.rules
    (rules === nothing || rect.height < 1 || rect.width < PAIR_W) && return buf
    pairs = bonded_colors(rules)
    if isempty(pairs)
        set_string!(buf, rect.x, rect.y, "no bonds", tstyle(:text_dim); max_x=right(rect))
        return buf
    end
    ncols = max(1, rect.width ÷ PAIR_W)
    capacity = ncols * rect.height
    shown = length(pairs) > capacity ? capacity - 1 : length(pairs)
    for i in 1:shown
        c1, c2 = pairs[i]
        selected = (c1, c2) == minmax(m.pair...)
        x = rect.x + PAIR_W * ((i - 1) % ncols)
        y = rect.y + (i - 1) ÷ ncols
        under = selected && m.focus === :rules
        set_char!(buf, x, y, colorlabel(c1), colorstyle(m, c1; bold=selected, underline=under))
        set_string!(buf, x + 1, y, " ─ ", bondstyle(m, c1, c2; bold=selected, underline=under); max_x=right(rect))
        set_char!(buf, x + 4, y, colorlabel(c2), colorstyle(m, c2; bold=selected, underline=under))
    end
    if shown < length(pairs)
        x = rect.x + PAIR_W * (shown % ncols)
        y = rect.y + shown ÷ ncols
        set_string!(buf, x, y, string("+", length(pairs) - shown), tstyle(:text_dim); max_x=right(rect))
    end
    return buf
end

### Entry point

# Documented on the stub in `Roly`, so that Documenter finds it without loading this extension.
function ruleeditor(spcs::ParticleSpecies; output::Symbol=:rules)
    output in (:rules, :bonds, :matrix) || throw(ArgumentError("output must be :rules, :bonds, or :matrix"))

    m = EditorModel(spcs)
    app(m)

    rules = m.rules
    (rules === nothing || nbonds(rules) == 0) && return nothing
    output == :bonds && return bondsfromrules(rules)
    output == :matrix && return interactionmatrix(rules)
    return rules
end

end # module
