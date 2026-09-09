module Roly

using LinearAlgebra, StaticArrays, SparseArrays, Rotations, Statistics, Random
using Base.Iterators, DataStructures
using Graphs, NautyGraphs
using ReverseSearch

# Geometry primitives
export Pose, dimension, numtype, posetype
export Rotation, Angle2d, RotXYZ, RotMatrix3, rotation_angle, rotation_axis, SVector

# Binding sites
export BindingSite, color

# Bodies and their symmetry groups
export Polyhedron
export Tetrahedron, Cube, Octahedron, Dodecahedron, Icosahedron, Pyramid, Prism, Antiprism
export RotationGroup, Cyclic, Dihedral, Tetrahedral, Octahedral, Icosahedral, grouporder

# Particle species
export ParticleSpecies, SpeciesAndPose
export nsites, bindingsite, bindingsites, graphrep, symmetrynumber

# Assembly system
export BindingRules, interactionmatrix
export ncolors, nspecies, nbonds, bonded_colors, bonded_sites, bonded_species, isinert, species

# Polyforms
export Polyform, nparticles, bindingrules, composition
export bonds, bondindex, interior_edges, exterior_edges
public canonbindingsite, canonbindingsites

# Enumeration
export ACCEPT, REJECT, BREAK
export RSStatus, Finished, MaxVerticesReached, MaxDepthReached, BreakTriggered
export polyenum, polygen, countpolyforms, PolyformCount

# Species
export PolygonParticleSpecies, UnitNgon, UnitTriangle, UnitSquare, UnitHexagon
export PolyhedronParticleSpecies
export UnitTetrahedron, UnitCube, UnitOctahedron, UnitDodecahedron, UnitIcosahedron
export UnitPyramid, UnitPrism, UnitAntiprism
export PatchyParticleSpecies, PatchyDisk, PatchySphere

# Public, but not exported: reach for these as `Roly.faces(p)`, or import them by name.
# They are stable API, but specific enough to a body, an encoding or a species that putting
# them in every user's namespace is not worth it.

# Polyhedron geometry
public corners, faces, facevertices, nfaces, nedges
public facecentroid, facecentroids, facenormal, facenormals, edgemidpoint
public inradius, minedgelength

# Graph encodings and the symmetry they record
public dartencoding, cycleencoding
public rotationgroup, faceorbits, facesym, siteorbits, stabilizerorders, sitelabel
public permutationgroup, check_encoding

# What a bond fixes about relative orientation
public contact_pairing, standard_twist, twistfreedom, twist

# The `ParticleSpecies` interface, implemented rather than called
public isconvex, bounding_radius, could_contact, overlap
public sat_overlap, edgenormals, polyhedron

include("utils.jl")
include("pose.jl")
include("bindingsite.jl")
include("particlespecies.jl")
include("encoding.jl")
include("bindingrules.jl")
include("particle.jl")
include("polyform.jl")
include("enumeration.jl")

include("species/polygonparticlespecies.jl")
include("species/polyhedronparticlespecies.jl")
include("species/patchyparticlespecies.jl")

export ruleeditor
export render, polyformplot, polyformplot!

"""
    ruleeditor(species::ParticleSpecies; output=:rules)

Open a terminal editor that builds binding rules geometrically, and return them once you
accept. Provided by the Tachikoma extension, so it needs `import Tachikoma`, which is preferable
to `using` because Tachikoma exports `render` and `Rect` as well.

Rules can be built two ways, and the editor keeps both on screen at once. Placing particles reads
bonds off the geometry: each step attaches a copy of a species to one of the structure's free
binding sites, choosing which free site to attach to, which site of the incoming particle meets
it, and in 3D which twist of the bond to use. The rules then come from every pair of placed
particles, so contacts a particle makes with neighbors it was not attached to count too, which
is how an arrangement reveals rules you did not set out to specify. Editing the interaction
matrix sets bonds directly, and an edit holds against whatever the geometry goes on to produce,
so the two can be used together.

The editor opens on two panes, the rules and the polyforms they enumerate, with a third for
building structures that `b` brings in. `tab` moves the focus between the panes on screen and the
arrow keys and enter act on whichever has it. The rules are always shown, since they are what the
editor produces.

The editor starts with one instance of `species`. In the rules pane `a` adds another, of the same
shape but distinct colors, and `d` drops the last; in the construction pane, which shows them as
a strip marking the active one, digits `1`-`9` pick which to place. Every instance gets a color range of its own in the rules whether or not it is
ever placed, and one carrying a bond or a placement is not dropped.

With the rules focused, the arrow keys move a cursor over the interaction matrix, rows with up
and down and columns with left and right, and enter turns that pair of colors on or off. The two
colors are picked out in the species drawings above the matrix, so moving the cursor is how two
sites are selected and enter bonds them. The matrix and the list of bonding pairs are both shown,
the matrix saying what is possible and the list what is set; the cursor drives the matrix, being
the editable one, and the list marks whichever pair it is on. A matrix too big to sit beside the
list moves to a full-height pane of its own.

A site is drawn in a shade of its species' own hue, so its color says which particle it belongs
to, and a matrix cell is split between the colors of the two sites it joins, the row's above and
the column's below. Each particle carries a small filled triangle at its center pointing away
from its first binding site, showing which way it is turned, left off wherever the particle is
drawn too small for it to read as a direction.

With the construction focused, the arrow keys walk the cursor around the structure's perimeter
one free site at a time, taking whichever of the current site's two neighbors better matches the
direction pressed, so no site is skipped and no key is dead; `,` and `.` step in perimeter order.
`r`/`R` turn the pending particle by changing which of its sites meets the cursor, `t`/`T` pick
the twist in 3D, enter attaches, backspace undoes, `n` starts a disconnected component, and `c`
clears the drawing while keeping the rules it produced.

Seating a particle against the chosen site can put it against several at once, which is how a
ring closes. Every contact the pending particle would make is labelled before it is committed
and the sidebar counts the ones beyond the site aimed at, so a closure shows in advance. A
placement that would overlap an existing particle is refused, and is drawn crossed out in gray
beforehand rather than looking legal until enter does nothing.

The enumeration pane draws the polyforms the rules allow. It runs on `e` rather than
automatically, since a permissive rules set can take long enough to be felt between keystrokes;
changing the rules marks the last run out of date rather than repeating it. `s`/`S` set the
largest structure to look for and `x`/`X` the number to stop after, which is what bounds the
work. Structures are numbered by their position in the run; with the pane focused the arrow keys
pick one, and the box beside the grid redraws it as large as the pane allows, titled with its
number, naming every site in its own color and writing both colors where two meet, with its
composition vector underneath.

`-` and `=` zoom the construction view and `0` refits it, those being the only pane with a camera.
`q` accepts from anywhere. A line along the foot reports what the editor last did, or why it
declined to, and counts the species, bonds and placed particles.

The construction view holds its zoom while you build, so attaching a particle leaves the rest of
the structure where it was, and it zooms out only when the structure would leave the pane.
Zooming by hand holds the cursor still rather than the middle of the structure, and switches the
automatic fitting off until `0` asks for it again.

Binding sites are drawn with the color the returned rules give them, so a label names the same
row the interaction matrix does. Where two sites meet, both colors are written side by side, the
pending bond included, and the pending particle's remaining sites are named where they sit so
that turning it turns the labels with it.

`output` controls the return value:
- `:rules` (default): a `BindingRules` object, ready for `polyenum`/`polygen`.
- `:bonds`: the `nx4` integer bonds matrix `[species1 site1 species2 site2; ...]`,
  suitable for pasting back into code as `BindingRules(bonds, species)`.
- `:matrix`: the `Symmetric{Bool}` color-indexed interaction matrix, suitable for
  pasting back as `BindingRules(intmat, species)`.

If no bond was formed, returns `nothing` regardless of `output`.

Only 2D species are supported; a 3D structure is better read with [`render`](@ref).
"""
function ruleeditor end

"""
    render(p; hidedecorations=true, kwargs...)

Draw a [`Polyform`](@ref) or [`ParticleSpecies`](@ref) into a fresh figure, picking a 2D or 3D
axis to match its dimension. Provided by the Makie extension, so it needs a backend loaded.

3D output needs a backend with a depth buffer, GLMakie or WGLMakie. CairoMakie sorts
primitives instead of depth-testing them, so 3D polyforms show artifacts where faces meet.

`bindingrules=rules` draws sites no bond can use as inert, which works for a bare species too.
"""
function render end

"""
    polyformplot(p; kwargs...)

Makie recipe behind [`render`](@ref), drawing a polyform or species into a new axis.
Provided by the Makie extension.
"""
function polyformplot end

"""
    polyformplot!(ax, p; kwargs...)

Draw a polyform or species onto an axis you already have, the mutating form of
[`polyformplot`](@ref). Provided by the Makie extension.
"""
function polyformplot! end

end
