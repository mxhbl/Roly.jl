# Workflow

```@meta
CurrentModule = Roly
```

## Defining binding rules

A `BindingRules` object holds a list of particle species and the bonds allowed between their binding sites.
A bond is written `[species_i site_i species_j site_j]`, meaning site `site_i` of species `species_i` may bind site `site_j` of species `species_j`.

```jldoctest workflow
julia> using Roly

julia> bonds = [1 3 2 3;
                2 2 3 2;
                2 1 4 1;
                3 1 4 1];

julia> rules = BindingRules(bonds, UnitTriangle)
2d BindingRules[n=4, k=4]
```

Pass a single species if all blocks have the same shape, otherwise a vector with one species per index used in `bonds`.

In 2D Roly ships [`UnitTriangle`](@ref), [`UnitSquare`](@ref), [`UnitHexagon`](@ref) and [`PolygonParticleSpecies`](@ref) for regular polygons, and [`PatchyDisk`](@ref) for a disk with sites on its rim.
In 3D it ships [`UnitTetrahedron`](@ref), [`UnitCube`](@ref), [`UnitOctahedron`](@ref), [`UnitDodecahedron`](@ref), [`UnitIcosahedron`](@ref), [`UnitPyramid`](@ref), [`UnitPrism`](@ref), [`UnitAntiprism`](@ref) and [`PolyhedronParticleSpecies`](@ref) for any convex polyhedron with one site per face, and [`PatchySphere`](@ref) for a sphere whose patches inherit a polyhedron's rotation group.
See [Custom particle species](custom_species.md) to define your own.

## Colors and symmetry

You describe a species by *coloring* its binding sites, and the bond table refers to those colors.
Which sites the particle cannot tell apart follows from the coloring and the geometry, and Roly derives it.

```jldoctest workflow
julia> symmetrynumber(PolyhedronParticleSpecies(Cube()))                  # every face distinct
1

julia> symmetrynumber(PolyhedronParticleSpecies(Cube(); colors=fill(1, 6)))  # all faces alike
24

julia> caps = [abs(Roly.facenormal(Cube(), i)[3]) > 0.5 ? 2 : 1 for i in 1:6];

julia> symmetrynumber(PolyhedronParticleSpecies(Cube(); colors=caps))     # caps apart from sides
8
```

Faces come in no particular order, so the third example picks the caps by their normals.
The answer 8 is `D_4`: telling two opposite faces apart leaves the 4-fold axis through them and the 2-fold axes across it.

Two more keywords say what a bond at a face *means*, rather than which bonds exist.
`locking` decides whether a site holds its partner in the orientation its frame names, and `twists` turns a site about its normal to pick which orientation that is.
See [Orientation and twists](orientation.md).

## Symmetry groups

Bodies are built by name: [`Cube`](@ref), [`Prism`](@ref), [`Antiprism`](@ref) and the rest all return a [`Polyhedron`](@ref).
[`rotationgroup`](@ref) lists the rotations a body has, and [`grouporder`](@ref) counts those of a named group.

```jldoctest workflow
julia> length(Roly.rotationgroup(Cube()))
24

julia> grouporder(Octahedral())
24

julia> grouporder(Dihedral(5))
10
```

The named groups are [`Cyclic`](@ref), [`Dihedral`](@ref), [`Tetrahedral`](@ref), [`Octahedral`](@ref) and [`Icosahedral`](@ref), the rotation-only point groups a rigid body can have.

## Sketching rules interactively

[`ruleeditor`](@ref) builds a bond table geometrically: grow a structure by attaching blocks face to face and it reads the rules off every pair of touching sites.
It is provided by a package extension, so Tachikoma has to be loaded.
Load it with `import` rather than `using`: Tachikoma exports `render` and `Rect` as well, and `using` both packages would leave those names ambiguous.

```julia
using Roly
import Tachikoma
rules = ruleeditor(UnitSquare)  # any species, in 2D or 3D
```

The editor opens on two panes, the rules and the polyforms they enumerate, with a third for building structures that `b` brings in.
`Tab` moves the focus between the panes on screen, and the arrow keys and `Enter` act on whichever has it; the sidebar keeps a key section per pane so that what is available does not change under you.
The rules are always shown, since they are what the editor produces, and a design can be made as a bond table alone.
A line along the foot reports what the editor last did, or why it declined to, and keeps a running count of species, bonds and placed particles.

### Editing the rules

With the rules pane focused, the arrow keys move a cursor over the interaction matrix, up and down along the rows and left and right along the columns, and `Enter` turns that pair of colors on or off.
The two colors the cursor names are picked out in the species drawings above the matrix, so a cell can be read as the sites it stands for: moving the cursor is how two sites are selected, and `Enter` bonds them.
Both tables are on screen at once, the matrix beside the list of bonding pairs: the matrix says what is possible and the list says what is set, which is what you want once the relation is sparse and the grid is mostly dots.
The cursor drives the matrix, since that is the editable one, and the list simply marks whichever pair it is on.
A matrix too big to sit beside the list moves to a full-height pane of its own, which holds a far larger one than a corner of the rules pane does; if it will not fit even there the rules pane says so rather than leaving a gap.

A site is drawn in a shade of its species' own hue, so its color says which particle it belongs to.
A bond cell in the matrix is split between the colors of the two sites it joins, the row's in the upper half and the column's in the lower.

Each particle carries a small filled triangle at its center pointing away from its first binding site, which is what tells two otherwise identical particles apart when they are turned differently.
A triangle rather than an arrow: a shaft with a head needs more dots than a particle at this scale has to give, while a triangle is a direction and a marker at once.
It is left off wherever the particle is drawn too small for it to read as a direction.

`a` adds a species and `d` drops the last one, each getting its own contiguous range of colors, so adding one widens the matrix by that species' worth of rows whether or not it is ever placed.
Only the last can be dropped, so that the colors of the others do not shift under bonds already set, and one carrying a bond or a placement is kept with a note in the status line rather than silently taking those with it.

### Building a structure

`b` brings in the construction pane, which is the other way to arrive at a bond table: place blocks and the rules are read off the contacts.
Each attachment puts a copy of a species on one of the structure's free binding sites, which takes two choices: where to attach, and how the incoming particle is turned.
The arrow keys walk the cursor around the perimeter, one free site at a time, each press taking whichever of the current site's two neighbors better matches the direction pressed, so a site is never skipped over and no key is ever dead; `,` and `.` step in perimeter order regardless of direction.
`r` and `R` turn the pending particle, which changes which of its sites meets the cursor.
The species appear as a strip along the top of this pane, marking the one the next placement would use, since that is the only thing choosing one decides; `Enter` attaches, backspace undoes, digits `1` to `9` pick the species to place, `n` starts a disconnected component, and `c` clears the drawing without clearing the rules it produced, so an arrangement can be built to discover a bond and then cleared away to build the next one.

Seating a particle against the chosen site can put it against several at once, which is how a ring closes.
Every contact the pending particle would make is labelled before it is committed, and the sidebar counts the ones beyond the site aimed at, so a closure is visible in advance rather than only after the fact.

A placement that would put the particle inside one already there is refused.
That is checked while the pending particle is being drawn as well as when it is committed, so it appears crossed out in a neutral gray and the sidebar names what is in the way, rather than looking legal until `Enter` does nothing.
Gray rather than red, since red is one of the species colors and a red ghost would read as a different species.

The view keeps its zoom as you build, so what is already drawn does not move, and it zooms out only when the structure would leave the pane.
`-` and `=` zoom by hand, holding the cursor still on screen rather than the middle of the structure, so zooming in follows the site being worked on; they belong to this pane, being the only one with a camera.
Zooming by hand also switches the automatic fitting off, since a view chosen deliberately should not be undone by the next attachment; `0` refits to the current contents and hands fitting back.

Editing does not detach the rules from the drawing.
A bond set by hand stays set as the structure grows, and one cleared by hand stays clear even where the geometry keeps producing that contact, so the two halves of the editor can be used together: place blocks to discover the bonds an arrangement implies, then adjust the table directly.

### Previewing the enumeration

The rightmost pane draws the polyforms the current rules allow, smallest first.
It runs only when `e` is pressed, since a permissive rules set can take long enough to be felt between keystrokes; changing the rules marks the last run out of date rather than repeating it, and the pane header says so.
`s` and `S` set the largest structure to look for and `x` and `X` the number to stop after, the latter being what actually bounds the work: the search gives up after that many however many the rules admit.

Each structure is numbered by its position in the run.
With the pane focused the arrow keys pick one, and the box beside the grid redraws the selection as large as the pane allows, titled with its number.
Every site is named there in its own color, with both colors written side by side wherever two meet, so a bond in an enumerated structure can be read straight off against the matrix.
Underneath it is the composition vector: how many particles of each species the structure has, then how many bonds of each bonding color pair.
That box is worth reaching for: a thumbnail is only a few braille dots across, so a small structure can lose edges shorter than a dot, and the same structure enlarged is drawn properly.

```
╭─ Editor ───────────╮╭─ Rules ──────────────────────────────────╮╭─ Construction ───────────── 4 particles ─╮╭─ Enumeration ───────── 200+ ≤ 6 ─╮
│── Any pane ────────││species 1                                 ││▶■1                                       ││1            2   ⢀⣀⣀⣀             │
│tab   next pane     ││                                          ││                                          ││   ⡏⠉⠉⡉⠉⠉⡇       ⢸  ⢸             │
│b     hide build    ││                  ⡏⠉⠉3⠉⠉⡇                 ││                                          ││   ⡇ ⣰⣷⡀ ⡇       ⢸⠤⠤⢼             │
│q     accept        ││                  2 ⣰⣷⡀ 4                 ││                                          ││   ⡇ ⠉⠉⠁ ⡇       ⢸  ⢸             │
│                    ││                  ⡇ ⠉⠉⠁ ⡇                 ││                                          ││────────────     ⠘⠒⠒⠚             │
│── Rules ───────────││                  ⠉⠉⠉1⠉⠉⠁                 ││                                          ││                                  │
│  2 ─ 4  bonded     ││──────────────────────────────────────────││                                          ││3 ⢀⣀⣀⣀⣀⣀⣀⣀   4 ⢀⣀⣀⣀⣀⣀⣀⣀           │
│↑↓←→  cell          ││  1 2 3 4                         1 ─ 2   ││                                          ││  ⢸   ⡇  ⢸     ⢸   ⡇  ⢸           │
│enter bond          ││1 · ▀ ▀ ▀                         1 ─ 3   ││                                          ││  ⠸⠤⠤⠤⡧⠤⠤⢼     ⢸⠤⠤⠤⡧⠤⠤⠼           │
│a / d add/drop      ││2 ▀ · · ▀                         1 ─ 4   ││                                          ││      ⡇  ⢸     ⢸   ⡇              │
│                    ││3 ▀ · · ·                         2 ─ 4   ││             ⢰⠒⠒3⠒⢲⠒⠒3⠒⢲                  ││      ⠓⠒⠒⠚     ⠘⠒⠒⠒⠃              │
│── Build ───────────││4 ▀ ▀ · ·                                 ││             2 ⢠⣧ 24⢠⣧ 4                  ││                                  │
│  at 4.2 site 1     ││                                          ││             ⢸⣀⣉14⣸⣀⣉13⣸                  ││5 ⢀⣀⣀⣀⣀⣀⣀⣀   6  ⢀⣀⣀⡀              │
│                    ││                                          ││             ⢸  ⣀⡄⢸  ⡄ ⢸                  ││  ⢸   ⡇  ⢸      ⢸  ⡇              │
│↑↓←→  site          ││                                          ││             3 ⠉⠛⠇21⠼⠿⠄4                  ││  ⢸⠤⠤⠤⡧⠤⠤⢼      ⢸⠉⠉⡏⠉⢹            │
│, .   step          ││                                          ││             ⢰⠒⠒21⢲⠒⠒1⠒⠚                  ││  ⢸   ⡇  ⢸      ⠈⠉⠉⡏⠉⢹            │
│r / R turn          ││                                          ││             4 ⠹⡿⠁2                       ││  ⠘⠒⠒⠒⠓⠒⠒⠚         ⠓⠒⠚            │
│enter attach        ││                                          ││             ⢸⣀⣀3⣀⣸                       ││                                  │
│bksp  undo          ││                                          ││                                          ││                                  │
│1-9   species       ││                                          ││                                          ││                                  │
│n / c part/clear    ││                                          ││                                          ││                                  │
│                    ││                                          ││                                          ││                                  │
╰────────────────────╯╰──────────────────────────────────────────╯╰──────────────────────────────────────────╯╰──────────────────────────────────╯
                                                                                                                      1 species  4 bonds  4 placed
```

The faint outline at the top is the pending attachment, drawn before it is committed.
The digits are binding site colors, numbered as the returned rules number them rather than per species, so a label names the same row the interaction matrix does.
Where two sites meet, both colors are written side by side, the pending bond at the cursor included, so the pair being created can be read off before it is committed.
The pending particle's other sites are named too, so turning it with `r` turns the labels with it.

Attachment decides only where a particle goes; the rules still come from every pair of particles in the structure.
So closing a ring reports bonds that no attachment named: the session above took three attachments to build a 2×2 block, and the last square turned out to touch a second neighbor, which is the fourth pair in the matrix.

### Three dimensions

A 3D species is drawn in isometric projection: an orthographic view from the direction (1, 1, 1), which is what makes depth easy to resolve, since nothing changes size with distance and a nearer convex particle simply covers what stands behind it.

Particles are drawn as braille outlines, the edges of the faces that turn toward the camera.
Hidden edges go in two steps.
Within a particle, dropping the faces turned away removes exactly the edges on its far side, the body being convex.
Between particles, the drawing runs farthest first and each particle erases the dots its silhouette covers from everything already drawn behind it.
The one exception is the box beside the enumeration grid, which fills the faces instead, in three shades of the species' hue picked by which way each face turns; a filled face is a whole cell, twice the width and four times the height of a braille dot, so it needs that much room to read.

`[` and `]` turn the camera about the vertical axis and `{` and `}` tilt it.
Most of the time neither is needed, because the camera follows the cursor: stepping onto a face on the far side of the structure turns the camera far enough to bring that face into view and no further, in the plane the view direction and the face normal span, so the structure stays recognizable across the move.
The cursor is what you drive and the camera is what follows it, which is why the arrow keys mean the same thing in 3D as in 2D.
`t` and `T` pick which twist of the bond the incoming particle takes, a choice a 2D bond does not leave open.

Only the sites on the near side are named, and a bond between two placed particles is left unlabelled: its two faces meet inside the solid, where a label would sit on whichever particle happens to stand in front of it.
The species drawings in the rules pane name every site regardless, those on the far side dimmed, since that is the reference drawing of the species and a color in the matrix has to be findable in it.

Species with no polyhedron behind them are drawn as the silhouette of their bounding sphere, a circle, in the same way a 2D species with no corners is.
[`PatchySphere`](@ref) therefore comes out as a circle with its near patches named, and occludes as a sphere.
What the drawing does assume is convexity, which every built-in 3D species has.

Below, three cubes bonded 1-2 and 3-4, with a fourth pending, and the enumeration run to size 4.
The selected structure in the box on the right is filled; here its face colors are shown as shading.

```
╭─ Editor ───────────╮╭─ Rules ──────────────────────────────────╮╭─ Construction ───────────── 3 particles ─╮╭─ Enumeration ───────────── 28 ≤ 4 ─╮
│── Any pane ────────││species 1                                 ││▶■1                                       ││1    ⣀⢄⡀     2  ⢀⡠⢄⡀                │
│tab   next pane     ││                    ⣀⢄⡀                   ││                                          ││  ⢠⣒⠉  ⠈⢑⣢     ⡮⢅⡀⢀⡠⠔⠤⣀             │
│b     hide build    ││                 ⢠⣒⠉ 4⠈⢑⣢                 ││                                          ││  ⢸ ⠉⠒⡔⠊⠁⢸     ⡇ ⢸⠓⢄⡀⡠⠔⡇            │
│q     accept        ││                 ⢸ 1⠒⡔2⠁⢸                 ││                                          ││  ⠸⣀  ⡇ ⢀⡸     ⠑⠢⢸  ⢸  ⡇            │
│                    ││                 ⠸⣀5 36⢀⡸                 ││                                          ││    ⠉⠒⠗⠊⠁         ⠉⠒⠼⠒⠉             │
│── Rules ───────────││                   ⠉⠒⠗⠊⠁                  ││                                          ││                                    │
│  3 ─ 4  bonded     ││──────────────────────────────────────────││                                          ││3    ⢀⢄⡀     4                      │
│↑↓←→  cell          ││  1 2 3 4 5 6                     1 ─ 2   ││                ⣀⢄⡀                       ││    ⡾⢥⣠⠼⡆     ⢠⣔⠊⠉⢀⢄⡀⠉⢒⣤            │
│enter bond          ││1 · ▀ · · · ·                     3 ─ 4   ││             ⣠⠔⠊ 6⠈⠑⢤⣀                    ││    ⢇ ⡇⣀⠇     ⢸ ⠉⡾⢥⣠⠼⡆⠁⢸            │
│a / d add/drop      ││2 ▀ · · · · ·                             ││             ⡏⠙⠲⢤⣀⠤⠒⠉4⠉⠒⠤⡀                ││    ⡇⠉⠋⠁⡇      ⠉⠒⢇⡀⣇⡠⠇⠊⠁            │
│                    ││3 · · · ▀ · ·                             ││             ⡇ 3 ⡏⠑⠦⣀⢀⡠⠖⠉⡇                ││    ⠈⠒⠗⠉          ⠈⠁                │
│── Build ───────────││4 · · ▀ · · ·                             ││             ⠑⠢⣀ ⡇ 5 ⡇ 6 ⡇                ││                                    │
│  at 3.2 site 4     ││5 · · · · · ·                             ││                ⠉⠣⢄⡀ ⡇⢀⣀⠤⠒⠤⣀              ││5   ⣀⠤⣀      6     ⢀⢄⡀              │
│                    ││6 · · · · · ·                             ││                 ⡇ ⠈⠑⡶⢍⡀ 2 ⣀⠭⡆            ││   ⡟⠦⣠⠔⠑⢢⡀     ⢀⡠⠔⡾⢥⣀⠬⢳             │
│↑↓←→  site          ││                                          ││                 ⢇⡀6 ⡇ 24⡔⠊⠁ ⡇            ││   ⠣⢄⡏⠙⡞⢉⡇     ⢸⠉⠒⢇ ⡇⢀⡸             │
│, .   step          ││                                          ││                  ⠈⠒⠤⡇ 6 ⡇ 3 ⡇            ││     ⠈⠑⠋⠁⡇     ⠘⠢⢄⡇⠉⠋⠁⢸             │
│r / R turn          ││                                          ││                     ⠈⠑⠤⣀⣇⡠⠒⠉             ││────────────      ⠈⠒⠗⠊⠁             │
│[ ]   turn          ││                                          ││                         ⠁                ││                                    │
│{ }   tilt          ││                                          ││                                          ││                                    │
│enter attach        ││                                          ││                                          ││                                    │
│bksp  undo          ││                                          ││                                          ││                                    │
╰────────────────────╯╰──────────────────────────────────────────╯╰──────────────────────────────────────────╯╰────────────────────────────────────╯
                                                                                                                        1 species  2 bonds  3 placed
```

Pass `output=:bonds` or `output=:matrix` for a copy-pasteable result instead of a `BindingRules`:

```julia
bonds = ruleeditor(UnitSquare; output=:bonds)  # n×4 integer matrix
rules   = BindingRules(bonds, UnitSquare)         # reproduces the same rules
```

## Enumerating polyforms

`polyenum` walks every polyform the rules allow.

```jldoctest workflow
julia> result = polyenum(rules; maxsize=20, maxstrs=100_000);

julia> result.nstructures
16

julia> result.largest_size
5

julia> result.status
Finished::RSStatus = 0
```

Cap `maxsize` (particles per polyform) or `maxstrs` (total polyforms) when the rules allow unbounded growth.
`status` says why the run stopped: `Finished`, `MaxDepthReached`, `MaxVerticesReached` or `BreakTriggered`.

## Storing polyforms

`polygen` returns the polyforms in a `Vector`, sorted by size.

```jldoctest workflow
julia> polys = polygen(rules; maxsize=20);

julia> length(polys)
16
```

## Counting polyforms

`countpolyforms` counts without storing anything, switching to an unbiased sampled estimate when exact enumeration gets too expensive.

```jldoctest workflow
julia> c = countpolyforms(rules);

julia> c.n
16.0

julia> c.exact
true

julia> c.uncertainty
0.0
```

It requires an explicit `maxsize` when the rules allow polyforms of unbounded size.

## Applying constraints

`polyenum` takes a callback that runs at each polyform, receiving it and its size and returning one of three signals:

- `ACCEPT` counts the polyform and keeps exploring it,
- `REJECT` skips the polyform and everything grown from it,
- `BREAK` stops the enumeration.

```jldoctest workflow
julia> constraint(s, _) = composition(s)[4] <= 1 ? ACCEPT : REJECT;

julia> polyenum(constraint, rules).nstructures
14
```

`REJECT` prunes a whole subtree, so the constraint must be monotone.
If a polyform violates it, everything grown from it must violate it too.

## Visualizing polyforms

Load any [Makie](https://docs.makie.org) backend to activate the plotting extension.

```julia
using GLMakie      # or CairoMakie, WGLMakie, ...
render(polys[end]) # display a single polyform
```

[`render`](@ref) picks a 2D or 3D axis to match the polyform.
**Use GLMakie or WGLMakie for 3D**, since CairoMakie sorts primitives instead of depth-testing them and shows seams where faces meet.
2D output is fine in any backend.

A species renders too, which is the quickest way to see how its faces are colored and which way its sites face.

```julia
render(PolyhedronParticleSpecies(Prism(3); colors=[1, 2, 2, 2, 1]))
render(UnitCube; bindingrules=rules)   # sites no bond can use are drawn inert
```

[`polyformplot!`](@ref)`(ax, poly)` draws onto an existing Makie axis.
