# Roly.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://goodrichgroup.github.io/Roly.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://goodrichgroup.github.io/Roly.jl/dev/)
[![Build Status](https://github.com/goodrichgroup/Roly.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/goodrichgroup/Roly.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/goodrichgroup/Roly.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/goodrichgroup/Roly.jl)


Roly.jl (_<ins>R</ins>everse-Search P<ins>oly</ins>form Enumerator_) is a Julia package for the enumeration of arbitrary polyforms via [reverse search](https://en.wikipedia.org/wiki/Reverse-search_algorithm). It makes it possible to exhaustively enumerate polyforms (aggregates formed by connecting arbitrarily shaped building blocks at their binding sites) in 2D or 3D and provides an interface to define your own building block geometries and binding rules. Roly.jl is under active development, and breaking changes can occur at any time. Because of its dependencies, Roly.jl currently does not work on Windows.

## Installation
To install Roly.jl directly from your Julia REPL, first press `]` to enter Pkg mode, and then run
```
pkg> add https://github.com/goodrichgroup/Roly.jl
```

## Basic Usage
Enumeration in Roly.jl starts from a `BindingRules` object, which is a list of building block geometries together with an interaction matrix that specifies which binding sites are allowed to bind to each other.
The allowed polyforms can then be enumerated with `polyenum`, generated and stored with `polygen`, or counted (exactly or approximately) with `countpolyforms`.

### Defining Binding Rules
To illustrate the basic process, let's construct a system consisting of four species of triangular building blocks. Binding rules are defined as a list of bonds, where every bond is specified in the form `[species_i site_i species_j site_j]`. For example, `[1 3 2 3]` indicates that site 3 of species 1 is allowed to bind to site 3 of species 2. Roly already comes with definitions for simple polygonal particle geometries (e.g. `UnitTriangle`, `UnitSquare`, `UnitHexagon`), convex polyhedra (e.g. `UnitCube`, `UnitPrism(n)`), as well as patchy particles (e.g. `PatchyDisk`, `PatchySphere`). The `BindingRules` constructor takes either a list of geometries or a single geometry if all building blocks are identically shaped.
```julia
using Roly

bonds = [1 3 2 3;
         2 2 3 2;
         2 1 4 1;
         3 1 4 1]
rules = BindingRules(bonds, UnitTriangle)
```

The [documentation](https://goodrichgroup.github.io/Roly.jl/dev/workflow/) lists every built-in geometry and shows how coloring a species' binding sites sets its symmetry. In 3D a bond also fixes how the two blocks are turned relative to one another, which the [orientation page](https://goodrichgroup.github.io/Roly.jl/dev/orientation/) explains. To implement your own particle species, see [custom particle species](https://goodrichgroup.github.io/Roly.jl/dev/custom_species/).

### Sketching Binding Rules interactively
Instead of writing the bonds matrix out by hand, you can build one geometrically with `editrules`, a terminal editor that grows a structure by attaching blocks face to face and reads the binding rules off every pair of touching sites.
It comes from a package extension, so Tachikoma has to be loaded. Load it with `import` rather than `using`: Tachikoma exports `render` and `Rect` too, and `using` both packages would make those names ambiguous.
```julia
using Roly
import Tachikoma
rules = editrules(UnitSquare)  # any species, in 2D or 3D
```

The editor opens on two panes, the rules and the polyforms they enumerate, with a third for building structures that `b` brings in. `Tab` moves the focus between the panes on screen and the arrow keys and `Enter` act on whichever has it, with the sidebar keeping a key section per pane. The rules are always shown, since they are what the editor produces, and a design can be made as a bond table alone. A line along the foot reports what the editor last did and keeps a running count of species, bonds and placed particles.

With the rules focused, the arrow keys move a cursor over the interaction matrix and `Enter` turns that pair of colors on or off. The two colors are picked out in the species drawings above, so moving the cursor is how two sites are selected and `Enter` bonds them. Both tables are shown at once, the matrix beside the list of bonding pairs; the cursor drives the matrix and the list marks whichever pair it is on. A matrix too big to sit beside the list moves to a full-height pane of its own. A site is drawn in a shade of its species' hue, and a matrix cell is split between the colors of the two sites it joins. Each particle carries a small filled triangle at its center pointing away from its first site, to show which way it is turned.

With the construction focused, the arrow keys walk a cursor around the perimeter one free site at a time, taking whichever of its two neighbors better matches the key, so nothing is skipped and no key is dead (`,` and `.` step it in order). The species appear as a strip along the top of the pane, marking the one the next placement would use. `r`/`R` turn the pending particle, `Enter` attaches, backspace undoes, and digits `1`–`9` pick the species to place. Structures are built in construction windows: `n` opens another and `w` cycles through them, and the rules are read off all of them at once, so an arrangement can be kept while the next is built beside it. `c` clears the window in front, taking with it the bonds only that window produced. The rules are read off every pair of placed particles, and every contact a pending particle would make is labelled before it is committed, so a closure that creates more than the one bond you aimed at shows in advance. A placement that would overlap an existing particle is refused, and is drawn crossed out in gray beforehand rather than looking legal until `Enter` does nothing. A bond set or cleared by hand holds against whatever the geometry goes on to produce, so the two ways of working combine.

The rightmost pane draws the polyforms the rules allow. It runs on `e` rather than automatically, since a permissive rules set can take long enough to be felt between keystrokes; `s`/`S` set the largest structure to look for and `x`/`X` the number to stop after. Focus the pane and the arrow keys pick a structure; the box beside the grid redraws it large, naming every site in its own color and writing both colors where two meet, with the composition vector underneath. A typical session looks like this:

```
╭─ Editor ───────────╮╭─ Rules ──────────────────────────────────╮╭─ Construction ───────────── 4 particles ─╮╭─ Enumeration ─────────── 200+ ≤ 6 ─╮
│── Any pane ────────││species 1                                 ││▶■1                                       ││1            2   ⢀⣀⣀⣀               │
│tab   next pane     ││                                          ││                                          ││   ⡏⠉⠉⡉⠉⠉⡇       ⢸  ⢸               │
│b     hide build    ││                  ⡏⠉⠉3⠉⠉⡇                 ││                                          ││   ⡇ ⣰⣷⡀ ⡇       ⢸⠤⠤⢼               │
│q     accept        ││                  2 ⣰⣷⡀ 4                 ││                                          ││   ⡇ ⠉⠉⠁ ⡇       ⢸  ⢸               │
│                    ││                  ⡇ ⠉⠉⠁ ⡇                 ││                                          ││────────────     ⠘⠒⠒⠚               │
│── Rules ───────────││                  ⠉⠉⠉1⠉⠉⠁                 ││                                          ││                                    │
│  1 ╌ 1             ││──────────────────────────────────────────││                                          ││3 ⢀⣀⣀⣀⣀⣀⣀⣀   4 ⢀⣀⣀⣀⣀⣀⣀⣀             │
│↑↓←→  cell          ││  1 2 3 4                         1 ─ 2   ││                                          ││  ⢸   ⡇  ⢸     ⢸   ⡇  ⢸             │
│enter bond          ││1 · ▀ ▀ ▀                         1 ─ 3   ││                                          ││  ⠸⠤⠤⠤⡧⠤⠤⢼     ⢸⠤⠤⠤⡧⠤⠤⠼             │
│a / d add/drop      ││2 ▀ · · ▀                         1 ─ 4   ││                                          ││      ⡇  ⢸     ⢸   ⡇                │
│                    ││3 ▀ · · ·                         2 ─ 4   ││                                          ││      ⠓⠒⠒⠚     ⠘⠒⠒⠒⠃                │
│── Build ───────────││4 ▀ ▀ · ·                                 ││             ⢰⠒⠒3⠒⢲⠒⠒3⠒⢲                  ││                                    │
│  at 4.2 site 1     ││                                          ││             2 ⢠⣧ 24⢠⣧ 4                  ││5 ⢀⣀⣀⣀⣀⣀⣀⣀   6  ⢀⣀⣀⡀                │
│                    ││                                          ││             ⢸⣀⣉14⣸⣀⣉13⣸                  ││  ⢸   ⡇  ⢸      ⢸  ⡇                │
│↑↓←→  site          ││                                          ││             ⢸  ⣀⡄⢸  ⡄ ⢸                  ││  ⢸⠤⠤⠤⡧⠤⠤⢼      ⢸⠉⠉⡏⠉⢹              │
│, .   step          ││                                          ││             3 ⠉⠛⠇21⠼⠿⠄4                  ││  ⢸   ⡇  ⢸      ⠈⠉⠉⡏⠉⢹              │
│r / R turn          ││                                          ││             ⢰⠒⠒21⢲⠒⠒1⠒⠚                  ││  ⠘⠒⠒⠒⠓⠒⠒⠚         ⠓⠒⠚              │
│enter attach        ││                                          ││             4 ⠹⡿⠁2                       ││                                    │
│bksp  undo          ││                                          ││             ⢸⣀⣀3⣀⣸                       ││7            8  ⢀⣀⣀⣀⣀⣀              │
│- = 0 zoom/fit      ││                                          ││                                          ││ ⢸⠉⠉⢹⠉⠉⢹        ⢸  ⡇ ⢸              │
│1-9   species       ││                                          ││                                          ││ ⠸⠤⠤⢼⠤⠤⢼⠤⠤⢤     ⠈⠉⠉⡏⠉⢹              │
│n / w new/next      ││                                          ││                                          ││    ⢸  ⢸  ⢸        ⡏⠉⢹              │
│c     clear/close   ││                                          ││                                          ││    ⠈⠉⠉⠉⠉⠉⠉        ⠓⠒⠚              │
│                    ││                                          ││                                          ││                                    │
╰────────────────────╯╰──────────────────────────────────────────╯╰──────────────────────────────────────────╯╰────────────────────────────────────╯
                                                                                                                        1 species  4 bonds  4 placed
```

The faint outline at the top is the pending attachment, drawn before you commit it. The digits are binding site colors, numbered the way the returned rules number them, so they name the same rows the matrix does. Where two sites meet, both colors are written side by side, including the pending bond at the cursor, so you can read off which pair you are about to create before committing it.

Because the rules are read off every pair of particles rather than only the pairs you attached, closing a ring reveals bonds you never asked for. In the session above three attachments built a 2×2 block, and the last square turned out to touch a second neighbor as well, which is why the matrix carries four color pairs.

A 3D species is drawn in isometric projection, as the edges of the faces that turn toward the camera, shaded by which way each face turns so that the drawing reads as a lit solid. Dropping the faces turned away removes a particle's own hidden edges, and drawing farthest first, each particle erasing what its silhouette covers of the ones behind it, removes the rest. `[`/`]` step the camera through eight viewpoints, shared by every pane, and it never moves on its own: what follows it is the cursor, since the sites offered to attach to are the ones the camera can see. A small set of labelled axes says which way you are looking. `t`/`T` pick which twist of a bond to use, a choice only 3D has.

Pass `output=:bonds` or `output=:matrix` to get a copy-pasteable representation instead of a `BindingRules`, useful for pinning a specific design in code:
```julia
bonds = editrules(UnitSquare; output=:bonds)  # n×4 integer matrix
rules = BindingRules(bonds, UnitSquare)           # reproduces the same rules
```

### Enumeration
Once you have defined a set of binding rules, use `polyenum` to enumerate all allowed polyforms:
```julia
result = polyenum(rules; maxsize=20, maxstrs=100_000)
result.nstructures   # number of polyforms found
result.largest_size  # size of the largest polyform found
result.status        # Finished, MaxDepthReached, MaxVerticesReached, or BreakTriggered
```
The simple system we have chosen here only allows 16 different polyforms to form. In general however, the number of polyforms might be unbounded, so it is advisable to always impose either a maximal size (`maxsize`) or a maximal count (`maxstrs`). To store all polyforms in memory for further processing, use `polygen`, which returns a list sorted by size:
```julia
strs = polygen(rules; maxsize=20, maxstrs=100_000)
```

### Counting
To count polyforms without storing them, or to estimate when full enumeration is too expensive, use `countpolyforms`:
```julia
c = countpolyforms(rules)
c.n            # count (exact or estimated mean)
c.exact        # true if the count is exact
c.uncertainty  # standard error of the estimate (0 if exact)
```
`countpolyforms` enumerates exactly up to a configurable budget and switches to importance-sampled estimation beyond it. Pass `maxsize` for systems that allow unbounded growth.

### Incorporating constraints
It is often desirable to impose additional constraints on generated polyforms. For example, to enumerate only polyforms with at most one particle of species 4:
```julia
constraint(s, n) = composition(s)[4] <= 1 ? ACCEPT : REJECT
polyenum(constraint, rules)
```
Warning: To ensure well-defined behavior, if a polyform `s` violates the constraint, all larger polyforms that can be generated by adding particles to `s` must also violate the constraint.

### Visualization
Roly.jl provides a [Makie](https://docs.makie.org) extension. Load any Makie backend to activate it:
```julia
using GLMakie  # or CairoMakie, WGLMakie, ...
render(s)      # display a single polyform, or a species
```
`render` picks a 2D or 3D axis to match. For 3D use GLMakie or WGLMakie, since CairoMakie sorts primitives rather than depth-testing them. `polyformplot!` can be used to draw onto an existing Makie axis.

## Citation
If you use Roly.jl in your work, please cite [our paper](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.134.058204) below:
```
@article{roly2025, 
         year = {2025}, 
         title = {{Accessing Semiaddressable Self-Assembly with Efficient Structure Enumeration}}, 
         author = {Hübl, Maximilian C. and Goodrich, Carl P.}, 
         journal = {Physical Review Letters}, 
         issn = {0031-9007}, 
         doi = {10.1103/physrevlett.134.058204}, 
         pages = {058204}, 
         number = {5}, 
         volume = {134}, 
         keywords = {}
}
```
