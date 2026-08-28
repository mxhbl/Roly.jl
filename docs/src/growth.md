# Bounded and unbounded rules

Some rule sets run out of structures; others grow forever.
[`isunbounded`](@ref) settles which, and when the answer is yes, [`growthwitness`](@ref) hands back the motion that repeats and [`tilings`](@ref) the repeat unit, if there is one.

The figures below use CairoMakie because the documentation builds headless.
For 3D of your own, prefer GLMakie or WGLMakie, which depth-test instead of sorting.

```@setup growth
using Roly, CairoMakie, LinearAlgebra, Rotations
CairoMakie.activate!(type = "png")
```

## A bounded system runs out

Squares whose site 1 binds site 2 can only turn the same way at every bond, so a chain of them folds back on itself.

```@example growth
ring = BindingRules([1 1 1 2], UnitSquare)
isunbounded(ring)
```

`false` here is a proof rather than a search giving up: the motion between two repeats of a chain state is a quarter turn, and every copy of a pure rotation stays the same distance from its axis, so infinitely many of them cannot fit.
The whole system is therefore finite and [`countpolyforms`](@ref) counts it exactly.

```@example growth
countpolyforms(ring)
```

Four structures, and the largest closes into a ring with no open site left.

```@example growth
polys = polygen(ring)
biggest = polys[argmax(nparticles.(polys))]
nparticles(biggest), nbonds(biggest), length(opensitelocs(biggest))
```

```@example growth
render(biggest)
```

## An unbounded system, and its repeat unit

Binding opposite sides of a square instead gives the square lattice.

```@example growth
square = BindingRules([1 2 1 4; 1 3 1 1], UnitSquare)
isunbounded(square)
```

[`growthwitness`](@ref) says *why*: a single particle repeats under a pure translation.
The angle comes back as `2π` rather than `0` because the motion accumulates its turn, which is why the test wraps before asking whether it turns at all.

```@example growth
w = growthwitness(square)
(period = w.period, translation = w.generator.x, angle = rotation_angle(w.generator.psi))
```

Here the repeat is also a *tiling*: one particle whose open sites all close onto translated copies of itself.
[`tilings`](@ref) finds the closures, and each one carries the lattice vectors, which of the rules' bonds it spends, and whether it leaves any site open.

```@example growth
cell = first(polygen(square; maxsize=1))
complete = filter(t -> iscomplete(t), tilings(cell))
first(complete)
```

```@example growth
tilelatticevectors(cell)
```

Two vectors, so the cell tiles the plane.
`tilings` also returns the partial closures — one vector only, which is an infinite strip rather than a tiling — so filter on `complete` when you want the real thing.

```@example growth
length(tilings(cell)), length(complete)
```

The cell itself is the repeat unit:

```@example growth
render(cell; bindingrules=square)
```

A cell that needs several copies is found by raising `maxorder`; with the default of `1` a system whose neighbour is rotated rather than translated reports no tiling at all.

```@example growth
turn = BindingRules([1 1 1 2; 1 3 1 4], UnitSquare)
mono = first(polygen(turn; maxsize=1))
length(tilings(mono; maxorder=1)), length(tilings(mono; maxorder=2))
```

## Counting what cannot be enumerated

An unbounded system has infinitely many structures, so [`countpolyforms`](@ref) needs a `maxsize`.
Up to that size it enumerates exactly while the count is small, and switches to sampling the reverse-search tree beyond `exact_budget`.

```@example growth
countpolyforms(square; maxsize=9)
```

The `≈` and the error bar mark an estimate rather than a count.
Raising `maxsize` costs time but not memory: `maxsize=14` gives about `1.0e7` structures, and the estimate stays within a couple of percent across trials.

## Unbounded without ever repeating

In the plane, growing forever and closing periodically are the same thing.
In space they are not, and regular tetrahedra are the standard example: glued face to face they form the Boerdijk–Coxeter helix, whose twist per tetrahedron is an irrational fraction of a turn.

```@example growth
tetrahedra = PolyhedronParticleSpecies(Tetrahedron(); colors=fill(1, 4))
helixrules = BindingRules([1 1 1 1], tetrahedra)
w = growthwitness(helixrules)
axis = rotation_axis(w.generator.psi)
(angle = rotation_angle(w.generator.psi), acos_minus_two_thirds = acos(-2/3), pitch = dot(w.generator.x, axis))
```

The motion is a screw: it turns by `acos(-2/3)` and climbs by `0.316` along its axis every tetrahedron.
Copies far enough apart along that axis can never touch, which is what makes the growth infinite — and no multiple of that angle is ever a full turn, so the chain never closes.
[`canchain`](@ref), which looks for a *periodic* closure, therefore finds nothing here, while `isunbounded` still says yes.

```@example growth
canchain(helixrules), isunbounded(helixrules)
```

This is the one place the two disagree, and the reason to reach for `isunbounded`.

The helix is visible in the structures themselves. Among the chains of eight tetrahedra, the longest reaches well short of eight straight steps, because it is curling.

```@example growth
function ischain(p)
    degree = zeros(Int, nparticles(p))
    for (a, b) in bonds(p)
        degree[a.particle] += 1
        degree[b.particle] += 1
    end
    return all(<=(2), degree)
end

chains = filter(ischain, polygen(helixrules; maxsize=8))
helix = chains[argmax(nparticles.(chains))]
xs = [p.pose.x for p in helix.particles]
(endtoend = norm(xs[end] - xs[1]), straight = (nparticles(helix) - 1) * norm(w.generator.x))
```

```@example growth
render(helix)
```
