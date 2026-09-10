# API reference

```@meta
CurrentModule = Roly
```

## Poses

```@docs
Pose
dimension
numtype
posetype
```

## Binding sites

```@docs
BindingSite
color
```

### Bonds and twists

How many distinct ways a partner may attach at a site, and which one a given attachment is.
See [Orientation and twists](orientation.md) for what the numbers mean.

```@docs
twistfreedom(::BindingSite)
twistfreedom(::BindingSite, ::BindingSite)
twist
standard_twist
contact_pairing
```

## Polyhedra and graph encodings

```@docs
Polyhedron
corners
faces
facevertices
nfaces
nedges
facecentroid
facenormal
edgemidpoint
minedgelength
inradius
rotationgroup(::Polyhedron)
rotationgroup(::ParticleSpecies)
faceorbits
facesym
siteorbits
sitelabel
permutationgroup
check_encoding
stabilizerorders(::ParticleSpecies)
stabilizerorders(::Any, ::Any, ::Any)
dartencoding
cycleencoding
```

### Rotation groups

```@docs
RotationGroup
Cyclic
Dihedral
Tetrahedral
Octahedral
Icosahedral
grouporder
```

### Solids

```@docs
Tetrahedron
Cube
Octahedron
Dodecahedron
Icosahedron
Pyramid
Prism
Antiprism
```

## Particle species

```@docs
ParticleSpecies
SpeciesAndPose
bindingsite
bindingsites
nsites
graphrep
isconvex
symmetrynumber
bounding_radius
could_contact
overlap
sat_overlap
edgenormals
```

### Built-in species

```@docs
PolygonParticleSpecies
UnitNgon
UnitTriangle
UnitSquare
UnitHexagon
PolyhedronParticleSpecies
polyhedron
UnitTetrahedron
UnitCube
UnitOctahedron
UnitDodecahedron
UnitIcosahedron
UnitPyramid
UnitPrism
UnitAntiprism
PatchyParticleSpecies
PatchyDisk
PatchySphere
```

## Binding rules

```@docs
BindingRules
interactionmatrix
nspecies
nbonds
ncolors
species
bonded_colors
bonded_sites
bonded_species
isinert
```

## Polyforms

```@docs
Polyform
nparticles
bindingrules
composition
canonbindingsite
canonbindingsites
rotationgroup(::Polyform)
permutationgroup(::Polyform)
bonds
bondindex
interior_edges
exterior_edges
```

## Enumeration

```@docs
polyenum
polygen
countpolyforms
PolyformCount
```

An enumeration reports why it stopped as an `RSStatus`: `Finished`, `MaxDepthReached`, `MaxVerticesReached` or `BreakTriggered`.
A callback returns `ACCEPT`, `REJECT` or `BREAK`; see [Applying constraints](workflow.md#Applying-constraints).

## Visualization

Provided by the Makie extension, so a backend has to be loaded.

```@docs
render
polyformplot
polyformplot!
```

## Interactive rule editor

Provided by the Tachikoma extension, so `import Tachikoma` is required.

```@docs
editrules
```
