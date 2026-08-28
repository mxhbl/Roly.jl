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
sitetype
particletype
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
Contact
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
setcolors!
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
SymmetricUnitTriangle
SymmetricUnitSquare
SymmetricUnitHexagon
PolyhedronParticleSpecies
polyhedron
UnitTetrahedron
UnitCube
UnitOctahedron
UnitDodecahedron
UnitIcosahedron
SymmetricUnitTetrahedron
SymmetricUnitCube
SymmetricUnitOctahedron
SymmetricUnitDodecahedron
SymmetricUnitIcosahedron
UnitPyramid
UnitPrism
UnitAntiprism
PatchyParticleSpecies
PatchyDisk
PatchySphere
```

### Meta-species

A [`Polyform`](@ref) wrapped as a species, so that assemblies can be built out of assemblies.

```@docs
MetaParticleSpecies
BindingRules(::MetaParticleSpecies)
BindingRules(::AbstractVector{<:MetaParticleSpecies})
polyform
inducedrules
recast
```

## Binding rules

```@docs
BindingRules
SpeciesSite
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
AbstractPolyform
Polyform
ParticleSite
nparticles
bindingrules
composition
canonbindingsite
canonbindingsites
exposedsites
opensites
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

## Tilings and unbounded growth

Whether a rule set admits arbitrarily large structures, and the repeat unit when it does.
See [Bounded and unbounded rules](growth.md).

```@docs
isunbounded
growthwitness
chainstatebound
canchain
tilings
tilingenum
Tiling
unitcell
latticevectors
bondtypes
iscomplete
tilingorder
isunitcell
tilelatticevectors
cantile
```

## Visualization

Provided by the Makie extension, so a backend has to be loaded.

```@docs
render
polyformplot
polyformplot!
```

## Interactive rule editor

```@docs
ruleeditor
```
