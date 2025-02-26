# Roly.jl

[![Build Status](https://github.com/mxhbl/Roly.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mxhbl/Roly.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mxhbl/Roly.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mxhbl/Roly.jl)


Roly.jl (_Reverse-Search Polyform Enumerator_) is a Julia package for the enumeration of arbitrary polyforms via [reverse search](https://en.wikipedia.org/wiki/Reverse-search_algorithm). It allows you to exhaustively enumerate structures composed out of a set of arbitrarily shaped building blocks in 2D or 3D and provides an interface to define your own building block geometries and binding rules. Roly.jl is under active development, and breaking changes can occur at any time. Because of its depencencies, Roly.jl currently requires a POSIX operating system (this requirement will be lifted in the future).

## Installation
To install Roly.jl and its dependencies directly from your Julia REPL, first press `]` to enter Pkg mode, and then run
```
pkg> add https://github.com/mxhbl/Roly.jl
```

## Basic Usage
Structure enumeration in Roly.jl starts from an `AssemblySystem`, which is a list of building block geometries together with an interaction matrix that specificies which binding site of which building block are allowed to bind.
The allowed structures can then be enumerated with the function `polyenum`, or generated and stored with `polygen`.

### Assembly System Definition
To illustrate the basic process, let's construct an assembly system consisting of four species of triangular building blocks. You can define the binding rules via a matrix, where every row defines an allowed bond in the form `[species_i site_i species_j site_j]`. For example, the row `[1 3 2 3]` indicates that site 3 of species 1 is allowed to bind to site 3 of species 2. Roly already comes with definitions for simple polygonal and polyhedral building block geometries, allowing us to specify e.g. triangular building blocks via a `UnitTriangleGeometry`. The `AssemblySystem` constructor takes either a list of geometries or a single geometry if all building blocks are identically shaped.
```
using Roly

bonds = [1 3 2 3;
         2 2 3 2;
         2 1 4 1;
         3 1 4 1]
asys = AssemblySystem(bonds, UnitTriangleGeometry)
```

### Enumeration of Structures
Once you have defined an assembly system, you can use `polyenum` to enumerate all possible structures:
```
n_strs, largest_strsize = polyenum(asys; max_size=20, max_strs=100_000)
```
The simple system we have chosen here only allows 16 different structures to form. In general however, the number of structures might be unbounded and it is advisable to always impose either a maximal structure size (`max_size`), or maximal number of structures (`max_strs`) to be enumerated. In addition, `polyenum` allows the user to pass functions for processing, selectively storing, or rejecting structures.

### Generation of Structures
If you want to store all possible structures in memory for further processing, you can either use `polyenum` with an aggregation function, or use `polygen`. `polygen` stores all structures in memory, making it less memory efficient than `polyenum`.
`polygen` uses the same interface as `polyenum`, but simply returns a list of structures:
```
strs = polygen(asys; max_size=20, max_strs=100_000)
```

### Visualization
To visualize 2D structures, use [RolyVis.jl](https://github.com/mxhbl/RolyVis.jl):
```
using RolyVis
draw_polyforms(strs, asys, "filename.pdf//png")
```

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
