# Roly.jl
Roly.jl (_Reverse-Search Polyform Enumerator_) is a Julia package for the enumeration of arbitrary polyforms via [reverse search](https://en.wikipedia.org/wiki/Reverse-search_algorithm). It allows you to exhaustively enumerate structures composed out of a set of arbitrarily shaped building blocks in 2D or 3D and provides an interface to define your own building block geometries and binding rules.
## Installation
To install Roly.jl and its dependencies directly from your Julia REPL, first press `]` to enter Pkg mode, and then run
```
pkg> add https://github.com/maxhbl/NautyGraphs.jl
pkg> add https://github.com/maxhbl/Roly.jl
```
Once Roly.jl is installed, it can be imported with
```
using Roly
```
## Basic Usage
Structure enumeration in Roly.jl starts from an `AssemblySystem`, which is a list of building block geometries together with an interaction matrix that specificies which binding site of which building block are allowed to bind.
The allowed structures can then be enumerated with the function `polyenum`, or generated and stored with `polygen`.

### Assembly System Definition
To illustrate the basic process, let's construct an assembly system consisting of four species of triangular building blocks. You can define the binding rules via a matrix, where every row defines an allowed bond in the form `[species_i site_i species_j site_j]`. Roly already comes with definitions for simple polygonal and polyhedral building block geometries, allowing us to specify triangular building blocks via a `UnitTriangleGeometry`. The `AssemblySystem` constructor takes either a list of geometries or a single geometry if all building blocks are identially shaped.
```
bonds = [1 3 2 3;
         2 2 3 2;
         2 1 4 1;
         3 1 4 1]
asys = AssemblySystem(bonds, UnitTriangleGeometry)
```

### Enumeration of Structures
Once you have defined an assembly system, you can use `polyenum` to enumerate all possible structures:
```
n_strs, largest_strsize = polyenum(asys; max_size=20, max_strs=100_000)`
```
The simple system we have chosen here only allows 16 different structures to form. In general however, the number of structures might be unbounded and it is advisable to always impose either a maximal structure size (`max_size`), or maximal number of structures (`max_strs`) to be enumerated. In addition, `polyenum` allows the user to pass functions for processing, selectively storing, or rejecting structures (see docs).

### Generation of Structures
If you want to store all possible structures in memory for further processing, you can either use `polyenum` with an aggregation function, or use `polygen`. `polygen` stores all structures in memory instead of using reverse search, making it slightly faster, but less memory efficient than `polyenum`.
`polygen` uses the same interface as `polyenum`, but simply returns a list of structures:
```
strs = polygen(asys; max_size=20, max_strs=100_000)`
```
The structures can then be further processed, or rendered using
```
draw_polyforms(strs, asys, "<<FILENAME>>.<<png / pdf>>")
```

## Advanced Usage
TODO: make a doc page that includes
- Custom Geometry Definition
- Function interface of `polyenum` / `polygen`
- Accessing `Polyform` attributes
- Utility tools, like cycle-free checking, unboundedness checking, extracting composition of structures
