# Roly.jl

```@meta
CurrentModule = Roly
```

Roly.jl (*Reverse-search Polyform enumerator*) enumerates polyforms: aggregates built by connecting arbitrarily shaped blocks at their binding sites, in 2D or 3D.

Given a set of building blocks and the bonds allowed between their sites, Roly.jl can:

- enumerate every geometrically valid polyform up to a chosen size,
- generate and store polyforms for further processing,
- estimate how many there are when exact enumeration is too expensive,
- build binding rules interactively by attaching particles face to face ([`editrules`](workflow.md#Sketching-rules-interactively)).

## Installation

From the Julia REPL, press `]` to enter Pkg mode, then:

```
pkg> add https://github.com/goodrichgroup/Roly.jl
```

Roly.jl currently does not support Windows.

## Where to next

- [Workflow](workflow.md) covers the basics: defining binding rules, coloring a species' sites, then enumerating, generating, counting and visualizing polyforms.
- [Custom particle species](custom_species.md) explains how to define your own shapes.
- [Orientation and twists](orientation.md) covers what a bond fixes about relative orientation and how to control it.
  Read it if a species does not enumerate what you expected.
- [API reference](api.md) lists every exported name.

## Citation

If you use Roly.jl in your work, please cite:

```
@article{roly2025,
    year = {2025},
    title = {{Accessing Semiaddressable Self-Assembly with Efficient Structure Enumeration}},
    author = {Hübl, Maximilian C. and Goodrich, Carl P.},
    journal = {Physical Review Letters},
    doi = {10.1103/physrevlett.134.058204},
    pages = {058204},
    number = {5},
    volume = {134}
}
```
