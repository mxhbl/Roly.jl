module Roly

using LinearAlgebra, StaticArrays, SparseArrays
using Base.Iterators, DataStructures
using Graphs, NautyGraphs

export polyenum, polygen
export AssemblySystem, Polyform, rhash, composition, compositions
export PolygonGeometry, UnitTriangleGeometry, UnitSquareGeometry, UnitPentagonGeometry, UnitHexagonGeometry, UnitCubeGeometry

include("utils.jl")
include("geometry_utils.jl")
include("geometry.jl")
include("polyform.jl")
include("assembly_systems.jl")
include("concatenation.jl")
include("enumeration.jl")
end