using Test
using Roly
using Rotations, StaticArrays, LinearAlgebra, Random
using NautyGraphs
using CairoMakie
using Tachikoma: Tachikoma

include("testutils.jl")

@testset verbose=true "Roly" begin
    include("pose.jl")
    include("bindingsite.jl")
    include("encoding.jl")
    include("particle.jl")
    include("bindingrules.jl")
    include("polyform.jl")
    include("enumeration.jl")
    include("utils.jl")
    include("plotting.jl")
    include("custom_species.jl")
    include("symmetry.jl")
    include("ruleeditor.jl")

    @testset verbose=true "species" begin
        include("species/polygonparticlespecies.jl")
        include("species/polyhedronparticlespecies.jl")
        include("species/patchyparticlespecies.jl")
    end
end;
