using Test
using Roly

@testset verbose=true "Roly" begin
    include("geometry.jl")
    include("structure.jl")
    include("enumeration.jl")
end;