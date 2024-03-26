using Test

@testset verbose=true "Crafts" begin
    include("geometry.jl")
    include("structure.jl")
    include("enumeration.jl")
end;