using Test
using Roly

@testset verbose=true "Roly" begin
    include("geometry.jl")
    include("polyform.jl")
    include("utils.jl")
    include("enumeration.jl")
end;