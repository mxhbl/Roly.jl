@testset "assembly_system" begin

    sys1 = AssemblySystem([1 1 2 1; 2 3 3 1], UnitSquareGeometry)
    sys2 = AssemblySystem([3 3 2 2; 2 4 1 1], UnitSquareGeometry)

    @test rhash(sys1) == rhash(sys2)

    sys3 = AssemblySystem([3 3 2 3; 2 4 1 1], UnitSquareGeometry)

    @test rhash(sys1) != rhash(sys3)
    
    sys4 = AssemblySystem([1 1 2 1; 2 3 3 1; 1 3 1 3], UnitSquareGeometry)
    sys5 = AssemblySystem([1 1 2 1; 2 3 3 1; 1 1 1 1], UnitSquareGeometry)
    sys6 = AssemblySystem([1 1 2 1; 2 3 3 1; 1 2 1 2], UnitSquareGeometry)
    sys7 = AssemblySystem([1 1 2 1; 2 3 3 1; 1 2 1 3], UnitSquareGeometry)

    @test rhash(sys1) != rhash(sys4)
    @test rhash(sys1) != rhash(sys5)
    @test rhash(sys1) != rhash(sys6)
    @test rhash(sys1) != rhash(sys7)

    @test rhash(sys4) != rhash(sys5)
    @test rhash(sys4) != rhash(sys6)
    @test rhash(sys4) != rhash(sys7)

    sys8 = AssemblySystem([1 1 1 1;], UnitSquareGeometry)
    sys9 = AssemblySystem([1 1 1 2;], UnitSquareGeometry)

    @test rhash(sys8) != rhash(sys9)

end