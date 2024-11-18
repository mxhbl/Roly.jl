@testset "assembly_system" begin

    sys1 = AssemblySystem([1 1 2 1; 2 3 3 1], UnitSquareGeometry)
    sys2 = AssemblySystem([3 3 2 2; 2 4 1 1], UnitSquareGeometry)

    @test asys_hash(sys1) == asys_hash(sys2)

    sys3 = AssemblySystem([3 3 2 3; 2 4 1 1], UnitSquareGeometry)

    @test asys_hash(sys1) != asys_hash(sys3)
    
    sys4 = AssemblySystem([1 1 2 1; 2 3 3 1; 1 3 1 3], UnitSquareGeometry)
    sys5 = AssemblySystem([1 1 2 1; 2 3 3 1; 1 1 1 1], UnitSquareGeometry)
    sys6 = AssemblySystem([1 1 2 1; 2 3 3 1; 1 2 1 2], UnitSquareGeometry)
    sys7 = AssemblySystem([1 1 2 1; 2 3 3 1; 1 2 1 3], UnitSquareGeometry)

    @test asys_hash(sys1) != asys_hash(sys4)
    @test asys_hash(sys1) != asys_hash(sys5)
    @test asys_hash(sys1) != asys_hash(sys6)
    @test asys_hash(sys1) != asys_hash(sys7)

    @test asys_hash(sys4) != asys_hash(sys5)
    @test asys_hash(sys4) != asys_hash(sys6)
    @test asys_hash(sys4) != asys_hash(sys7)

    sys8 = AssemblySystem([1 1 1 1;], UnitSquareGeometry)
    sys9 = AssemblySystem([1 1 1 2;], UnitSquareGeometry)

    @test asys_hash(sys8) != asys_hash(sys9)

end