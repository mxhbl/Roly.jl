@testset "enumeration" begin
    I16 = AssemblySystem(
        [1 3 2 3;
         2 1 4 1;
         2 2 3 2;
         3 1 4 1],
        UnitTriangleGeometry)


    I137 = AssemblySystem(
        [1 2 2 1;
         1 3 3 1;
         1 3 3 2;
         1 3 4 3;
         2 1 3 3;
         3 2 4 2;
         3 1 4 2;
         1 2 5 1], 
        UnitTriangleGeometry)

    Icyc = AssemblySystem(
        [1 3 2 3;
         2 2 3 2;
         3 3 1 1;
         3 3 4 3;
         4 2 1 1],
        UnitTriangleGeometry)

    nstrs_16_enum = polyenum(I16)[1]
    nstrs_137_enum = polyenum(I137)[1]
    nstrs_cyc_enum = polyenum(Icyc)[1]

    @test nstrs_16_enum == 16
    @test nstrs_137_enum == 137
    @test nstrs_cyc_enum == 283

    nstrs_16_gen = polygen(I16) |> length
    nstrs_137_gen = polygen(I137) |> length
    nstrs_cyc_gen = polygen(Icyc) |> length

    @test nstrs_16_gen == 16
    @test nstrs_137_gen == 137
    @test nstrs_cyc_gen == 283

    I_polymino = AssemblySystem(
        [1 1 1 1;
         1 2 1 2;
         1 3 1 3; 
         1 4 1 4], 
         UnitSquareGeometry, 
         ones(4))

    # Number of one-sided polyminoes (https://oeis.org/A000988) [starting from 1]
    n_polyminoes = [1, 1, 2, 7, 18, 60, 196, 704, 2500]
    n_polyminoes_cumulative = cumsum(n_polyminoes)

    I_polyiamond = AssemblySystem(
            [1 1 1 1;
             1 2 1 2;
             1 3 1 3], 
             UnitTriangleGeometry, 
             ones(3))
    
    # Number of one-sided polyiamonds (https://oeis.org/A006534)
    n_polyiamonds = [1, 1, 1, 4, 6, 19, 43, 120, 307, 866]
    n_polyiamonds_cumulative = cumsum(n_polyiamonds)

    I_polycube = AssemblySystem(
        [1 1 1 1; 
         1 2 1 2;
         1 3 1 3;
         1 4 1 4;
         1 5 1 5;
         1 6 1 6], UnitCubeGeometry, ones(Int, 24))

    # Number of one-sided polycubes (https://oeis.org/A006534)
    n_polycubes = [1, 1, 2, 8, 29, 166, 1023]
    n_polycubes_cumulative = cumsum(n_polycubes)

    nstrs_polymino_enum = [polyenum(I_polymino, max_size=i)[1] for i in 1:length(n_polyminoes_cumulative)]
    nstrs_polyiamond_enum = [polyenum(I_polyiamond, max_size=i)[1] for i in 1:length(n_polyiamonds_cumulative)]
    nstrs_polycube_enum = [polyenum(I_polycube, max_size=i)[1] for i in 1:length(n_polycubes_cumulative)]

    @test nstrs_polymino_enum == n_polyminoes_cumulative
    @test nstrs_polyiamond_enum == n_polyiamonds_cumulative
    @test nstrs_polycube_enum == n_polycubes_cumulative

    nstrs_polymino_gen = [length(polygen(I_polymino, max_size=i)) for i in 1:length(n_polyminoes_cumulative)]
    nstrs_polyiamond_gen = [length(polygen(I_polyiamond, max_size=i)) for i in 1:length(n_polyiamonds_cumulative)]
    nstrs_polycube_gen = [length(polygen(I_polycube, max_size=i)) for i in 1:length(n_polycubes_cumulative)]

    @test nstrs_polymino_gen == n_polyminoes_cumulative
    @test nstrs_polyiamond_gen == n_polyiamonds_cumulative
    @test nstrs_polycube_gen == n_polycubes_cumulative
end;