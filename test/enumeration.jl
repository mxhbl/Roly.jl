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

    nstrs_16 = polyenum(I16, max_size=Inf)[1]
    nstrs_137 = polyenum(I137, max_size=Inf)[1]
    nstrs_cyc = polyenum(Icyc, max_size=Inf)[1]

    @test nstrs_16 == 16
    @test nstrs_137 == 137
    @test nstrs_cyc == 283


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

    nstrs_polymino = [polyenum(I_polymino, max_size=i)[1] for i in 1:length(n_polyminoes_cumulative)]
    nstrs_polyiamond = [polyenum(I_polyiamond, max_size=i)[1] for i in 1:length(n_polyiamonds_cumulative)]

    @test nstrs_polymino == n_polyminoes_cumulative
    @test nstrs_polyiamond == n_polyiamonds_cumulative

    # if length(workers()) > 1
    #     rmprocs(workers())
    # end
    # @add_enumworkers 3

    # @test nstrs_16 == polyenum_distributed(I16, max_size=Inf, verbose=false)[1]
    # @test nstrs_137 == polyenum_distributed(I137, max_size=Inf, verbose=false)[1]
    # @test nstrs_cyc == polyenum_distributed(Icyc, max_size=Inf, verbose=false)[1]
    # @test nstrs_polymino == [polyenum_distributed(I_polymino, max_size=i, verbose=false)[1] for i in 1:length(n_polyminoes_cumulative)]
    # @test nstrs_polyiamond == [polyenum_distributed(I_polyiamond, max_size=i, verbose=false)[1] for i in 1:length(n_polyiamonds_cumulative)]
end;