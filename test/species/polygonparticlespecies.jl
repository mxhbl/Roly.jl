using Roly: PolygonParticleSpecies, UnitTriangle, UnitSquare, UnitHexagon,
            nsites, dimension, isconvex, numtype, bindingsites, graphrep, setcolors!, color,
            could_contact, overlap, symmetrynumber

@testset "PolygonParticleSpecies" begin
    for n in (3, 4, 6)
        ps = PolygonParticleSpecies(n)
        @test nsites(ps) == n
        @test dimension(ps) == 2
        @test isconvex(ps)
        @test numtype(ps) == Float64
        @test nv(graphrep(ps)) == n

        sites = collect(bindingsites(ps))
        @test length(sites) == n
        r_in = 0.5 * cot(π / n)
        for s in sites
            @test norm(s.pose.x) ≈ r_in atol=1e-10
        end
    end

    ps32 = PolygonParticleSpecies(3, Float32(1.0))
    @test numtype(ps32) == Float32

    ps_labels = PolygonParticleSpecies(4; colors=[1, 1, 1, 1])
    @test nsites(ps_labels) == 4

    ps_colors = PolygonParticleSpecies(4; colors=[1, 2, 1, 2])
    @test color(bindingsite(ps_colors, 1)) == 1
    @test color(bindingsite(ps_colors, 2)) == 2
    @test color(bindingsite(ps_colors, 3)) == 1

    ps = PolygonParticleSpecies(4)
    ps2 = copy(ps)
    @test nsites(ps2) == nsites(ps)
    setcolors!(ps2, [10, 20, 30, 40])
    @test color(bindingsite(ps2, 1)) == 10
    @test color(bindingsite(ps, 1)) != 10

    io = IOBuffer()
    show(io, ps)
    @test contains(String(take!(io)), "PolygonParticleSpecies")

    id = Pose{2}()
    far = Pose(SVector(100.0, 0.0), Angle2d(0.0))
    adj = Pose(SVector(1.0, 0.0), Angle2d(0.0))

    @test could_contact(UnitSquare => id, UnitSquare => adj)
    @test !could_contact(UnitSquare => id, UnitSquare => far)
    @test overlap(UnitSquare => id, UnitSquare => id)
    @test !overlap(UnitSquare => id, UnitSquare => far)

    # Two squares d > rmin+rmin=1.0 but overlapping due to π/4 rotation -- requires SAT.
    @test overlap(UnitSquare => id, UnitSquare => Pose(SVector(1.05, 0.0), Angle2d(π/4)))

    pent = PolygonParticleSpecies(5)
    @test overlap(pent => id, pent => Pose(SVector(0.3, 0.3), Angle2d(π/5)))
    @test !overlap(pent => id, pent => Pose(SVector(10.0, 5.0), Angle2d(π/3)))

    # Two aligned unit squares diagonally offset
    at(x, y, θ=0.0) = Pose(SVector(x, y), Angle2d(θ))
    for (dx, dy) in ((0.9, 0.9), (0.99, 0.99), (0.6, 0.6), (0.9, 0.0), (0.0, 0.0))
        @test overlap(UnitSquare => id, UnitSquare => at(dx, dy))
    end
    for (dx, dy) in ((1.0, 1.0), (1.1, 1.1), (1.01, 0.0), (2.0, 0.5))
        @test !overlap(UnitSquare => id, UnitSquare => at(dx, dy))
    end
    
    # The same for the other two regular tilings, at their own lattice-adjacent orientations.
    @test overlap(UnitHexagon => id, UnitHexagon => at(0.9, 0.9))
    @test overlap(UnitTriangle => id, UnitTriangle => at(0.1, 0.1, π))

    @test UnitTriangle isa PolygonParticleSpecies
    @test nsites(UnitTriangle) == 3
    @test UnitSquare isa PolygonParticleSpecies
    @test nsites(UnitSquare) == 4
    @test UnitHexagon isa PolygonParticleSpecies
    @test nsites(UnitHexagon) == 6

    @test symmetrynumber(UnitSquare) == 1
    @test symmetrynumber(PolygonParticleSpecies(4; colors=[1,1,1,1])) == 4

    @test symmetrynumber(UnitTriangle) == 1
    @test symmetrynumber(PolygonParticleSpecies(3; colors=[1,1,1])) == 3

    @test symmetrynumber(UnitHexagon) == 1
    @test symmetrynumber(PolygonParticleSpecies(6; colors=[1,1,1,1,1,1])) == 6

    # the symmetric counterparts differ from the defaults in coloring alone: one color throughout,
    # so the body keeps its rotation group where a color per edge leaves it none
    for (plain, sym, n) in ((UnitTriangle, SymmetricUnitTriangle, 3),
                            (UnitSquare, SymmetricUnitSquare, 4),
                            (UnitHexagon, SymmetricUnitHexagon, 6))
        @test [color(bindingsite(plain, i)) for i in 1:n] == 1:n
        @test [color(bindingsite(sym, i)) for i in 1:n] == fill(1, n)
        @test [bindingsite(plain, i).pose for i in 1:n] == [bindingsite(sym, i).pose for i in 1:n]
        @test symmetrynumber(plain) == 1
        @test symmetrynumber(sym) == n
    end
end
