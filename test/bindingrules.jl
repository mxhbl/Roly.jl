using Roly: nspecies, nbonds, nsites, dimension, species,
            interactionmatrix, bonded_colors, bonded_sites, bonded_species,
            sitesofcolor, speciesofcolor, isinert, possible_attachments

using Roly: PolygonParticleSpecies, PolyhedronParticleSpecies, PatchyDisk, Cube, Prism,
            Tetrahedron, Polyhedron, facenormal, nfaces

using Roly: possible_attachments, distinct_attachments, collect_attachments,
            raise!, Polyform, graphrep, PolyhedronParticleSpecies, PolygonParticleSpecies,
            Cube, Prism, nfaces, color, bindingsites

using NautyGraphs: NautyDiGraph
using StaticArrays: SVector

@testset "bindingrules" begin
    rules = BindingRules([1 1 1 3; 1 2 1 4], UnitSquare)

    @test nspecies(rules) == 1
    @test nbonds(rules) == 2
    @test nsites(rules) == 4
    @test dimension(rules) == 2
    @test length(species(rules)) == 1
    @test species(rules, 1) === species(rules)[1]

    imat = interactionmatrix(rules)
    @test size(imat) == (4, 4)
    @test issymmetric(imat)

    bc = bonded_colors(rules)
    @test length(bc) == 2
    @test all(c isa NTuple{2,Int} for c in bc)

    bs = bonded_sites(rules)
    @test length(bs) == 2

    bsp = bonded_species(rules)
    @test length(bsp) == 2
    @test all(==((1, 1)), bsp)

    c1 = color(rules, SpeciesSiteLoc(1, 1))
    c3 = color(rules, SpeciesSiteLoc(1, 3))
    @test imat[c1, c3]
    @test sitesofcolor(rules, c1) == [SpeciesSiteLoc(1, 1)]
    @test speciesofcolor(rules, c1) == 1

    @test !isinert(rules, c1)
    @test !isinert(rules, SpeciesSiteLoc(1, 1))

    @test possible_attachments(rules, c1) == [SpeciesSiteLoc(1, 3)]
    @test possible_attachments(rules, c3) == [SpeciesSiteLoc(1, 1)]

    sys1bond = BindingRules([1 1 1 3], UnitSquare)
    c2 = color(sys1bond, SpeciesSiteLoc(1, 2))
    @test isinert(sys1bond, c2)
    @test isinert(sys1bond, SpeciesSiteLoc(1, 2))

    io = IOBuffer()
    show(io, rules)
    @test contains(String(take!(io)), "BindingRules")

    sys1 = BindingRules([1 1 2 1], UnitTriangle)
    sys2 = BindingRules([1 3 2 2], UnitTriangle)
    sys3 = BindingRules([1 3 3 2], UnitTriangle)

    g1 = graphrep(sys1)
    g2 = graphrep(sys2)
    g3 = graphrep(sys3)

    @test g1 ≃ g2
    @test !(g1 ≃ g3)

    # on-lattice check
    poly(n, a=1.0) = PolygonParticleSpecies(n, a; colors=fill(1, n))
    # Reach past the constructor to a polygon that is not regular.
    stretch(ps, f) = typeof(ps)(copy(ps.g), copy(ps.sites),
                                [SVector(f * c[1], c[2]) for c in ps.corners],
                                ps.rmin, ps.rmax, ps.skin)
    sides(p) = [i for i in 1:nfaces(p) if abs(facenormal(p, i)[3]) < 1e-8]
    faced(p, sticky; kw...) =
        PolyhedronParticleSpecies(p; colors=[i in sticky ? 1 : 2 for i in 1:nfaces(p)], kw...)

    # A system is on-lattice when every species tiles space at one size and no bond can leave
    # the tiling.
    for (name, rules) in [
        ("squares",             BindingRules([1 1 1 1], poly(4))),
        ("triangles",           BindingRules([1 1 1 1], poly(3))),
        ("hexagons",            BindingRules([1 1 1 1], poly(6))),
        ("squares, 2 species",  BindingRules([1 1 2 1], [poly(4), poly(4)])),
    ]
        @test rules._onlattice
    end

    # test for non-tiling
    for (name, rules) in [
        ("pentagons",           BindingRules([1 1 2 1], [poly(5), poly(5)])),
        ("squares, two sizes",  BindingRules([1 1 2 1], [poly(4), poly(4, 2.0)])),
        # Tile together, but five triangles and a square close a ring at 390 degrees, so two
        # cells overlap with their centers apart.
        ("squares + triangles", BindingRules([1 1 2 1], [poly(4), poly(3)])),
        # Equal bounding radius, different edge: comparing radii would call this one lattice.
        ("square + big triangle", BindingRules([1 1 2 1],
                                      [poly(4), poly(3, sqrt(2) / sqrt(3) * 2)])),
        # `corners` is a field, so a caller can build a species the constructor never would.
        # Four sites and the right inradius are not enough; the shape has to be regular.
        ("stretched square",    BindingRules([1 1 1 1], stretch(poly(4), 1.4))),
        # 3D does not opt in.
        ("polycubes",           BindingRules([1 1 1 1], faced(Cube(), 1:6))),
        ("square prisms",       BindingRules([1 1 1 1], faced(Prism(4, 1.0; h=2.0),
                                                             sides(Prism(4, 1.0; h=2.0))))),
        ("tetrahedra",          BindingRules([1 1 1 1],
                                    PolyhedronParticleSpecies(Tetrahedron(); colors=fill(1, 4)))),
        ("patchy disks",        BindingRules([1 1 1 1],
                                    PatchyDisk([0.0, 2π/3, 4π/3]; colors=fill(1, 3)))),
    ]
        @test !rules._onlattice
    end

    # compare with un-shortcutted versions
    for (n, want) in ((3, [1, 1, 1, 4, 6, 19, 43, 120]),    # https://oeis.org/A006534
                      (4, [1, 1, 2, 7, 18, 60, 196]),       # https://oeis.org/A000988
                      (6, [1, 1, 3, 10, 33, 147]))          # https://oeis.org/A006535
        rules = BindingRules([1 1 1 1], poly(n))
        @test rules._onlattice
        counts = [polyenum(rules; maxsize=i)[1] for i in eachindex(want)]
        @test counts == cumsum(want)
        slow = withoutlattice(rules)
        @test !slow._onlattice
        @test [polyenum(slow; maxsize=i)[1] for i in eachindex(want)] == counts
    end
    
    # orbit testing

    # A cube with every face alike has one orbit of six sites, so six compatible mates collapse
    # to one representative. With every face distinct there is nothing to collapse.
    alike = BindingRules([1 1 1 1], PolyhedronParticleSpecies(Cube(); colors=fill(1, 6)))
    @test length(possible_attachments(alike, 1)) == 6
    @test length(distinct_attachments(alike, 1)) == 1

    distinct = BindingRules(reduce(vcat, [[1 i 1 j] for i in 1:6 for j in i:6]),
                            PolyhedronParticleSpecies(Cube()))
    for c in 1:Roly.ncolors(distinct)
        @test length(distinct_attachments(distinct, c)) == length(possible_attachments(distinct, c))
    end

    # check that representatives are enough, and everything else is just duplicates
    function children(poly, sitesof)
        rules = Roly.bindingrules(poly)
        out = Set{NautyDiGraph}()
        for orig_v in poly.canon2orig
            part = Roly.particle_from_leadingvertex(poly, orig_v)
            isnothing(part) && continue
            for k in 1:Roly.nsites(part, rules)
                site = bindingsite(part, rules, k)
                Roly._isbound_vertex(poly, part, first(site.vertices); canonidxs=false) && continue
                Roly.isinert(rules, color(site)) && continue
                for loc in sitesof(rules, color(site))
                    mate = bindingsite(Roly.species(rules, loc.species), loc.site)
                    for r in 0:(Roly._ndistincttwists(site, mate) - 1)
                        trial = copy(poly)
                        ismissing(raise!(trial, site, loc, r)) && continue
                        push!(out, copy(graphrep(trial)))
                    end
                end
            end
        end
        return out
    end

    p3 = Prism(3, 1.0; h=2.0)
    systems = [
        ("cubes",     alike),
        ("prisms",    BindingRules([1 first(sides(p3)) 1 first(sides(p3))],
                          PolyhedronParticleSpecies(p3;
                              colors=[i in sides(p3) ? 1 : 2 for i in 1:nfaces(p3)]))),
        ("hexagons",  BindingRules([1 1 1 1], PolygonParticleSpecies(6, 1.0; colors=fill(1, 6)))),
        ("triangles", BindingRules([1 1 1 1], PolygonParticleSpecies(3, 1.0; colors=fill(1, 3)))),
    ]
    for (name, rules) in systems
        poly = Polyform(rules, 1)
        for _ in 1:3     # monomer, then grow, so the host is asymmetric in later rounds too
            @test children(poly, possible_attachments) == children(poly, distinct_attachments)
            nxt = nothing
            for (site, loc, r) in collect_attachments(poly)
                trial = copy(poly)
                ismissing(raise!(trial, site, loc, r)) && continue
                nxt = trial; break
            end
            isnothing(nxt) && break
            poly = nxt
        end
    end

    # Constructor from an interaction matrix (round-trip via intmat).
    ref = BindingRules([1 1 1 3; 1 2 1 4], UnitSquare)
    im  = Matrix{Bool}(interactionmatrix(ref))
    reconstructed = BindingRules(im, UnitSquare)
    @test interactionmatrix(reconstructed) == interactionmatrix(ref)
    @test nbonds(reconstructed) == nbonds(ref)
    @test bonded_colors(reconstructed) == bonded_colors(ref)

    # Multi-species intmat round-trip.
    ref2 = BindingRules([1 1 2 1], UnitTriangle)  # 2 species x 3 sites = 6 colors
    im2  = Matrix{Bool}(interactionmatrix(ref2))
    reconstructed2 = BindingRules(im2, [copy(UnitTriangle) for _ in 1:nspecies(ref2)])
    @test interactionmatrix(reconstructed2) == interactionmatrix(ref2)
    @test bonded_colors(reconstructed2) == bonded_colors(ref2)

    # Error cases.
    @test_throws ArgumentError BindingRules(falses(3, 4), UnitSquare)          # not square
    @test_throws ArgumentError BindingRules(falses(3, 3), UnitSquare)          # wrong ncolors
    bad_asymmetric = Bool[false true false false;
                          false false false false;
                          false false false false;
                          false false false false]
    @test_throws ArgumentError BindingRules(bad_asymmetric, UnitSquare)        # not symmetric
end
