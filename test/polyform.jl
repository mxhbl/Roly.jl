using Roly: Polyform, nparticles, nsites, bindingrules, symmetrynumber, dimension,
            bonds, bondindex, composition, interior_edges, exterior_edges, tocanon, toorig,
            BindingRules, UnitSquare, nbonds, raise!, lower!, bindingsites, graphrep,
            opensitelocs, collect_attachments, particle_from_leadingvertex,
            PolygonParticleSpecies, species, polygen

@testset "polyform" begin
    # A bond joins two sites, and `contact_pairing` pairs gcd(k1, k2) of their vertices, so a
    # dart-encoded bond reaches the graph as several edges. Counting edges would report one
    # cube-to-cube bond four times, and `composition` would say the same.
    let cubes = BindingRules([1 1 1 1], PolyhedronParticleSpecies(Cube(); colors=fill(1, 6)))
        dimer = first(p for p in polygen(cubes; maxsize=2) if nparticles(p) == 2)
        @test length(collect(Roly.exterior_edges(dimer))) == 4
        @test nbonds(dimer) == 1
        @test length(collect(bonds(dimer))) == 1
        @test composition(dimer) == [2, 1]
        # every bond is listed once, and joins two sites nothing else uses
        for p in polygen(cubes; maxsize=4)
            bs = collect(bonds(p))
            @test length(bs) == nbonds(p)
            @test allunique(bs)
            @test allunique(reduce(vcat, [[b.first, b.second] for b in bs]; init=[]))
        end
    end

    rules = BindingRules([1 1 1 3; 1 2 1 4], UnitSquare)

    empty_poly = Polyform(rules)
    @test nparticles(empty_poly) == 0
    @test nsites(empty_poly) == 0

    mono = Polyform(rules, 1)
    @test nparticles(mono) == 1
    @test nsites(mono) == 4
    @test dimension(mono) == 2
    @test bindingrules(mono) === rules
    @test symmetrynumber(mono) == 1

    @test mono == Polyform(rules, 1)
    @test mono != empty_poly

    mono2 = copy(mono)
    @test mono == mono2
    copy!(empty_poly, mono)
    @test mono == empty_poly

    io = IOBuffer()
    show(io, mono)
    @test contains(String(take!(io)), "Polyform")

    @test !isempty(collect(interior_edges(mono)))
    @test isempty(collect(exterior_edges(mono)))
    @test isempty(collect(bonds(mono)))

    comp = composition(mono)
    @test comp[1] == 1
    @test sum(comp[2:end]) == 0

    @test length(collect(bindingsites(mono))) == nsites(mono)

    open_sites = opensitelocs(mono)
    @test length(open_sites) == 4

    for v in eachindex(mono.canon2orig)
        @test tocanon(mono, toorig(mono, v)) == v
    end

    di = copy(mono)
    site, loc = first(collect_attachments(di))
    @test !isnothing(raise!(di, site, loc))
    @test nparticles(di) == 2
    @test nsites(di) == 8

    @test !isempty(collect(exterior_edges(di)))
    @test length(collect(bonds(di))) == 1

    for e in exterior_edges(di)
        idx = bondindex(di, e.src, e.dst)
        @test !isnothing(idx)
        @test 1 <= idx <= nbonds(rules)
    end

    comp_di = composition(di)
    @test comp_di[1] == 2
    @test sum(comp_di[2:end]) == 1

    lower!(di)
    @test nparticles(di) == 1
    lower!(di)
    @test nparticles(di) == 0
    @test isnothing(lower!(di))

    ### Check raise / remove equality
    rules = BindingRules([1 1 1 4; 1 4 2 2; 1 1 3 2; 1 1 3 4; 2 1 3 3; 2 1 3 1], UnitSquare)
    poly = polygen(rules)[end]

    poly_raise = copy(poly)
    lower!(poly_raise)

    attachments = collect_attachments(poly_raise)
    site, loc = first(attachments)
    i = 1
    while ismissing(raise!(poly_raise, attachments[i]...))
        i += 1
        siteandloc = attachments[i]
    end

    @test poly == poly_raise
    @test poly_raise.canon2orig == poly.canon2orig
    @test poly_raise.orig2canon == poly.orig2canon

    using Graphs: has_edge, nv
    using Roly: graphvertices, PolyhedronParticleSpecies, UnitCube, Cube

    # A monomer's canon2orig is the permutation `canonize!` applied to the species graph,
    # not the identity.
    for spcs in (PolygonParticleSpecies(6; colors=[1, 2, 1, 2, 1, 2]),
                 PolygonParticleSpecies(4; colors=[1, 1, 1, 1]))
        s = BindingRules([1 1 1 3], spcs)
        mono = Polyform(s, 1)
        for e in edges(graphrep(species(s, 1)))
            @test has_edge(graphrep(mono), tocanon(mono, e.src), tocanon(mono, e.dst))
        end
        @test nv(graphrep(mono)) == nv(graphrep(species(s, 1)))
    end

    # Same for a 3D species, whose 24-vertex graph really does get permuted.
    s3 = BindingRules([1 1 1 2], dartspecies(Cube()))
    mono3 = Polyform(s3, 1)
    for e in edges(graphrep(species(s3, 1)))
        @test has_edge(graphrep(mono3), tocanon(mono3, e.src), tocanon(mono3, e.dst))
    end

    # Removing a particle compacts the graph's vertex numbering, so the surviving
    # particles' vertex blocks have to be adjusted
    sys_pm = BindingRules([1 1 1 3; 1 2 1 4], PolygonParticleSpecies(4, 1.0; colors=[1, 1, 1, 1]))
    for p in polygen(sys_pm; maxsize=6)
        nparticles(p) < 2 && continue
        q = copy(p)
        lower!(q)
        blocks = sort(reduce(vcat, [collect(graphvertices(pt, sys_pm)) for pt in q.particles]))
        @test blocks == 1:nv(graphrep(q))
    end
end
