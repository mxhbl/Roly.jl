@testset "environments" begin
    using Roly: subpolyform

    # 1-2 dimer: one active site per species
    dimerrules = BindingRules([1 1 2 1], UnitSquare)
    # 1-2-3 chain of squares
    chainrules = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare)
    # threefold-symmetric central triangle capped by up to three copies of a second species
    central = PolygonParticleSpecies(3; colors=[1, 1, 1])
    outer = PolygonParticleSpecies(3)
    trianglerules = BindingRules([1 1 2 1; 1 2 2 1; 1 3 2 1], [central, outer])
    # self-binding square: site 1 binds site 1 of another copy
    selfrules = BindingRules([1 1 1 1], UnitSquare)

    chainstrs = polygen(chainrules; maxsize=3)
    trimer = chainstrs[findfirst(s -> nparticles(s) == 3, chainstrs)]
    speciesof(poly, p) = poly.particles[p].speciesindex
    mid = findfirst(p -> speciesof(trimer, p) == 2, 1:3)
    end1 = findfirst(p -> speciesof(trimer, p) == 1, 1:3)

    # subpolyform: composition, canonical identity of the full subset, recomputed symmetry
    sub = subpolyform(trimer, [p for p in 1:3 if p != end1])
    @test nparticles(sub) == 2
    @test composition(sub) ==
          composition(trimer) .- composition(Polyform(chainrules, 1)) .- [zeros(Int, 3); 1; 0]

    whole = subpolyform(trimer, 1:3)
    @test graphrep(whole) == graphrep(trimer)
    @test symmetrynumber(whole) == symmetrynumber(trimer)

    # the fully capped triangle is threefold symmetric, and subpolyform must recompute sigma
    tristrs = polygen(trianglerules; maxsize=4)
    full = tristrs[findfirst(s -> nparticles(s) == 4, tristrs)]
    @test symmetrynumber(subpolyform(full, 1:4)) == symmetrynumber(full) == 3

    # extraction: the middle particle sees the whole trimer, an end only the dimer
    envmid = ParticleEnvironment(trimer, mid; depth=1)
    @test nv(envmid.graph) == 12
    envend = ParticleEnvironment(trimer, end1; depth=1)
    @test nv(envend.graph) == 8
    @test nv(ParticleEnvironment(trimer, end1; depth=2).graph) == 12

    # the same ball reached from a different polyform gives the same class: in the trimer the ball
    # around the species-1 end is the 1-2 dimer with the onward site unresolved, which is the same
    # marked graph as the standalone dimer with that site unbound
    dimer12 = chainstrs[findfirst(s -> composition(s)[1:3] == [1, 1, 0], chainstrs)]
    root12 = findfirst(p -> speciesof(dimer12, p) == 1, 1:2)
    @test ParticleEnvironment(dimer12, root12; depth=1) == envend

    # symmetric-root pinning: the central triangle's classes are one per cap count (nv 3, 6, 9, 12)
    # regardless of which site carries a cap, plus the outer species unbound/bound (nv 3, 6)
    envs = particleenvironments(trianglerules; depth=1)
    @test sort(nv.(getfield.(envs, :graph))) == [3, 3, 6, 6, 9, 12]

    # enumerator counts, by hand
    @test length(particleenvironments(dimerrules; depth=1)) == 4
    @test length(particleenvironments(trianglerules; depth=1)) == 6
    @test length(particleenvironments(trianglerules; depth=1, maxsize=2)) == 4

    # enumerator vs extraction: the reverse-search route and per-root extraction agree exactly
    for (rules, maxsize) in ((dimerrules, 4), (trianglerules, 6), (chainrules, 5))
        enumerated = Set(particleenvironments(rules; depth=1))
        extracted = Set{eltype(enumerated)}()
        for s in polygen(rules; maxsize)
            union!(extracted, particleenvironments(s; depth=1))
        end
        @test enumerated == extracted
    end

    # streaming emits every class exactly once (the reverse-search parent-map test)
    seen = []
    result = particleenvironments(trianglerules; depth=1) do env, n
        push!(seen, env)
        ACCEPT
    end
    @test allunique(seen)
    @test length(seen) == 6
    @test result.status == Finished

    # batch extraction agrees with single-root extraction
    batch = particleenvironments(trimer; depth=1)
    @test length(batch) == 3
    @test all(batch[p] == ParticleEnvironment(trimer, p; depth=1) for p in 1:3)

    # bond environments: reversal is an involution, distinct bonds give distinct classes
    bes = bondenvironments(envmid)
    @test length(bes) == 2
    @test all(reverse(reverse(b)) == b for b in bes)
    @test bes[1] != bes[2]   # the two chain bonds join different species pairs

    # crop of the full polyform agrees with the crop of the containing particle environment
    b = first(Iterators.filter(b -> mid in (b.first.particle, b.second.particle), bonds(trimer)))
    direct = BondEnvironment(trimer, b; depth=0)
    @test direct in bes || reverse(direct) in bes

    # a bond between two copies of the same species is symmetric
    selfdimer = polygen(selfrules; maxsize=2)[end]
    @test nparticles(selfdimer) == 2
    se = BondEnvironment(selfdimer, first(bonds(selfdimer)); depth=1)
    @test reverse(se) == se

    @test_throws ArgumentError bondenvironments(ParticleEnvironment(trimer, mid; depth=0))

    # hierarchy closure: every polyform of the dimer system has radius <= 1, so depth 2 adds nothing
    e1 = Set(particleenvironments(dimerrules; depth=1))
    e2 = particleenvironments(dimerrules; depth=2)
    @test length(e2) == 4
    @test Set(PolyformEnvironment(e.graph, e.rootvertices, 1, e.rules) for e in e2) == e1

    # crop reproduces direct extraction at the smaller radius, on every root of every polyform
    for s in polygen(chainrules; maxsize=5), p in 1:nparticles(s)
        env2 = ParticleEnvironment(s, p; depth=2)
        @test crop(env2, 1) == ParticleEnvironment(s, p; depth=1)
        @test crop(env2, 2) === env2
        @test nv(crop(env2, 0).graph) == nv(graphrep(Roly.species(chainrules, 1)))
    end
    # enumerated depth-2 classes crop onto valid depth-1 classes
    e1chain = Set(particleenvironments(chainrules; depth=1))
    @test all(crop(e, 1) in e1chain for e in particleenvironments(chainrules; depth=2))
    @test_throws ArgumentError crop(first(e1chain), 2)

    # bond-environment crops and root extraction agree with direct computation on every bond
    for s in polygen(chainrules; maxsize=5), b in bonds(s)
        be = BondEnvironment(s, b; depth=1)
        @test crop(be, 0) == BondEnvironment(s, b; depth=0)
        @test crop(be, 1) === be
        @test rootenvironment(be, 1, 1) == ParticleEnvironment(s, b.first.particle; depth=1)
        @test rootenvironment(be, 2, 1) == ParticleEnvironment(s, b.second.particle; depth=1)
        @test rootenvironment(be, 1, 0) == ParticleEnvironment(s, b.first.particle; depth=0)
    end
end
