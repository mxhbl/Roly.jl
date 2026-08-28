@testset "API conventions" begin
    using Roly: SpeciesSiteLoc, ParticleSiteLoc

    exported = [n for n in names(Roly) if Base.isexported(Roly, n)]
    # names Roly defines, as opposed to the handful it re-exports from Rotations and StaticArrays
    ours = filter(exported) do n
        v = getproperty(Roly, n)
        return (v isa Function || v isa Type) && parentmodule(v) === Roly
    end
    @test length(ours) > 60

    # A name reaching every user's namespace has to be one word. Underscores are for internals,
    # where the two halves are usually a qualifier and a noun that only make sense together.
    @test isempty(filter(n -> occursin('_', string(n)) && !endswith(string(n), "!"), ours))

    # Nothing exported takes or yields a bare graph vertex. Which numbering a vertex is in,
    # canonical or original, is an invariant callers inside the package keep; it must not be
    # something a user of the package has to know about.
    for n in (:bondindex, :interior_edges, :exterior_edges, :subpolyform)
        @test Base.ispublic(Roly, n)
        @test !Base.isexported(Roly, n)
    end

    # The functions that do take one demand to be told which numbering, with no default to fall
    # into. `UndefKeywordError` rather than a wrong answer.
    rules = BindingRules([1 1 1 3; 1 2 1 4], UnitSquare)
    dimer = first(p for p in polygen(rules; maxsize=2) if nparticles(p) == 2)
    part = dimer.particles[1]
    @test_throws UndefKeywordError Roly._same_particle(dimer, 1, 2)
    @test_throws UndefKeywordError Roly._isbound_vertex(dimer, part, 1)
    @test_throws UndefKeywordError Roly._vertex_to_particle_site(dimer, 1)

    # The two site addresses are different types, so neither can stand in for the other.
    @test SpeciesSiteLoc(1, 2) != ParticleSiteLoc(1, 2)
    @test !isa(SpeciesSiteLoc(1, 2), ParticleSiteLoc)
    @test_throws MethodError color(rules, ParticleSiteLoc(1, 1))
    @test_throws MethodError raise!(copy(dimer), bindingsite(dimer, first(opensitelocs(dimer))), ParticleSiteLoc(1, 1))
    # but each still destructures like the pair it replaced
    @test (SpeciesSiteLoc(3, 4)...,) == (3, 4)
    @test (ParticleSiteLoc(3, 4)...,) == (3, 4)

    # each address indexes the thing it names
    @test bindingsite(dimer, ParticleSiteLoc(1, 1)) == bindingsite(dimer.particles[1], rules, 1)
    @test bindingsite(rules, SpeciesSiteLoc(1, 3)) == bindingsite(Roly.species(rules, 1), 3)
    @test_throws MethodError bindingsite(dimer, SpeciesSiteLoc(1, 1))
    @test_throws MethodError bindingsite(rules, ParticleSiteLoc(1, 1))

    # `bonds` speaks ParticleSiteLoc, the rules speak SpeciesSiteLoc
    @test eltype(collect(bonds(dimer))) == Pair{ParticleSiteLoc,ParticleSiteLoc}
    # the site accessors return addresses; the site itself is one index away
    @test eltype(opensitelocs(dimer)) == ParticleSiteLoc
    @test eltype(exposedsitelocs(dimer)) == ParticleSiteLoc
    @test opensitelocs(dimer) ⊆ exposedsitelocs(dimer)
    @test all(l -> bindingsite(dimer, l) isa BindingSite, exposedsitelocs(dimer))
    @test eltype(Roly.possible_attachments(rules, 1)) == SpeciesSiteLoc
end
