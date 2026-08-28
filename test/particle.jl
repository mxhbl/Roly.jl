using Roly: Particle, graphvertices, leadingvertex, speciesindex, overlap, could_contact

@testset "particle" begin
    rules = BindingRules([1 1 1 3; 1 2 1 4], UnitSquare)

    part = Particle(rules, 1; leadingvertex=1)
    @test leadingvertex(part) == 1
    @test speciesindex(part) == 1

    part2 = Particle(rules, 1; leadingvertex=1)
    gvs = graphvertices(part2, rules)
    @test first(gvs) == 1

    part_shifted = Particle(rules, 1; leadingvertex=length(gvs) + 1)
    gvs2 = graphvertices(part_shifted, rules)
    @test first(gvs2) == length(gvs) + 1
    @test length(gvs2) == length(gvs)

    @test nsites(part, rules) == 4

    site1 = bindingsite(part, rules, 1)
    @test first(site1.vertices) == 1
    site1_shifted = bindingsite(part_shifted, rules, 1)
    @test first(site1_shifted.vertices) == length(gvs) + 1

    @test length(collect(bindingsites(part, rules))) == nsites(part, rules)

    v = SVector(3.0, 0.0)
    part_far = translate(part, v)
    @test part_far.pose.x ≈ v

    rot = Angle2d(π/2)
    @test (rot * part).pose.psi ≈ rot * part.pose.psi
    @test (part * rot).pose.psi ≈ part.pose.psi * rot

    @test overlap(part, part2, rules)
    @test !overlap(part, part_far, rules)
    @test !could_contact(part, part_far, rules)

    io = IOBuffer()
    show(io, part)
    @test contains(String(take!(io)), "Particle")
    show(io, MIME"text/plain"(), part)
    @test contains(String(take!(io)), "speciesindex")
end
