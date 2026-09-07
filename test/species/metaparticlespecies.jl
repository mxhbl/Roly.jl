@testset "MetaParticleSpecies" begin
    using Roly: setcolors!, opensitelocs, overlap, polyform
    using Graphs: ne

    # one species chaining through opposite sites, so a dimer keeps one open end of each kind
    chainlike = BindingRules([1 1 1 3], UnitSquare)
    chainstrs = polygen(chainlike; maxsize=4)
    dimer = chainstrs[findfirst(s -> nparticles(s) == 2, chainstrs)]
    open = [bindingsite(dimer, l) for l in opensitelocs(dimer)]
    @test length(open) == 2

    mp = MetaParticleSpecies(dimer)
    @test dimension(mp) == 2
    @test nsites(mp) == length(open)
    @test polyform(mp) == dimer
    @test [color(bindingsite(mp, i)) for i in 1:nsites(mp)] == [color(s) for s in open]
    @test [bindingsite(mp, i).pose for i in 1:nsites(mp)] == [s.pose for s in open]
    @test bounding_radius(mp) >= maximum(norm(p.pose.x) for p in dimer.particles)
    # `_ndistincttwists` needs the twist freedoms, which are 1 in 2D
    @test stabilizerorders(mp) == ones(Int, nsites(mp))
    @test _ndistincttwists(bindingsite(mp, 1), bindingsite(mp, 2)) == 1

    @testset "the encoding is the polyform's own graph" begin
        # bound sites stay as interior vertices, so the graph is larger than the site count
        @test nv(graphrep(mp)) == nv(graphrep(dimer))
        @test nv(graphrep(mp)) > nsites(mp)
        @test ne(graphrep(mp)) == ne(graphrep(dimer))

        # a symmetric polyform: two squares joined site1-site1, swapping them is an automorphism
        selfrules = BindingRules([1 1 1 1; 1 2 1 2], UnitSquare)
        selfstrs = polygen(selfrules; maxsize=2)
        sym = selfstrs[findfirst(s -> nparticles(s) == 2, selfstrs)]
        @test symmetrynumber(sym) == 2
        nopen = length(opensitelocs(sym))
        @test nopen == 2

        # the species inherits the polyform's symmetry exactly -- an encoding over the sites alone
        # could not have known it
        @test symmetrynumber(MetaParticleSpecies(sym)) == symmetrynumber(sym)
        @test symmetrynumber(MetaParticleSpecies(sym; colors=fill(9, nopen))) == 2
        # colours a bond can tell apart break it
        @test symmetrynumber(MetaParticleSpecies(sym; colors=1:nopen)) == 1
    end

    @testset "a meta-assembly is the assembly its particles would have formed" begin
        metasys = BindingRules([1 1 1 2], MetaParticleSpecies(dimer))
        mstrs = polygen(metasys; maxsize=2)
        pair = mstrs[findfirst(s -> nparticles(s) == 2, mstrs)]
        quad = chainstrs[findfirst(s -> nparticles(s) == 4, chainstrs)]

        gm, gp = graphrep(pair), graphrep(quad)
        @test nv(gm) == nv(gp) && ne(gm) == ne(gp)
        # identical up to the labeling, which differs by construction: open-site vertices carry
        # their colors and every species' labels are shifted into a global range
        flatten(g) = (h = copy(g); setlabels!(h, fill(Cint(1), nv(h))); nauty(h; canonize=true); h)
        @test flatten(gm) == flatten(gp)
    end

    @testset "choosing which sites to expose" begin
        # every unbound site is offered, bound ones never are
        @test length(exposedsitelocs(dimer)) == 2 * nsites(UnitSquare) - 2
        @test all(((p, k),) -> 1 <= p <= nparticles(dimer), exposedsitelocs(dimer))
        # the open ones are the non-inert ones, and are what the default exposes, in this order
        @test length(opensitelocs(dimer)) == length(open)
        @test opensitelocs(dimer) ==
              filter(l -> !Roly.isinert(chainlike, color(bindingsite(dimer, l))), exposedsitelocs(dimer))
        @test [color(bindingsite(mp, i)) for i in 1:nsites(mp)] ==
              [color(bindingsite(dimer, l)) for l in opensitelocs(dimer)]
        # the two views agree site for site
        @test [bindingsite(dimer, l) for l in exposedsitelocs(dimer)] ==
              [bindingsite(dimer.particles[p], chainlike, k) for (p, k) in exposedsitelocs(dimer)]

        # a site that is inert inside the polyform becomes usable by being named and given a
        # colour the new rules speak
        inert = setdiff(exposedsitelocs(dimer), opensitelocs(dimer))
        @test !isempty(inert)
        activated = MetaParticleSpecies(dimer, inert[1:2]; colors=[1, 2])
        @test nsites(activated) == 2
        @test [color(bindingsite(activated, i)) for i in 1:2] == [1, 2]
        # and those sites really do bind under rules written over the new colors
        # `exposeinert` names the same set without spelling it out
        @test nsites(MetaParticleSpecies(dimer; exposeinert=true)) == length(exposedsitelocs(dimer))
        @test nsites(MetaParticleSpecies(dimer)) == length(opensitelocs(dimer))

        activerules = BindingRules([1 1 1 2], activated)
        @test [polyenum(activerules; maxsize=k)[1] for k in 1:3] == [1, 2, 3]
        # the squares behind those sites cannot meet under `chainlike`, so recasting onto it
        # fails; the rules those meta bonds induce bond the pair they stand for, and it works
        activepair = polygen(activerules; maxsize=2)[end]
        @test_throws ArgumentError recast(activepair, chainlike)
        projected = inducedrules(activerules)
        @test nbonds(projected) == nbonds(chainlike) + 1
        @test nparticles(recast(activepair, projected)) == 2nparticles(dimer)

        # exposure order is the caller's
        pair = opensitelocs(dimer)
        @test [color(bindingsite(MetaParticleSpecies(dimer, reverse(pair)), i)) for i in 1:2] ==
              reverse([color(bindingsite(dimer, l)) for l in opensitelocs(dimer)])

        bound = [ParticleSiteLoc(p, k) for p in 1:nparticles(dimer) for k in 1:nsites(UnitSquare)
                 if ParticleSiteLoc(p, k) ∉ Set(exposedsitelocs(dimer))]
        @test length(bound) == 2
        @test_throws ArgumentError MetaParticleSpecies(dimer, bound[1:1])
        @test_throws ArgumentError MetaParticleSpecies(dimer, [(1, 99)])
        @test_throws ArgumentError MetaParticleSpecies(dimer, [(99, 1)])
        @test_throws ArgumentError MetaParticleSpecies(dimer, [pair[1], pair[1]])
        @test_throws ArgumentError MetaParticleSpecies(dimer, Tuple{Int,Int}[])
    end

    @testset "exposure decides the symmetry" begin
        # the symmetric stack again: exposing both equivalent ends keeps the swap, exposing one
        # of them cannot
        selfrules = BindingRules([1 1 1 1; 1 2 1 2], UnitSquare)
        sym = let strs = polygen(selfrules; maxsize=2)
            strs[findfirst(s -> nparticles(s) == 2, strs)]
        end
        live = opensitelocs(sym)
        @test length(live) == 2
        @test symmetrynumber(MetaParticleSpecies(sym, live)) == 2
        @test symmetrynumber(MetaParticleSpecies(sym, live[1:1])) == 1
    end

    # a meta-particle whose sites sit on its bounding sphere -- one wrapping a patchy particle,
    # whose radius is the patch radius -- meets a neighbour at exactly the sum of the two radii.
    # The generic `could_contact` asks for strictly less than that, so without a method of its own
    # the pair is never looked at and the meta-species bonds to nothing
    let patchy = PatchyDisk([0.0, π/2, π, 3π/2]; colors=[1, 2, 1, 2])
        prules = BindingRules([1 1 1 2], patchy)
        pmono = first(polygen(prules; maxsize=1))
        base = [p for p in polygen(prules; maxsize=2) if nparticles(p) == 2]
        meta = [p for p in polygen(BindingRules(MetaParticleSpecies(pmono)); maxsize=2)
                if nparticles(p) == 2]
        @test !isempty(base)
        @test length(meta) == length(base)
        @test all(nbonds(m) == 1 for m in meta)
    end

    @testset "keywords and rejections" begin
        recolored = MetaParticleSpecies(dimer; colors=[4, 9])
        @test [color(bindingsite(recolored, i)) for i in 1:2] == [4, 9]
        @test_throws DimensionMismatch MetaParticleSpecies(dimer; colors=[1])

        # a saturated polyform exposes nothing: the 1-2-3 chain's ends are inert
        chain = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare)
        trimer = let strs = polygen(chain; maxsize=3)
            strs[findfirst(s -> nparticles(s) == 3, strs)]
        end
        @test isempty(opensitelocs(trimer))
        @test_throws ArgumentError MetaParticleSpecies(trimer)

    end

    @testset "3D: symmetry and twist freedom come from the polyform" begin
        # every face alike, so a cube is 24-fold and each face 4-fold about its normal
        cube = PolyhedronParticleSpecies(Cube(); colors=fill(1, 6))
        @test symmetrynumber(cube) == 24
        @test stabilizerorders(cube) == fill(4, 6)

        stack = BindingRules([1 1 1 1], cube)
        cubedimer = let strs = polygen(stack; maxsize=2)
            strs[findfirst(s -> nparticles(s) == 2, strs)]
        end
        # 4 turns about the stacking axis, times swapping the two cubes
        @test symmetrynumber(cubedimer) == 8

        cmp = MetaParticleSpecies(cubedimer)
        @test dimension(cmp) == 3
        @test symmetrynumber(cmp) == symmetrynumber(cubedimer)
        @test nsites(cmp) == 10
        # a dart-encoded face spans four vertices, which stay contiguous and in order
        @test all(length(bindingsite(cmp, i).vertices) == 4 for i in 1:nsites(cmp))
        @test nv(graphrep(cmp)) == nv(graphrep(cubedimer))

        # the two faces opposite the bond keep the whole 4-fold twist, since turning about them
        # carries the stack onto itself; the eight side faces keep none of it
        stabs = stabilizerorders(cmp)
        @test count(==(4), stabs) == 2
        @test count(==(1), stabs) == 8
        @test all(bindingsite(cmp, i).sitesym == 4 for i in 1:nsites(cmp))

        # a cube with distinctly coloured faces has no symmetry to inherit
        plaindimer = let strs = polygen(BindingRules([1 1 1 2], UnitCube); maxsize=2)
            strs[findfirst(s -> nparticles(s) == 2, strs)]
        end
        @test stabilizerorders(MetaParticleSpecies(plaindimer)) == [1, 1]

        # and blocks assemble out of blocks in 3D too
        metasys = BindingRules([1 1 1 1], cmp)
        @test polyenum(metasys; maxsize=2)[1] >= 2
    end

    @testset "setcolors! relabels the sites and leaves the interior alone" begin
        ps = MetaParticleSpecies(dimer)
        sitevs = Set(v for i in 1:nsites(ps) for v in bindingsite(ps, i).vertices)
        interior = [v for v in 1:nv(graphrep(ps)) if v ∉ sitevs]
        before = labels(graphrep(ps))[interior]

        setcolors!(ps, [7, 8])
        @test [color(bindingsite(ps, i)) for i in 1:2] == [7, 8]
        @test labels(graphrep(ps))[interior] == before
        # site labels sit above every interior one, so a color can never alias interior structure
        @test all(labels(graphrep(ps))[v] > maximum(before) for v in sitevs)
        @test stabilizerorders(ps) == [1, 1]
        @test_throws ArgumentError setcolors!(ps, [1, 2, 3])

        # equal colors on the two ends make them interchangeable only if the polyform agrees;
        # this chain dimer is not symmetric, so it stays rigid
        setcolors!(ps, [5, 5])
        @test symmetrynumber(ps) == symmetrynumber(dimer)
    end

    @testset "overlap delegates to the constituents" begin
        P = Roly.posetype(mp)
        @test overlap(mp => one(P), mp => one(P))
        @test !overlap(mp => one(P), mp => Pose(SVector(50.0, 0.0), Angle2d(0.0)))
        @test !overlap(mp => one(P), mp => Pose(SVector(2.0, 0.0), Angle2d(0.0)))
        @test overlap(mp => one(P), mp => Pose(SVector(0.5, 0.0), Angle2d(0.0)))
    end

    @testset "assembling blocks out of blocks" begin
        metasys = BindingRules([1 1 1 2], MetaParticleSpecies(dimer))
        # the two ends are distinguishable, so blocks chain in exactly one way per length
        @test [polyenum(metasys; maxsize=k)[1] for k in 1:4] == [1, 2, 3, 4]
        blocks = polygen(metasys; maxsize=3)
        @test sort(nparticles.(blocks)) == [1, 2, 3]
        triple = blocks[findfirst(s -> nparticles(s) == 3, blocks)]
        @test nparticles(triple) * nparticles(dimer) == 6
    end

    @testset "a block assembles like the shape it makes" begin
        using Roly: nbonds, nfaces, facenormal

        persize(rules, K) = [count(p -> nparticles(p) == k, polygen(rules; maxsize=K)) for k in 1:K]
        outward(s) = s.pose.psi[:, 1]
        # exposed sites that face opposite ways. Unique when the block has one site per direction,
        # which is what lets a bond through it fix the neighbour completely
        function facing(sites)
            pairs = NTuple{2,Int}[]
            for i in eachindex(sites), j in eachindex(sites)
                i < j && isapprox(outward(sites[j]), -outward(sites[i]); atol=1e-8) || continue
                push!(pairs, (i, j))
            end
            return pairs
        end
        # the site directly across the block, with no sideways offset between the two. Needed
        # where a side carries several sites, so that facing alone leaves a choice
        function across(sites)
            return filter(facing(sites)) do (i, j)
                n = outward(sites[i])
                d = sites[i].pose.x - sites[j].pose.x
                return norm(d - dot(d, n) * n) < 1e-8
            end
        end
        # a distinct color per exposed site, bonding only the given pairs
        function keyed(poly, pairing)
            sites = [bindingsite(poly, l) for l in opensitelocs(poly)]
            ps = MetaParticleSpecies(poly; colors=1:length(sites))
            pairs = pairing(sites)
            @test sort(collect(Iterators.flatten(pairs))) == 1:length(sites)
            return BindingRules(reduce(vcat, [[1 a 1 b] for (a, b) in pairs]), ps)
        end

        ### several meta-species lifted together, bonding across each other as well as within
        mono = MetaParticleSpecies(Polyform(chainlike, 1))
        di = MetaParticleSpecies(dimer)
        both = BindingRules([mono, di])
        @test nspecies(both) == 2
        metas = polygen(both; maxsize=2)
        # a chain of squares either way, and a meta-pair holds as many squares as its copies do
        chains = Set(graphrep(p) for p in polygen(chainlike; maxsize=4))
        flat = [recast(m, chainlike) for m in metas]
        @test all(graphrep(f) in chains for f in flat)
        @test all(zip(metas, flat)) do (m, f)
            nparticles(f) == sum(nparticles(polyform(species(both, p.speciesindex))) for p in m.particles)
        end
        # 3 squares can only be a monomer bonded to a dimer, so the lift crossed the two species
        @test sort(unique(nparticles.(flat))) == [1, 2, 3, 4]
        # a lifted system induces the rules it was lifted from, and nothing more
        @test interactionmatrix(inducedrules(both)) == interactionmatrix(chainlike)

        # there is nothing to lift when the polyforms were built under different rules
        @test_throws ArgumentError BindingRules([mono, MetaParticleSpecies(Polyform(BindingRules([1 1 1 1],
                                                                                                UnitSquare), 1))])
        # nor when a recoloring has taken the sites out of the underlying rules' vocabulary
        @test_throws ArgumentError BindingRules(MetaParticleSpecies(dimer; colors=[99, 100]))

        ### wrapping a single particle changes nothing
        for (rules, K) in ((BindingRules([1 1 1 3; 1 2 1 4], UnitSquare), 5),
                           (BindingRules([1 1 1 1], PolygonParticleSpecies(3; colors=fill(1, 3))), 5),
                           (BindingRules([1 1 1 1],
                                         PolyhedronParticleSpecies(Cube(); colors=fill(1, 6))), 4))
            wrapped = BindingRules(MetaParticleSpecies(Polyform(rules, 1)))
            @test persize(wrapped, K) == persize(rules, K)
        end

        # Which lattice-animal count a system gives is set by its coloring and its dimension. A
        # distinct color on every site leaves the tile no symmetry, so structures are counted up
        # to translation alone -- the *fixed* sequences. One color throughout leaves the tile its
        # own symmetry; in 2D the encoding is a digraph, so a reflection is not an automorphism
        # and the counts are *one-sided*, while in 3D flipping a flat arrangement over is a proper
        # rotation, which makes them *free*. All four below extend correctly by one more term.

        ### two triangles make a rhombus, which tiles the plane the way a square does
        iamonds = BindingRules([1 1 1 1], PolygonParticleSpecies(3; colors=fill(1, 3)))
        rhombus = only(p for p in polygen(iamonds; maxsize=2) if nparticles(p) == 2)
        @test length(opensitelocs(rhombus)) == 4
        # fixed polyominoes, https://oeis.org/A001168
        @test persize(keyed(rhombus, facing), 5) ==
              persize(BindingRules([1 1 1 3; 1 2 1 4], UnitSquare), 5) == [1, 2, 6, 19, 63]

        ### six triangles make a hexagon, which tiles like one under either coloring
        ring = only(p for p in polygen(iamonds; maxsize=6)
                    if nparticles(p) == 6 && length(opensitelocs(p)) == 6)
        # a distinct color per site leaves the block no symmetry: fixed polyhexes,
        # https://oeis.org/A001207
        @test persize(keyed(ring, facing), 4) ==
              persize(BindingRules([1 1 1 4; 1 2 1 5; 1 3 1 6], UnitHexagon), 4) == [1, 3, 11, 44]
        # the colors it inherits are all one, so it keeps the hexagon's full 6-fold symmetry, and
        # the counts are one-sided: https://oeis.org/A006535
        @test symmetrynumber(MetaParticleSpecies(ring)) == 6
        @test persize(BindingRules(MetaParticleSpecies(ring)), 5) ==
              persize(BindingRules([1 1 1 1], PolygonParticleSpecies(6; colors=fill(1, 6))), 5) ==
              [1, 1, 3, 10, 33]

        ### the same construction in 3D: six triangular prisms make a hexagonal prism, again
        ### fully symmetric. A reflection of a flat arrangement is a rotation about an in-plane
        ### axis here, so the very same coloring gives the free counts rather than the one-sided
        sides(p) = [i for i in 1:nfaces(p) if abs(facenormal(p, i)[3]) < 1e-8]
        sticky(shp, s) = PolyhedronParticleSpecies(shp; colors=[i in s ? 1 : 2 for i in 1:nfaces(shp)])
        tri3, hex3 = Prism(3, 1.0; h=2.0), Prism(6, 1.0; h=2.0)
        prismrules(shp) = (s = sides(shp); BindingRules([1 first(s) 1 first(s)], sticky(shp, s)))
        ring3 = only(p for p in polygen(prismrules(tri3); maxsize=6)
                     if nparticles(p) == 6 && length(opensitelocs(p)) == 6)
        mp3 = MetaParticleSpecies(ring3)
        @test nsites(mp3) == 6
        # the ring is as symmetric as the hexagonal prism it makes, D_6 of order 12
        @test symmetrynumber(mp3) == symmetrynumber(ring3) == 12
        # free polyhexes, https://oeis.org/A000228
        @test persize(BindingRules(mp3), 4) == persize(prismrules(hex3), 4) == [1, 1, 3, 7]

        # the same split without a block, on plain squares: one-sided in 2D, free in 3D
        # https://oeis.org/A000988 and https://oeis.org/A000105
        @test persize(BindingRules([1 1 1 1], PolygonParticleSpecies(4; colors=fill(1, 4))), 5) ==
              [1, 1, 2, 7, 18]
        @test persize(prismrules(Prism(4, 1.0; h=2.0)), 5) == [1, 1, 2, 5, 12]

        ### a 2x2 block of squares, whose sides carry two sites each
        sqrules = BindingRules([1 1 1 3; 1 2 1 4], UnitSquare)
        block = only(p for p in polygen(sqrules; maxsize=4)
                     if nparticles(p) == 4 && length(opensitelocs(p)) == 8)
        # keeping the colors it inherits lets a block meet its neighbour half a block over,
        # sharing one edge instead of two, which no single square can do
        loose = polygen(BindingRules(MetaParticleSpecies(block)); maxsize=2)
        @test count(m -> nparticles(m) == 2, loose) > 1
        @test any(m -> nbonds(m) == 1, loose)
        # every one of them is still an assembly of squares
        direct = Set(graphrep(p) for p in polygen(sqrules; maxsize=8))
        flatloose = [recast(m, sqrules) for m in loose]
        @test all(graphrep(f) in direct for f in flatloose)
        @test all(nparticles(f) == 4nparticles(m) for (m, f) in zip(loose, flatloose))
        @test all(nbonds(f) == 4nparticles(m) + nbonds(m) for (m, f) in zip(loose, flatloose))

        ### the numbering: a meta-polyform's own vertices are already the recast polyform's, which
        ### is what lets the tiling cells name sites without renumbering anything
        for (m, f) in zip(loose, flatloose)
            @test nv(graphrep(f)) == nv(graphrep(m))
            # every vertex range a meta-site occupies is the range of a site of the recast polyform
            recastranges = Set(bindingsite(f, ParticleSiteLoc(p, k)).vertices
                               for p in 1:nparticles(f) for k in 1:nsites(f.particles[p], sqrules))
            @test all(bindingsite(m, i).vertices in recastranges for i in 1:nsites(m))
            # each copy's first particle starts where that copy's own vertex block starts
            underlying = Roly._underlying_particles(m)
            @test issubset([q.leadingvertex for q in m.particles], [q.leadingvertex for q in underlying])
            @test [q.leadingvertex for q in underlying] == [q.leadingvertex for q in f.particles]
            # and the bonds between copies are bonds of the recast polyform, by the same ranges
            recastbonds = Set{NTuple{2,UnitRange{Int}}}()
            for (a, b) in bonds(f)
                r1, r2 = bindingsite(f, a).vertices, bindingsite(f, b).vertices
                push!(recastbonds, (r1, r2))
                push!(recastbonds, (r2, r1))
            end
            @test all((bindingsite(m, a).vertices, bindingsite(m, b).vertices) in recastbonds
                      for (a, b) in bonds(m))
        end

        # giving the two sites of a side different colors leaves the offset nothing to bond to,
        # and the block tiles exactly as a square does
        @test persize(keyed(block, across), 5) == [1, 2, 6, 19, 63]

        # meta rules need not stand for bonds of the underlying rules. these bond a west site to
        # a north one, which no two squares can do under `sqrules`
        bent = MetaParticleSpecies(block)
        west, north = findfirst(i -> color(bindingsite(bent, i)) == 2, 1:nsites(bent)),
        findfirst(i -> color(bindingsite(bent, i)) == 3, 1:nsites(bent))
        pair = first(p for p in polygen(BindingRules([1 west 1 north], bent); maxsize=2)
                     if nparticles(p) == 2)
        @test_throws ArgumentError recast(pair, sqrules)
        # the induced rules say what the meta bond claims about squares, and then it recasts
        bentrules = inducedrules(bindingrules(pair))
        @test nbonds(bentrules) == nbonds(sqrules) + 1
        @test nparticles(recast(pair, bentrules)) == 8

        ### recasting is a change of lens, not something meta-particles own: a polyform reads
        ### under any rules that permit it, and refuses under rules that do not
        chain3 = only(p for p in polygen(chainlike; maxsize=3) if nparticles(p) == 3)
        wider = recast(chain3, sqrules)
        @test nparticles(wider) == 3
        @test graphrep(wider) in Set(graphrep(p) for p in polygen(sqrules; maxsize=3))
        @test nbonds(wider) == nbonds(chain3) == 2
        # the block bonds two squares side by side, which `chainlike` leaves inert
        @test_throws ArgumentError recast(block, chainlike)
        # and a species standing for itself has to be the species the target has in its place:
        # squares do not become triangles just because both systems have a species 1
        @test_throws ArgumentError recast(chain3, BindingRules([1 1 1 2], PolygonParticleSpecies(3)))

        ### a substitution replaces a species by a polyform, meta-species or not
        @test nparticles(recast(Polyform(chainlike, 1), chainlike;
                                 substitutions=Dict(1 => dimer))) == 2
        @test graphrep(recast(Polyform(chainlike, 1), chainlike;
                               substitutions=Dict(1 => dimer))) == graphrep(dimer)

        # a substitution written by hand has to be a polyform of the target rules, whatever
        # species it happens to wear
        twospecies = BindingRules([1 1 2 3; 1 2 2 4], [UnitSquare, UnitSquare])
        @test_throws ArgumentError recast(Polyform(chainlike, 1), chainlike;
                                          substitutions=Dict(1 => Polyform(twospecies, 1)))

        ### two species, which must not be conflated when the block is recast
        two = BindingRules([1 1 2 3; 1 2 2 4], [UnitSquare, UnitSquare])
        direct2 = Set(graphrep(p) for p in polygen(two; maxsize=6))
        for seed in (p for p in polygen(two; maxsize=2) if nparticles(p) == 2)
            metas = polygen(BindingRules(MetaParticleSpecies(seed)); maxsize=3)
            @test !isempty(metas)
            @test all(graphrep(recast(m, two)) in direct2 for m in metas)
            @test all(composition(recast(m, two))[1:2] == nparticles(m) .* composition(seed)[1:2]
                      for m in metas)
        end

        ### 3D polycubes, where a block keeps the symmetry of the polyform it wraps
        cuberules = BindingRules([1 1 1 1], PolyhedronParticleSpecies(Cube(); colors=fill(1, 6)))
        direct3 = Set(graphrep(p) for p in polygen(cuberules; maxsize=6))
        for seed in (p for p in polygen(cuberules; maxsize=3) if nparticles(p) > 1)
            mp = MetaParticleSpecies(seed)
            @test symmetrynumber(mp) == symmetrynumber(seed)
            metas = polygen(BindingRules(mp); maxsize=nparticles(seed) == 2 ? 3 : 2)
            @test all(graphrep(recast(m, cuberules)) in direct3 for m in metas)
        end
    end

    @testset "copy is independent" begin
        ps = MetaParticleSpecies(dimer)
        before = copy(labels(graphrep(ps)))
        cp = copy(ps)
        setcolors!(cp, [11, 12])
        @test [color(bindingsite(ps, i)) for i in 1:2] != [11, 12]
        @test [color(bindingsite(cp, i)) for i in 1:2] == [11, 12]
        @test labels(graphrep(ps)) == before
    end

    @testset "the encoding is checked like any other species" begin
        # `check_encoding` runs in the constructor now, against the polyform's own rotations
        # rather than a group derived from the sites, which could not see the polyform behind them
        @test Roly.check_encoding(mp) === mp
        @test length(Roly.permutationgroup(mp)) == symmetrynumber(mp)

        selfrules = BindingRules([1 1 1 1; 1 2 1 2], UnitSquare)
        sym = polygen(selfrules; maxsize=2)[end]
        @test length(Roly.rotationgroup(sym)) == 2

        # the polyform's rotations are only the candidates. Coloring the two exposed sites apart
        # costs the meta-species a symmetry the polyform keeps, so its group is a subgroup
        for (cols, order) in (([9, 9], 2), ([4, 5], 1))
            ps = MetaParticleSpecies(sym; colors=cols)
            @test Roly.check_encoding(ps) === ps
            @test length(Roly.rotationgroup(ps)) == order
            @test length(Roly.permutationgroup(ps)) == symmetrynumber(ps) == order
        end

        # and so does exposing only one of the two
        one = MetaParticleSpecies(sym, opensitelocs(sym)[1:1])
        @test length(Roly.rotationgroup(one)) == symmetrynumber(one) == 1
    end

    @testset "recoloring tracks the orbits, not the colors" begin
        # this dimer has no symmetry relating its two open sites, so they sit in different orbits
        # whatever they are colored, and recoloring cannot move the labeling
        ps = MetaParticleSpecies(dimer)
        before = copy(labels(graphrep(ps)))
        setcolors!(ps, [11, 12])
        @test labels(graphrep(ps)) == before
        setcolors!(ps, [7, 7])
        @test labels(graphrep(ps)) == before

        # where a symmetry does relate the two sites, they share an orbit and a label until a
        # coloring tells them apart
        selfrules = BindingRules([1 1 1 1; 1 2 1 2], UnitSquare)
        sym = polygen(selfrules; maxsize=2)[end]
        sitelabels(q) = [labels(graphrep(q))[first(bindingsite(q, i).vertices)] for i in 1:nsites(q)]
        together = MetaParticleSpecies(sym; colors=[9, 9])
        @test length(unique(sitelabels(together))) == 1
        setcolors!(together, [4, 5])
        @test length(unique(sitelabels(together))) == 2
        @test symmetrynumber(together) == 1
    end
end
