@testset "plotting" begin
    # every species must have a `plot_particlespecies!`
    # method reachable through `render`, and a polyform of them must draw as well.
    ext = Base.get_extension(Roly, :MakieExt)
    @test !isnothing(ext)

    tri = PolygonParticleSpecies(3)
    disk = PatchyDisk([0.0, 2π / 3, 4π / 3])
    # Distinct colors per face, so a rules table naming one bond leaves the other four inert
    # and the inert-site rendering below has something to grey out.
    cube = PolyhedronParticleSpecies(Cube())
    sphere = PatchySphere(Cube(), 1.0)

    for spcs in (tri, disk, cube, sphere, UnitTetrahedron, UnitDodecahedron, UnitPrism(5))
        @test render(spcs) isa Figure
    end

    # A 3D polyform, drawn particle by particle at its bonded pose.
    opposite(p, i) = findfirst(
        j -> isapprox(dot(facenormal(p, i), facenormal(p, j)), -1; atol=1e-8), 1:nfaces(p)
    )
    pairs = unique([minmax(i, opposite(Cube(), i)) for i in 1:6])
    rules = BindingRules(reduce(vcat, [[1 a 1 b] for (a, b) in pairs]), cube)
    polys = polygen(rules; maxsize=3)
    @test render(polys[end]) isa Figure

    # A 2D polyform, so the shared color resolution is covered in both dimensions.
    sys2d = BindingRules([1 1 1 1], tri)
    polys2d = polygen(sys2d; maxsize=3)
    @test render(polys2d[end]) isa Figure

    ### a meta-particle draws as the polyform it wraps
    metaseed = first(p for p in polys if nparticles(p) == 2)
    mp = MetaParticleSpecies(metaseed)
    @test render(mp) isa Figure
    @test render(mp; bindingrules=BindingRules(mp)) isa Figure
    # and a whole meta-assembly draws too
    @test render(polygen(BindingRules(mp); maxsize=2)[end]) isa Figure

    # its exposed sites keep their colors; everything else is one wash of the species color
    si, exposed = ext._metasites(mp, nothing, nothing, nothing)
    @test length(exposed) == nsites(mp)
    @test Set(keys(exposed)) == Set(Roly.opensitelocs(metaseed))
    parts = ext._metaparts(mp, Roly.Pose{3,Float64}(), nothing, nothing, nothing, nothing)
    @test length(parts) == nparticles(metaseed)
    wash = ext.RGBf(ext.INERT_COLOR)
    for (p, (_, _, _, sitecolor)) in enumerate(parts), k in 1:nsites(cube)
        want = get(exposed, ParticleSiteLoc(p, k), wash)
        @test sitecolor(0, k) == want
    end
    # every exposed color stays clear of the interior one. A species palette ramps from pale to
    # dark, so a tinted species color would land on its pale end and the two would not read apart
    chan(c) = (ext.RGBf(c).r, ext.RGBf(c).g, ext.RGBf(c).b)
    apart(a, b) = sum(abs, chan(a) .- chan(b))
    @test minimum(apart(c, wash) for c in values(exposed)) > 0.25

    # the sites a bond inside the cluster consumes are exactly the washed ones
    @test count(((p, k),) -> get(exposed, ParticleSiteLoc(p, k), wash) == wash,
                [(p, k) for p in 1:nparticles(metaseed) for k in 1:nsites(cube)]) ==
          nparticles(metaseed) * nsites(cube) - nsites(mp)

    # a site the *new* rules leave inert is drawn inert, whatever it was inside the polyform.
    # The default colors repeat across the two copies, so naming one site makes every site
    # sharing its color live -- colors are the whole interface to the rules.
    onebond = BindingRules([1 1 1 2], mp)
    _, inertexposed = ext._metasites(mp, 1, onebond, nothing)
    inertcount = count(i -> Roly.isinert(onebond, SpeciesSiteLoc(1, i)), 1:nsites(mp))
    @test 0 < inertcount < nsites(mp)
    @test count(==(ext.INERT_COLOR), values(inertexposed)) == inertcount

    # the merged mesh is the constituents' meshes, and nothing else
    geom = ext.particlemesh(mp)
    one = ext.particlemesh(cube)
    @test length(geom[1]) == nparticles(metaseed) * length(one[1])
    @test length(geom[2]) == nparticles(metaseed) * length(one[2])

    # a tiling draws as a patch of itself: its cell, repeated along each lattice vector
    lattice = BindingRules([1 2 1 4; 1 3 1 1], UnitSquare)
    tiled = first(filter(iscomplete, tilings(first(polygen(lattice; maxsize=1)))))
    @test render(tiled) isa Figure
    @test render(tiled; repeats=1) isa Figure
    @test render(tiled; ghost=1, showlattice=false) isa Figure
    # the arrow is drawn rather than asked for: shaft plus two barbs, six points per segment pair
    @test length(ext._openarrow(SVector(0.0, 0.0), SVector(1.0, 0.0))) == 6
    @test length(ext._openarrow(SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 2.0))) == 6

    # An explicitly passed `speciesindex` picks the palette instead of being overridden by
    # the default of 1.
    _, _, colors1, bonding1 = ext._resolve_colors(cube, 1, nothing, nothing)
    _, _, colors2, _ = ext._resolve_colors(cube, 2, nothing, nothing)
    @test length(colors1) == length(colors2) == nsites(cube)
    @test colors1 != colors2
    # With no rules in scope every site counts as bonding, so nothing is tinted away.
    @test all(bonding1)

    # By default inert site are greyed out, which is what markers do.
    rules = BindingRules([1 1 1 6], cube)
    _, pal, colors, bonding = ext._resolve_colors(cube, 1, rules, nothing)
    @test count(!, bonding) == 4
    for i in 1:nsites(cube)
        @test (colors[i] == ext.INERT_COLOR) == !bonding[i]
    end
    # A polyhedron keeps the species hue throughout and separates by tint instead.
    _, pal, colors, _ = ext._resolve_colors(
        cube, 1, rules, nothing; bond_tint=ext.FACE_TINT, inert_color=nothing
    )
    for i in 1:nsites(cube)
        want = ext._tint(pal[i], bonding[i] ? ext.FACE_TINT : ext.BODY_TINT)
        @test colors[i] == want
    end
    # Tinting is skipped entirely when a `sitecolor` callback is given.
    _, _, colors, _ = ext._resolve_colors(cube, 1, rules, (_, i) -> ext.INERT_COLOR)
    @test all(==(ext.INERT_COLOR), colors)
end
