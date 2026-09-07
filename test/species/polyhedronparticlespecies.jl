using Roly
using Roly: PolyhedronParticleSpecies, UnitTetrahedron, UnitCube, UnitOctahedron,
            UnitDodecahedron, UnitIcosahedron, UnitPyramid, UnitPrism, UnitAntiprism,
            Polyhedron, Tetrahedron, Cube, Octahedron, Dodecahedron, Icosahedron,
            Pyramid, Prism, Antiprism, polyhedron, corners, nfaces, facedegree, facecentroid,
            faceorbits, rotationgroup, inradius, edgemidpoint,
            nsites, dimension, isconvex, numtype, bindingsites, graphrep, setcolors!, color,
            could_contact, overlap, symmetrynumber, nparticles, raise!, lower!,
            collect_attachments, tocanon, toorig, BindingRules, Polyform, nbonds,
            permutationgroup, PatchySphere

using Roly: PolyhedronParticleSpecies, Prism, Tetrahedron, BindingRules, Polyform,
            raise!, collect_attachments, symmetrynumber, graphrep, nfaces
using Rotations: RotMatrix3, rotation_angle
using LinearAlgebra: inv
using Graphs, NautyGraphs, LinearAlgebra, StaticArrays, Rotations, Random

@testset "PolyhedronParticleSpecies" begin
    nedges(p) = length(collect(Roly.exterior_edges(p)))
    solids = [
        ("UnitTetrahedron", UnitTetrahedron, Tetrahedron(), 4, 12),
        ("UnitCube", UnitCube, Cube(), 6, 24),
        ("UnitOctahedron", UnitOctahedron, Octahedron(), 8, 24),
        ("UnitDodecahedron", UnitDodecahedron, Dodecahedron(), 12, 60),
        ("UnitIcosahedron", UnitIcosahedron, Icosahedron(), 20, 60),
        ("UnitPyramid(5)", UnitPyramid(5), Pyramid(5), 6, 5),
        ("UnitPrism(5)", UnitPrism(5), Prism(5), 7, 10),
        ("UnitAntiprism(4)", UnitAntiprism(4), Antiprism(4), 10, 8),
    ]

    for (name, ps, shp, nf, order) in solids
        @test dimension(ps) == 3
        @test numtype(ps) === Float64
        @test isconvex(ps)
        @test nsites(ps) == nf
        @test polyhedron(ps) === shp || nfaces(polyhedron(ps)) == nf

        # Default labels are all distinct, so the sparse encoding is used.
        @test nv(graphrep(ps)) == nf
        @test symmetrynumber(ps) == 1

        # One site per face, sitting at the face centroid.
        for i in 1:nf
            @test isapprox(bindingsite(ps, i).pose.x, facecentroid(shp, i); atol=1e-10)
            @test length(bindingsite(ps, i).vertices) == 1
        end
        # The closest site is at the inradius, and everything is inside the bound.
        @test isapprox(minimum(norm(bindingsite(ps, i).pose.x) for i in 1:nf), inradius(shp))
        @test all(norm(c) <= Roly.bounding_radius(ps) + 1e-10 for c in corners(shp))

        @test occursin("PolyhedronParticleSpecies", sprint(show, ps))
        @test occursin("$nf sites", sprint(show, ps))
    end

    # The package convention: outward normal on local x. Additionally fix local z
    # at the midpoint of the face's first edge. Asked of the constructor rather than of the
    # `Unit*` constants, which turn some of their faces off it -- see `matingtwists`
    for (name, _, shp, nf, _) in solids
        plain = PolyhedronParticleSpecies(shp)
        for i in 1:nf
            psi = bindingsite(plain, i).pose.psi
            @test isapprox(psi[:, 1], Roly.facenormal(shp, i); atol=1e-10)
            v = edgemidpoint(shp, i, 1) - facecentroid(shp, i)
            @test isapprox(psi[:, 3], normalize(v); atol=1e-10)
            @test isapprox(det(psi), 1; atol=1e-10)      # proper rotation
        end
    end

    # The twist references have to agree up to `stab` under every label-preserving rotation:
    # that is what `_propagate_faces` establishes, and what makes bonds between symmetry-related
    # faces equivalent. `_canonical_faces` on its own only gets them to agree up to `sitesym`,
    # which is strictly weaker wherever a face is more symmetric than the body around it.
    # Tetrahedron and octahedron are here for the cases with no translation-mated faces at all:
    # a tetrahedron has no antiparallel pair, and an octahedron's are related by inversion,
    # which is not a proper rotation. Neither has an orientation-free convention to find, and
    # neither needs one; propagation asks only for consistency along the rotations there are.
    for (name, shp) in [("Cube", Cube()), ("Prism(4,h=2)", Prism(4, 1.0; h=2.0)),
                        ("Prism(6)", Prism(6)), ("Prism(3)", Prism(3)),
                        ("Antiprism(4)", Antiprism(4)), ("Dodecahedron", Dodecahedron()),
                        ("Tetrahedron", Tetrahedron()), ("Octahedron", Octahedron())]
        ps = PolyhedronParticleSpecies(shp; colors=faceorbits(shp))
        sitelabel(i) = Roly.sitelabel(ps, i)
        centroids = [facecentroid(shp, i) - sum(corners(shp)) / length(corners(shp))
                     for i in 1:nfaces(shp)]
        faceat(x) = findfirst(i -> isapprox(centroids[i], x; atol=1e-8), 1:nfaces(shp))

        group = filter(Roly.rotationgroup(shp)) do Q
            all(i -> (j = faceat(Q * centroids[i]); !isnothing(j) && sitelabel(j) == sitelabel(i)),
                1:nfaces(shp))
        end
        @test length(group) == symmetrynumber(ps)

        for Q in group, i in 1:nfaces(shp)
            j = faceat(Q * centroids[i])
            bi, bj = Roly.bindingsite(ps, i), Roly.bindingsite(ps, j)
            @test bj.stab == bi.stab
            # Q carries site i's frame onto site j's, up to a turn about j's normal lying in
            # j's stabilizer -- not merely in its sitesym.
            @test any(0:(bj.stab - 1)) do m
                isapprox(Q * bi.pose.psi, bj.pose.psi * RotX(2π * m / bj.stab); atol=1e-8)
            end
        end
    end

    # Bonding two faces that a symmetry relates gives congruent dimers, so a
    # space-filling solid assembles into one lattice rather than several.
    for (name, shp) in [("Cube", Cube()), ("Prism(6)", Prism(6)), ("Prism(3)", Prism(3))]
        sides = [i for i in 1:nfaces(shp) if abs(Roly.facenormal(shp, i)[3]) < 1e-8]
        ps = PolyhedronParticleSpecies(shp; colors=[i in sides ? 1 : 2 for i in 1:nfaces(shp)])
        rules = BindingRules([1 first(sides) 1 first(sides)], ps)
        step = 2 * norm(facecentroid(shp, first(sides)))
        poly = Polyform(rules, 1)
        grown = 0
        for (site, loc, r) in collect_attachments(poly)
            trial = copy(poly)
            ismissing(raise!(trial, site, loc, r)) && continue
            a, b = trial.particles
            @test isapprox(norm(a.pose.x - b.pose.x), step; atol=1e-8)
            # Every neighbour lands in the plane the side faces span.
            @test isapprox((b.pose.x - a.pose.x)[3], 0; atol=1e-8)
            grown += 1
        end
        # One candidate per sticky face: the mate sites collapse to one representative, while
        # the host sites do not, so canonization is what removes the duplicates.
        @test grown == length(sides)
    end

    for (name, _, shp, nf, order) in solids
        # Faces grouped by geometric orbit recover the solid's rotation group.
        geo = PolyhedronParticleSpecies(shp; colors=faceorbits(shp))
        @test symmetrynumber(geo) == order == length(rotationgroup(shp))
        # Distinct labels give 1, regardless of encoding.
        @test symmetrynumber(PolyhedronParticleSpecies(shp)) == 1
        @test symmetrynumber(dartspecies(shp)) == 1
    end

    # encoding

    # The graph symmetry number must equal the rotational symmetry of its site arrangement.
    # The constructors enforce this, so every
    # species built above already satisfies it; check it explicitly across the library.
    for (_, ps, shp, _, order) in solids
        # Distinct labels: both encodings describe the arrangement, and both give 1.
        for build in (PolyhedronParticleSpecies, dartspecies, cyclespecies)
            distinct = build(shp)
            @test symmetrynumber(distinct) == length(permutationgroup(distinct)) == 1
        end
        # Repeated labels: only the dart encoding carries the rotation group, and forcing the
        # sparse one is rejected rather than silently reporting the cyclic order.
        geo = dartspecies(shp; colors=faceorbits(shp))
        @test symmetrynumber(geo) == length(permutationgroup(geo)) == order
        order == nfaces(shp) ||
            @test_throws ArgumentError cyclespecies(shp; colors=faceorbits(shp))
    end
    for (shp, order) in [(Tetrahedron(), 12), (Cube(), 24), (Dodecahedron(), 60), (Prism(5), 10)]
        sphere = PatchySphere(shp, 1.0; colors=faceorbits(shp))
        @test symmetrynumber(sphere) == length(permutationgroup(sphere)) == order
    end

    # Site stabilizers
    # How much of a site's own symmetry the whole particle keeps. At most its sitesym, and the
    # ratio is how many distinct ways a partner can attach there: turns in the stabilizer put
    # the same body in the same place with only its sites permuted.
    sitesyms(ps) = [bindingsite(ps, i).sitesym for i in 1:nsites(ps)]
    ntwists(ps) = sitesyms(ps) .÷ Roly.stabilizerorders(ps)

    # A cube keeps all four turns about a face normal, so a face-to-face bond has one
    # twist and nothing changes for polycubes.
    cube = PolyhedronParticleSpecies(Cube(); colors=fill(1, 6))
    @test sitesyms(cube) == fill(4, 6)
    @test Roly.stabilizerorders(cube) == fill(4, 6)
    @test ntwists(cube) == fill(1, 6)

    # Distinguishing the caps costs the side faces two of those turns, since a quarter turn
    # about a side normal carries the other sides onto caps.
    caps = [abs(n[3]) > 0.5 ? 2 : 1 for n in Roly.facenormals(Cube())]
    capped = PolyhedronParticleSpecies(Cube(); colors=caps)
    @test Roly.stabilizerorders(capped) == [c == 2 ? 4 : 2 for c in caps]
    @test ntwists(capped) == [c == 2 ? 1 : 2 for c in caps]

    # A triangular prism's side faces are squares
    # when h == a (sitesym 4) but the prism is only 2-fold about them, so a partner can attach
    # two ways: in the plane, or tipped out of it.
    tri = PolyhedronParticleSpecies(Prism(3); colors=faceorbits(Prism(3)))
    @test sitesyms(tri) == [3, 4, 4, 4, 3]
    @test Roly.stabilizerorders(tri) == [3, 2, 2, 2, 3]
    @test ntwists(tri) == [1, 2, 2, 2, 1]

    # Make the prism taller, so that faces become rectangles: sitesym and stabilizer agree at 2,
    # leaving a single twist.
    tall = Prism(3, 1.0; h=2.0)
    tallps = PolyhedronParticleSpecies(tall; colors=faceorbits(tall))
    @test sitesyms(tallps) == [3, 2, 2, 2, 3]
    @test ntwists(tallps) == fill(1, 5)

    # A stabilizer always divides the sitesym, and always divides the symmetry number.
    for ps in (cube, capped, tri, tallps, UnitDodecahedron, UnitAntiprism(4))
        stabs = Roly.stabilizerorders(ps)
        @test all(sitesyms(ps) .% stabs .== 0)
        @test all(symmetrynumber(ps) .% stabs .== 0)
    end

    # Deriving the labeling from the coloring cannot claim a symmetry that is not there: one
    # color on all four faces of a triangular pyramid does not make its base equivalent to its
    # sides, and the derivation splits them, leaving the 3-fold axis rather than a tetrahedral
    # 12.
    pyramid = PolyhedronParticleSpecies(Pyramid(3); colors=fill(1, 4))
    @test length(unique(Roly.labels(graphrep(pyramid)))) == 2
    @test symmetrynumber(pyramid) == length(permutationgroup(pyramid)) == 3
    # The sparse encoding imposes a cyclic order, which is not the symmetry of a tetrahedral
    # or octahedral patch arrangement: it claims n where the truth is |G|. This is the failure
    # that using `cycleencoding`/`dartencoding` does not rule out on its own.
    @test_throws ArgumentError cyclesphere(Tetrahedron(), 1.0; colors=fill(1, 4))
    @test_throws ArgumentError cyclesphere(Cube(), 1.0; colors=fill(1, 6))
    # The dart encoding is the one that describes them
    @test symmetrynumber(dartsphere(Tetrahedron(), 1.0; colors=fill(1, 4))) == 12

    # A cube with the two caps distinguished from the four sides keeps the 4-fold axis
    # and the 2-fold axes through it: D_4, of order 8.
    caps = [abs(n[3]) > 0.5 ? 2 : 1 for n in Roly.facenormals(Cube())]
    @test symmetrynumber(PolyhedronParticleSpecies(Cube(); colors=caps)) == 8
    # Singling out one face leaves only the rotations about its normal.
    @test symmetrynumber(PolyhedronParticleSpecies(Cube(); colors=[2, 1, 1, 1, 1, 1])) == 4

    shp = Cube()
    # The constructor takes the sparse encoding when it is provably equivalent, and the dart
    # encoding when a repeated label means the graph has to carry the rotation group.
    @test nv(graphrep(PolyhedronParticleSpecies(shp))) == 6
    @test nv(graphrep(PolyhedronParticleSpecies(shp; colors=fill(1, 6)))) == 24
    @test nv(graphrep(dartspecies(shp))) == 24
    @test nv(graphrep(cyclespecies(shp))) == 6
    @test all(length(bindingsite(dartspecies(shp), i).vertices) == 4 for i in 1:6)

    @test_throws ArgumentError PolyhedronParticleSpecies(shp; colors=1:5)

    ps = PolyhedronParticleSpecies(Cube(); colors=[3, 1, 4, 1, 5, 9])
    @test color(bindingsite(ps, 1)) == 3
    @test color(bindingsite(ps, 5)) == 5

    cp = copy(ps)
    # A recoloring that groups the sites the same way needs no new graph, so it applies in place.
    setcolors!(cp, [30, 10, 40, 10, 50, 90])
    @test color(bindingsite(cp, 1)) == 30
    @test color(bindingsite(ps, 1)) == 3          # the copy is independent
    @test_throws ArgumentError setcolors!(cp, [1, 2])

    # A recoloring that groups sites differently does need a new graph. This cube's colors are nearly all
    # distinct, so it was built as a bare 6-cycle; making every face alike would ask that cycle
    # to report the full 24, which it cannot, so it errors.
    @test symmetrynumber(cp) == symmetrynumber(ps) == 1
    @test_throws ArgumentError setcolors!(cp, fill(10, 6))
    @test symmetrynumber(cp) == 1                  # and the failed call leaves it untouched
    @test [color(bindingsite(cp, i)) for i in 1:6] == [30, 10, 40, 10, 50, 90]
    # Built with those colors from the start, it is the dart encoding and the symmetry is there.
    @test symmetrynumber(PolyhedronParticleSpecies(Cube(); colors=fill(10, 6))) == 24

    # A coloring is the whole statement, so recoloring has to carry the labeling and the
    # stabilizers with it. Given a graph roomy enough to hold the result -- here the dart
    # encoding, asked for explicitly -- a prism whose faces all start out distinct can be
    # recolored into its full D_3 and the derived quantities follow.
    prism = dartspecies(Prism(3, 1.0; h=2.0); colors=1:5)
    @test symmetrynumber(prism) == length(permutationgroup(prism)) == 1
    @test Roly.stabilizerorders(prism) == fill(1, 5)

    caps = [i for i in 1:5 if abs(Roly.facenormal(Prism(3, 1.0; h=2.0), i)[3]) > 1e-8]
    setcolors!(prism, [i in caps ? 7 : 8 for i in 1:5])
    @test symmetrynumber(prism) == length(permutationgroup(prism)) == 6
    # Caps are 3-fold about their normals and the prism is 3-fold about them; the rectangular
    # sides are 2-fold and so is the prism about those.
    @test [Roly.bindingsite(prism, i).stab for i in 1:5] == [i in caps ? 3 : 2 for i in 1:5]
    @test length(unique(Roly.labels(graphrep(prism)))) == 2

    id = Pose{3,Float64,RotMatrix3{Float64}}(SVector(0.0, 0.0, 0.0), one(RotMatrix3{Float64}))
    shifted(d, R=one(RotMatrix3{Float64})) =
        Pose{3,Float64,RotMatrix3{Float64}}(SVector(d, 0.0, 0.0), R)

    @test overlap(UnitCube => id, UnitCube => id)
    @test !overlap(UnitCube => id, UnitCube => shifted(5.0))
    @test !could_contact(UnitCube => id, UnitCube => shifted(5.0))
    @test could_contact(UnitCube => id, UnitCube => shifted(1.0))

    # Face to face at unit separation: touching, not overlapping.
    @test !overlap(UnitCube => id, UnitCube => shifted(1.0))
    @test overlap(UnitCube => id, UnitCube => shifted(0.9))

    # Rotated by 45 degrees the bounding spheres still reach and the inscribed spheres do
    # not, so these go through the separating-axis test rather than either fast path.
    rot45 = RotMatrix3{Float64}(RotZ(π / 4))
    @test overlap(UnitCube => id, UnitCube => shifted(1.05, rot45))
    @test !overlap(UnitCube => id, UnitCube => shifted(1.3, rot45))

    # An edge-on-edge configuration, which no face normal separates: this is what the
    # edge cross-product axes are for
    edge = RotMatrix3{Float64}(RotZ(π / 4) * RotX(π / 4))
    @test !overlap(UnitCube => id, UnitCube => shifted(1.5, edge))
    @test overlap(UnitCube => id, UnitCube => shifted(0.95, edge))

    for (name, ps, _, _, _) in solids
        @test overlap(ps => id, ps => id)
        @test !overlap(ps => id, ps => shifted(10.0))
    end

    # check SAT overlap against in and out radii
    function fullsat(s1, pose1, s2, pose2)
        axes = Iterators.flatten((
            (pose1.psi * n for n in s1.normals), (pose2.psi * n for n in s2.normals),
            (cross(pose1.psi * e1, pose2.psi * e2)
             for e1 in s1.edgedirections, e2 in s2.edgedirections)))
        return Roly.sat_overlap(axes, corners(s1), pose1, corners(s2), pose2, s1.skin + s2.skin)
    end
    P(x, R=one(RotMatrix3{Float64})) = Pose{3,Float64,RotMatrix3{Float64}}(SVector{3}(x), R)

    Random.seed!(20260812)
    for shp in (Cube(), Prism(6), Prism(3), Tetrahedron(), Dodecahedron(), Octahedron(),
                Antiprism(4), Polyhedron([SVector(x, y, z) for x in (-1.0, 1.0)
                                          for y in (-2.0, 2.0) for z in (-3.0, 3.0)]))
        ps = PolyhedronParticleSpecies(shp)
        rmin, rmax = Roly.inradius(shp), Roly.bounding_radius(shp)
        disagreements = 0
        overlapping = 0
        for _ in 1:4000
            u = normalize(SVector{3}(randn(3)))
            t = u * (rand() < 0.6 ? 2rmin + 2(rmax - rmin) * rand() : 4rmax * rand())
            # Mostly aligned, which is the case the fast path exists for, but not only.
            pose2 = P(t, rand() < 0.7 ? one(RotMatrix3{Float64}) :
                         RotMatrix3(rand(RotMatrix{3,Float64})))
            got = overlap(ps => P(zeros(3)), ps => pose2)
            got == fullsat(ps, P(zeros(3)), ps, pose2) || (disagreements += 1)
            got && (overlapping += 1)
        end
        @test disagreements == 0
        # The sampling has to straddle the boundary, or the agreement above is vacuous.
        @test 0.1 < overlapping / 4000 < 0.9
    end

    # Number of proper rigid motions mapping a polyform onto itself, computed from the
    # particle poses alone.
    function geometric_symmetry(poly)
        xs = [p.pose.x for p in poly.particles]
        Rs = [p.pose.psi for p in poly.particles]
        n = length(xs)
        return count(1:n) do j
            Q = Rs[j] * inv(Rs[1])
            t = xs[j] - Q * xs[1]
            all(1:n) do i
                xi, Ri = Q * xs[i] + t, Q * Rs[i]
                any(k -> isapprox(xs[k], xi; atol=1e-8) && isapprox(Rs[k], Ri; atol=1e-8), 1:n)
            end
        end
    end

    squarefaces(p) = [i for i in 1:nfaces(p) if facedegree(p, i) == 4]
    systems = [
        ("tetra", Tetrahedron(), [1 1 1 2; 1 3 1 4], 5),
        ("cube/opposite", Cube(), [1 1 1 6; 1 2 1 5; 1 3 1 4], 5),
        ("cube/full", Cube(), reduce(vcat, [[1 i 1 j] for i in 1:6 for j in i:6]), 3),
        ("prism5/sides", Prism(5),
         reduce(vcat, [[1 i 1 j] for i in squarefaces(Prism(5)) for j in squarefaces(Prism(5)) if j >= i]), 4),
        ("dodeca", Dodecahedron(), [1 1 1 12; 1 2 1 11; 1 3 1 10], 4),
    ]

    for (name, shp, bonds, maxn) in systems
        results = map((cyclespecies, dartspecies)) do build
            rules = BindingRules(bonds, build(shp))
            polys = polygen(rules; maxsize=maxn)
            (counts=[count(p -> nparticles(p) == k, polys) for k in 1:maxn],
             sigmas=sort([(nparticles(p), symmetrynumber(p)) for p in polys]),
             polys=polys)
        end
        cyc, dart = results

        # With all face labels distinct the sparse encoding is provably equivalent to
        # the dart encoding: same structures, same symmetry numbers.
        @test cyc.counts == dart.counts
        @test cyc.sigmas == dart.sigmas
        @test first(cyc.counts) == 1

        # The graph symmetry number must equal the assembly's actual rotational
        # symmetry. A pairing that lost the handedness of a face-to-face bond would
        # show up here as a graph automorphism the geometry does not have.
        for polys in (cyc.polys, dart.polys)
            for p in polys
                @test symmetrynumber(p) == geometric_symmetry(p)
            end
        end
    end

    rules = BindingRules([1 1 1 6; 1 2 1 5; 1 3 1 4], UnitCube)
    polys = polygen(rules; maxsize=4)

    for p in polys
        nparticles(p) < 2 && continue
        # Vertex bookkeeping stays consistent.
        for v in eachindex(p.canon2orig)
            @test tocanon(p, toorig(p, v)) == v
        end
        # Particle vertex blocks cover the graph exactly; a stale leading vertex after a
        # removal would leave a block hanging off the end.
        blocks = sort(reduce(vcat, [collect(Roly.graphvertices(pt, rules)) for pt in p.particles]))
        @test blocks == 1:nv(graphrep(p))
        # A connected polyform has at least n-1 bonds, and this cube is cycle-encoded, so each
        # bond reaches the graph as exactly one edge.
        @test nbonds(p) >= nparticles(p) - 1
        @test nedges(p) == nbonds(p)
    end

    # Raising then lowering returns the original structure.
    big = polys[end]
    rebuilt = copy(big)
    lower!(rebuilt)
    @test nparticles(rebuilt) == nparticles(big) - 1
    for (site, loc) in collect_attachments(rebuilt)
        trial = copy(rebuilt)
        ismissing(raise!(trial, site, loc)) && continue
        trial == big && (@test trial.canon2orig == big.canon2orig; break)
    end

    # Sites of different sizes may bond: a prism's 4-vertex square faces to its 5-vertex
    # pentagonal ones. `contact_pairing` joins them at gcd(4, 5) = 1 vertex, so the bond
    # carries one edge instead of four.
    prism = dartspecies(Prism(5))
    squares = [i for i in 1:7 if facedegree(Prism(5), i) == 4]
    pentagons = [i for i in 1:7 if facedegree(Prism(5), i) == 5]

    mixed = BindingRules([1 squares[1] 1 pentagons[1]], prism)
    @test mixed isa BindingRules
    same = BindingRules([1 squares[1] 1 squares[2]], prism)
    @test same isa BindingRules
    # Four dart pairs for a square-to-square bond, one reference pair for the mixed one, and
    # one bond either way.
    @test nedges(polygen(same; maxsize=2)[end]) == 4
    @test nedges(polygen(mixed; maxsize=2)[end]) == 1
    @test nbonds(polygen(same; maxsize=2)[end]) == nbonds(polygen(mixed; maxsize=2)[end]) == 1

    # twists

    shp = Prism(3)
    sides = [i for i in 1:nfaces(shp) if abs(Roly.facenormal(shp, i)[3]) < 1e-8]
    colors = [i in sides ? 1 : 2 for i in 1:nfaces(shp)]
    counts(rules) = [polyenum(rules; maxsize=i)[1] for i in 1:5]

    # Turning a whole orbit by the same amount turns the partner by twice that, since the offset
    # lands on both sides of the face-to-face flip. Here that is 2*90 = 180 degrees, which is a
    # symmetry of the prism about a side face (stab = 2), so the lattice is untouched however
    # far the sides are turned
    for t in (0, π/2, π, 3π/2)
        ps = PolyhedronParticleSpecies(shp; colors,
                                       twists=[i in sides ? t : 0.0 for i in 1:nfaces(shp)])
        @test symmetrynumber(ps) == 6
        @test counts(BindingRules([1 first(sides) 1 first(sides)], ps)) == [1, 2, 3, 6, 10]
    end

    # A side of this prism is square, so a quarter turn carries its frame onto an equivalent one
    # and marks nothing, however alone it stands. A whole 2pi turn least of all: it is the
    # identity, and reading it as a mark would cost the body symmetry it plainly has.
    for t in (π/2, π, 3π/2, 2π)
        whole = PolyhedronParticleSpecies(shp; colors,
                                          twists=[i == first(sides) ? t : 0.0 for i in 1:nfaces(shp)])
        @test symmetrynumber(whole) == 6
        @test length(unique(Roly.labels(graphrep(whole)))) == 2
    end

    # What the face's own symmetry cannot absorb does mark it, and then turning one face of an
    # orbit differently is a deliberate break that the labeling has to record, or the graph would
    # keep claiming a symmetry the frames no longer have. The twist is folded into the key
    # `siteorbits` groups by, so the orbit splits.
    split = PolyhedronParticleSpecies(shp; colors,
                                      twists=[i == first(sides) ? π/4 : 0.0 for i in 1:nfaces(shp)])
    @test symmetrynumber(split) == length(permutationgroup(split)) == 2
    @test length(unique(Roly.labels(graphrep(split)))) == 3
    @test symmetrynumber(PolyhedronParticleSpecies(shp; colors)) == 6

    # On a particle with no symmetry to hide behind, the twist angle is the bond twist outright.
    # Turning only the mate's face turns the partner by exactly one dart step, not two.
    function dimer(t)
        ps = dartspecies(Tetrahedron(); twists=[0.0, t, 0.0, 0.0])
        @test symmetrynumber(ps) == 1
        poly = Polyform(BindingRules([1 1 1 2], ps), 1)
        for (site, loc, r) in collect_attachments(poly)
            trial = copy(poly)
            ismissing(raise!(trial, site, loc, r)) && continue
            return trial
        end
        return nothing
    end
    a, b = dimer(0.0), dimer(2π/3)
    @test !isnothing(a) && !isnothing(b)
    @test isapprox(rotation_angle(RotMatrix3(a.particles[2].pose.psi * inv(b.particles[2].pose.psi))),
                   2π / 3; atol=1e-8)
    # A triangular face has three darts, so three steps is the identity.
    @test isapprox(dimer(2π).particles[2].pose.psi, a.particles[2].pose.psi; atol=1e-8)

    # the encoding picks it up
    @test graphrep(a) != graphrep(b)

    # Any angle is allowed, not just whole dart steps. A tetrahedron face has three darts, so 2π/3 lands
    # back on a dart and 0.4 does not; both are allowed and both turn the partner by exactly
    # the angle asked for.
    for θ in (0.4, 1.0, 2π/3 + 0.1)
        a, b = dimer(0.0), dimer(θ)
        @test isapprox(rotation_angle(RotMatrix3(a.particles[2].pose.psi *
                                                 inv(b.particles[2].pose.psi))), θ; atol=1e-8)
    end
    # A twist shared across an orbit leaves the symmetry alone, since turns about a site's own
    # normal commute with its stabilizer; turning one face differently splits the orbit.
    uniform = PolyhedronParticleSpecies(shp; colors, twists=[i in sides ? 0.37 : 0.0 for i in 1:nfaces(shp)])
    @test symmetrynumber(uniform) == length(permutationgroup(uniform)) == 6
    partial = PolyhedronParticleSpecies(shp; colors, twists=[i == first(sides) ? 0.37 : 0.0 for i in 1:nfaces(shp)])
    @test symmetrynumber(partial) == length(permutationgroup(partial)) == 2

    @test_throws ArgumentError PolyhedronParticleSpecies(Cube(); twists=[0.0, 1.0])
    @test_throws ArgumentError PolyhedronParticleSpecies(Cube(); locking=[true, false])

    # bodies that are not regular: opposing faces need not be congruent, and a face need not be
    # symmetric at all. Nothing about squaring opposing faces up assumes otherwise -- a face whose
    # only symmetry is the identity has no whole step but zero, so it is never turned
    box = Roly.Polyhedron(SVector{3,Float64}[(x, 2y, 3z) for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)])
    frustum = Roly.Polyhedron(SVector{3,Float64}[(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
                                                 (-0.4, -0.4, 1), (0.4, -0.4, 1), (0.4, 0.4, 1),
                                                 (-0.4, 0.4, 1)])
    # a box of three different side lengths: every face a rectangle, twofold rather than fourfold
    @test [Roly.facesym(box, i) for i in 1:nfaces(box)] == fill(2, 6)
    @test symmetrynumber(PolyhedronParticleSpecies(box; colors=fill(1, 6))) == 4
    # a frustum: its two square caps oppose each other but are not the same size, and its sides
    # are trapezoids with no symmetry at all
    @test sort(unique(Roly.facesym(frustum, i) for i in 1:nfaces(frustum))) == [1, 4]
    @test symmetrynumber(PolyhedronParticleSpecies(frustum; colors=fill(1, 6))) == 4
    @test PolyhedronParticleSpecies(frustum) isa PolyhedronParticleSpecies

    # the symmetric counterparts differ from the defaults in coloring alone: one color throughout,
    # so the body keeps its rotation group where a color per face leaves it none
    for (plain, sym, sigma) in ((UnitTetrahedron, SymmetricUnitTetrahedron, 12),
                                (UnitCube, SymmetricUnitCube, 24),
                                (UnitOctahedron, SymmetricUnitOctahedron, 24),
                                (UnitDodecahedron, SymmetricUnitDodecahedron, 60),
                                (UnitIcosahedron, SymmetricUnitIcosahedron, 60))
        n = nsites(plain)
        @test nsites(sym) == n
        @test [color(bindingsite(plain, i)) for i in 1:n] == 1:n
        @test [color(bindingsite(sym, i)) for i in 1:n] == fill(1, n)
        # the sites sit in the same places; only their twist references and colors differ
        @test [bindingsite(plain, i).pose.x for i in 1:n] == [bindingsite(sym, i).pose.x for i in 1:n]
        @test symmetrynumber(plain) == 1
        @test symmetrynumber(sym) == sigma
    end

    # `UnitCube` carries the twists that make opposite faces mate square-on, so a rule set bonding
    # them tiles space. `SymmetricUnitCube` gets there the other way: one color leaves every face
    # a fourfold stabilizer, so a bond admits all four twists and the translation is among them
    @test isunitcell(first(polygen(BindingRules([1 1 1 6; 1 2 1 5; 1 3 1 4], UnitCube); maxsize=1)))
    @test isunitcell(first(polygen(BindingRules([1 1 1 1], SymmetricUnitCube); maxsize=1)))
end
