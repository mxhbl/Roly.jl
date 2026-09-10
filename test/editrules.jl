const RE = Base.get_extension(Roly, :TachikomaExt)

# Attach a particle to the free site `(placement, site)` using site `incoming` of the active
# species, the way pressing enter in the editor does.
function _attach!(m, placement, site, incoming)
    i = findfirst(==((placement, site)), m.free)
    isnothing(i) && error("site ($placement, $site) is not free")
    m.focus = :construction  # enter means attach only in the construction pane
    m.anchor = i
    m.incoming = incoming
    RE.update!(m, Tachikoma.KeyEvent(:enter))
    return m
end

# Run the enumeration, which the editor does only when asked.
_enumerate!(m) = RE.update!(m, Tachikoma.KeyEvent('e'))

function _frame(tb)
    area = Tachikoma.Rect(1, 1, tb.width, tb.height)
    return Tachikoma.Frame(tb.buf, area, Tachikoma.GraphicsRegion[], [])
end

function _render(m, w=170, h=32)
    tb = Tachikoma.TestBackend(w, h)
    RE.view(m, _frame(tb))
    return tb
end

# Bearing of each free site about the structure's centroid, which is the order the anchor
# cycles in.
function _bearings(m)
    sites = RE.absolutesites(RE.placements(m), m.species)
    c = sum(pose.x for (_, pose) in RE.placements(m)) / length(RE.placements(m))
    return [(d=sites[i][k].pose.x - c; atan(d[2], d[1])) for (i, k) in m.free]
end

# A structure of `n` species, each a copy of `spcs` attached to a different site of the first.
function _manyspecies(n, spcs=UnitHexagon)
    m = RE.EditorModel(spcs)
    m.focus = :construction  # digits pick the active species, which only the build pane has
    for d in 2:n
        RE.update!(m, Tachikoma.KeyEvent(Char('0' + d)))
        _attach!(m, 1, d, 1)
    end
    return m
end

# `n` species, none of them placed, which is what the rules pane alone can produce.
function _addspecies(n, spcs=UnitHexagon)
    m = RE.EditorModel(spcs)
    for _ in 2:n
        RE.update!(m, Tachikoma.KeyEvent('a'))
    end
    return m
end

# The distinct face colors covering a region of a buffer, and how many cells each covers. A 3D
# particle is drawn as cell backgrounds, which `find_text` cannot see, the faces being spaces.
function _fills(buf, rect)
    seen = Dict{UInt8,Int}()
    for x in rect.x:(rect.x + rect.width - 1), y in rect.y:(rect.y + rect.height - 1)
        bg = buf.content[Tachikoma.buf_index(buf, x, y)].style.bg
        bg isa Tachikoma.Color256 && (seen[bg.code] = get(seen, bg.code, 0) + 1)
    end
    return seen
end

# Attach at whichever free site the cursor is on. In 3D the free list holds only the sites the
# camera shows, so naming one by index does not carry across viewpoints.
function _attachhere!(m, incoming)
    m.focus = :construction
    m.incoming = incoming
    RE.update!(m, Tachikoma.KeyEvent(:enter))
    return m
end

# The characters covering a region, which is what a braille drawing puts there.
function _chars(buf, rect)
    return [
        buf.content[Tachikoma.buf_index(buf, x, y)].char for
        x in rect.x:(rect.x + rect.width - 1), y in rect.y:(rect.y + rect.height - 1)
    ]
end

# Cells holding one of the axis indicator's letters, which it draws bold. The sidebar's own text
# has letters in it too, so the style is what tells them apart.
function _axisletters(tb)
    n = 0
    for x in 1:tb.width, y in 1:tb.height
        cell = tb.buf.content[Tachikoma.buf_index(tb.buf, x, y)]
        n += (cell.char in ('x', 'y', 'z') && cell.style.bold)
    end
    return n
end

# How squarely the site the cursor is on faces the camera. Negative means it is on the far side.
function _anchorfacing(m)
    s = RE.anchorsite(m)
    return dot(s.pose.psi * SVector(1.0, 0.0, 0.0), RE.camera(m).view)
end

@testset "editrules" begin
    @test !isnothing(RE)

    # A fresh model holds one particle, every site of it free, and no rule yet.
    m = RE.EditorModel(UnitSquare)
    @test length(RE.placements(m)) == 1
    @test length(m.free) == nsites(UnitSquare)
    @test nbonds(m.rules) == 0
    @test RE.ntwists(m) == 1  # a 2D bond fixes the partner's orientation outright

    # An attachment bonds the anchor's color to the incoming site's color, and consumes both.
    _attach!(m, 1, 2, 4)
    @test RE.inferred_bonds(RE.placements(m), m.species) == [(2, 4)]
    @test length(m.free) == 2 * nsites(UnitSquare) - 2
    @test RE.bonds_matrix(RE.placements(m), m.species) == [1 2 1 4]
    @test !isnothing(m.rules)

    m = RE.EditorModel(UnitSquare)
    _attach!(m, 1, 1, 3)
    @test RE.inferred_bonds(RE.placements(m), m.species) == [(1, 3)]

    m = RE.EditorModel(UnitTriangle)
    _attach!(m, 1, 3, 3)
    @test RE.inferred_bonds(RE.placements(m), m.species) == [(3, 3)]

    m = RE.EditorModel(UnitHexagon)
    _attach!(m, 1, 3, 6)
    @test RE.inferred_bonds(RE.placements(m), m.species) == [(3, 6)]

    # Rules are read off every pair of particles, not just the pairs that were attached. Here
    # the fourth square closes a 2x2 block: it bonds (1, 2) to the particle it was attached to,
    # and incidentally (1, 4) to the one across the block, which no attachment named.
    m = RE.EditorModel(UnitSquare)
    _attach!(m, 1, 1, 3)
    _attach!(m, 1, 2, 4)
    @test RE.inferred_bonds(RE.placements(m), m.species) == [(1, 3), (2, 4)]
    _attach!(m, 2, 2, 1)
    @test RE.inferred_bonds(RE.placements(m), m.species) == [(1, 2), (1, 3), (1, 4), (2, 4)]
    @test size(RE.bonds_matrix(RE.placements(m), m.species), 1) == 4

    # An empty grid yields nothing to infer.
    m = RE.EditorModel(UnitSquare)
    @test isempty(RE.inferred_bonds(RE.placements(m), m.species))
    @test size(RE.bonds_matrix(RE.placements(m), m.species)) == (0, 4)
    @test isnothing(RE.buildrules(RE.placements(m), m.species))

    # A pentagon chain curls back on itself, and the attachment that would overlap is refused.
    m = RE.EditorModel(UnitNgon(5))
    _attach!(m, 1, 4, 5)
    _attach!(m, 2, 4, 5)
    @test length(RE.placements(m)) == 3
    _attach!(m, 3, 4, 5)
    @test length(RE.placements(m)) == 3
    @test startswith(m.message, "overlaps")

    # `.` walks the perimeter in order and reaches every free site; `,` goes back.
    m = RE.EditorModel(UnitSquare)
    m.focus = :construction
    visited = Int[]
    for _ in eachindex(m.free)
        push!(visited, m.anchor)
        RE.update!(m, Tachikoma.KeyEvent('.'))
    end
    @test sort(visited) == collect(eachindex(m.free))
    @test m.anchor == 1  # a full lap returns to the start
    RE.update!(m, Tachikoma.KeyEvent(','))
    @test m.anchor == length(m.free)

    # The editor opens on the rules, with the construction pane put away.
    m = RE.EditorModel(UnitSquare)
    @test m.focus === :rules
    @test !m.showconstruction

    # Tab cycles the visible panes; the construction pane joins the cycle when `b` shows it.
    RE.update!(m, Tachikoma.KeyEvent(:tab))
    @test m.focus === :enumeration
    RE.update!(m, Tachikoma.KeyEvent(:tab))
    @test m.focus === :rules
    RE.update!(m, Tachikoma.KeyEvent('b'))
    @test m.showconstruction
    @test m.focus === :construction  # showing it also moves there, since that is why you pressed it
    RE.update!(m, Tachikoma.KeyEvent(:tab))
    @test m.focus === :enumeration
    RE.update!(m, Tachikoma.KeyEvent(:tab))
    @test m.focus === :rules
    RE.update!(m, Tachikoma.KeyEvent('b'))
    @test !m.showconstruction

    # Arrow keys step one site along the perimeter, so nothing is skipped, and always to one of
    # the two neighbors, so no key is ever dead. Which neighbor is decided by bearing alone.
    m = RE.EditorModel(UnitSquare)
    _attach!(m, 1, 1, 3)
    _attach!(m, 1, 2, 4)
    _attach!(m, 2, 2, 1)
    _sitepos(n) = RE.absolutesites(RE.placements(m), m.species)[m.free[n][1]][m.free[n][2]].pose.x
    _cos(from, to, d) = (v=_sitepos(to) - _sitepos(from); (v[1] * d[1] + v[2] * d[2]) / hypot(v...))
    n = length(m.free)
    for start in 1:n, dir in (:up, :down, :left, :right)
        m.anchor = start
        RE.update!(m, Tachikoma.KeyEvent(dir))
        nxt, prv = mod1(start + 1, n), mod1(start - 1, n)
        @test m.anchor in (nxt, prv)          # a neighbor, never a jump
        @test m.anchor != start               # never dead
        d = getfield(RE.ARROWS, dir)
        @test _cos(start, m.anchor, d) >= _cos(start, m.anchor == nxt ? prv : nxt, d) - 1e-12
    end

    # `r` turns the pending particle by moving which of its sites meets the anchor.
    m = RE.EditorModel(UnitSquare)
    m.focus = :construction
    RE.update!(m, Tachikoma.KeyEvent('r'))
    @test m.incoming == 2
    RE.update!(m, Tachikoma.KeyEvent('R'))
    @test m.incoming == 1
    RE.update!(m, Tachikoma.KeyEvent('R'))
    @test m.incoming == nsites(UnitSquare)

    # The turn shows in the drawing, because the ghost's sites are named where they sit.
    _cons = (RE.SIDEBAR_W + RE.PANE_W + 1):(RE.SIDEBAR_W + 2 * RE.PANE_W)
    _canvas(tb) = [collect(Tachikoma.row_text(tb, y))[_cons] for y in 1:32]
    m = RE.EditorModel(UnitSquare)
    m.showconstruction = true
    m.focus = :construction
    before = _canvas(_render(m))
    RE.update!(m, Tachikoma.KeyEvent('r'))
    @test _canvas(_render(m)) != before

    # A contact carries both of the colors it joins, written side by side.
    m = RE.EditorModel(UnitSquare)
    m.showconstruction = true
    _attach!(m, 1, 1, 3)
    @test RE.contacts(RE.placements(m), m.species) ==
        [(RE.absolutesites(RE.placements(m), m.species)[1][1].pose.x, 1, 1, 1, 3)]
    @test !isnothing(Tachikoma.find_text(_render(m), "13"))

    # Colors are the ones the rules assign, so each species gets its own range rather than every
    # species being labelled 1, 2, 3 and disagreeing with the interaction matrix.
    m = RE.EditorModel(UnitSquare)
    m.showconstruction = true
    _attach!(m, 1, 1, 3)
    RE.update!(m, Tachikoma.KeyEvent('2'))  # `_attach!` left the focus on the build pane
    _attach!(m, 1, 2, 4)
    gc = RE.sitecolors(m)
    @test [gc[(1, k)] for k in 1:4] == 1:4
    @test [gc[(2, k)] for k in 1:4] == 5:8
    @test ncolors(RE.buildrules(RE.placements(m), m.species)) == 8
    # The bond between species 1 site 2 and species 2 site 4 is (2, 8), and reads that way.
    @test (2, 8) in RE.bonded_colors(RE.buildrules(RE.placements(m), m.species))
    @test !isnothing(Tachikoma.find_text(_render(m), "28"))

    # The pending bond is labelled the same way, anchor colour then the ghost site meeting it.
    m = RE.EditorModel(UnitSquare)
    m.showconstruction = true
    m.focus = :construction
    m.anchor = findfirst(==((1, 2)), m.free)
    m.incoming = 4
    @test !isnothing(Tachikoma.find_text(_render(m), "24"))
    RE.update!(m, Tachikoma.KeyEvent('r'))  # turn the ghost, the pairing changes with it
    @test !isnothing(Tachikoma.find_text(_render(m), "21"))

    # Digits add species instances lazily and make them active.
    m = RE.EditorModel(UnitSquare)
    m.focus = :construction
    @test length(m.species) == 1
    RE.update!(m, Tachikoma.KeyEvent('3'))
    @test m.active_species == 3
    @test length(m.species) == 3
    RE.update!(m, Tachikoma.KeyEvent('1'))
    @test m.active_species == 1
    @test length(m.species) == 3  # does not shrink; the sidebar decides what to show
    RE.update!(m, Tachikoma.KeyEvent('0'))
    @test m.active_species == 1   # digits outside 1-9 are ignored

    # Backspace undoes the last attachment, and never removes the seed particle.
    m = RE.EditorModel(UnitSquare)
    _attach!(m, 1, 1, 3)
    RE.update!(m, Tachikoma.KeyEvent(:backspace))
    @test length(RE.placements(m)) == 1
    RE.update!(m, Tachikoma.KeyEvent(:backspace))
    @test length(RE.placements(m)) == 1

    # `n` opens another construction window, `w` cycles between them, and the rules are read off
    # all of them at once, so an arrangement can be kept while the next is built beside it.
    m = RE.EditorModel(UnitSquare)
    m.focus = :construction
    _attach!(m, 1, 1, 3)
    @test length(m.windows) == 1
    @test interactionmatrix(m.rules)[1, 3]
    RE.update!(m, Tachikoma.KeyEvent('n'))
    @test length(m.windows) == 2
    @test m.window == 2
    @test length(RE.placements(m)) == 1              # the new window starts empty
    @test interactionmatrix(m.rules)[1, 3]           # the first window's bond stands
    _attach!(m, 1, 2, 4)
    @test interactionmatrix(m.rules)[1, 3]           # and both windows contribute
    @test interactionmatrix(m.rules)[2, 4]
    RE.update!(m, Tachikoma.KeyEvent('w'))
    @test m.window == 1
    @test length(RE.placements(m)) == 2              # back to the first arrangement
    RE.update!(m, Tachikoma.KeyEvent('W'))
    @test m.window == 2

    # Clearing a window takes with it the bonds only that window produced, and leaves the rest.
    RE.update!(m, Tachikoma.KeyEvent('c'))
    @test length(RE.placements(m)) == 1
    @test interactionmatrix(m.rules)[1, 3]           # still the first window's
    @test !interactionmatrix(m.rules)[2, 4]          # this window's is gone with it
    # Clearing an already-clear window closes it, so the cycle does not fill up with empties.
    RE.update!(m, Tachikoma.KeyEvent('c'))
    @test length(m.windows) == 1
    @test m.window == 1
    @test interactionmatrix(m.rules)[1, 3]
    # The last window is cleared rather than closed, there always being one to build in.
    RE.update!(m, Tachikoma.KeyEvent('c'))
    RE.update!(m, Tachikoma.KeyEvent('c'))
    @test length(m.windows) == 1
    @test nbonds(m.rules) == 0

    # A bond two windows both produce survives either of them being cleared.
    m = RE.EditorModel(UnitSquare)
    m.focus = :construction
    _attach!(m, 1, 1, 3)
    RE.update!(m, Tachikoma.KeyEvent('n'))
    _attach!(m, 1, 1, 3)
    @test interactionmatrix(m.rules)[1, 3]
    RE.update!(m, Tachikoma.KeyEvent('c'))
    @test interactionmatrix(m.rules)[1, 3]

    # A bond set by hand is not the geometry's to take away.
    m = RE.EditorModel(UnitSquare)
    m.focus = :construction
    _attach!(m, 1, 1, 3)
    m.pair = (1, 3)
    RE.togglebond!(m, m.pair)                        # off
    RE.togglebond!(m, m.pair)                        # and on again, now by hand
    RE.update!(m, Tachikoma.KeyEvent('c'))
    @test interactionmatrix(m.rules)[1, 3]

    # Two species: the returned rules keep both, renumbered from one.
    m = RE.EditorModel(UnitSquare)
    m.focus = :construction
    RE.update!(m, Tachikoma.KeyEvent('2'))
    _attach!(m, 1, 1, 3)
    rules = RE.buildrules(RE.placements(m), m.species)
    @test nspecies(rules) == 2

    # `q` and escape accept; nothing else does.
    m = RE.EditorModel(UnitSquare)
    RE.update!(m, Tachikoma.KeyEvent(:up))
    @test !m.quit
    RE.update!(m, Tachikoma.KeyEvent('q'))
    @test m.quit
    m = RE.EditorModel(UnitSquare)
    RE.update!(m, Tachikoma.KeyEvent(:escape))
    @test m.quit

    # End to end: an L trimer of squares gives rules that enumerate.
    m = RE.EditorModel(UnitSquare)
    _attach!(m, 1, 1, 3)
    _attach!(m, 1, 2, 4)
    rules = RE.buildrules(RE.placements(m), m.species)
    @test rules isa BindingRules
    @test polyenum(rules; maxsize=3, maxstrs=100).nstructures > 0

    # Rendering: every pane is drawn, and the rules pane shows both tables at once.
    m = RE.EditorModel(UnitSquare)
    m.showconstruction = true
    _attach!(m, 1, 1, 3)
    tb = _render(m)
    @test !isnothing(Tachikoma.find_text(tb, "╭─ Enumeration"))
    @test !isnothing(Tachikoma.find_text(tb, "Editor"))
    @test !isnothing(Tachikoma.find_text(tb, "╭─ Construction"))
    @test !isnothing(Tachikoma.find_text(tb, "╭─ Rules"))
    @test !isnothing(Tachikoma.find_text(tb, "species 1"))
    @test !isnothing(Tachikoma.find_text(tb, string(RE.bondglyph())))  # the matrix marks the bond
    @test !isnothing(Tachikoma.find_text(tb, "1 ─ 3"))                 # and the list names it

    # `b` hides the construction pane; the rules stay, since they are what the editor produces.
    RE.update!(m, Tachikoma.KeyEvent('b'))
    @test isnothing(Tachikoma.find_text(_render(m), "╭─ Construction"))
    @test !isnothing(Tachikoma.find_text(_render(m), "╭─ Rules"))
    RE.update!(m, Tachikoma.KeyEvent('b'))
    @test !isnothing(Tachikoma.find_text(_render(m), "╭─ Construction"))

    # A narrow terminal drops the enumeration first, then the construction.
    @test isnothing(Tachikoma.find_text(_render(m, RE.SIDEBAR_W + 2 * RE.PANE_W, 32), "╭─ Enumeration"))
    @test !isnothing(Tachikoma.find_text(_render(m, RE.SIDEBAR_W + 2 * RE.PANE_W, 32), "╭─ Construction"))
    @test isnothing(Tachikoma.find_text(_render(m, RE.SIDEBAR_W + RE.PANE_W, 32), "╭─ Construction"))
    @test !isnothing(Tachikoma.find_text(_render(m, RE.SIDEBAR_W + RE.PANE_W, 32), "╭─ Rules"))

    # The enumeration waits to be asked, so it opens saying so rather than having run.
    m = RE.EditorModel(UnitSquare)
    @test m.stale
    @test !isnothing(Tachikoma.find_text(_render(m), "press e to enumerate"))

    # A bond can be set from the rules pane, without placing anything.
    @test m.focus === :rules
    @test m.pair == (1, 1)
    RE.update!(m, Tachikoma.KeyEvent(:down))
    RE.update!(m, Tachikoma.KeyEvent(:right))
    @test m.pair == (2, 2)
    @test !interactionmatrix(m.rules)[2, 2]
    RE.update!(m, Tachikoma.KeyEvent(:enter))
    @test interactionmatrix(m.rules)[2, 2]
    @test nbonds(m.rules) == 1
    RE.update!(m, Tachikoma.KeyEvent(:enter))  # enter toggles, so it comes off again
    @test !interactionmatrix(m.rules)[2, 2]

    # The cursor wraps, and a hand-set bond survives further attachments.
    m = RE.EditorModel(UnitSquare)
    for _ in 1:nsites(UnitSquare)
        RE.update!(m, Tachikoma.KeyEvent(:right))
    end
    @test m.pair == (1, 1)
    m.pair = (2, 2)
    RE.update!(m, Tachikoma.KeyEvent(:enter))
    m.focus = :construction
    _attach!(m, 1, 1, 3)
    @test interactionmatrix(m.rules)[2, 2]   # kept
    @test interactionmatrix(m.rules)[1, 3]   # and the geometry still contributes

    # Clearing a bond the geometry keeps producing holds it cleared.
    m = RE.EditorModel(UnitSquare)
    _attach!(m, 1, 1, 3)
    @test interactionmatrix(m.rules)[1, 3]
    RE.togglebond!(m, (1, 3))
    @test !interactionmatrix(m.rules)[1, 3]
    _attach!(m, 2, 2, 4)                     # a further contact does not bring it back
    @test !interactionmatrix(m.rules)[1, 3]
    @test interactionmatrix(m.rules)[2, 4]

    # The pair list is a display beside the matrix, not a second cursor: it marks whichever pair
    # the matrix cursor is on, and the arrows always move the matrix.
    m = RE.EditorModel(UnitSquare)
    _attach!(m, 1, 1, 3)
    _attach!(m, 1, 2, 4)
    m.focus = :rules
    @test length(bonded_colors(m.rules)) == 2
    m.pair = (1, 1)
    RE.update!(m, Tachikoma.KeyEvent(:right))
    @test m.pair == (1, 2)  # a column step, not a jump to the next listed bond

    # The enumeration runs only when asked, and changing the rules marks the last run stale.
    m = RE.EditorModel(UnitSquare)
    @test isempty(m.polyforms)
    _attach!(m, 1, 1, 3)
    @test m.stale
    @test isempty(m.polyforms)   # attaching does not pay for an enumeration
    _enumerate!(m)
    @test !m.stale
    @test !isempty(m.polyforms)
    @test all(p -> nparticles(p) > 0, m.polyforms)
    @test all(p -> nparticles(p) <= m.maxsize, m.polyforms)
    @test length(m.polyforms) <= m.maxstrs
    chain = length(m.polyforms)
    @test !isnothing(Tachikoma.find_text(_render(m), "╭─ Enumeration"))

    # A rule set that also bonds 2 to 4 allows strictly more structures.
    RE.togglebond!(m, (2, 4))
    @test m.stale
    _enumerate!(m)
    @test length(m.polyforms) > chain

    # Taking every bond away empties it again.
    RE.togglebond!(m, (2, 4))
    RE.togglebond!(m, (1, 3))
    _enumerate!(m)
    @test isempty(m.polyforms)
    @test m.enuminfo == "no bonds"

    # The bounds are adjustable, and moving either marks the run stale rather than re-running.
    m = RE.EditorModel(UnitSquare)
    _attach!(m, 1, 1, 3)
    _enumerate!(m)
    @test !m.stale
    m.focus = :enumeration
    size0, strs0 = m.maxsize, m.maxstrs
    RE.update!(m, Tachikoma.KeyEvent('S'))
    @test m.maxsize == size0 + 1
    @test m.stale
    RE.update!(m, Tachikoma.KeyEvent('s'))
    @test m.maxsize == size0
    RE.update!(m, Tachikoma.KeyEvent('X'))
    @test m.maxstrs > strs0
    RE.update!(m, Tachikoma.KeyEvent('x'))
    @test m.maxstrs == strs0
    for _ in 1:20
        RE.update!(m, Tachikoma.KeyEvent('s'))
    end
    @test m.maxsize == first(RE.MAXSIZE_RANGE)   # clamped
    for _ in 1:20
        RE.update!(m, Tachikoma.KeyEvent('X'))
    end
    @test m.maxstrs == last(RE.MAXSTRS_RANGE)

    # A smaller cap really does cut the run short.
    m = RE.EditorModel(UnitSquare)
    _attach!(m, 1, 1, 3)
    RE.togglebond!(m, (2, 4))
    m.maxsize = 3
    _enumerate!(m)
    @test all(p -> nparticles(p) <= 3, m.polyforms)

    # The enumeration pane takes a selection, which the detail box draws on its own.
    m = RE.EditorModel(UnitSquare)
    _attach!(m, 1, 1, 3)
    _enumerate!(m)
    m.focus = :enumeration
    @test m.selected == 1
    RE.update!(m, Tachikoma.KeyEvent(:right))
    @test m.selected == 2
    RE.update!(m, Tachikoma.KeyEvent(:left))
    @test m.selected == 1
    RE.update!(m, Tachikoma.KeyEvent(:left))
    @test m.selected == 1  # clamped rather than wrapping, so the ends stay put
    @test !isnothing(Tachikoma.find_text(_render(m), "1 particle"))
    RE.update!(m, Tachikoma.KeyEvent(:right))
    @test !isnothing(Tachikoma.find_text(_render(m), "2 particles"))

    # The detail box reports the selection's composition vector: species counts, then bond counts.
    tb = _render(m)
    @test !isnothing(Tachikoma.find_text(tb, "spc"))
    @test !isnothing(Tachikoma.find_text(tb, "bnd"))

    # It also names every site, writing both colors where two meet, so a bond can be read off
    # against the matrix. The dimer bonds colors 1 and 3, so that pair appears.
    @test !isnothing(Tachikoma.find_text(tb, "13"))

    # A particle carries a stroke away from its first site, so its orientation is visible. It is
    # drawn onto the same canvas as the outline, so turning a particle changes the drawing even
    # where no label moves.
    cv = Tachikoma.Canvas(20, 10)
    v = RE.fitworld([(1, one(Roly.posetype(UnitSquare)))], [UnitSquare], 20, 10)
    pts = RE.outline(UnitSquare, one(Roly.posetype(UnitSquare)))
    @test all(iszero, cv.dots)
    RE.drawheading!(cv, v, one(Roly.posetype(UnitSquare)), UnitSquare, pts)
    @test any(!iszero, cv.dots)

    # Too small to read as a direction, so nothing is drawn.
    small = Tachikoma.Canvas(3, 2)
    vsmall = RE.fitworld([(1, one(Roly.posetype(UnitSquare)))], [UnitSquare], 3, 2)
    RE.drawheading!(small, vsmall, one(Roly.posetype(UnitSquare)), UnitSquare, pts)
    @test all(iszero, small.dots)

    # Clearing the only window empties the rules with it, the geometry being all that produced
    # them; keeping an arrangement means giving the next one a window of its own.
    m = RE.EditorModel(UnitSquare)
    m.showconstruction = true
    _attach!(m, 1, 1, 3)
    @test interactionmatrix(m.rules)[1, 3]
    RE.update!(m, Tachikoma.KeyEvent('c'))
    @test length(RE.placements(m)) == 1
    @test isempty(RE.inferred_bonds(RE.placements(m), m.species))
    @test nbonds(m.rules) == 0

    # A placement can meet several particles at once, and every bond it would make is reported
    # before it is committed, not only the one aimed at.
    m = RE.EditorModel(UnitSquare)
    @test length(RE.pendingcontacts(m, RE.previewpose(m))) == 1  # nothing to close yet
    _attach!(m, 1, 1, 3)
    _attach!(m, 1, 2, 4)
    m.anchor = findfirst(==((2, 2)), m.free)
    m.incoming = 1
    pending = RE.pendingcontacts(m, RE.previewpose(m))
    @test length(pending) == 2  # the site aimed at, and the one across the block
    @test (2, 2) in [(i, k) for (_, i, k, _) in pending]
    @test (3, 1) in [(i, k) for (_, i, k, _) in pending]

    # Committing it really does make both.
    before = length(RE.inferred_bonds(RE.placements(m), m.species))
    RE.update!(m, Tachikoma.KeyEvent(:enter))
    @test length(RE.inferred_bonds(RE.placements(m), m.species)) == before + 2

    # A placement that would overlap is refused, and looks refused before enter is pressed: the
    # ghost is drawn in the error color and the sidebar names what is in the way.
    disk = PatchyDisk([0.0, π / 6])   # patches close enough that two neighbours intersect
    m = RE.EditorModel(disk)
    m.showconstruction = true
    _attach!(m, 1, 1, 2)
    m.anchor = findfirst(==((1, 2)), m.free)
    m.incoming = 2
    @test RE.blockedby(m, RE.previewpose(m)) == 2
    @test !isnothing(Tachikoma.find_text(_render(m), "blocked by 2"))
    RE.update!(m, Tachikoma.KeyEvent(:enter))
    @test length(RE.placements(m)) == 2               # refused
    @test m.message == "overlaps particle 2"
    @test m.messagekind === :warning

    # A refused placement is crossed out rather than given a heading, and drawn in a neutral gray:
    # red is one of the species colors, so a red ghost reads as another species.
    pose = one(Roly.posetype(UnitSquare))
    pts = RE.outline(UnitSquare, pose)
    v = RE.fitworld([(1, pose)], [UnitSquare], 20, 10)
    cv = Tachikoma.Canvas(20, 10)
    RE.drawreject!(cv, v, pts)
    cross = copy(cv.dots)
    @test any(!iszero, cross)
    fill!(cv.dots, 0x00)
    RE.drawheading!(cv, v, pose, UnitSquare, pts)
    @test cross != cv.dots

    # Patches far enough apart seat cleanly and nothing is blocked.
    m = RE.EditorModel(PatchyDisk([0.0, π / 3]))
    m.showconstruction = true
    _attach!(m, 1, 1, 2)
    m.anchor = findfirst(==((1, 2)), m.free)
    m.incoming = 2
    @test isnothing(RE.blockedby(m, RE.previewpose(m)))
    RE.update!(m, Tachikoma.KeyEvent(:enter))
    @test length(RE.placements(m)) == 3

    # Species are added and dropped from the rules pane, which is where they matter when there
    # is no construction to place them in.
    m = RE.EditorModel(UnitSquare)
    @test m.focus === :rules
    RE.update!(m, Tachikoma.KeyEvent('a'))
    RE.update!(m, Tachikoma.KeyEvent('a'))
    @test length(m.species) == 3
    @test m.active_species == 3        # the new one becomes active
    @test ncolors(m.rules) == 12

    # A species carrying a bond is kept, with a note rather than a silent loss.
    m.pair = (1, 9)
    RE.update!(m, Tachikoma.KeyEvent(:enter))
    @test nbonds(m.rules) == 1
    RE.update!(m, Tachikoma.KeyEvent('d'))
    @test length(m.species) == 3
    @test m.message == "species 3 has bonds"
    @test m.messagekind === :warning

    # Clear the bond and it goes.
    RE.update!(m, Tachikoma.KeyEvent(:enter))
    RE.update!(m, Tachikoma.KeyEvent('d'))
    @test length(m.species) == 2
    @test ncolors(m.rules) == 8
    @test m.messagekind === :info

    # A placed species is kept too, and the last one is never dropped.
    m = RE.EditorModel(UnitSquare)
    RE.update!(m, Tachikoma.KeyEvent('a'))
    m.focus = :construction
    RE.update!(m, Tachikoma.KeyEvent('2'))
    _attach!(m, 1, 1, 3)
    m.focus = :rules
    RE.update!(m, Tachikoma.KeyEvent('d'))
    @test length(m.species) == 2
    @test m.message == "species 2 is placed"
    m = RE.EditorModel(UnitSquare)
    RE.update!(m, Tachikoma.KeyEvent('d'))
    @test length(m.species) == 1
    @test m.messagekind === :warning

    # Digits pick the active species, and only the construction pane has an active species.
    m = RE.EditorModel(UnitSquare)
    RE.update!(m, Tachikoma.KeyEvent('3'))
    @test length(m.species) == 1   # ignored in the rules pane
    m.focus = :construction
    RE.update!(m, Tachikoma.KeyEvent('3'))
    @test length(m.species) == 3
    @test m.active_species == 3

    # The species live in the construction pane, as a strip marking the one the next placement
    # would use, so nothing about them is on screen while there is nothing to place them in.
    m = RE.EditorModel(UnitSquare)
    RE.update!(m, Tachikoma.KeyEvent('a'))
    @test isnothing(Tachikoma.find_text(_render(m), "▶"))
    RE.update!(m, Tachikoma.KeyEvent('b'))
    @test m.showconstruction
    tb = _render(m)
    @test !isnothing(Tachikoma.find_text(tb, "▶"))
    # The strip is inside the construction pane, not the sidebar.
    @test Tachikoma.find_text(tb, "▶").x > RE.SIDEBAR_W

    # Zoom belongs to the construction pane, whose drawing it scales, and does nothing elsewhere.
    m = RE.EditorModel(UnitSquare)
    m.showconstruction = true
    _attach!(m, 1, 1, 3)
    _render(m)
    scale = m.scale
    m.focus = :rules
    RE.update!(m, Tachikoma.KeyEvent('='))
    @test m.scale == scale
    @test !m.manualzoom
    m.focus = :construction
    RE.update!(m, Tachikoma.KeyEvent('='))
    @test m.scale > scale

    # Each enumerated structure is numbered by its position, and the detail box says which it is
    # showing.
    m = RE.EditorModel(UnitSquare)
    _attach!(m, 1, 1, 3)
    RE.togglebond!(m, (2, 4))
    _enumerate!(m)
    @test length(m.polyforms) > 3
    m.focus = :enumeration
    RE.update!(m, Tachikoma.KeyEvent(:right))
    RE.update!(m, Tachikoma.KeyEvent(:right))
    @test m.selected == 3
    @test !isnothing(Tachikoma.find_text(_render(m), "#3"))

    # The status line carries the last message and a running count.
    m = RE.EditorModel(UnitSquare)
    RE.update!(m, Tachikoma.KeyEvent('a'))
    tb = _render(m)
    @test !isnothing(Tachikoma.find_text(tb, "added species 2"))
    @test !isnothing(Tachikoma.find_text(tb, "2 species"))

    # A patchy disk works the same way as a polygon: it is drawn as its bounding circle with its
    # patches named on the rim, and the heading triangle points away from patch 1, which is the
    # only cue a disk gives about its orientation at all.
    disk = PatchyDisk([0.0, 2π / 3, 4π / 3])
    m = RE.EditorModel(disk)
    m.showconstruction = true
    @test length(m.free) == nsites(disk)
    _attach!(m, 1, 1, 2)
    @test RE.inferred_bonds(RE.placements(m), m.species) == [(1, 2)]
    _enumerate!(m)
    @test !isempty(m.polyforms)
    tb = _render(m)
    @test !isnothing(Tachikoma.find_text(tb, "╭─ Construction"))
    cv = Tachikoma.Canvas(20, 10)
    v = RE.fitworld([(1, one(Roly.posetype(disk)))], [disk], 20, 10)
    RE.drawheading!(cv, v, one(Roly.posetype(disk)), disk, RE.outline(disk, one(Roly.posetype(disk))))
    @test any(!iszero, cv.dots)

    # Species created from the rules pane appear in the gallery, not only placed ones.
    m = RE.EditorModel(UnitSquare)
    RE.update!(m, Tachikoma.KeyEvent('a'))
    @test isempty(RE.usedspecies(RE.placements(m))) == false
    @test RE.usedspecies(RE.placements(m)) == [1]   # species 2 has never been placed
    tb = _render(m)
    @test !isnothing(Tachikoma.find_text(tb, "species 1"))
    @test !isnothing(Tachikoma.find_text(tb, "species 2"))

    # Selecting species 2 and then 3 without placing either leaves both in the rules, so the
    # matrix grows rather than swapping one for the other.
    m = RE.EditorModel(UnitSquare)
    m.focus = :construction
    @test ncolors(m.rules) == 4
    RE.update!(m, Tachikoma.KeyEvent('2'))
    @test ncolors(m.rules) == 8
    RE.update!(m, Tachikoma.KeyEvent('3'))
    @test ncolors(m.rules) == 12
    RE.update!(m, Tachikoma.KeyEvent('4'))
    @test ncolors(m.rules) == 16
    RE.update!(m, Tachikoma.KeyEvent('1'))
    @test ncolors(m.rules) == 16  # selecting a lower one does not drop the others

    # The camera holds still while the structure fits, so attaching does not move what is
    # already drawn.
    m = RE.EditorModel(UnitSquare)
    m.showconstruction = true
    _render(m)
    scale0 = m.scale
    _attach!(m, 1, 1, 3)
    _render(m)
    @test m.scale == scale0
    @test (m.cx, m.cy) == (0.0, 0.0)

    # A chain long enough to leave the view zooms out, and never zooms back in on its own.
    for _ in 1:12
        _attach!(m, length(RE.placements(m)), 1, 3)
        _render(m)
    end
    @test m.scale < scale0
    zoomed = m.scale

    # Trimming the chain back down leaves the camera where it was, so the remaining particles
    # do not jump, until `0` asks for a refit.
    for _ in 1:6
        RE.update!(m, Tachikoma.KeyEvent(:backspace))
        _render(m)
    end
    @test m.scale == zoomed
    RE.update!(m, Tachikoma.KeyEvent('0'))
    _render(m)
    @test m.scale > zoomed

    # `-` and `=` zoom, and hold the cursor still on screen while doing it.
    m = RE.EditorModel(UnitSquare)
    m.showconstruction = true
    _attach!(m, 1, 1, 3)
    _render(m)
    s0 = m.scale
    p = RE.anchorposition(m)
    offset(m) = ((p[1] - m.cx) * m.scale, (p[2] - m.cy) * m.scale)
    was = offset(m)
    RE.update!(m, Tachikoma.KeyEvent('='))
    @test m.scale > s0
    @test m.manualzoom
    @test all(isapprox.(offset(m), was))
    RE.update!(m, Tachikoma.KeyEvent('-'))
    @test m.scale ≈ s0
    @test all(isapprox.(offset(m), was))
    RE.update!(m, Tachikoma.KeyEvent('+'))  # a synonym for `=`
    @test m.scale > s0

    # A view chosen by hand is not undone by the next attachment, until `0` asks for a refit.
    zoomed = m.scale
    _attach!(m, 1, 2, 4)
    _render(m)
    @test m.scale == zoomed
    RE.update!(m, Tachikoma.KeyEvent('0'))
    _render(m)
    @test !m.manualzoom

    # Either direction is clamped to a fixed range about the opening view.
    for _ in 1:100
        RE.update!(m, Tachikoma.KeyEvent('='))
    end
    @test m.scale ≈ m.basescale * RE.ZOOM_RANGE
    for _ in 1:400
        RE.update!(m, Tachikoma.KeyEvent('-'))
    end
    @test m.scale ≈ m.basescale / RE.ZOOM_RANGE

    # Free sites are ordered clockwise around the structure, so stepping the anchor walks the
    # perimeter rather than jumping between particles in the order they were placed.
    m = RE.EditorModel(UnitSquare)
    @test issorted(_bearings(m); rev=true)
    _attach!(m, 1, 1, 3)
    _attach!(m, 1, 2, 4)
    _attach!(m, 2, 2, 1)
    _attach!(m, 4, 2, 4)
    @test issorted(_bearings(m); rev=true)
    @test length(unique(first.(m.free))) > 1  # the walk does interleave particles

    # Attaching leaves the cursor next to the site it just used, rather than at whatever now
    # holds the old index.
    m = RE.EditorModel(UnitSquare)
    m.anchor = 3
    before = RE.anchorsite(m).pose.x
    m.incoming = 3
    RE.update!(m, Tachikoma.KeyEvent(:enter))
    moved = sqrt(sum(abs2, RE.anchorsite(m).pose.x - before))
    @test moved < 2 * Roly.bounding_radius(UnitSquare)

    # Gallery scaling: the width is filled first, so species are laid out across rather than in
    # one tall column, and the rows take what height is left.
    @test RE.gallerylayout(2, 4, 38, 14) == (2, 1, RE.GALLERY_ROW)
    @test RE.gallerylayout(6, 4, 38, 21)[1] == 3
    @test isnothing(RE.gallerylayout(6, 4, 38, 4))   # rows too short to read
    @test isnothing(RE.gallerylayout(1, 4, 8, 20))   # too narrow for any drawing

    # When it declines, a one-line-per-species summary takes over, which the gallery never
    # draws, so finding it means the fallback ran. It uses columns too.
    m = _addspecies(7)
    @test !isnothing(Tachikoma.find_text(_render(m, 104, 12), "3 ■"))

    # A matrix too big for the rules pane moves to a full-height pane of its own rather than
    # being dropped for the pair list.
    m = _manyspecies(6)
    @test ncolors(m.rules) == 36
    tb = _render(m, 150, 26)
    @test !isnothing(Tachikoma.find_text(tb, "╭─ Matrix"))
    @test isnothing(Tachikoma.find_text(tb, "wider window"))

    # Too big even for that, and the rules pane says so rather than leaving a gap.
    @test !isnothing(Tachikoma.find_text(_render(m, 104, 26), "wider window"))
    @test isnothing(Tachikoma.find_text(_render(m, 104, 26), "╭─ Matrix"))

    # Small enough and it stays inline, beside the pair list.
    m = _manyspecies(4, UnitSquare)
    tb = _render(m, 170, 32)
    @test isnothing(Tachikoma.find_text(tb, "╭─ Matrix"))
    @test isnothing(Tachikoma.find_text(tb, "wider window"))

    # A patchy disk works the same way as a polygon: it is drawn as its bounding circle with its
    # patches named on the rim, and the heading triangle points away from patch 1, which is the
    # only cue a disk gives about its orientation at all.
    disk = PatchyDisk([0.0, 2π / 3, 4π / 3])
    m = RE.EditorModel(disk)
    m.showconstruction = true
    @test length(m.free) == nsites(disk)
    _attach!(m, 1, 1, 2)
    @test RE.inferred_bonds(RE.placements(m), m.species) == [(1, 2)]
    _enumerate!(m)
    @test !isempty(m.polyforms)
    tb = _render(m)
    @test !isnothing(Tachikoma.find_text(tb, "╭─ Construction"))
    cv = Tachikoma.Canvas(20, 10)
    v = RE.fitworld([(1, one(Roly.posetype(disk)))], [disk], 20, 10)
    RE.drawheading!(cv, v, one(Roly.posetype(disk)), disk, RE.outline(disk, one(Roly.posetype(disk))))
    @test any(!iszero, cv.dots)

    # Species created from the rules pane appear in the gallery, not only placed ones.
    m = RE.EditorModel(UnitSquare)
    RE.update!(m, Tachikoma.KeyEvent('a'))
    @test isempty(RE.usedspecies(RE.placements(m))) == false
    @test RE.usedspecies(RE.placements(m)) == [1]   # species 2 has never been placed
    tb = _render(m)
    @test !isnothing(Tachikoma.find_text(tb, "species 1"))
    @test !isnothing(Tachikoma.find_text(tb, "species 2"))

    # Selecting species 2 and then 3 without placing either leaves both in the rules, so the
    # matrix grows rather than swapping one for the other.
    m = RE.EditorModel(UnitSquare)
    m.focus = :construction
    @test ncolors(m.rules) == 4
    RE.update!(m, Tachikoma.KeyEvent('2'))
    @test ncolors(m.rules) == 8
    RE.update!(m, Tachikoma.KeyEvent('3'))
    @test ncolors(m.rules) == 12
    RE.update!(m, Tachikoma.KeyEvent('4'))
    @test ncolors(m.rules) == 16
    RE.update!(m, Tachikoma.KeyEvent('1'))
    @test ncolors(m.rules) == 16  # selecting a lower one does not drop the others

    # The camera holds still while the structure fits, so attaching does not move what is
    # already drawn.
    m = RE.EditorModel(UnitSquare)
    m.showconstruction = true
    _render(m)
    scale0 = m.scale
    _attach!(m, 1, 1, 3)
    _render(m)
    @test m.scale == scale0
    @test (m.cx, m.cy) == (0.0, 0.0)

    # A chain long enough to leave the view zooms out, and never zooms back in on its own.
    for _ in 1:12
        _attach!(m, length(RE.placements(m)), 1, 3)
        _render(m)
    end
    @test m.scale < scale0
    zoomed = m.scale

    # Trimming the chain back down leaves the camera where it was, so the remaining particles
    # do not jump, until `0` asks for a refit.
    for _ in 1:6
        RE.update!(m, Tachikoma.KeyEvent(:backspace))
        _render(m)
    end
    @test m.scale == zoomed
    RE.update!(m, Tachikoma.KeyEvent('0'))
    _render(m)
    @test m.scale > zoomed

    # `-` and `=` zoom, and hold the cursor still on screen while doing it.
    m = RE.EditorModel(UnitSquare)
    m.showconstruction = true
    _attach!(m, 1, 1, 3)
    _render(m)
    s0 = m.scale
    p = RE.anchorposition(m)
    offset(m) = ((p[1] - m.cx) * m.scale, (p[2] - m.cy) * m.scale)
    was = offset(m)
    RE.update!(m, Tachikoma.KeyEvent('='))
    @test m.scale > s0
    @test m.manualzoom
    @test all(isapprox.(offset(m), was))
    RE.update!(m, Tachikoma.KeyEvent('-'))
    @test m.scale ≈ s0
    @test all(isapprox.(offset(m), was))
    RE.update!(m, Tachikoma.KeyEvent('+'))  # a synonym for `=`
    @test m.scale > s0

    # A view chosen by hand is not undone by the next attachment, until `0` asks for a refit.
    zoomed = m.scale
    _attach!(m, 1, 2, 4)
    _render(m)
    @test m.scale == zoomed
    RE.update!(m, Tachikoma.KeyEvent('0'))
    _render(m)
    @test !m.manualzoom

    # Either direction is clamped to a fixed range about the opening view.
    for _ in 1:100
        RE.update!(m, Tachikoma.KeyEvent('='))
    end
    @test m.scale ≈ m.basescale * RE.ZOOM_RANGE
    for _ in 1:400
        RE.update!(m, Tachikoma.KeyEvent('-'))
    end
    @test m.scale ≈ m.basescale / RE.ZOOM_RANGE

    # Free sites are ordered clockwise around the structure, so stepping the anchor walks the
    # perimeter rather than jumping between particles in the order they were placed.
    m = RE.EditorModel(UnitSquare)
    @test issorted(_bearings(m); rev=true)
    _attach!(m, 1, 1, 3)
    _attach!(m, 1, 2, 4)
    _attach!(m, 2, 2, 1)
    _attach!(m, 4, 2, 4)
    @test issorted(_bearings(m); rev=true)
    @test length(unique(first.(m.free))) > 1  # the walk does interleave particles

    # Attaching leaves the cursor next to the site it just used, rather than at whatever now
    # holds the old index.
    m = RE.EditorModel(UnitSquare)
    m.anchor = 3
    before = RE.anchorsite(m).pose.x
    m.incoming = 3
    RE.update!(m, Tachikoma.KeyEvent(:enter))
    moved = sqrt(sum(abs2, RE.anchorsite(m).pose.x - before))
    @test moved < 2 * Roly.bounding_radius(UnitSquare)
end

@testset "editrules 3d" begin
    # The isometric camera looks along (1, 1, 1), with an orthonormal pair of screen axes and
    # the world's z axis up on the screen rather than down.
    cam = RE.ISOCAM
    @test cam.view ≈ normalize(SVector(1.0, 1.0, 1.0))
    @test abs(dot(cam.view, cam.right)) < 1e-12
    @test abs(dot(cam.view, cam.up)) < 1e-12
    @test abs(dot(cam.right, cam.up)) < 1e-12
    @test norm(cam.right) ≈ 1
    @test norm(cam.up) ≈ 1
    @test dot(cam.up, SVector(0.0, 0.0, 1.0)) > 0
    # Right-handed, so a right-handed world frame stays right-handed on screen: from (1, 1, 1)
    # the x, y and z axes run counterclockwise, not clockwise.
    @test cross(cam.right, cam.up) ≈ cam.view
    @test all(c -> cross(c.right, c.up) ≈ c.view, RE.VIEWPOINTS)

    # A 2D point is already in the projection plane, so the camera leaves it alone and every 2D
    # drawing is unchanged by the projection being there at all.
    @test RE.plane(cam, SVector(3.0, -4.0)) == (3.0, -4.0)

    # The projection is orthographic, so a displacement perpendicular to the view keeps its
    # length. This is what lets the camera be fitted to a bounding radius in 3D as in 2D.
    d = normalize(cross(cam.view, SVector(0.0, 0.3, 1.0)))
    u, v = RE.plane(cam, d)
    @test hypot(u, v) ≈ 1

    # Turning the camera moves the projection but not the world.
    turned = RE.Camera(RE.ISO_AZIMUTH + π / 2, RE.ISO_ELEVATION)
    @test RE.plane(turned, SVector(1.0, 0.0, 0.0)) != RE.plane(cam, SVector(1.0, 0.0, 0.0))
    @test abs(dot(turned.view, turned.right)) < 1e-12

    # A convex particle shows half its faces, and the faces it shows are the ones whose outward
    # normal points at the camera.
    pose = one(Roly.posetype(UnitCube))
    faces = RE.visiblefaces(UnitCube, pose, cam)
    @test length(faces) == 3
    @test all(length(first(f)) == 4 for f in faces)
    @test sort(last.(faces)) == [1, 2, 3]  # one per shade

    # A 3D species with no polyhedron behind it falls back to the silhouette of its bounding
    # sphere, the way a 2D one falls back to a circle.
    sphere = PatchySphere(Cube(), 0.9)
    fallback = RE.visiblefaces(sphere, one(Roly.posetype(sphere)), cam)
    @test length(fallback) == 1
    @test length(fallback[1][1]) == 24

    # Exactly the sites on the near half are drawn, in 3D. In 2D nothing is ever hidden.
    @test count(k -> RE.facing(cam, Roly.bindingsite(UnitCube, k)), 1:nsites(UnitCube)) == 3
    @test all(RE.facing(cam, Roly.bindingsite(UnitSquare, k)) for k in 1:nsites(UnitSquare))

    # A bond between two placed particles is drawn where their faces meet, which in 3D is inside
    # the solid and so not drawn at all.
    @test RE.interiorbonds(UnitSquare)
    @test !RE.interiorbonds(UnitCube)

    # The scanline fill paints the grid points inside the polygon, and nothing for a polygon with
    # no interior.
    painted = Tuple{Int,Int}[]
    RE.fillpolygon!((x, y) -> push!(painted, (x, y)), [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)])
    @test length(painted) == 16
    @test all(0 <= p[1] <= 3 && 0 <= p[2] <= 3 for p in painted)
    empty!(painted)
    RE.fillpolygon!((x, y) -> push!(painted, (x, y)), [(0.0, 0.0), (4.0, 4.0)])
    @test isempty(painted)

    # The eight viewpoints are all the same isometric view, one per octant, and between them
    # they show every face of a convex particle.
    @test length(RE.VIEWPOINTS) == 8
    @test all(c -> all(≈(1 / sqrt(3)), abs.(c.view)), RE.VIEWPOINTS)
    @test length(unique(sign.(c.view) for c in RE.VIEWPOINTS)) == 8
    @test all(k -> any(c -> RE.facing(c, Roly.bindingsite(UnitCube, k)), RE.VIEWPOINTS), 1:nsites(UnitCube))
    # Consecutive viewpoints differ by a quarter turn or a flip, never both, so stepping the
    # camera never jumps to the far side of the structure.
    @test all(n -> dot(RE.VIEWPOINTS[n].view, RE.VIEWPOINTS[mod1(n + 1, 8)].view) ≈ 1 / 3, eachindex(RE.VIEWPOINTS))

    # The editor opens on the first of them, with the cursor already on a site it shows.
    m = RE.EditorModel(UnitCube)
    @test m.viewpoint == 1
    @test RE.camera(m) === RE.ISOCAM
    @test _anchorfacing(m) > 0
    @test length(RE.placements(m)) == 1
    # Only the sites the camera shows are offered, which for a convex particle is half of them.
    @test length(m.free) == nsites(UnitCube) ÷ 2

    # Attaching reads the bond off the geometry, exactly as in 2D.
    m.focus = :construction
    _attachhere!(m, 2)
    @test length(RE.placements(m)) == 2
    @test length(RE.inferred_bonds(RE.placements(m), m.species)) == 1
    @test nbonds(m.rules) == 1

    # The camera never moves on its own: the sites offered to attach to are the ones it shows,
    # so the cursor stays on the near side and the view stays where the user put it.
    m = RE.EditorModel(UnitCube)
    m.focus = :construction
    @test all(((i, k),) -> RE.facing(RE.camera(m), RE.absolutesites(RE.placements(m), m.species)[i][k]), m.free)
    was = m.viewpoint
    for key in (:up, :down, :left, :right, ',', '.')
        for _ in 1:6
            RE.update!(m, Tachikoma.KeyEvent(key))
            @test m.viewpoint == was
            @test _anchorfacing(m) > 0
        end
    end

    # Walking the cursor of a grown structure stays on the near side too.
    m = RE.EditorModel(UnitCube)
    m.focus = :construction
    _attachhere!(m, 2)
    for _ in 1:12
        RE.update!(m, Tachikoma.KeyEvent(:right))
        @test _anchorfacing(m) > 0
    end
    for c in (',', '.')
        for _ in 1:8
            RE.update!(m, Tachikoma.KeyEvent(c))
            @test _anchorfacing(m) > 0
        end
    end

    # Turning the camera offers the sites the new viewpoint shows, and every site of the
    # structure is reachable across the eight of them.
    reachable = Set{Tuple{Int,Int}}()
    for n in 1:length(RE.VIEWPOINTS)
        RE.lookfrom!(m, n)
        union!(reachable, m.free)
        @test all(((i, k),) -> RE.facing(RE.camera(m), RE.absolutesites(RE.placements(m), m.species)[i][k]), m.free)
    end
    RE.lookfrom!(m, 1)
    @test length(reachable) == 2 * nsites(UnitCube) - 2  # every free site of the two cubes

    # The cursor comes with the camera, landing on whichever offered site is nearest to where it
    # was rather than on whatever now holds the old index.
    before = RE.anchorsite(m).pose.x
    RE.turn!(m, 1)
    @test issetequal(m.free, RE.freesites(RE.placements(m), m.species, RE.camera(m)))
    moved = sqrt(sum(abs2, RE.anchorsite(m).pose.x - before))
    @test moved <= 2 * Roly.bounding_radius(UnitCube)

    # The two camera keys step along the list and wrap round it, from any pane.
    m = RE.EditorModel(UnitCube)
    for pane in (:rules, :construction, :enumeration)
        m.focus = pane
        m.viewpoint = 1
        RE.update!(m, Tachikoma.KeyEvent(']'))
        @test m.viewpoint == 2
        RE.update!(m, Tachikoma.KeyEvent('['))
        @test m.viewpoint == 1
        RE.update!(m, Tachikoma.KeyEvent('['))
        @test m.viewpoint == length(RE.VIEWPOINTS)
    end
    m.viewpoint = 1
    for _ in 1:length(RE.VIEWPOINTS)
        RE.update!(m, Tachikoma.KeyEvent(']'))
    end
    @test m.viewpoint == 1

    # Turning is a 2D no-op, there being nothing to turn, and 2D offers every free site.
    flat = RE.EditorModel(UnitSquare)
    RE.turn!(flat, 3)
    @test flat.viewpoint == 1
    @test length(flat.free) == nsites(UnitSquare)

    # `c` clears the drawing and puts the camera back where it opened.
    m = RE.EditorModel(UnitCube)
    m.focus = :construction
    _attachhere!(m, 2)
    RE.turn!(m, 3)
    RE.update!(m, Tachikoma.KeyEvent('c'))
    @test length(RE.placements(m)) == 1
    @test m.viewpoint == 1

    # The camera fits the projection of the structure, not its world coordinates: a bond along
    # the view axis still takes room on screen.
    m = RE.EditorModel(UnitCube)
    box = RE.worldbox(RE.placements(m), m.species, cam)
    r = Roly.bounding_radius(UnitCube)
    @test box[1] ≈ -r && box[2] ≈ r && box[3] ≈ -r && box[4] ≈ r

    # A face's edges take the shade its orientation gives it, the same three a filled drawing
    # would use, so the wireframe still reads as a lit solid.
    ramp = RE.faceramp(RE.speciesrgb(1))
    @test length(unique(c.code for c in ramp)) == 3
    @test all(length(unique(c.code for c in RE.faceramp(RE.speciesrgb(i)))) == 3 for i in 1:8)

    # A 3D scene is drawn as braille, inside the pane it was given and nowhere else, with no
    # cell filled: the drawing is dots, so it reads at four times the vertical resolution a fill
    # would.
    buf = Tachikoma.Buffer(Tachikoma.Rect(1, 1, 60, 24))
    rect = Tachikoma.Rect(10, 5, 30, 12)
    v = RE.fitworld(RE.placements(m), m.species, rect.width, rect.height)
    RE.drawparticles!(buf, rect, v, RE.placements(m), m.species)
    @test isempty(_fills(buf, rect))
    @test count(!=(Tachikoma.EMPTY_CHAR), _chars(buf, rect)) > 5
    @test all(==(Tachikoma.EMPTY_CHAR), _chars(buf, Tachikoma.Rect(1, 1, 9, 24)))
    @test all(==(Tachikoma.EMPTY_CHAR), _chars(buf, Tachikoma.Rect(41, 1, 20, 24)))

    # A particle standing behind another is hidden by it, the drawing being composited in depth
    # order rather than every particle drawing through whatever stands in front of it.
    P = Roly.posetype(UnitCube)
    small = PatchySphere(Cube(), 0.4)  # small enough to sit inside the cube's silhouette
    front = one(P)
    fixed = RE.World(20.0, 0.0, 0.0, rect.width, rect.height, cam)
    function _wire(ps)
        b = Tachikoma.Buffer(Tachikoma.Rect(1, 1, 60, 24))
        RE.drawparticles!(b, rect, fixed, ps, [UnitCube, small])
        return _chars(b, rect)
    end
    alone = _wire([(1, front)])
    back = P(front.x - 3 * cam.view, front.psi)
    @test count(!=(Tachikoma.EMPTY_CHAR), _wire([(2, back)])) > 5  # it does draw, on its own
    @test _wire([(1, front), (2, back)]) == alone
    @test _wire([(1, front), (2, P(front.x + 3 * cam.view, front.psi))]) != alone

    # The whole editor draws a 3D species, and the sidebar offers the camera keys that only 3D
    # has, naming the viewpoint they step through.
    m = RE.EditorModel(UnitCube)
    m.showconstruction = true
    m.focus = :construction
    _attachhere!(m, 2)
    tb = _render(m)
    @test !isempty(Tachikoma.find_text(tb, "Construction"))
    @test !isempty(Tachikoma.find_text(tb, "2 particles"))
    @test !isempty(Tachikoma.find_text(tb, "view 1/8"))
    @test isnothing(Tachikoma.find_text(_render(RE.EditorModel(UnitSquare)), "view 1/8"))

    # The enumeration runs on 3D rules and draws the structures it finds.
    _enumerate!(m)
    @test !m.stale
    @test length(m.polyforms) > 1
    @test all(Roly.dimension(Roly.species(bindingrules(p))[1]) == 3 for p in m.polyforms)
    m.focus = :enumeration
    tb = _render(m)
    @test !isempty(Tachikoma.find_text(tb, "Enumeration"))

    # Zooming holds the cursor still in the projection plane, the same as it does in 2D.
    m = RE.EditorModel(UnitCube)
    m.focus = :construction
    _render(m)
    s0 = m.scale
    RE.update!(m, Tachikoma.KeyEvent('='))
    @test m.scale > s0
    @test m.manualzoom
    RE.update!(m, Tachikoma.KeyEvent('-'))
    @test m.scale ≈ s0

    # The axis indicator says which way the camera is looking, and turns with it.
    axesbuf = Tachikoma.Buffer(Tachikoma.Rect(1, 1, 60, 24))
    axesrect = Tachikoma.Rect(5, 5, 20, 10)
    RE.drawaxes!(axesbuf, axesrect, RE.ISOCAM)
    letters = filter(c -> c in ('x', 'y', 'z'), vec(_chars(axesbuf, axesrect)))
    @test sort(letters) == ['x', 'y', 'z']
    # Too small a pane leaves it off rather than drawing over the structure.
    tiny = Tachikoma.Buffer(Tachikoma.Rect(1, 1, 60, 24))
    RE.drawaxes!(tiny, Tachikoma.Rect(5, 5, RE.AXES_W, RE.AXES_H), RE.ISOCAM)
    @test all(==(Tachikoma.EMPTY_CHAR), _chars(tiny, Tachikoma.Rect(1, 1, 60, 24)))
    # It appears in the construction pane and the inspector of a 3D editor, and in neither of a
    # 2D one, where there is no orientation to report.
    m = RE.EditorModel(UnitCube)
    m.showconstruction = true
    @test _axisletters(_render(m)) == 3
    m.pair = (4, 5)
    RE.togglebond!(m, m.pair)
    _enumerate!(m)
    @test !isempty(m.polyforms)
    @test _axisletters(_render(m)) == 6  # the inspector carries one too
    flat2 = RE.EditorModel(UnitSquare)
    flat2.showconstruction = true
    flat2.pair = (1, 3)
    RE.togglebond!(flat2, flat2.pair)
    _enumerate!(flat2)
    @test !isempty(flat2.polyforms)
    @test _axisletters(_render(flat2)) == 0

    # Windows work the same in 3D: a second one keeps the first one's bonds on the books.
    for sp in (UnitCube, UnitIcosahedron, UnitPrism(3))
        m = RE.EditorModel(sp)
        m.focus = :construction
        _attachhere!(m, 2)
        was = RE.effectivematrix(m)
        @test any(was)
        RE.update!(m, Tachikoma.KeyEvent('n'))
        @test length(m.windows) == 2
        @test length(RE.placements(m)) == 1
        @test RE.effectivematrix(m) == was
        RE.update!(m, Tachikoma.KeyEvent('c'))  # closes the empty window
        @test length(m.windows) == 1
        @test RE.effectivematrix(m) == was
    end

    # Color labels are one character while they can be, `1`-`9` then `a`-`z` then `A`-`Z`, and
    # never repeat inside that range: two sites sharing a label would name different things the
    # same way.
    @test RE.labelwidth(1) == 1
    @test RE.labelwidth(RE.NARROW_COLORS) == 1
    @test RE.labelwidth(RE.NARROW_COLORS + 1) == 2
    @test length(unique(RE.colorlabel(c) for c in 1:RE.NARROW_COLORS)) == RE.NARROW_COLORS
    @test RE.colorlabel(9) == '9'
    @test RE.colorlabel(10) == 'a'
    @test RE.colorlabel(36) == 'A'
    @test RE.colorlabel(RE.NARROW_COLORS) == 'Z'
    # Past that they are two characters, in decimal, which reads without counting the alphabet.
    @test RE.colorlabel(7, 2) == "07"
    @test RE.colorlabel(62, 2) == "62"
    @test length(unique(RE.colorlabel(c, 2) for c in 1:RE.MAX_COLORS)) == RE.MAX_COLORS

    # A species is refused once its colors would run past what can be labelled, rather than the
    # labels wrapping round.
    big = RE.EditorModel(UnitIcosahedron)
    while RE.roomforspecies(big)
        n = length(big.species)
        RE.update!(big, Tachikoma.KeyEvent('a'))
        @test length(big.species) == n + 1
    end
    @test RE.totalcolors(big) + nsites(UnitIcosahedron) > RE.MAX_COLORS
    n = length(big.species)
    RE.update!(big, Tachikoma.KeyEvent('a'))
    @test length(big.species) == n
    @test occursin("no room", big.message)
    big.focus = :construction
    RE.update!(big, Tachikoma.KeyEvent('9'))  # a digit naming a species there is no room for
    @test length(big.species) == n
    @test occursin("no room", big.message)

    # The label width follows the color count, and the matrix widens to match.
    small = RE.EditorModel(UnitCube)
    @test small.labelw == 1
    wide = RE.EditorModel(UnitIcosahedron)
    RE.update!(wide, Tachikoma.KeyEvent('a'))
    RE.update!(wide, Tachikoma.KeyEvent('a'))
    RE.update!(wide, Tachikoma.KeyEvent('a'))
    @test RE.totalcolors(wide) > RE.NARROW_COLORS
    @test wide.labelw == 2
    @test RE.matrixwidth(wide) > RE.matrixwidth(small)
    @test RE.matrixwidth(wide) == 3 * (ncolors(wide.rules) + 1)  # label plus separator per color
    @test RE.label(wide, 62) == "62"
    @test RE.label(small, 6) == "6"
    @test !isnothing(_render(wide))  # and the whole editor still draws at that width
end
