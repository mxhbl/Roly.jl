const INERT_COLOR = colorant"#E7E7E7"
const DEFAULT_COLORS = [
    [colorant"#B2D0EA", colorant"#54A3E4", colorant"#1A78C6", colorant"#00549A"],
    [colorant"#F39F9D", colorant"#ED7E7C", colorant"#E75451", colorant"#BE2E2C"],
    [colorant"#FEE3B5", colorant"#FFC563", colorant"#FFB12F", colorant"#FFA000"],
    [colorant"#DBB9E4", colorant"#D99DE8", colorant"#B571C4", colorant"#8E4E9C"],
    [colorant"#AAF5FF", colorant"#86E4F1", colorant"#43CFE2", colorant"#20BCD2"],
    [colorant"#BCF3E1", colorant"#3DD4A4", colorant"#19B684", colorant"#008D60"],
    [colorant"#FFD4B6", colorant"#F7AF7D", colorant"#FF8835", colorant"#ED6304"],
    [colorant"#C5C9FE", colorant"#A9AFFF", colorant"#7C85FF", colorant"#4854FD"],
    [colorant"#F4FFC4", colorant"#E0F290", colorant"#B9D63E", colorant"#859F13"],
    [colorant"#FFB2C5", colorant"#FF8EA9", colorant"#FF5C83", colorant"#E32654"],
    [colorant"#95FABA", colorant"#4DD980", colorant"#21AB53", colorant"#008832"]
]

# How far a color is blended towards white. Surface that cannot bond is drawn as a pale ghost
# of the species color, so an assembly reads as a whole while the sites that can actually bond
# keep the visual weight. A polyhedron's whole surface is tinted, since it is bulk rather than
# a marker; a patchy particle's body is fainter still, since its patches sit on top of it.
const BODY_TINT = 0.72
const FACE_TINT = 0.35
const PATCHY_BODY_TINT = 0.68

function species_palette(spcs_index::Int, n::Int)
    base = DEFAULT_COLORS[mod1(spcs_index, length(DEFAULT_COLORS))]
    return cgrad(base, n; categorical=true)
end

"""
    species_basecolor(spcs_index)

The color that identifies a species irrespective of its site count. Taken from the middle of
the ramp rather than from `species_palette`, whose first entry is the palest color there is
and washes out to nothing once tinted.
"""
species_basecolor(spcs_index::Int) = DEFAULT_COLORS[mod1(spcs_index, length(DEFAULT_COLORS))][3]

"""
    _tint(c, t)

Blend `c` towards white by the fraction `t`, keeping its hue.
"""
_tint(c, t::Real) = Makie.Colors.weighted_color_mean(1 - t, RGBf(c), colorant"white")

"""
    _shade(c, t)

Blend `c` towards black by the fraction `t`, keeping its hue.
"""
_shade(c, t::Real) = Makie.Colors.weighted_color_mean(1 - t, RGBf(c), colorant"black")

# A fixed world-space key light, and how much of the color is ambient versus directional.
const LIGHT_DIRECTION = normalize(Vec3f(0.35, -0.6, 0.72))
const LIGHT_AMBIENT = 0.62
const LIGHT_DIFFUSE = 0.38

"""
    _lambert(c, n)

Darken `c` by Lambert's cosine law for a surface with unit normal `n`, preserving its alpha.

Curved surfaces need shading to read as curved, but Makie applies shading per plot, and the
whole point of merging particles into one mesh is that everything lives in one plot. Baking
the shading into the vertex colors keeps a single flat-shaded mesh, and keeps control of how
far the lighting is allowed to darken a deliberately faint color.
"""
function _lambert(c, n)
    lit = LIGHT_AMBIENT + LIGHT_DIFFUSE * max(0, dot(n, LIGHT_DIRECTION))
    return RGBAf(_shade(c, 1 - lit), Makie.Colors.alpha(RGBAf(c)))
end

"""
    _resolve_colors(spcs, speciesindex, rules, sitecolor)

Work out which palette a species is drawn in, which of its sites can bond, and what color
each site gets.

`speciesindex` selects the palette; when it is `nothing` it is looked up in `rules`, falling
back to `1`. Sites no bond in `rules` can use are drawn in `inert_color` if one is given, and
otherwise as a tint of their own color by `inert_tint`; bonding sites are tinted by
`bond_tint`. Without a `rules` there is no notion of a bond, so every site counts as bonding.
`sitecolor` is an optional callback `(speciesindex, siteindex) -> color` that overrides all
of that, tinting included.

Returns `(speciesindex, palette, sitecolors, bonding)`.
"""
function _resolve_colors(
    spcs, speciesindex, rules, sitecolor; bond_tint=0, inert_tint=BODY_TINT, inert_color=INERT_COLOR
)
    n = nsites(spcs)
    if isnothing(speciesindex)
        speciesindex = isnothing(rules) ? 1 : something(findfirst(==(spcs), species(rules)), 1)
    end
    pal = species_palette(speciesindex, n)
    bonding = isnothing(rules) ? trues(n) : [!isinert(rules, SpeciesSiteLoc(speciesindex, i)) for i in 1:n]

    inert(i) = isnothing(inert_color) ? _tint(pal[i], inert_tint) : RGBf(inert_color)
    colors = if isnothing(sitecolor)
        [bonding[i] ? _tint(pal[i], bond_tint) : inert(i) for i in 1:n]
    else
        [sitecolor(speciesindex, i) for i in 1:n]
    end
    return speciesindex, pal, colors, bonding
end
