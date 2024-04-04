import Luxor
using Plots: @colorant_str, palette, color
using Statistics: mean

const default_colors = palette([colorant"#7E9EC4", colorant"#E64B48", colorant"#FFC200",
                          colorant"#C87CDA", colorant"#78F8F0", colorant"#34DA3A", colorant"#FF629A"])

function draw_ngon(n, a; c="grey90", shades=nothing)
    if isnothing(shades)
        shades = [1 - i / n for i in 0:(n - 1)]
    end
    
    colors = color(c) .* shades

    a_out = a
    a_in = a * 5 / 8
    corner_r = a / 30 * n

    angles = [i * 2π / n for i in 0:(n - 1)]
    orientation = 2π / n * (floor(n // 2) - 1) - π / 2
    if iseven(n)
        orientation -= π / n
    end

    interior_angle = (n - 2) / n * π

    Luxor.setline(a / 30)

    Luxor.@layer begin
        for (c, θ) in zip(colors, angles)
            Luxor.@layer begin
                Luxor.rotate(θ)
                ngon_outer = Luxor.ngonside(Luxor.O, a_out, n, orientation; vertices=true)
                clip_points = 1.2 * [ngon_outer[2], Luxor.Point(0, 0), ngon_outer[1]]
                clip_poly = Luxor.poly(clip_points, :clip)

                Luxor.sethue(c)
                Luxor.polysmooth(ngon_outer, corner_r, :fill)
            end
        end

        # Draw the outline
        Luxor.sethue("black")
        ngon_outer = Luxor.ngonside(Luxor.O, a_out, n, orientation; vertices=true)
        Luxor.polysmooth(ngon_outer, corner_r, :stroke)

        # Draw the hole in the center
        ngon_inner = Luxor.ngonside(Luxor.O, a_in, n, orientation; vertices=true)

        Luxor.sethue("white")
        Luxor.polysmooth(ngon_inner, corner_r, :fill)

        Luxor.sethue("black")
        Luxor.polysmooth(ngon_inner, corner_r, :stroke)

        # Draw lines delimiting the binding sides
        offset = 1 / sin(interior_angle / 2) - 1
        multiplyer = n == 3 ? 2 : 1 # don't ask me why this is needed
        for (ngi, ngo) in zip(ngon_inner, ngon_outer)
            Luxor.line(ngi * (1 - multiplyer * corner_r / a_in * offset),
                       ngo * (1 - multiplyer * corner_r / a_out * offset), :stroke)
        end
    end
end

function draw_polyform(pform, assembly_system; r=200, species_colors=nothing)
    spes = species(pform)

    if isnothing(species_colors)
        species_colors = Dict(i => "grey90" for i in values(spes))
    end

    d = r
    xs = [d * SVector(1, -1) .* x for x in pform.xs]
    ψs = [ψ.θ for ψ in pform.ψs]

    polygons = [length(assembly_system.geometries[k])
                for k in species(pform)]

    # Make sure the center of mass is centered
    x_com = mean(xs)
    xs .-= Ref(x_com)

    # Draw all the triangles at the correct positions
    for (particle, (x, ψ)) in enumerate(zip(xs, ψs))
        Luxor.@layer begin
            Luxor.translate(x...)
            Luxor.rotate(-ψ * π)
            draw_ngon(polygons[particle], r; c=species_colors[spes[particle]])
        end
    end

    return
end

function draw_polyform_grid(structures, assembly_system; r=50, box_size=5, ncols=4, species_colors=nothing,
                         structure_names=nothing)
    nstructures = length(structures)
    nrows = nstructures ÷ ncols
    if nstructures % ncols != 0
        nrows += 1
    end

    structures = [structures[i] for i in sort(collect(keys(structures)))]

    tiles = Luxor.Tiler(box_size * r * ncols, box_size * r * nrows, nrows, ncols)

    if isnothing(structure_names)
        i = 1
        for (structure, (pos, n)) in zip(structures, tiles)
            Luxor.@layer begin
                Luxor.translate(pos)
                draw_polyform(structure, assembly_system; r=r, species_colors=species_colors)
            end
            i += 1
        end
    else
        i = 1
        for (structure, (pos, n), name) in zip(structures, tiles, structure_names)
            Luxor.@layer begin
                Luxor.translate(pos)
                draw_polyform(structure, assembly_system; r=r, species_colors=species_colors)
                Luxor.fontsize(r / 2)
                Luxor.text(string(name), Luxor.Point(-2 * r / 2, -3 * r / 2))
            end
            i += 1
        end
    end
end

function draw_polyforms(polyforms::AbstractVector{<:Polyform}, assembly_system::AssemblySystem, filename::String; r=50, ncols=4, box_size=5, species_colors=default_colors, size=(1000, 1000))
    extension = splitext(filename)[end]
    if extension == ".png"
        return Luxor.@png(draw_polyform_grid(polyforms, assembly_system, r=r, box_size=box_size, ncols=ncols, species_colors=species_colors), size[1], size[2], filename)
    elseif extension == ".pdf"
        return Luxor.@pdf(draw_polyform_grid(polyforms, assembly_system, r=r, box_size=box_size, ncols=ncols, species_colors=species_colors), size[1], size[2], filename)
    else
        error("Can only render png or pdf.")
    end
end