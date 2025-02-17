
function open_bond(p::Polyform, assembly_system::AssemblySystem, j::Integer)
    for v in Iterators.reverse(p.canonical_order)
        !iszero(p.bond_partners[v]) && continue
        label = vertex2label(p, assembly_system, v)
        partner_label, Δj = @views find_nth(!iszero, assembly_system.intmat[:, label], j)
        j -= Δj
        if isnothing(partner_label)
            continue
        else
            return v, partner_label
        end
    end

    return 0, 0
end

function get_sitepos(p::Polyform, geometries, v)
    particle, side = vertex2particle(p, assembly_system, v)
    spc = species(p)[particle]
    geom = geometries[spc]
    return p.xs[particle] + rotate(geom.xs[side], p.ψs[particle])
end

function attach_monomer!(p::Polyform{D,T,F}, v::Integer, partner_label::Integer, assembly_system::AssemblySystem, canonize::Bool=false) where {D,T,F}
    p.bond_partners[v] != 0 && return false

    nvp = nvertices(p)
    spcs2, site2 = label2spcssite(partner_label, assembly_system)
    bblock = buildingblocks(assembly_system)[spcs2]
    v2 = particle2vertex(bblock, assembly_system, 1, site2)

    xs2, ψs2 = get_attached_coordinates(p, bblock, assembly_system, v, v2)
    success, new_edges = find_bonds(p, bblock, assembly_system; xs2=xs2, ψs2=ψs2)
    !success && return false

    append!(p.bond_partners, zeros(T, nvertices(bblock)))
    NautyGraphs.blockdiag!(p.anatomy, bblock.anatomy)
    for edge in new_edges
        vi, vj = edge.src, edge.dst
        vj += nvp

        add_edge!(p.anatomy, vi, vj)
        add_edge!(p.anatomy, vj, vi)
        p.bond_partners[vi] = vj
        p.bond_partners[vj] = vi
    end

    append!(p.xs, xs2)
    append!(p.ψs, ψs2)
    append!(p.species, spcs2)
    resize!(p.canonical_order, nvp + nvertices(bblock))

    if canonize
        canonize!(p)
    else
        p.σ = 0
        p.canonical_order .= 1:nvertices(p)
    end
    return true
end

function get_attached_coordinates(p1::Polyform, p2::Polyform, assembly_system::AssemblySystem, v1, v2)
    geoms = geometries(assembly_system)
    spcs1 = species(p1)
    spcs2 = species(p2)

    part1, site1 = vertex2particle(p1, assembly_system, v1)
    part2, site2 = vertex2particle(p2, assembly_system, v2)

    Δx, Δψ = attachment_offset(site1, site2, geoms[spcs1[part1]], geoms[spcs2[part2]])
    Δx = p1.xs[part1] + rotate(Δx, p1.ψs[part1])
    Δψ = Δψ * p1.ψs[part1]

    xs2, ψs2 = copy(p2.xs), copy(p2.ψs)
    grab!(xs2, ψs2, part2, Δx, Δψ)
    return xs2, ψs2
end

function find_bonds(p1::Polyform, p2::Polyform, assembly_system::AssemblySystem; xs1=nothing, ψs1=nothing, xs2=nothing, ψs2=nothing, ignore_parts1=nothing, ignore_parts2=nothing)
    geoms = geometries(assembly_system)
    intmat = interaction_matrix(assembly_system)
    spcs1 = species(p1)
    spcs2 = species(p2)
    
    xs1 = isnothing(xs1) ? p1.xs : xs1
    ψs1 = isnothing(ψs1) ? p1.ψs : ψs1
    xs2 = isnothing(xs2) ? p2.xs : xs2
    ψs2 = isnothing(ψs2) ? p2.ψs : ψs2

    new_edges = Edge{Cint}[]
    for (i, (xi, ψi)) in enumerate(zip(xs1, ψs1)), (j, (xj, ψj)) in enumerate(zip(xs2, ψs2))
        if !isnothing(ignore_parts1) && i in ignore_parts1 continue end
        if !isnothing(ignore_parts2) && j in ignore_parts2 continue end
        cstat, bound_sites = contact_status(xj - xi, ψi, ψj, geoms[spcs1[i]], geoms[spcs2[j]])
        # Return if particles overlap
        !cstat && return false, new_edges

        for sites in bound_sites
            si, sj = sites

            vs_i = particle2multivertex(p1, assembly_system, i, si)
            vs_j = particle2multivertex(p2, assembly_system, j, sj)
            label_i = p1.anatomy.labels[first(vs_i)]
            label_j = p2.anatomy.labels[first(vs_j)]

            # Return if a disallowed bond is found
            !intmat[label_i, label_j] && return false, new_edges

            bond_to_edge!(new_edges, vs_i, vs_j)
        end
    end
    return true, new_edges
end

function bond_to_edge!(edges::AbstractVector{<:Edge}, vs_i::AbstractVector{<:Integer}, vs_j::AbstractVector{<:Integer})
    ni, nj = length(vs_i), length(vs_j)
    shared_symmetry = max(ni, nj) // min(ni, nj)

    if (ni == nj == 1) || !isinteger(shared_symmetry)
        edge = eltype(edges)(first(vs_i), first(vs_j))
        push!(edges, edge)
        return
    end

    sym = numerator(shared_symmetry)
    if ni <= nj
        vs1, vs2 = vs_i, reverse(vs_j[1:sym:end])
    else
        vs1, vs2 = vs_i[1:sym:end], reverse(vs_j)
    end

    for (v1, v2) in zip(vs1, vs2)
        edge = eltype(edges)(v1, v2)
        push!(edges, edge)
    end
    return
end

function raise!(p::Polyform{D,T,F}, k::Integer, assembly_system::AssemblySystem) where {D,T,F}
    v, partner_label = open_bond(p, assembly_system, k)
    iszero(v) && return false, k

    success = attach_monomer!(p, v, partner_label, assembly_system)
    !success && return raise!(p, k+1, assembly_system)

    canonize!(p)
    return true, k
end
    

function lower!(p::Polyform, assembly_system::AssemblySystem)
    n = size(p)
    n == 0 && return false

    k = 1
    if n > 1
        for l in Iterators.reverse(p.canonical_order)
            part, _ = vertex2particle(p, assembly_system, l)
            vertices = particle2multivertex(p, assembly_system, part)
            if !are_bridge(p.anatomy, vertices)
                k = l
                break
            end
        end
    end

    del_part, _ = vertex2particle(p, assembly_system, k)
    del_vs = particle2multivertex(p, assembly_system, del_part)
    sort!(del_vs)

    deleteat!(p.xs, del_part)
    deleteat!(p.ψs, del_part)
    deleteat!(p.species, del_part)

    # open up bonds that are now free
    partners = p.bond_partners[del_vs]
    for part in partners
        if part == 0
            continue
        end
        p.bond_partners[part] = 0
    end
    deleteat!(p.bond_partners, del_vs)
    
    vertex_shift(v) = sum(x -> x <= v, del_vs)
    @views p.bond_partners .-= vertex_shift.(p.bond_partners)

    rem_vertices!(p.anatomy, del_vs)

    deleteat!(p.canonical_order, del_vs)
    canonize!(p)
    return true
end