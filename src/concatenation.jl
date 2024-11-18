
function open_bond(p::Polyform, assembly_system::AssemblySystem, j::Integer)
    for v in 1:nvertices(p)
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
    particle, side = vertex2particle(p, v)
    spc = species(p)[particle]
    geom = geometries[spc]
    return p.xs[particle] + rotate(geom.xs[side], p.ψs[particle])
end

function tile_check(p::Polyform{D}, assembly_system::AssemblySystem, tile_edges, v1, v2) where {D}
    @assert D == 2  # TODO only works in 2d

    if maximum(p.cyclic) >= D
        return false
    end
    v = findfirst(x->x==1, p.cyclic)
    w = p.bond_partners[v]

    xs1, ψs1 = get_attached_coordinates(p, p, assembly_system, v, w)
    xs2, ψs2 = get_attached_coordinates(p, p, assembly_system, v1, v2)

    can_tile, new_edges = find_bonds(p, p, assembly_system, xs1=xs1, ψs1=ψs1, xs2=xs2, ψs2=ψs2)
    !can_tile && return false

    for e in new_edges
        if (e ∉ tile_edges || reverse(e) ∉ tile_edges) ||
            p.cyclic[e.src] == 0
            return false
        end
    end

    xs2, ψs2 = get_attached_coordinates(p, p, assembly_system, v2, v1)

    can_tile, new_edges = find_bonds(p, p, assembly_system, xs1=xs1, ψs1=ψs1, xs2=xs2, ψs2=ψs2)
    !can_tile && return false

    for e in new_edges
        if (e ∉ tile_edges || reverse(e) ∉ tile_edges) ||
            p.cyclic[e.src] == 0
            return false
        end
    end

    return true
end

function attach_monomer!(p::Polyform{D,T,F}, v::Integer, partner_label::Integer, assembly_system::AssemblySystem, fillhash::Bool=false; findcycles=false) where {D,T,F}
    p.bond_partners[v] != 0 && return false

    spcs2, site2 = label2spcssite(partner_label, assembly_system)
    bblock = buildingblocks(assembly_system)[spcs2]
    v2 = particle2vertex(bblock, 1, site2)

    xs2, ψs2 = get_attached_coordinates(p, bblock, assembly_system, v, v2)
    duplicate_species = (i for (i, s) in enumerate(species(p)) if s == spcs2 && p.ψs[i] ≈ only(ψs2))

    success, new_edges = find_bonds(p, bblock, assembly_system; xs2=xs2, ψs2=ψs2)
    !success && return false

    istiling = false
    already_tiling = false
    iscap = false
    if findcycles
        tile_edges = eltype(new_edges)[]
        for dup in duplicate_species
            vdup = particle2vertex(p, dup, site2)

            vpartner = p.bond_partners[vdup]
            xs_tile, ψs_tile = get_attached_coordinates(p, p, assembly_system, v, vdup)
            if vpartner == 0  # vdup is open
                # Test whether tiling is sterically possible
                if istiling
                    already_tiling = istiling
                end
                istiling, tile_edges = find_bonds(p, p, assembly_system, xs2=xs_tile, ψs2=ψs_tile)
                (istiling && already_tiling) && return false
                ## TODO if there are parallel bonds, we need to do extra steps
                ## If parallel bonds form, we also need to construct the non-cylic str (i.e. just do normal attachment)
            else
                # TODO: optimize this
                if vpartner != 0
                    blocking_particle, _ = vertex2particle(p, vpartner)
                    blocking_vertices = particle2multivertex(p, blocking_particle)
                    disconnected_vertices = [k for k in vertices(p.anatomy) if !has_path(p.anatomy, v, k, exclude_vertices=blocking_vertices)]
                    vdup in disconnected_vertices && continue
                end

                # Test whether tiling would be sterically possible if blocking particle was not there
                # TODO need to test every path individually
                pathitr = PathIterator(p.anatomy, v, vdup)
                for path in pathitr
                    forbidden_particles = unique([vertex2particle(p, k)[1] for k in findall(iszero, path)])
                    could_tile, _ = find_bonds(p, p, assembly_system, xs2=xs_tile, ψs2=ψs_tile, ignore_parts1=forbidden_particles, ignore_parts2=forbidden_particles)
                    if could_tile
                        iscap = true
                        break
                    end
                end
            end
        end
    end

    if iscap && !istiling
        return false
    end

    if !istiling && any(!iszero, p.cyclic)
        return false
    end
    if istiling && any(!iszero, p.cyclic)
        v, w = first(tile_edges).src, first(tile_edges).dst
        check = tile_check(p, assembly_system, tile_edges, v, w)
        !check && return false
    end
    # TODO: do a tile check if there are multiple cycles
    # only allow d cycles if the structure tiles in d dimensions (or the cyclies are "parallel")
    # Perform tile check by extending all cycles simultaneously and
    # check if they are compatible!
    # > d cycles must be impossible, but make sure we dont overreject here
    # if there are multiple choices
    
    new_edges = !istiling ? new_edges : tile_edges
    Δnv = !istiling ? nvertices(p) : 0
    if !istiling
        append!(p.xs, xs2)
        append!(p.ψs, ψs2)
        append!(p.species, spcs2)
        p.encoder = concatenate(p.encoder, bblock.encoder)
        append!(p.bond_partners, zeros(T, nvertices(bblock)))
        append!(p.cyclic, zeros(Bool, nvertices(bblock)))
        NautyGraphs.blockdiag!(p.anatomy, bblock.anatomy)
    end

    for edge in new_edges
        vi, vj = edge.src, edge.dst
        vj += Δnv

        add_edge!(p.anatomy, vi, vj)
        add_edge!(p.anatomy, vj, vi)
        p.bond_partners[vi] = vj
        p.bond_partners[vj] = vi

        tile_class = maximum(p.cyclic) + 1
        if istiling 
            p.cyclic[vi] = tile_class 
            p.cyclic[vj] = tile_class 
        end
    end

    if fillhash
        n, _, _, hashval = NautyGraphs._nautyhash(p.anatomy)
        p.σ = n
        # TODO: clean this up
        p.anatomy.hashval = hashval
    else
        p.σ = 0
    end
    return true
end

function get_attached_coordinates(p1::Polyform, p2::Polyform, assembly_system::AssemblySystem, v1, v2)
    geoms = geometries(assembly_system)
    spcs1 = species(p1)
    spcs2 = species(p2)

    part1, site1 = vertex2particle(p1, v1)
    part2, site2 = vertex2particle(p2, v2)

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

            vs_i = particle2multivertex(p1, i, si)
            vs_j = particle2multivertex(p2, j, sj)
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
    

function lower!(p::Polyform)
    # if the structure has a cyclic bond, we should
    # remove the bond instead of a particle
    # I guess we always remove the bond with the highest label?
    n = size(p)
    n == 0 && return false

    k = 1
    if n > 1
        while true
            part, _ = vertex2particle(p, k)
            vertices = particle2multivertex(p, part)
            !are_bridge(p.anatomy, vertices) && break
            k += 1
        end
    end

    del_part, _ = vertex2particle(p, k)
    del_vs = particle2multivertex(p, del_part)
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
    deleteat!(p.encoder, del_part)

    canonize!(p)
    return true
end