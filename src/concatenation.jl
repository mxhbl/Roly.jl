function open_bond(p::Polyform, assembly_system::AssemblySystem, j::Integer)
    for v in 1:length(p.bond_partners)
        p.bond_partners[v] != 0 && continue

        label = vertex2label(p, assembly_system, v)
        
        partner_label, Δj = @views find_nth(!iszero, assembly_system.intmat[:, label], j)
        j -= Δj
        if isnothing(partner_label)
            continue
        else
            return v, partner_label
        end
    end

    return nothing, nothing
end

function get_sitepos(p::Polyform, geometries, v)
    particle, side = vertex2particle(p, v)
    spc = species(p)[particle]
    geom = geometries[spc]
    return p.xs[particle] + Roly.rotate(geom.xs[side], p.ψs[particle])
end


function tile_check(p::Polyform, assembly_system, vert_i, vert_j)
    geoms = geometries(assembly_system)
    spcs = species(p)
    x0_i = get_sitepos(p, geoms, vert_i)
    x0_j = get_sitepos(p, geoms, vert_j)
    xs1 = p.xs
    xs2 = p.xs .+ Ref(x0_j - x0_i)
    ψs1 = ψs2 = p.ψs

    anatomy_edges = Edge[]

    for (i, (xi, ψi)) in enumerate(zip(xs1, ψs1)), (j, (xj, ψj)) in enumerate(zip(xs2, ψs2))
    cstat, pairs = contact_status(xj - xi, ψi, ψj, geoms[spcs[i]], geoms[spcs[j]])
        if !cstat
            # println("overlap")
            return false, anatomy_edges
        end
        for pair in pairs
            face_i, face_j = pair

            vs_i = p.encoder.fwd[i][face_i]
            vs_j = p.encoder.fwd[j][face_j]
            lbl_i = p.anatomy.labels[first(vs_i)]
            lbl_j = p.anatomy.labels[first(vs_j)]

            if !interaction_matrix(assembly_system)[lbl_i, lbl_j]
                # println("blocked bond")
                return false, anatomy_edges
            end
        
            aedges = add_bondedges(convert(Vector{Cint}, vs_i), convert(Vector{Cint}, vs_j))
            append!(anatomy_edges, aedges)
        end
    end

    return true, anatomy_edges
end

function attach_monomer!(p::Polyform{D,T,F}, v::Integer, partner_label::Integer, assembly_system::AssemblySystem, fillhash::Bool=false) where {D,T,F}
    # v:: vertex of polyform p
    # blabel:: site label of building block tb attached
    p.bond_partners[v] != 0 && return false

    n = size(p)
    nvert = nvertices(p)
    geoms = geometries(assembly_system)
    spcs = species(p)

    part_i, site_i = vertex2particle(p, v)
    spcs_i = spcs[part_i]

    spcs_j, site_j = label2spcssite(partner_label, assembly_system)
    bblock = buildingblocks(assembly_system)[spcs_j]

    geom_i = geoms[spcs_i]
    geom_j = geoms[spcs_j]

    Δx, Δψ = attachment_offset(site_i, site_j, geom_i, geom_j)
    ψj = Δψ * p.ψs[part_i]
    xj = p.xs[part_i] + rotate(Δx, p.ψs[part_i])

    new_edges = Edge{Cint}[]
    imat = interaction_matrix(assembly_system)
    for (i, (xi, ψi)) in enumerate(zip(p.xs, p.ψs))
        cstat, bound_sites = contact_status(xj - xi, ψi, ψj, geoms[spcs[i]], geom_j)
        # Return if particles overlap
        !cstat && return false

        for sites in bound_sites
            si, sj = sites

            vs_i = particle2vertices(p, i, si)
            vs_j = particle2vertices(bblock, 1, sj)
            label_i = p.anatomy.labels[first(vs_i)]
            label_j = bblock.anatomy.labels[first(vs_j)]

            # Return if a disallowed bond forms
            !imat[label_i, label_j] && return false

            bond_to_edge!(new_edges, vs_i, vs_j .+ nvert)
        end
    end

    push!(p.xs, xj)
    push!(p.ψs, ψj)
    push!(p.species, spcs_j)
    p.encoder = concatenate(p.encoder, bblock.encoder)
    append!(p.bond_partners, zeros(T, nvertices(bblock)))

    NautyGraphs.blockdiag!(p.anatomy, bblock.anatomy)
    for edge in new_edges
        vi, vj = edge.src, edge.dst

        add_edge!(p.anatomy, vi, vj)
        add_edge!(p.anatomy, vj, vi)
        p.bond_partners[vi] = vj
        p.bond_partners[vj] = vi
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

function bond_to_edge!(edges::AbstractVector{<:Edge}, vs_i::AbstractVector{<:Integer}, vs_j::AbstractVector{<:Integer})
    ni, nj = length(vs_i), length(vs_j)
    shared_symmetry = max(ni, nj) // min(ni, nj)

    if (ni == nj == 1) || !isinteger(shared_symmetry)
        edge = eltype(edges)(first(vs_i), first(vs_j))
        push!(edges, edge)
        return
    end
   
    if ni <= nj
        vs1, vs2 = vs_i, vs_j
    else
        vs1, vs2 = vs_j, vs_i
    end
    
    for (v1, v2) in zip(vs1, reverse(vs2[1:numerator(shared_symmetry):end]))
        edge = eltype(edges)(v1, v2)
        push!(edges, edge)
    end
    return 
end

function grow!(p::Polyform{D,T,F}, k::Integer, assembly_system::AssemblySystem) where {D,T,F}
    v, partner_label = open_bond(p, assembly_system, k)

    if isnothing(v)
        return false, k
    end

    success = attach_monomer!(p, v, partner_label, assembly_system)
    !success && return grow!(p, k+1, assembly_system)

    canonize!(p)
    return true, k
end
    

function shrink!(p::Polyform)
    n = size(p)
    n == 0 && return false

    # TODO: optimize
    k = 0
    if n > 1
        while true
            idx = p.encoder.bwd[end-k][1]
            idxs = [i for part in p.encoder.fwd[idx] for i in part]
            if !are_bridge(p.anatomy, idxs)
                break
            end
            k += 1
        end
    end

    idx = p.encoder.bwd[end-k][1]
    del_vs = sort([i for part in p.encoder.fwd[idx] for i in part])

    deleteat!(p.xs, idx)
    deleteat!(p.ψs, idx)
    deleteat!(p.species, idx)

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
    deleteat!(p.encoder, idx)

    canonize!(p)
    return true
end