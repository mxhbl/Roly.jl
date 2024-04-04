function open_bond(p::Polyform, j::Integer, assembly_sytem::AssemblySystem)
    for v in 1:length(p.bond_partners)
        if p.bond_partners[v] > 0
            continue
        end
        #TODO :FIX THIS INTERFACE AND DECIDE WHAT LABEL SHOULD DO
        spcs = species(p)[p.encoder.bwd[v][1]]
        offset = spcs > 1 ? assembly_sytem._sites_sum[spcs-1] : 0
        lbl = p.encoder.bwd[v][2] + offset
        for (k, fk) in enumerate(@view assembly_sytem.intmat[:, lbl])
            if fk > 0
                j -= 1
            end
            if j == 0
                return (v, k)
            end
        end
    end

    return nothing
end

function all_open_bonds(p::Polyform, assembly_sytem::AssemblySystem)
    bonds = Tuple{Int,Int}[]
    j = 1
    b = open_bond(p, j, assembly_sytem)
    while !isnothing(b)
        push!(bonds, b)
        j += 1
        b = open_bond(p, j, assembly_sytem)
    end
    return bonds
end

function attach_monomer!(p::Polyform{D,T,F}, bond::Tuple{<:Integer,<:Integer}, assembly_system::AssemblySystem, fillhash::Bool=false) where {D,T,F}
    n = size(p)
    nvert = nvertices(p)
    geoms = geometries(assembly_system)

    v, bblock_side_label = bond
    i, face_i, _ = p.encoder.bwd[v]
    spcs_j, face_j = irg_unflatten(bblock_side_label, assembly_system._sites_sum)
    bblock = buildingblocks(assembly_system)[spcs_j]

    spcs = species(p)
    spcs_i = spcs[i]
    geom_i = geoms[spcs_i]
    geom_j = geoms[spcs_j]

    Δx, Δψ = attachment_offset(T(face_i), T(face_j), geom_i, geom_j)
    ψj = Δψ * p.ψs[i]
    xj = p.xs[i] + rotate(Δx, p.ψs[i])

    anatomy_edges = Edge{Cint}[]
    for (i, (xi, ψi)) in enumerate(zip(p.xs, p.ψs))
        cstat, pairs = contact_status(xj - xi, ψi, ψj, geoms[spcs[i]], geom_j)
        if !cstat
            # println("overlap")
            return false
        end
        for pair in pairs
            face_i, face_j = pair

            vs_i = p.encoder.fwd[i][face_i]
            vs_j = bblock.encoder.fwd[1][face_j]
            lbl_i = p.anatomy.labels[first(vs_i)]
            lbl_j = bblock.anatomy.labels[first(vs_j)]

            if !interaction_matrix(assembly_system)[lbl_i, lbl_j]
                # println("blocked bond")
                return false
            end
        
            aedges = contract_faces(convert(Vector{Cint}, vs_i), vs_j .+ nvert)
            append!(anatomy_edges, aedges)
        end
    end

    NautyGraphs.blockdiag!(p.anatomy, bblock.anatomy)

    push!(p.xs, xj)
    push!(p.ψs, ψj)

    push!(p.species, spcs_j)
    p.encoder = concatenate(p.encoder, bblock.encoder)
    append!(p.bond_partners, zeros(T, nvertices(bblock.encoder))) #TODO: make n_faces a higher level interface
    for ae in anatomy_edges
        add_edge!(p.anatomy, ae)
        add_edge!(p.anatomy, reverse(ae))

        vi, vj = ae.src, ae.dst
        p.bond_partners[vi] = vj
        p.bond_partners[vj] = vi
    end

    #TODO MAKE THIS PRETTY
    if fillhash
        _, n = NautyGraphs._fill_hash!(p.anatomy)
        p.σ = n
    else
        p.σ = 0
    end
    return true
end

function contract_faces(vs_i::AbstractVector{T}, vs_j::AbstractVector{T}) where {T}
    #TODO: handle non-identical compatible and incompatible symmetry
    # symmetry_order = gcd(length(vs_i), length(vs_j))
    if length(vs_i) == 1 && length(vs_j) == 1
        return [Edge{T}(only(vs_i), only(vs_j))]
    else
        #TODO optimize
        idxs_j = [1; length(vs_j):-1:2]
    end
    return [Edge{T}(vi, vj) for (vi, vj) in zip(vs_i, vs_j[idxs_j])]
end

function grow!(p::Polyform{D,T,F}, k::Integer, assembly_system::AssemblySystem) where {D,T,F}
    bond = open_bond(p, k, assembly_system)
    if isnothing(bond)
        return false, k
    end

    success = attach_monomer!(p, bond, assembly_system)
    if !success
        return grow!(p, k+1, assembly_system)
    end

    canonize!(p)
    return true, k
end
    

function shrink!(p::Polyform)
    n = size(p)
    if n == 0
        return false
    end

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