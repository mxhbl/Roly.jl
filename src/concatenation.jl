function open_bond(p::Polyform, j::Integer, interaction_matrix::AbstractMatrix)
    fs = faces(p)
    for v in 1:length(p.bond_partners)
        if p.bond_partners[v] > 0
            continue
        end

        fi = fs[v]
        for (k, fk) in enumerate(@view interaction_matrix[:, fi])
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

function all_open_bonds(p::Polyform, interaction_matrix::AbstractMatrix)
    bonds = Tuple{Int,Int}[]
    j = 1
    b = open_bond(p, j, interaction_matrix)
    while !isnothing(b)
        push!(bonds, b)
        j += 1
        b = open_bond(p, j, interaction_matrix)
    end
    return bonds
end

function attach_monomer!(p::Polyform{D,T,F}, bond::Tuple{<:Integer,<:Integer}, assembly_system::AssemblySystem, fillhash::Bool=false) where {D,T,F}
    n, _, nf = size(p)
    geoms = geometries(assembly_system)

    v, bblock_side_label = bond
    i, face_i, _ = @view p.encoder.bwd[:, v]
    spcs_j, face_j = irg_unflatten(bblock_side_label, assembly_system._sides_sum)
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

            aidx_i = first(p.encoder.fwd[i][face_i])
            aidx_j = face_j
            lbl_i = p.anatomy.labels[aidx_i]
            lbl_j = bblock.anatomy.labels[aidx_j]

            if !interaction_matrix(assembly_system)[lbl_i, lbl_j]
                # println("blocked bond")
                return false
            end
            
            anatomy_edge = Edge(Cint(aidx_i), Cint(aidx_j+nf))
            push!(anatomy_edges, anatomy_edge)
            push!(anatomy_edges, reverse(anatomy_edge))
        end
    end

    NautyGraphs.blockdiag!(p.anatomy, bblock.anatomy)
    for ae in anatomy_edges
        add_edge!(p.anatomy, ae)
    end

    push!(p.xs, xj)
    push!(p.ψs, ψj)

    push!(p.species, spcs_j)
    p.encoder = concatenate(p.encoder, bblock.encoder)
    append!(p.bond_partners, zeros(T, bblock.encoder.n_vertices)) #TODO: make n_faces a higher level interface
    for e in anatomy_edges
        vi, vj = e.src, e.dst
        p.bond_partners[vi] = vj
        p.bond_partners[vj] = vi
    end

    #TODO MAKE THIS PRETTY, THIS WILL CANONIZE A STRUCTURE TWICE IF GROW! is CALLED
    if fillhash
        _, n = NautyGraphs._fill_hash!(p.anatomy)
        p.σ = n
    else
        p.σ = 0
    end
    return true
end

function grow!(p::Polyform{D,T,F}, k::Integer, assembly_system::AssemblySystem) where {D,T,F}
    bond = open_bond(p, k, interaction_matrix(assembly_system))
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
    n, _, _ = size(p)
    if n == 0
        return false
    end

    # TODO: optimize
    k = 0
    if n > 1
        while true
            idx = p.encoder.bwd[1, end-k]
            idxs = [i for part in p.encoder.fwd[idx] for i in part]
            if !are_bridge(p.anatomy, idxs)
                break
            end
            k += 1
        end
    end

    idx = p.encoder.bwd[1, end-k]
    del_vs = sort([i for part in p.encoder.fwd[idx] for i in part])

    deleteat!(p.xs, idx)
    deleteat!(p.ψs, idx)
    deleteat!(p.species, idx)

    # open up bonds that are now free
    # TODO: this should only be relevant when we shrink a monomer, could be handled better
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