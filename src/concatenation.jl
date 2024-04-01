function open_bond(p::Polyform, j::Integer, interaction_matrix::AbstractMatrix)
    n = nv(p.anatomy)
    fs = faces(p)
    l = j
    oneighs = zeros(Int, n)
    for ai in 1:n
        nneighs = NautyGraphs.outneighbors!(oneighs, p.anatomy, ai)
        bound = false
        for i_neigh in 1:nneighs
            if has_edge(p.anatomy, oneighs[i_neigh], ai)
                bound = true
                break
            end
        end
        if bound
            continue
        end

        fi = fs[ai]
        for (k, fk) in enumerate(@view interaction_matrix[:, fi])
            if fk > 0
                l -= 1
            end
            if l == 0
                return (ai, k)
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
    n, nk, nf = size(p)
    geoms = geometries(assembly_system)
    bblocks = buildingblocks(assembly_system)

    ai, j = bond
    i, face_i = @view p.translator.bwd[:, ai]
    species_j, face_j = irg_unflatten(j, assembly_system._sides_sum)

    spec = species(p)
    sp_i = species(p)[i]
    geom_i = geoms[sp_i]
    geom_j = geoms[species_j]

    Δx, Δψ = attachment_offset(T(face_i), T(face_j), geom_i, geom_j)
    ψj = Δψ * p.ψs[i]
    xj = p.xs[i] + rotate(Δx, p.ψs[i])

    anatomy_edges = Edge{Cint}[]
    for (i, (xi, ψi)) in enumerate(zip(p.xs, p.ψs))
        cstat, pairs = contact_status(xj - xi, ψi, ψj, geoms[spec[i]], geom_j)
        # If the new structure is invalid, call grow! again with j->j+1
        if !cstat
            # println("overlap")
            return false # return the new index as well
        end
        for pair in pairs
            fi, fj = pair

            ai = p.translator.fwd[i][fi]
            aj = fj
            li = p.anatomy.labels[ai]
            lj = bblocks[species_j].anatomy.labels[aj]

            if !interaction_matrix(assembly_system)[li, lj]
                # println("blocked bond")
                return false
            end
            
            anatomy_edge = Edge(Cint(ai), Cint(aj+nf))
            push!(anatomy_edges, anatomy_edge)
            push!(anatomy_edges, reverse(anatomy_edge))
        end
    end

    NautyGraphs.blockdiag!(p.anatomy, bblocks[species_j].anatomy)
    for ae in anatomy_edges
        add_edge!(p.anatomy, ae)
    end

    push!(p.xs, xj)
    push!(p.ψs, ψj)

    push!(p.species, species_j)
    p.translator = vcat(p.translator, bblocks[species_j].translator)
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
            idx = p.translator.bwd[1, end-k]
            idxs = p.translator.fwd[idx]
            if !are_bridge(p.anatomy, idxs)
                break
            end
            k += 1
        end
    end

    idx = p.translator.bwd[1, end-k]

    deleteat!(p.xs, idx)
    deleteat!(p.ψs, idx)
    deleteat!(p.species, idx)
    rem_vertices!(p.anatomy, p.translator.fwd[idx])

    deleteat!(p.translator, idx)

    canonize!(p)
    return true
end