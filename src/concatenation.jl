function concatenate_structures(si::Structure{T,F}, sj::Structure{T,F}, loc_i::SideLoc{T},
    loc_j::SideLoc{T}, assembly_system::AssemblySystem) where {T,F}

    ni, _, nfi = convert(Tuple{T,T,T}, size(si))
    geometries = assembly_system.geometries
    interaction_matrix = assembly_system.intmat

    i, face_i = loc_i
    j, face_j = loc_j

    species_i = si.species
    species_j = sj.species

    Δx, Δψ = attachment_offset(face_i, face_j, geometries[species_i[i]], geometries[species_j[j]])
    qi = si.positions

    qj = deepcopy(sj.positions)
    grab_at!(qj, j)
    rotate!(qj, qi.ψs[i] + Δψ)
    shift!(qj, qi.xs[i] + rotate(Δx, qi.ψs[i]))

    anatomy = blockdiag(si.anatomy, sj.anatomy)

    for (i, (xi, ψi)) in enumerate(qi)
        for (j, (xj, ψj)) in enumerate(qj)
            cstat, pairs = contact_status(xj - xi, ψi, ψj, geometries[species_i[i]], geometries[species_j[j]])
            if !cstat
                return nothing
            end
            for pair in pairs
                fi, fj = pair

                ai = si.translator.fwd[i][fi]
                aj = sj.translator.fwd[j][fj]

                li = si.anatomy.labels[ai]
                lj = sj.anatomy.labels[aj]

                if !interaction_matrix[li, lj]
                    return nothing
                end
                
                add_edge!(anatomy, Edge(ai, aj+nfi))
                add_edge!(anatomy, Edge(aj+nfi, ai)) # Reverse
            end
        end
    end

    positions = vcat(qi, qj)
    species = vcat(si.species, sj.species)
    translator = vcat(si.translator, sj.translator)

    return Structure{T,F}(anatomy, translator, species, positions)
end

function open_bond(s::Structure, j::Integer, interaction_matrix::AbstractMatrix)
    n = nv(s.anatomy)
    fs = faces(s)
    l = j
    oneighs = zeros(Int, n)
    for ai in 1:n
        nneighs = NautyGraphs.outneighbors!(oneighs, s.anatomy, ai)
        bound = false
        for i_neigh in 1:nneighs
            if has_edge(s.anatomy, oneighs[i_neigh], ai)
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

function all_open_bonds(s::Structure, interaction_matrix::AbstractMatrix)
    bonds = Tuple{Int,Int}[]
    j = 1
    b = open_bond(s, j, interaction_matrix)
    while !isnothing(b)
        push!(bonds, b)
        j += 1
        b = open_bond(s, j, interaction_matrix)
    end
    return bonds
end

function attach_monomer!(s::Structure{T,F}, bond::Tuple{<:Integer,<:Integer}, assembly_system::AssemblySystem, fillhash::Bool=false) where {T,F}
    n, nk, nf = size(s)
    geometries = assembly_system.geometries

    ai, j = bond
    i, face_i = @view s.translator.bwd[:, ai]
    species_j, face_j = irg_unflatten(j, assembly_system._sides_sum)

    spec = species(s)
    sp_i = species(s)[i]
    geom_i = geometries[sp_i]
    geom_j = geometries[species_j]

    Δx, Δψ = attachment_offset(T(face_i), T(face_j), geom_i, geom_j)
    ψj = Δψ + s.positions.ψs[i]
    # ψj = quaternion_multiply(Δψ, s.positions.ψs[i])

    xj = s.positions.xs[i] + rotate(Δx, s.positions.ψs[i])

    anatomy_edges = Edge{Cint}[]
    for (i, (xi, ψi)) in enumerate(s.positions)
        cstat, pairs = contact_status(xj - xi, ψi, ψj, geometries[spec[i]], geom_j)
        # If the new structure is invalid, call grow! again with j->j+1
        if !cstat
            # println("overlap")
            return false # return the new index as well
        end
        for pair in pairs
            fi, fj = pair

            ai = s.translator.fwd[i][fi]
            aj = fj
            li = s.anatomy.labels[ai]
            lj = assembly_system.monomers[species_j].anatomy.labels[aj]

            if !assembly_system.intmat[li, lj]
                # println("blocked bond")
                return false
            end
            
            anatomy_edge = Edge(Cint(ai), Cint(aj+nf))
            push!(anatomy_edges, anatomy_edge)
            push!(anatomy_edges, reverse(anatomy_edge))
        end
    end

    NautyGraphs.blockdiag!(s.anatomy, assembly_system.monomers[species_j].anatomy)
    for ae in anatomy_edges
        add_edge!(s.anatomy, ae)
    end

    push!(s.positions.xs, xj)
    push!(s.positions.ψs, ψj)

    push!(s.species, species_j)
    s.translator = vcat(s.translator, assembly_system.monomers[species_j].translator)
    #TODO MAKE THIS PRETTY, THIS WILL CANONIZE A STRUCTURE TWICE IF GROW! is CALLED
    if fillhash
        _, n = NautyGraphs._fill_hash!(s.anatomy)
        s.σ = n
    else
        s.σ = 0
    end
    return true
end

function grow!(s::Structure{T,F}, k::Integer, assembly_system::AssemblySystem) where {T,F}
    bond = open_bond(s, k, assembly_system.intmat)
    if isnothing(bond)
        return false, k
    end

    success = attach_monomer!(s, bond, assembly_system)
    if !success
        return grow!(s, k+1, assembly_system)
    end

    canonize!(s)
    return true, k
end
    

function shrink!(s::Structure)
    n, _, _ = size(s)
    if n == 0
        return false
    end

    # TODO: optimize
    k = 0
    if n > 1
        while true
            idx = s.translator.bwd[1, end-k]
            idxs = s.translator.fwd[idx]
            if !are_bridge(s.anatomy, idxs)
                break
            end
            k += 1
        end
    end

    idx = s.translator.bwd[1, end-k]

    deleteat!(s.positions.xs, idx)
    deleteat!(s.positions.ψs, idx)
    deleteat!(s.species, idx)
    rem_vertices!(s.anatomy, s.translator.fwd[idx])

    deleteat!(s.translator, idx)

    canonize!(s)
    return true
end