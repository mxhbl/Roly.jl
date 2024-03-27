function concatenate_structures(poly_i::Polyform{T,F}, poly_j::Polyform{T,F}, loc_i::SideLoc{T},
    loc_j::SideLoc{T}, assembly_system::AssemblySystem) where {T,F}

    ni, _, nfi = convert(Tuple{T,T,T}, size(poly_i))
    geometries = assembly_system.geometries
    interaction_matrix = assembly_system.intmat

    i, face_i = loc_i
    j, face_j = loc_j

    species_i = poly_i.species
    species_j = poly_j.species

    Δx, Δψ = attachment_offset(face_i, face_j, geometries[species_i[i]], geometries[species_j[j]])
    qi = poly_i.positions

    qj = deepcopy(poly_j.positions)
    grab_at!(qj, j)
    rotate!(qj, qi.ψs[i] + Δψ)
    shift!(qj, qi.xs[i] + rotate(Δx, qi.ψs[i]))

    anatomy = blockdiag(poly_i.anatomy, poly_j.anatomy)

    for (i, (xi, ψi)) in enumerate(qi)
        for (j, (xj, ψj)) in enumerate(qj)
            cstat, pairs = contact_status(xj - xi, ψi, ψj, geometries[species_i[i]], geometries[species_j[j]])
            if !cstat
                return nothing
            end
            for pair in pairs
                fi, fj = pair

                ai = poly_i.translator.fwd[i][fi]
                aj = poly_j.translator.fwd[j][fj]

                li = poly_i.anatomy.labels[ai]
                lj = poly_j.anatomy.labels[aj]

                if !interaction_matrix[li, lj]
                    return nothing
                end
                
                add_edge!(anatomy, Edge(ai, aj+nfi))
                add_edge!(anatomy, Edge(aj+nfi, ai)) # Reverse
            end
        end
    end

    positions = vcat(qi, qj)
    species = vcat(poly_i.species, poly_j.species)
    translator = vcat(poly_i.translator, poly_j.translator)

    return Polyform{T,F}(anatomy, translator, species, positions)
end

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

function attach_monomer!(p::Polyform{T,F}, bond::Tuple{<:Integer,<:Integer}, assembly_system::AssemblySystem, fillhash::Bool=false) where {T,F}
    n, nk, nf = size(p)
    geometries = assembly_system.geometries

    ai, j = bond
    i, face_i = @view p.translator.bwd[:, ai]
    species_j, face_j = irg_unflatten(j, assembly_system._sides_sum)

    spec = species(p)
    sp_i = species(p)[i]
    geom_i = geometries[sp_i]
    geom_j = geometries[species_j]

    Δx, Δψ = attachment_offset(T(face_i), T(face_j), geom_i, geom_j)
    ψj = Δψ + p.positions.ψs[i]
    # ψj = quaternion_multiply(Δψ, s.positions.ψs[i])

    xj = p.positions.xs[i] + rotate(Δx, p.positions.ψs[i])

    anatomy_edges = Edge{Cint}[]
    for (i, (xi, ψi)) in enumerate(p.positions)
        cstat, pairs = contact_status(xj - xi, ψi, ψj, geometries[spec[i]], geom_j)
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

    NautyGraphs.blockdiag!(p.anatomy, assembly_system.monomers[species_j].anatomy)
    for ae in anatomy_edges
        add_edge!(p.anatomy, ae)
    end

    push!(p.positions.xs, xj)
    push!(p.positions.ψs, ψj)

    push!(p.species, species_j)
    p.translator = vcat(p.translator, assembly_system.monomers[species_j].translator)
    #TODO MAKE THIS PRETTY, THIS WILL CANONIZE A STRUCTURE TWICE IF GROW! is CALLED
    if fillhash
        _, n = NautyGraphs._fill_hash!(p.anatomy)
        p.σ = n
    else
        p.σ = 0
    end
    return true
end

function grow!(p::Polyform{T,F}, k::Integer, assembly_system::AssemblySystem) where {T,F}
    bond = open_bond(p, k, assembly_system.intmat)
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

    deleteat!(p.positions.xs, idx)
    deleteat!(p.positions.ψs, idx)
    deleteat!(p.species, idx)
    rem_vertices!(p.anatomy, p.translator.fwd[idx])

    deleteat!(p.translator, idx)

    canonize!(p)
    return true
end