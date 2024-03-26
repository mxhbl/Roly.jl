using LinearAlgebra
using StaticArrays
using Graphs
using NautyGraphs

mutable struct EncTranslator{T}
    fwd::Vector{Vector{T}}
    bwd::Matrix{T}
    n::T
    m::T

    function EncTranslator(fwd::AbstractVector{<:AbstractVector{T}}) where {T}
        n = length(fwd)
        m = sum(length(f) for f in fwd; init=0)

        bwd = zeros(2, m)
        for i in 1:m
            for (j, f) in enumerate(fwd)
                if i ∈ f
                    bwd[:, i] .= [T(j), T(findfirst(x->x==i, f))]
                    break
                end
            end
        end

        return new{T}(fwd, bwd, n, m)
    end
    function EncTranslator(fwd::AbstractVector{<:AbstractVector{T}}, bwd::AbstractMatrix{T}) where {T}
        return new{T}(fwd, bwd, length(fwd), size(bwd, 2))
    end
end
EncTranslator{T}() where {T} = EncTranslator(Vector{T}[], zeros(T, (0, 0)))

function Base.permute!(t::EncTranslator, p::AbstractVector{<:Integer})
    @assert length(p) == t.m

    t.bwd .= @view t.bwd[:, p]

    inv_p = invperm(p)
    t.fwd .= [inv_p[f] for f in t.fwd]
    return
end


function Base.vcat(t::EncTranslator{T}, h::EncTranslator) where {T}
    fwd = vcat(t.fwd, [hf .+ t.m for hf in h.fwd]) #TODO: reduce allocation
    bwd = hcat(t.bwd, h.bwd .+ t.n*[T(1), T(0)])
    return EncTranslator(fwd, bwd)
end
function Base.deleteat!(t::EncTranslator, i::Integer)
    del_ai = t.fwd[i]

    ai_shifter(ai) = ai - sum(x->x < ai, del_ai) #TODO:check gt or gte
    k_shifter(k) = k - (k > i)

    deleteat!(t.fwd, i)
    for j in eachindex(t.fwd)
        t.fwd[j] .= ai_shifter.(t.fwd[j])
    end
    t.bwd = t.bwd[:, setdiff(1:end, del_ai)]
    t.bwd[1, :] .= k_shifter.(t.bwd[1, :])

    t.n -= 1
    t.m -= length(del_ai)
    return
end
Base.copy(t::EncTranslator) = EncTranslator([copy(f) for f in t.fwd], copy(t.bwd))
function Base.copy!(dest::EncTranslator{T}, src::EncTranslator{T}) where {T}
    # TODO: make this non-allocating
    # copy!(dest.fwd, src.fwd)
    resize!(dest.fwd, src.n)
    dest.fwd .= copy.(src.fwd)
    dest.bwd = copy(src.bwd)
    dest.n = src.n
    dest.m = src.m
    return dest
end

mutable struct Structure{T<:Integer,F<:AbstractFloat}
    anatomy::DirectedDenseNautyGraph{Cint}
    translator::EncTranslator{Int}
    species::Vector{T}
    positions::Coordinates2D{F}
    σ::T

    function Structure{T,F}(anatomy, translator, species, positions) where {T,F}
        s = new{T,F}(anatomy, translator, species, positions, 0)
        canonize!(s)
        return s
    end

    function Structure{T,F}(anatomy, translator, species, positions, σ) where {T,F}
        return new{T,F}(anatomy, translator, species, positions, σ)
    end
end
Structure(args...) = Structure{DefInt,DefFloat}(args...)
Structure{T,F}() where {T,F} = Structure{T,F}(DirectedDenseNautyGraph(0), EncTranslator{Int}(), zeros(T, 0), Coordinates2D{F}([], []), 0)

function canonize!(s::Structure)
    canon_perm, n = NautyGraphs.canonize!(s.anatomy)
    permute!(s.translator, canon_perm)
    s.σ = n
    return
end

faces(::Type{T}, s::Structure) where {T} = convert(Vector{T}, s.anatomy.labels)
faces(s::Structure) = faces(Int, s)
species(s::Structure) = s.species

function Base.size(s::Structure)
    e = ne(s.anatomy)
    v = nv(s.anatomy)
    return length(s.positions), (e - v) ÷ 2, v
end #  TODO: nbonds only valid in 2D
function Base.show(io::IO, s::Structure{T,F}) where {T,F}
    let (n, m, _) = size(s)
        println(io, "Structure{$T, $F}[n=$n, k=$m]")
    end
end

function Base.hash(s::Structure)
    return hash(s.anatomy)
end
function Base.hash(s::Structure, h::UInt64)
    return hash(s.anatomy, h)
end
function Base.:(==)(si::Structure, sj::Structure)
    return hash(si.anatomy) == hash(sj.anatomy)
end

function Base.copy(s::Structure{T,F}) where {T,F}
    return Structure{T,F}(copy(s.anatomy), copy(s.translator), copy(s.species), copy(s.positions), s.σ)
end
function Base.copy!(dest::Structure{T,F}, src::Structure{T,F}) where {T,F}
    copy!(dest.anatomy, src.anatomy)
    copy!(dest.translator, src.translator)
    copy!(dest.species, src.species)
    copy!(dest.positions, src.positions)
    dest.σ = src.σ
    return dest
end

function create_monomer(geometry::AbstractGeometry{T,F},
                        species_idx::T,
                        face_labels::AbstractVector{<:Integer}) where {T<:Integer,F<:Real}

    return Structure{T,F}(DirectedDenseNautyGraph(geometry.anatomy, face_labels),
                          EncTranslator([collect(1:length(geometry))]),
                          [species_idx],
                          Coordinates2D{F}([SA[0.0, 0.0]], Vector([0.0])))
end

struct AssemblySystem{T<:Integer, F<:AbstractFloat}
    interactions::InteractionDiagram{T}
    interaction_matrix::BitMatrix
    monomers::Vector{Structure{T,F}}
    geometries::Vector{PolygonGeometry{T,F}}
    face_translator::Vector{Tuple{T,T}}
end

function AssemblySystem(int_diag::InteractionDiagram, geoms::Vector{<:AbstractGeometry{T,F}}, face_labels=nothing) where {T,F}
    int_diag = convert(InteractionDiagram{T}, int_diag)
    n_species, _ = size(int_diag)
    n_sides = sum(length(geom) for geom in geoms; init=0)

    monomers = Structure{T,F}[]
    face_translator = Tuple{T,T}[]
    last_label = 0

    for i in 1:n_species
        if isnothing(face_labels)
            face_labels = collect(1:length(geoms[i]))
        end

        fl = convert(Vector{T}, face_labels .+ last_label)
        m = create_monomer(geoms[i], T(i), fl)
        push!(monomers, m)

        for j in 1:length(face_labels)
            push!(face_translator, (i, j))
        end

        last_label = fl[end]
    end

    interaction_matrix = falses(n_sides, n_sides)
    for i in 1:n_sides, j in i:n_sides
        if InteractionEdge{T}(face_translator[i], face_translator[j], 0) ∈ int_diag
            interaction_matrix[i, j] = 1
            interaction_matrix[j, i] = 1
        end
    end

    return AssemblySystem{T,F}(int_diag, interaction_matrix, monomers, geoms, face_translator)
end
function AssemblySystem(int_diag::InteractionDiagram, geom::AbstractGeometry{T,F}, face_labels=nothing) where {T,F}
    n_species, _ = size(int_diag)
    geoms = [geom for _ in 1:n_species]
    return AssemblySystem(int_diag, geoms, face_labels)
end
function AssemblySystem(A::Matrix, geom::Union{AbstractGeometry{T,F}, Vector{<:AbstractGeometry{T,F}}}, face_labels=nothing) where {T,F}
    A = convert(Matrix{T}, A)
    int_diag = InteractionDiagram(A)
    return AssemblySystem(int_diag, geom, face_labels)
end


function concatenate_structures(si::Structure{T,F}, sj::Structure{T,F}, loc_i::SideLoc{T},
    loc_j::SideLoc{T}, assembly_system::AssemblySystem) where {T,F}

    ni, _, nfi = convert(Tuple{T,T,T}, size(si))
    geometries = assembly_system.geometries
    interaction_matrix = assembly_system.interaction_matrix

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
    species_j, face_j = assembly_system.face_translator[j]

    spec = species(s)
    sp_i = species(s)[i]
    geom_i = geometries[sp_i]
    geom_j = geometries[species_j]

    Δx, Δψ = attachment_offset(T(face_i), T(face_j), geom_i, geom_j)
    ψj = Δψ + s.positions.ψs[i]
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

            if !assembly_system.interaction_matrix[li, lj]
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
    bond = open_bond(s, k, assembly_system.interaction_matrix)
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

function face_orientation(s::Structure, i::Integer, assembly_system::AssemblySystem)
    si, fi = s.translator.bwd[:, i]
    spec = s.species[si]
    ψ = Rational(s.positions.ψs[si])
    θ = Rational(assembly_system.geometries[spec].θs[fi])
    ϕ = Rational(assembly_system.geometries[spec].ϕs[fi])

    α = (θ + ϕ + ψ) % 2
    return α
end

function analyze_cycle(s::Structure, i::Integer, j::Integer)
    path = Graphs.a_star(s.anatomy, i, j)
    if length(path) == 0
        return 1 // 1
    end

    angles = [1//2, 0//1, -1//2]

    total_angle = 0 // 1
    angle = 0 // 1
    counter = 0
    for edge in path
        if !has_edge(s.anatomy, reverse(edge))
            counter += 1
            # i = s.species[s.translator.bwd[1, edge.src]]
            # n_sides = length(assembly_system.geometries[i])
            # angle += 2 // n_sides
        else
            angle = angles[counter]
            # if angle > 1 // 1
            #     angle -= 2 // 1
            # end
            total_angle += angle
            counter = 0
            angle = 0 // 1
        end
    end

    # if angle > 1 // 1
    #     angle -= 2 // 1
    # end
    angle = angles[counter]

    total_angle += angle

    return total_angle
end

function is_cyclic(s::Structure{T,F}, assembly_system::AssemblySystem) where {T,F}
    #WORKS ONLY IN 2D
    g = s.anatomy
    A = assembly_system.interaction_matrix

    open_sites = filter(v->outdegree(g, v)==1, vertices(g))
    filter!(s->sum(@view A[:, g.labels[s]]) > 0, open_sites)

    for (i, oi) in enumerate(open_sites), oj in @view open_sites[i:end]
        ti, tj = g.labels[oi], g.labels[oj]

        if A[ti, tj] > 0
            αi = face_orientation(s, oi, assembly_system)
            αj = face_orientation(s, oj, assembly_system)

            could_be_cyclic = abs(αi - αj) == 1
            if could_be_cyclic
                loc_i = Tuple{T,T}(s.translator.bwd[:, oi])
                loc_j = Tuple{T,T}(s.translator.bwd[:, oj])
                s_concat = concatenate_structures(s, s, loc_i, loc_j, assembly_system)
                !isnothing(s_concat) && return true
            end
        end
    end
    return false
end