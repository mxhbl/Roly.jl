struct AssemblySystem{D, T<:Integer, F<:AbstractFloat, G<:AbstractGeometry{T,F}}
    intmat::BitMatrix
    buildingblocks::Vector{Polyform{D,T,F}}
    geometries::Vector{G}
    n_species::Integer
    n_edges::Integer
    _sites_sum::Vector{T}
end
interaction_matrix(sys::AssemblySystem) = sys.intmat
intmat(sys::AssemblySystem) = interaction_matrix(sys)
buildingblocks(sys::AssemblySystem) = sys.buildingblocks
geometries(sys::AssemblySystem) = sys.geometries
Base.size(sys::AssemblySystem) = sys.n_species, sys.n_edges
Base.show(io::Core.IO, A::AssemblySystem{D,T,F}) where {D,T,F} = print(io, "AssemblySytem{$D,$T,$F}[n=$(A.n_species), k=$(A.n_edges)]")

function AssemblySystem(interactions::AbstractMatrix, geometries::Vector{<:AbstractGeometry{T,F}}, face_labels=nothing) where {T,F}
    ds = [dimension(g) for g in geometries]
    D = first(ds)
    @assert all(ds .== D)

    n_species = length(geometries)
    sites = [nsites(geom) for geom in geometries]

    monomers = Polyform{D,T,F}[]
    sites_sum = cumsum(sites)
    last_label = 0

    for i in 1:n_species
        if isnothing(face_labels)
            face_labels = collect(1:sites[i])
        end

        fl = convert(Vector{T}, face_labels .+ last_label)
        m = create_monomer(geometries[i], T(i), fl)
        push!(monomers, m)

        last_label = fl[end]
    end

    interaction_matrix = instantiate_interactionmatrix(interactions, sites)
    n_edges = sum(interaction_matrix) ÷ 2
    return AssemblySystem{D,T,F,eltype(geometries)}(interaction_matrix, monomers, geometries, n_species, n_edges, sites_sum)
end
function AssemblySystem(interactions::AbstractMatrix{<:Integer}, geometry::AbstractGeometry{T,F}, face_labels=nothing) where {T,F}
    n_species = maximum(interactions[:, [1, 3]])
    geometries = [geometry for _ in 1:n_species]
    return AssemblySystem(interactions, geometries, face_labels)
end

function spcs_site_to_siteidx(spcs::Integer, site::Integer, assembly_system::AssemblySystem)
    return irg_flatten(spcs, site, assembly_system._sites_sum)
end
function siteidx_to_spcs_site(site_index::Integer, assembly_system::AssemblySystem)
    return irg_unflatten(site_index, assembly_system._sites_sum)
end

function instantiate_interactionmatrix(interactions::AbstractMatrix{<:Integer}, sites)
    n_sites = sum(sites; init=0)
    sites_sum = cumsum(sites)

    interaction_matrix = falses(n_sites, n_sites)
    for edge in eachrow(interactions)
        a, b, c, d = edge
        i, j = irg_flatten(a, b, sites_sum), irg_flatten(c, d, sites_sum)
        interaction_matrix[i, j] = interaction_matrix[j, i] = true
    end
    return interaction_matrix
end
function instantiate_interactionmatrix(interactions::BitMatrix, args...)
    @assert interactions == interactions'
    interaction_matrix = BitMatrix(interactions)
    return interaction_matrix
end

function anatomy(asys::AssemblySystem)
    imat = interaction_matrix(asys)
    n_sites = size(imat, 1)
    n_edges = sum(imat) ÷ 2

    A = falses(n_sites + n_edges, n_sites + n_edges)
    edge_counter = 1
    for i in axes(imat, 1), j in i+1:size(imat, 2)
        if imat[i, j] == 0
            continue
        end

        A[i, n_sites+edge_counter] = true
        A[j, n_sites+edge_counter] = true
        A[n_sites+edge_counter, i] = true
        A[n_sites+edge_counter, j] = true
        edge_counter += 1
    end
    
    i = 1
    for geom in asys.geometries
        a = adjacency_matrix(geom.anatomy)
        n = size(a, 1)
        A[i:i+n-1, i:i+n-1] .+= a
        i += n
    end

    vertex_colors = [zeros(Cint,n_sites); ones(Cint, n_edges)]

    anatomy = NautyDiGraph(A, vertex_colors)
    return anatomy
end

function Base.hash(asys::AssemblySystem)
    a = anatomy(asys)
    a_prime = NautyDiGraph(adjacency_matrix(a)', a.labels)
    return hash(sort([hash(a), hash(a_prime)]))
end