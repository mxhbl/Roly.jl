using LinearAlgebra
using StaticArrays
using Graphs
using Base.Iterators
using DataStructures
using Distributed
using Serialization
using Base.Threads

function f!(k::Structure{T,F}, s::Structure{T,F}) where {T,F}
    copy!(k, s)
    shrink!(k)
    return k
end

function adj!(u::Structure{T,F}, v::Structure{T,F}, j::Integer, hashes::Vector{HashType},
              assembly_system::AssemblySystem) where {T,F}
    ### TODO: We should exploit that canonical labels are guaranteed (see Nauty User Guide p.4) to be in order of color
    ### So, if we know the color of the particle that would be removed, we can filter the possible offspring
    if size(v)[1] == 0
        if j > length(assembly_system.monomers)
            return false, j + 1
        end

        copy!(u, assembly_system.monomers[j])
        return true, j + 1
    end

    copy!(u, v)
    success, j = grow!(u, j, assembly_system)

    if success && hash(u.anatomy) ∈ hashes
        return adj!(u, v, j + 1, hashes, assembly_system)
    end

    return success, j + 1
end

FnorNothing = Union{<:Function,Nothing}

function polyrs(v₀::Structure{T,F},
                assembly_system::AssemblySystem{T,F},
                reducer::FnorNothing,
                reduce_op::FnorNothing,
                aggregator::FnorNothing,
                rejector::FnorNothing,
                unexplored_channel::Union{<:AbstractChannel,<:RemoteChannel,Nothing}=nothing,
                output_channel::Union{<:AbstractChannel,<:RemoteChannel,Nothing}=nothing,
                max_depth::Int=0,
                max_strs::Int=0) where {T,F}
    reduction = !isnothing(reducer)
    if reduction
        reduce_val = reducer(v₀)
    else
        reduce_val = nothing
    end
    aggregation = !isnothing(aggregator)
    if aggregation
        aggregate_val = typeof(aggregator(v₀)[2])[]
    else
        aggregate_val = nothing
    end
    rejecting = !isnothing(rejector)
    return_unexplored = !isnothing(unexplored_channel)

    d0 = size(v₀)[1]
    depth = 0
    n_strs = 0

    v = Structure{T,F}()
    u = Structure{T,F}()
    k = Structure{T,F}()

    copy!(v, v₀)

    js = [1]
    hashes = [HashType[]]

    if max_depth == 0
        max_depth = typemax(max_depth)
    end
    if max_strs == 0
        max_strs = typemax(max_strs)
    end

    lowest_depth = 0
    n_unexplored = 0
    break_triggered = false

    while true
        unexplored = false
        success, j_new = adj!(u, v, js[end], hashes[end], assembly_system)
        next = u
        if success
            js[end] = j_new
            push!(hashes[end], hash(next.anatomy))

            f!(k, next)

            if k == v
                if rejecting
                    reject_val = rejector(next, assembly_system)
                    # 0: dont reject, 1: reject post, 2: reject pre, 3: break immediately
                    reject_val === 2 && continue
                end
                if reduction
                    reduce_val = reduce_op(reduce_val, reducer(next))
                end
                if aggregation
                    aggregate, agr_val = aggregator(next)
                    aggregate && push!(aggregate_val, agr_val)
                end

                n_strs += 1
                depth += 1
                if depth > lowest_depth
                    lowest_depth = depth
                end

                if !isnothing(output_channel)
                    out = lowest_depth + d0, (reduce_val, aggregate_val, break_triggered)
                    try
                        put!(output_channel, out)
                    catch
                        return
                    end
                end

                if rejecting
                    if reject_val === 1
                        depth -= 1
                        continue
                    elseif reject_val === 3
                        break_triggered = true
                        break
                    end
                end

                if !return_unexplored && n_strs >= max_strs
                    break_triggered = true
                    break
                end

                if depth == max_depth || n_strs >= max_strs
                    unexplored = true
                end

                if unexplored
                    depth -= 1
                    return_unexplored && put!(unexplored_channel, copy(next))
                    n_unexplored += 1
                    continue
                end

                copy!(v, next)
                push!(js, 1)
                push!(hashes, HashType[])
            end

            continue
        end

        if depth > 0
            f!(k, v)    # TODO speed this up by storing structures
            copy!(v, k)

            depth -= 1

            pop!(js)
            pop!(hashes)
        else
            break
        end
    end

    if isnothing(output_channel)
        return (n_strs, lowest_depth), (reduce_val, aggregate_val, break_triggered)
    end
end

function polyenum(assembly_system::AssemblySystem{T,F};
                  max_size::Number=Inf, max_strs::Number=Inf,
                  reducer::FnorNothing=nothing,
                  reduce_op::FnorNothing=nothing,
                  aggregator::FnorNothing=nothing,
                  rejector::FnorNothing=nothing) where {T,F}
    v₀ = Structure{T,F}()
    max_size = isinf(max_size) ? 0 : convert(Int, max_size)
    max_strs = isinf(max_strs) ? 0 : convert(Int, max_strs)

    if !isnothing(reducer) && isnothing(reduce_op)
        reduce_op = Base.:+
    end

    out_base, out_vals = polyrs(v₀, assembly_system, reducer, reduce_op, aggregator,
                                rejector, nothing, nothing, max_size, max_strs)
    if isnothing(reducer) && isnothing(aggregator) && isnothing(rejector)
        return out_base
    else
        return out_base, out_vals
    end
end


function polygenerate(callback::Function, assembly_system::AssemblySystem{T,F};
                      max_size=Inf, max_strs=Inf) where {T,F}
    values = [callback(monomer) for monomer in assembly_system.monomers]
    hashes = Set{HashType}()
    queue = Queue{Structure{T,F}}()
    open_bonds = Dict{HashType,Vector{Tuple{Int,Int}}}()

    for monomer in assembly_system.monomers
        hashval = hash(monomer)

        enqueue!(queue, monomer)
        push!(hashes, hashval)
        open_bonds[hashval] = all_open_bonds(monomer, assembly_system.intmat)
    end

    u = Structure{T,F}()
    n_strs = length(assembly_system.monomers)

    while !isempty(queue) && n_strs < max_strs
        v = dequeue!(queue)
        n, _, nfv = size(v)
        if n >= max_size
            continue
        end

        open_bonds_v = open_bonds[hash(v)]
        bonds = []
        for bond in open_bonds_v
            ai, j = bond
            copy!(u, v)

            success = attach_monomer!(u, bond, assembly_system, true)
            if !success || (hash(u) ∈ hashes)
                continue
            end

            next = copy(u)
            species_j, aj =  irg_unflatten(j, assembly_system._sides_sum)
            hashval = hash(next)

            monomer_opens = open_bonds[hash(assembly_system.monomers[species_j])]
            new_opens = [b .+ (nfv, 0) for b in monomer_opens if b[1] != aj]
            open_bonds[hashval] = filter(x -> x ∉ bonds && x[1] != ai, open_bonds_v)
            append!(open_bonds[hashval], new_opens)

            push!(values, callback(next))
            push!(hashes, hashval)
            push!(bonds, bond)
            enqueue!(queue, next)
            n_strs += 1
            n_strs == max_strs && break
        end
    end

    return values
end
function polygenerate(assembly_system::AssemblySystem; kwargs...)
    return polygenerate(identity, assembly_system; kwargs...)
end
