function f!(k::Polyform{T,F}, s::Polyform{T,F}, assembly_system) where {T,F}
    copy!(k, s)
    lower!(k, assembly_system)
    return k
end

function adj!(u::Polyform{T,F}, v::Polyform{T,F}, j::Integer, hashes::Vector{HashType},
              assembly_system::AssemblySystem) where {T,F}
    ### TODO: We should exploit that canonical labels are guaranteed (see Nauty User Guide p.4) to be in order of color
    ### So, if we know the color of the particle that would be removed, we can filter the possible offspring
    if size(v) == 0
        bblocks = buildingblocks(assembly_system)
        if j > length(bblocks)
            return false, j + 1
        end

        copy!(u, bblocks[j])
        return true, j + 1
    end

    copy!(u, v)
    success, j = raise!(u, j, assembly_system)

    if success && rhash(u) ∈ hashes
        return adj!(u, v, j + 1, hashes, assembly_system)
    end

    return success, j + 1
end

FnorNothing = Union{<:Function,Nothing}

function polyrs(v₀::Polyform{D,T,F},
                assembly_system::AssemblySystem{D,T,F},
                reducer::FnorNothing,
                reduce_op::FnorNothing,
                aggregator::FnorNothing,
                rejector::FnorNothing,
                max_depth::Int=0,
                max_strs::Int=0) where {D,T,F}
    reducing = !isnothing(reducer)
    if reducing
        reduce_val = reducer(v₀, assembly_system)
    else
        reduce_val = nothing
    end
    aggregating = !isnothing(aggregator)
    if aggregating
        aggregate_val = typeof(aggregator(v₀, assembly_system)[2])[]
    else
        aggregate_val = nothing
    end
    rejecting = !isnothing(rejector)

    depth = 0
    n_strs = 0

    v = Polyform{D,T,F}()
    u = Polyform{D,T,F}()
    k = Polyform{D,T,F}()

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
    break_triggered = false

    while true
        success, j_new = adj!(u, v, js[end], hashes[end], assembly_system)
        next = u
        if success
            js[end] = j_new
            push!(hashes[end], rhash(next))

            f!(k, next, assembly_system)

            if k ≃ v
                if rejecting
                    reject_val = rejector(next, assembly_system)
                    # 0: dont reject, 1: reject post, 2: reject pre, 3: break immediately
                    reject_val == 2 && continue
                end
                if reducing
                    reduce_val = reduce_op(reduce_val, reducer(next, assembly_system))
                end
                if aggregating
                    aggregate, agr_val = aggregator(next, assembly_system)
                    aggregate && push!(aggregate_val, agr_val)
                end

                n_strs += 1
                depth += 1
                if depth > lowest_depth
                    lowest_depth = depth
                end

                if rejecting
                    if reject_val == 1
                        depth -= 1
                        continue
                    elseif reject_val == 3
                        break_triggered = true
                        break
                    end
                end

                if n_strs >= max_strs
                    break_triggered = true
                    break
                end

                if depth == max_depth
                    depth -= 1
                    continue
                end

                copy!(v, next)
                push!(js, 1)
                push!(hashes, HashType[])
            end

            continue
        end

        if depth > 0
            f!(k, v, assembly_system)    # TODO could speed this up by storing structures
            copy!(v, k)

            depth -= 1

            pop!(js)
            pop!(hashes)
        else
            break
        end
    end

    return (n_strs, lowest_depth), (reduce_val, aggregate_val, break_triggered)
end

function polyenum(assembly_system::AssemblySystem{D,T,F,G};
                  max_size=nothing, max_strs=nothing,
                  reducer::FnorNothing=nothing,
                  reduce_op::FnorNothing=Base.:+,
                  aggregator::FnorNothing=nothing,
                  rejector::FnorNothing=nothing) where {D,T,F,G}
    v₀ = Polyform{D,T,F}()
    max_size = isnothing(max_size) ? 0 : max_size
    max_strs = isnothing(max_strs) ? 0 : max_strs

    out_base, out_vals = polyrs(v₀, assembly_system, reducer, reduce_op, aggregator, rejector, max_size, max_strs)
    if isnothing(reducer) && isnothing(aggregator) && isnothing(rejector)
        return out_base
    else
        return out_base, out_vals
    end
end


function polygen(callback::Function, assembly_system::AssemblySystem{D,T,F,G};
                      max_size=Inf, max_strs=Inf) where {D,T,F,G}
    bblocks = buildingblocks(assembly_system)

    values = [callback(bblock) for bblock in bblocks]
    hashes = Set{HashType}()
    queue = Queue{Polyform{D,T,F}}()

    for bblock in bblocks
        hashval = rhash(bblock)
        enqueue!(queue, bblock)
        push!(hashes, hashval)
    end

    u = Polyform{D,T,F}()
    n_strs = length(bblocks)

    while !isempty(queue) && n_strs < max_strs
        v = dequeue!(queue)
        n = size(v)
        if n >= max_size
            continue
        end

        j = 1
        vert, partner_label = open_bond(v, assembly_system, j)
        while !iszero(vert)
            copy!(u, v)

            success = attach_monomer!(u, vert, partner_label, assembly_system, true)
            hashval = rhash(u)

            if success && (hashval ∉ hashes)
                next = copy(u)
            
                push!(values, callback(next))
                push!(hashes, hashval)
                enqueue!(queue, next)
                
                n_strs += 1
                n_strs == max_strs && break
            end

            j += 1
            vert, partner_label = open_bond(v, assembly_system, j)
        end
    end

    return values
end
function polygen(assembly_system::AssemblySystem; kwargs...)
    return polygen(identity, assembly_system; kwargs...)
end
