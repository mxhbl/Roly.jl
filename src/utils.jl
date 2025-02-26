HashType = UInt
DefInt, DefFloat = Int16, Float32

function are_bridge(g::AbstractNautyGraph, vs::AbstractVector{<:Integer})
    neighs = zeros(Bool, nv(g))
    buffer = zeros(Int, nv(g))

    for v in vs
        k = NautyGraphs.outneighbors!(buffer, g, v)
        for neigh in @view buffer[1:k]
            neighs[neigh] = true
        end
    end
    neighs[vs] .= false

    if length(neighs) < 2
        return false
    end

    forbidden = zeros(Bool, nv(g))
    forbidden[vs] .= true

    return !vertices_connected(g, findfirst(neighs), findall(neighs), forbidden)
end
function vertices_connected(g::NautyDiGraph, v0::Integer, 
    targets::AbstractVector{<:Integer}, forbidden::AbstractVector{<:Integer})

    n = nv(g)
    neighs = zeros(Int, n)
    explored = zeros(Bool, n)
    
    queue = zeros(Cint, n)
    queue[1] = v0
    q_start = 1
    q_end = 2

    explored[v0] = true
    while q_end > q_start
        v = queue[q_start]
        q_start += 1
        n_neighs = NautyGraphs.outneighbors!(neighs, g, v)
        for neigh in @view neighs[1:n_neighs]
            if forbidden[neigh]
                continue
            end

            if !explored[neigh]
                explored[neigh] = true

                if all(@view explored[targets])
                    return true
                end

                queue[q_end] = neigh
                q_end += 1                
            end
        end
    end

    return false
end


function irg_flatten(a::Integer, b::Integer, intervals::AbstractVector{<:Integer})
    if a == 1
        return b
    else
        return intervals[a-1] + b
    end
end
function irg_unflatten(i::Integer, intervals::AbstractVector{<:Integer})
    if i <= intervals[1]
        return (1, i)
    else
        a = searchsortedlast(intervals, i, lt=≤) + 1
        b = i - intervals[a - 1]
        return (a, b)
    end
end

function find_nth(f::Function, A, n)
    # Return the index of the nth true value of f.(A), as well as how many matches were found
    # TODO this is a weird function, make it better
    j = 0
    for (k, a) in enumerate(A)
        if f(a) j += 1 end
        j == n && return k, j
    end
    return nothing, j
end
find_nth(A, n) = find_nth(identity, A, n)


mutable struct PathIterator{G,T}
    g::G
    v::T
    w::T
end
PathIterator(g, v, w) = PathIterator(g, promote(v, w)...)
Base.IteratorSize(::PathIterator) = Base.SizeUnknown()

function Base.iterate(pathitr::PathIterator)
    path = zeros(Int, nv(pathitr.g))
    path[pathitr.v] = 1
    path = complete_path!(path, pathitr.g, pathitr.w)    
    return copy(path), path #TODO: optimize this
end

function Base.iterate(pathitr::PathIterator, state)
    path = state
    neighs = zeros(Int, nv(pathitr.g))

    path_length = maximum(path)

    u = pathitr.w
    for _ in 1:path_length-1
        k = path[u]
        path[u] = 0

        parent = findfirst(x->x==k-1, path)
        n_neighs = NautyGraphs.outneighbors!(neighs, pathitr.g, parent)
        idx = @views searchsortedfirst(neighs[1:n_neighs], u)

        if idx < n_neighs && path[neighs[idx + 1]] == 0
            path[neighs[idx + 1]] = k
            path = complete_path!(path, pathitr.g, pathitr.w, neighs)
            return copy(path), path
        else
            k -= 1
            u = parent
        end
    end

    return nothing
end

function complete_path!(path, g, w, neighs=nothing)
    ## Completes the path by depth first traversal
    ## Neighbors are picked in order of vertex number (NOT label)
    neighs = isnothing(neighs) ? zeros(Int, length(path)) : neighs

    v = argmax(path)
    v_last = 0
    k = sum(!iszero, path) + 1

    while path[w] == 0
        n_neighs = NautyGraphs.outneighbors!(neighs, g, v)
        n0 = let i=findfirst(x->x==v_last, @view neighs[1:n_neighs]) # use searchsortedfirst
            isnothing(i) ? 1 : i+1
        end
        v_last = 0

        success = false
        for neigh in @view neighs[n0:n_neighs] # Assume neighbors are sorted
            if path[neigh] == 0
                v = neigh
                path[v] = k
                k += 1
                success = true
                break
            end
        end

        if !success
            path[v] = 0
            k -= 1
            v_last = v
            v = argmax(path)
        end
    end

    return path
end
