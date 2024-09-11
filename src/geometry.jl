using LinearAlgebra
using SparseArrays
using NautyGraphs
using Graphs: cycle_digraph

abstract type AbstractGeometry{T<:Integer,F<:AbstractFloat} end
nsites(geom::AbstractGeometry) = length(geom.xs)

struct PolygonGeometry{T<:Integer,F<:AbstractFloat} <: AbstractGeometry{T,F}
    xs::Vector{Point{2,F}}                  # Displacement vectors from the center to the binding sites
    θs_ref::Vector{Angle{F}}                # Rotation necessary to rotate side i into reference position
    Δθs::Vector{Angle{F}}                   # Rotation to be applied to the added particle (assumed to be in reference orientation) once attached to side i 
    corners::Vector{Point{2,F}}             # Corners of the Polygon used for overlap checking
    anatomy::NautyDiGraph
    encoder::PolyEncoder{T}
    R_min::F
    R_max::F

    function PolygonGeometry(n::T, a::F) where {T,F}
        r_in = convert(F, 0.5a * cot(π / n))
        r_out = convert(F, 0.5a * csc(π / n))

        θs_ref = [Angle{F}(2/n * i) for i in 0:(n-1)]
        Δθs = [Angle{F}(1 - 2/n*i) for i in 0:(n-1)]

        rs = r_in * ones(F, n)
        corners = [convert.(F, r_out * [cos(π * (θ.θ + 1/n)), sin(π * (θ.θ + 1/n))]) for θ in θs_ref] #TODO: clean up the θ.θ

        xs = [SVector{2,F}(pol2cart(r, -θ.θ - 1/2)) for (r, θ) in zip(rs, θs_ref)]
        anatomy = NautyDiGraph(cycle_digraph(n))
        encoder = PolyEncoder([[[T(i)] for i in 1:length(xs)]])

        R_min = minimum(norm.(xs))
        R_max = maximum(norm.(corners))

        return new{T,F}(xs, θs_ref, Δθs, corners, anatomy, encoder, R_min, R_max)
    end
end

const UnitTriangleGeometry = PolygonGeometry(DefInt(3), DefFloat(1.0))
const UnitSquareGeometry = PolygonGeometry(DefInt(4), DefFloat(1.0))
const UnitPentagonGeometry = PolygonGeometry(DefInt(5), DefFloat(1.0))
const UnitHexagonGeometry = PolygonGeometry(DefInt(6), DefFloat(1.0))


struct PolyhedronGeometry{T,F<:AbstractFloat} <: AbstractGeometry{T,F}
    xs::Vector{Point{3,F}}                # Displacement vectors from the center to the binding sites
    θs_ref::Vector{Quaternion{F}}         # Rotation necessary to rotate side i into reference position   
    Δθs::Vector{Quaternion{F}}            # Rotation to be applied to the added particle (assumed to be in reference orientation) once attached to side i 
    corners::Vector{Point{3,F}}
    anatomy::NautyDiGraph                 # Oriented graph of the dual polyhedron associated with the bb symmetry group
    encoder::PolyEncoder{T}
    R_min::F
    R_max::F

    function PolyhedronGeometry{T}(shape::Symbol, a::F) where {T,F}
        if shape == :cube
            xs = a * [Point{3,F}(1, 0, 0),
                      Point{3,F}(0, 0, -1),
                      Point{3,F}(0, 1, 0),
                      Point{3,F}(0, 0, 1),
                      Point{3,F}(0, -1, 0),
                      Point{3,F}(-1, 0, 0)]

            corners = vec([a * SVector(i, j, k) for i in F[-1, 1], j in F[-1, 1], k in F[-1, 1]])
            θs_ref = [Quaternion{F}(1, 0, 0, 0),
                      inv(√2) * Quaternion{F}(1, 0, -1, 0),
                      inv(√2) * Quaternion{F}(1, 0, 0, -1),
                      inv(√2) * Quaternion{F}(1, 0, 1, 0),
                      inv(√2) * Quaternion{F}(1, 0, 0, 1),
                      Quaternion{F}(0, 0, 0, 1)]
            Δθs = [Quaternion{F}(0, 0, 0, 1),
                   Quaternion{F}(0, 0, 0, 1) * (inv(√2) * Quaternion{F}(1, 0, -1, 0)),
                   inv(√2) * Quaternion{F}(1, 0, 0, -1),
                   Quaternion{F}(0, 0, 0, 1) * (inv(√2) * Quaternion{F}(1, 0, 1, 0)),
                   inv(√2) * Quaternion{F}(1, 0, 0, 1),
                   Quaternion{F}(1, 0, 0, 0)]

            g0 = sparse([0 1 0 0;
                         0 0 1 0;
                         0 0 0 1;
                         1 0 0 0])
            G = blockdiag(g0, g0, g0, g0, g0, g0)
            edges = [1, 7], [2, 20], [3, 13], [4, 10], [5, 21], [6, 17], 
            [8, 9], [11, 16], [12, 22], [14, 19], [15, 23], [18, 24]
            for e in edges
                i, j = e
                G[i, j] = G[j, i] = 1
            end
            anatomy = NautyDiGraph(G)
            encoder = PolyEncoder([[collect(T, i:i+3) for i in 1:4:24]])

            R_min = minimum(norm.(xs))
            R_max = maximum(norm.(corners))
        else
            error("Shape $shape not supported.")
        end
        
        return new{T,F}(xs, θs_ref, Δθs, corners, anatomy, encoder, R_min, R_max)
    end
end

const UnitCubeGeometry = PolyhedronGeometry{DefInt}(:cube, DefFloat(1.))
#TODO Implement basic 3d shapes: Platonic solids and polygon extrusions

dimension(::PolygonGeometry) = 2
dimension(::PolyhedronGeometry) = 3

function attachment_offset(site_i::Integer, site_j::Integer,
                           geom_i::AbstractGeometry,
                           geom_j::AbstractGeometry)
    x_i, Δθ = geom_i.xs[site_i], geom_i.Δθs[site_i]
    x_j, θ_ref = geom_j.xs[site_j], geom_j.θs_ref[site_j]

    Δψ = Δθ * θ_ref
    Δx = x_i - rotate(x_j, Δψ)
    return Δx, Δψ
end

# TODO: this function is susceptible to roundoff errors, be very careful with thresh
function face_pairs(Δx::Point{D,F}, ψi::RotationOperator{F}, ψj::RotationOperator{F}, geom_i::AbstractGeometry, geom_j::AbstractGeometry, convex=false, thresh=1e-2) where {D,F}
    pairs = Tuple{Int,Int}[]

    for (i, z_i) in enumerate(geom_i.xs), (j, z_j) in enumerate(geom_j.xs)
        z_i = rotate(z_i, ψi)
        z_j = rotate(z_j, ψj)

        matching = maximum(abs, z_i - z_j - Δx) < F(thresh)
        if matching
            push!(pairs, (i, j))
            if convex
                break
            end
        end
    end
    return pairs
end

function overlap(Δx, ψ_i, ψ_j, geom_i, geom_j; buffer=0.1)
    # Assumes lattice!!
    return norm(Δx) < 2geom_i.R_min - buffer
    # return __check_overlap(geom_i.mesh, transform(geom_j.mesh, Δx, Δψ))
end

function contact_status(Δx::Point, ψi::RotationOperator, ψj::RotationOperator, geom_i::AbstractGeometry, geom_j::AbstractGeometry)
    if norm(Δx) > geom_i.R_max + geom_j.R_max
        return true, Tuple{Int,Int}[]
    end

    if overlap(Δx, ψi, ψj, geom_i, geom_j)
        return false, Tuple{Int,Int}[]
    end

    pairs = face_pairs(Δx, ψi, ψj, geom_i, geom_j)
    return true, pairs
end
