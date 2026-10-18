# Volume of the convex hull of a set of points in d dimensions, used for
# functional richness. The R package calls Qhull through
# `geometry::convhulln(options = "FA Qt")`; here an incremental
# (beneath-beyond) algorithm with simplicial facets is used.

struct _Facet
    verts::Vector{Int}
    normal::Vector{Float64}
    offset::Float64
end

function _hyperplane(P::AbstractMatrix, verts::Vector{Int}, interior::Vector{Float64})
    d = size(P, 2)
    p1 = P[verts[1], :]
    A = Matrix{Float64}(undef, d - 1, d)
    for k in 2:d
        A[k-1, :] = P[verts[k], :] .- p1
    end
    N = nullspace(A)
    size(N, 2) == 1 || return nothing
    nrm = N[:, 1]
    nrm ./= norm(nrm)
    off = dot(nrm, p1)
    if dot(nrm, interior) > off
        nrm = -nrm
        off = -off
    end
    return _Facet(verts, nrm, off)
end

# d + 1 affinely independent points, chosen greedily (farthest from the
# affine span of those already chosen); `nothing` if the points are flat
function _initial_simplex(P::AbstractMatrix, tol::Float64)
    n, d = size(P)
    chosen = [argmin(view(P, :, 1))]
    basis = Vector{Vector{Float64}}()
    for _ in 1:d
        o = P[chosen[1], :]
        best, bi = -1.0, 0
        for i in 1:n
            i in chosen && continue
            v = P[i, :] .- o
            for b in basis
                v .-= dot(v, b) .* b
            end
            nv = norm(v)
            nv > best && ((best, bi) = (nv, i))
        end
        (bi == 0 || best <= tol) && return nothing
        v = P[bi, :] .- o
        for b in basis
            v .-= dot(v, b) .* b
        end
        push!(basis, v ./ norm(v))
        push!(chosen, bi)
    end
    return chosen
end

"""
    convex_hull_volume(P::AbstractMatrix) -> Float64

Volume of the convex hull of the rows of `P` (n points in d dimensions;
area for d = 2, length for d = 1). Returns 0 when the points do not span d
dimensions (Qhull, used by R, raises an error in that case).
"""
function convex_hull_volume(P::AbstractMatrix{<:Real})
    P = Matrix{Float64}(P)
    n, d = size(P)
    d == 0 && return 0.0
    d == 1 && return maximum(P) - minimum(P)
    n < d + 1 && return 0.0
    scale = maximum(maximum(P; dims = 1) .- minimum(P; dims = 1))
    scale == 0 && return 0.0
    tol = 1e-10 * scale
    simplex = _initial_simplex(P, tol)
    simplex === nothing && return 0.0
    interior = vec(sum(P[simplex, :]; dims = 1)) ./ (d + 1)
    facets = _Facet[]
    for k in 1:d+1
        f = _hyperplane(P, simplex[[i for i in 1:d+1 if i != k]], interior)
        f === nothing && return 0.0
        push!(facets, f)
    end
    eps_vis = 1e-12 * scale
    inset = Set(simplex)
    # process far points first: fewer facets are created and later removed
    order = sortperm([norm(P[i, :] .- interior) for i in 1:n]; rev = true)
    for p in order
        p in inset && continue
        x = P[p, :]
        visible = [dot(f.normal, x) - f.offset > eps_vis for f in facets]
        any(visible) || continue
        ridgecount = Dict{Vector{Int},Int}()
        for (f, v) in zip(facets, visible)
            v || continue
            for k in 1:d
                r = sort!(f.verts[[i for i in 1:d if i != k]])
                ridgecount[r] = get(ridgecount, r, 0) + 1
            end
        end
        newf = _Facet[]
        for (r, c) in ridgecount
            c == 1 || continue
            f = _hyperplane(P, [r; p], interior)
            f === nothing && continue
            push!(newf, f)
        end
        facets = [f for (f, v) in zip(facets, visible) if !v]
        append!(facets, newf)
    end
    vol = 0.0
    M = Matrix{Float64}(undef, d, d)
    for f in facets
        for (k, v) in enumerate(f.verts)
            M[k, :] = P[v, :] .- interior
        end
        vol += abs(det(M))
    end
    return vol / factorial(d)
end
