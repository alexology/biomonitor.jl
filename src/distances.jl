# Distances, ordination and related tools (native versions of the ade4,
# ape and stats functions used by biomonitoR).

const DEFAULT_TACHET_BLOCKS = [8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8]

_D(d::DistMatrix) = d.d
_D(d::AbstractMatrix) = Matrix{Float64}(d)
_relabel(d::DistMatrix, D) = DistMatrix(d.labels, D)
_relabel(d::AbstractMatrix, D) = D

"double centring with weights (ade4 `bicenter.wt`)"
function bicenter(X::AbstractMatrix; row_wt = ones(size(X, 1)), col_wt = ones(size(X, 2)))
    rw = row_wt ./ sum(row_wt)
    cw = col_wt ./ sum(col_wt)
    rowmean = vec(sum(rw .* X; dims = 1))
    colmean = vec(sum(cw' .* X; dims = 2))
    colmean .-= sum(rowmean .* cw)
    return X .- rowmean' .- colmean
end

# eigen-decomposition with values in decreasing order (as R `eigen(symmetric = TRUE)`)
function _eigen_desc(M::AbstractMatrix)
    E = eigen(Symmetric(Matrix{Float64}((M + M') / 2)))
    o = sortperm(E.values; rev = true)
    return E.values[o], E.vectors[:, o]
end

_gower_delta(D) = -0.5 .* bicenter(D .* D)

# ---------------------------------------------------------------------------
# Fuzzy coding and Gower-like distance for fuzzy traits (ade4 prep.fuzzy and
# dist.ktab(type = "F"))
# ---------------------------------------------------------------------------

function _prep_fuzzy(V::AbstractMatrix, col_blocks::AbstractVector{<:Integer})
    sum(col_blocks) == size(V, 2) || throw(ArgumentError("non convenient data in col.blocks"))
    keepcol = [!all(isnan, view(V, :, j)) for j in axes(V, 2)]
    blockid = reduce(vcat, [fill(k, b) for (k, b) in enumerate(col_blocks)]; init = Int[])
    V = V[:, keepcol]
    blockid = blockid[keepcol]
    blocks = [count(==(k), blockid) for k in unique(blockid)]
    out = similar(V, Float64)
    k2 = 0
    for b in blocks
        k1 = k2 + 1
        k2 += b
        for i in axes(V, 1)
            s = sum(view(V, i, k1:k2))
            if isnan(s) || s == 0
                out[i, k1:k2] .= NaN
            else
                out[i, k1:k2] = V[i, k1:k2] ./ s
            end
        end
    end
    return out, blocks, keepcol
end

"""
    prep_fuzzy(t::TraitTable, col_blocks) -> TraitTable

Fuzzy-coded traits rescaled so that the modalities of each trait sum to 1
(ade4 `prep.fuzzy`). `col_blocks` gives the number of modalities of each
trait; rows whose modalities of a trait are all zero become `NaN` for that
trait.
"""
function prep_fuzzy(t::TraitTable, col_blocks::AbstractVector{<:Integer})
    V, _, keep = _prep_fuzzy(t.values, col_blocks)
    return TraitTable(t.taxa, t.traits[keep], V)
end

"""
    gower_fuzzy(V::AbstractMatrix, col_blocks; tol = 1e-8) -> Matrix
    gower_fuzzy(t::TraitTable, col_blocks) -> DistMatrix

Distance between taxa described by fuzzy-coded traits, as computed by
`ade4::dist.ktab(ktab.list.df(list(prep.fuzzy(...))), type = "F")`: for each
trait the overlap distance ``1 - \\sum p_i q_i / \\sqrt{\\sum p_i^2 \\sum q_i^2}``
(Manly), averaged over the traits and square-rooted. Traits with missing
values are skipped pair by pair. The input must already be prepared with
[`prep_fuzzy`](@ref) (rows of each block sum to 1 or are `NaN`).
"""
function gower_fuzzy(V::AbstractMatrix, col_blocks::AbstractVector{<:Integer}; tol::Real = 1e-8)
    n = size(V, 1)
    napres = any(isnan, V)
    total = zeros(n, n)
    ntvar = zeros(n, n)
    k2 = 0
    for b in col_blocks
        k1 = k2 + 1
        k2 += b
        B = V[:, k1:k2]
        ok = [!any(isnan, view(B, i, :)) for i in 1:n]
        for j in 1:n, i in 1:n
            if ok[i] && ok[j]
                if i != j
                    p = view(B, i, :)
                    q = view(B, j, :)
                    sp = sum(p); sq = sum(q)
                    total[i, j] += 1 - dot(p ./ sp, q ./ sq) / sqrt(sum(abs2, p ./ sp)) / sqrt(sum(abs2, q ./ sq))
                end
                ntvar[i, j] += 1
            end
        end
    end
    for i in 1:n, j in 1:n
        total[i, j] < tol && (total[i, j] = 0.0)
    end
    D = zeros(n, n)
    for j in 1:n, i in 1:n
        i == j && continue
        D[i, j] = napres ? sqrt(total[i, j] / ntvar[i, j]) : sqrt(total[i, j] / length(col_blocks))
    end
    return D
end
gower_fuzzy(t::TraitTable, col_blocks::AbstractVector{<:Integer}; kw...) =
    DistMatrix(t.taxa, gower_fuzzy(t.values, col_blocks; kw...))

"""
    euclidean_dist(V::AbstractMatrix; scale = true) -> Matrix
    euclidean_dist(t::TraitTable; scale = true) -> DistMatrix

Euclidean distance between rows, after centring and scaling the columns to
unit variance when `scale = true` (R `dist(scale(x))`).
"""
function euclidean_dist(V::AbstractMatrix; scale::Bool = true)
    X = Matrix{Float64}(V)
    if scale
        X = (X .- mean(X; dims = 1)) ./ std(X; dims = 1)
    end
    n = size(X, 1)
    return [sqrt(sum(abs2, view(X, i, :) .- view(X, j, :))) for i in 1:n, j in 1:n]
end
euclidean_dist(t::TraitTable; scale::Bool = true) = DistMatrix(t.taxa, euclidean_dist(t.values; scale))

# ---------------------------------------------------------------------------
# Euclidean property and corrections
# ---------------------------------------------------------------------------

"""
    is_euclid(d; tol = 1e-7) -> Bool

`true` if the distance matrix is Euclidean, i.e. its Gower-centred matrix
has no negative eigenvalue (relative tolerance `tol`, ade4 `is.euclid`).
"""
function is_euclid(d; tol::Real = 1e-7)
    λ, _ = _eigen_desc(_gower_delta(_D(d)))
    return λ[end] / λ[1] > -tol
end

"""
    cailliez(d; tol = 1e-7, cor_zero = true)

Cailliez (1983) correction: the smallest constant added to the distances
that makes them Euclidean. With `cor_zero = true` zero distances stay zero.
"""
function cailliez(d; tol::Real = 1e-7, cor_zero::Bool = true)
    if is_euclid(d)
        @warn "Euclidean distance found : no correction need"
        return d
    end
    D = _D(d)
    n = size(D, 1)
    M = [zeros(n, n) -bicenter(D .* D); -Matrix{Float64}(I, n, n) 2 .* bicenter(D)]
    λ = eigvals(M)
    c = maximum(real(l) for l in λ if imag(l) < tol)
    R = copy(D)
    for j in 1:n, i in 1:n
        i == j && continue
        if cor_zero
            R[i, j] > tol && (R[i, j] += c)
        else
            R[i, j] += c
        end
    end
    return _relabel(d, R)
end

"""
    lingoes(d; tol = 1e-7, cor_zero = true)

Lingoes (1971) correction: ``\\sqrt{d^2 + 2c}``, with c the absolute value
of the most negative eigenvalue.
"""
function lingoes(d; tol::Real = 1e-7, cor_zero::Bool = true)
    if is_euclid(d)
        @warn "Euclidean distance found : no correction need"
        return d
    end
    D = _D(d)
    λ, _ = _eigen_desc(_gower_delta(D))
    c = abs(λ[end])
    R = D .* D
    n = size(D, 1)
    for j in 1:n, i in 1:n
        i == j && continue
        if cor_zero
            R[i, j] > tol && (R[i, j] = sqrt(R[i, j] + 2c))
        else
            R[i, j] = sqrt(R[i, j] + 2c)
        end
    end
    return _relabel(d, R)
end

"""
    quasieuclid(d)

Euclidean distance on the principal coordinates with positive eigenvalues
(ade4 `quasieuclid`).
"""
function quasieuclid(d)
    if is_euclid(d)
        @warn "Euclidean distance found : no correction need"
        return d
    end
    D = _D(d)
    λ, U = _eigen_desc(_gower_delta(D))
    k = count(>(0), λ)
    X = U[:, 1:k] .* sqrt.(λ[1:k])'
    n = size(D, 1)
    return _relabel(d, [sqrt(sum(abs2, view(X, i, :) .- view(X, j, :))) for i in 1:n, j in 1:n])
end

_sqrt_dist(d) = _relabel(d, sqrt.(_D(d)))

function _correct(d, correction; tol = 1e-7, cor_zero = true)
    c = string(correction)
    c == "none" && return d
    c == "cailliez" && return cailliez(d; tol, cor_zero)
    c == "lingoes" && return lingoes(d; tol, cor_zero)
    c == "sqrt" && return _sqrt_dist(d)
    c == "quasi" && return quasieuclid(d)
    throw(ArgumentError("correction must be one of none, cailliez, lingoes, sqrt or quasi"))
end

# ---------------------------------------------------------------------------
# Principal coordinates
# ---------------------------------------------------------------------------

"""
    pcoa(d; nf = 2, tol = 1e-7) -> (coords, eig, rank)

Principal coordinate analysis with uniform weights (ade4 `dudi.pco`):
`coords` are the coordinates on the first `nf` axes (at most the number of
positive eigenvalues), `eig` the positive eigenvalues divided by n.
"""
function pcoa(d; nf::Integer = 2, tol::Real = 1e-7)
    D = _D(d)
    n = size(D, 1)
    w = fill(1 / n, n)
    delta = -0.5 .* bicenter(D .* D; row_wt = w, col_wt = w)
    delta = delta .* sqrt.(w) .* sqrt.(w)'
    λ, U = _eigen_desc(delta)
    λ[end] / λ[1] < -tol && @debug "Non euclidean distance"
    r = count(>(λ[1] * tol), λ)
    nf <= 0 && (nf = 2)
    nf = min(nf, r)
    coords = U[:, 1:r] .* sqrt.(λ[1:r])' ./ sqrt.(w)
    return (coords = coords[:, 1:nf], eig = λ[1:r], rank = r)
end

# classical scaling eigenvalues (R `cmdscale(d, k = n - 1, eig = TRUE)`)
function _cmdscale(D::AbstractMatrix)
    λ, U = _eigen_desc(_gower_delta(D))
    n = size(D, 1)
    k = n - 1
    ev = λ[1:k]
    pos = ev .> 0
    points = U[:, 1:k][:, pos] .* sqrt.(ev[pos])'
    return (points = points, eig = λ)
end

# ---------------------------------------------------------------------------
# Rao's quadratic entropy (ade4 divc and divcmax)
# ---------------------------------------------------------------------------

"""
    divcmax(d; epsilon = 1e-8) -> (value, vectors)

Maximal value of Rao's quadratic entropy for the Euclidean distance `d`
(ade4 `divcmax`, Pavoine et al. 2005).
"""
function divcmax(d; epsilon::Real = 1e-8)
    is_euclid(d) || throw(ArgumentError("Euclidean property is expected for dis"))
    D2 = _D(d) .^ 2 ./ 2
    n = size(D2, 1)
    x0 = vec(sum(D2; dims = 2)) ./ sum(D2)
    xk = copy(x0)
    relax = 0
    sature = falses(n)
    while true
        while true
            maxi = dot(xk, D2 * xk)
            deltaf = -2 .* (D2 * xk)
            sature = abs.(xk) .< epsilon
            if relax != 0
                sature[relax] = false
                relax = 0
            end
            yk = -deltaf
            yk[sature] .= 0
            ns = .!sature
            yk[ns] .-= mean(yk[ns])
            maximum(abs, yk) < epsilon && break
            neg = yk .< 0
            alpha_max = any(neg) ? minimum(-xk[neg] ./ yk[neg]) : Inf
            alpha_opt = -dot(xk, D2 * yk) / dot(yk, D2 * yk)
            alpha = (alpha_opt > alpha_max || alpha_opt < 0) ? alpha_max : alpha_opt
            abs(maxi - dot(xk .+ alpha .* yk, D2 * (xk .+ alpha .* yk))) < epsilon && break
            xk = xk .+ alpha .* yk
        end
        all(.!sature) && break
        vectD2 = D2 * xk
        u = 2 .* (mean(vectD2[.!sature]) .- vectD2[sature])
        minimum(u) >= 0 && break
        satu = findall(sature)
        relax = satu[findfirst(==(minimum(u)), u)]
    end
    num = copy(xk)
    num[num .< epsilon] .= 0
    value = _divc_one(num, _D(d))
    return (value = value, weights = num)
end

_divc_one(p, D) = sum(p) < 1e-16 ? 0.0 : dot(p, (D .^ 2) * p) / 2 / sum(p)^2

"""
    divc(A::AbstractMatrix, d = nothing; scale = false) -> Vector

Rao's quadratic entropy of each community (column of the taxa × communities
matrix `A`) for the distance `d` (ade4 `divc`). Without `d` all distances
are ``\\sqrt 2`` and the index is the Gini-Simpson diversity. With
`scale = true` values are divided by [`divcmax`](@ref).
"""
function divc(A::AbstractMatrix, d = nothing; scale::Bool = false)
    any(<(0), A) && throw(ArgumentError("Negative value in df"))
    n = size(A, 1)
    D = d === nothing ? (ones(n, n) - I) .* sqrt(2) : _D(d)
    size(D, 1) == n || throw(DimensionMismatch("Non convenient df"))
    res = [_divc_one(view(A, :, j), D) for j in axes(A, 2)]
    scale && (res ./= divcmax(d === nothing ? D : d).value)
    return res
end

# ---------------------------------------------------------------------------
# Mantel test, community distances, minimum spanning tree
# ---------------------------------------------------------------------------

"""
    mantel_test(d1, d2; nrepet = 999, rng = Random.default_rng()) -> (obs, pvalue)

Mantel test between two distance matrices: Pearson correlation of the
distances and permutation p-value (alternative "greater"), as
`ade4::mantel.randtest`.
"""
function mantel_test(d1, d2; nrepet::Integer = 999, rng::AbstractRNG = Random.default_rng())
    A = _D(d1)
    B = _D(d2)
    b = lowertri(B)
    obs = cor(lowertri(A), b)
    n = size(A, 1)
    cnt = 0
    for _ in 1:nrepet
        p = randperm(rng, n)
        cor(lowertri(A[p, p]), b) >= obs && (cnt += 1)
    end
    return (obs = obs, pvalue = (cnt + 1) / (nrepet + 1))
end

"""
    bray_curtis(M::AbstractMatrix) -> Matrix

Bray-Curtis dissimilarity between the rows of `M` (e.g. samples × taxa).
"""
function bray_curtis(M::AbstractMatrix)
    n = size(M, 1)
    s = vec(sum(M; dims = 2))
    return [i == j ? 0.0 : 1 - 2 * sum(min.(view(M, i, :), view(M, j, :))) / (s[i] + s[j]) for i in 1:n, j in 1:n]
end

"""
    jaccard_dist(M::AbstractMatrix) -> Matrix

``\\sqrt{1 - J}`` with J the Jaccard similarity of the rows of `M` (values
> 0 are presences), as `ade4::dist.binary(method = 1)`.
"""
function jaccard_dist(M::AbstractMatrix)
    P = M .> 0
    n = size(P, 1)
    D = zeros(n, n)
    for j in 1:n, i in 1:n
        i == j && continue
        a = count(view(P, i, :) .& view(P, j, :))
        o = count(view(P, i, :) .| view(P, j, :))
        D[i, j] = sqrt(1 - a / o)
    end
    return D
end

# minimum spanning tree adjacency matrix, translation of `ape::mst`
function _mst(D::AbstractMatrix)
    X = Matrix{Float64}(D)
    n = size(X, 1)
    N = zeros(Int, n, n)
    large = maximum(X) + 1
    for i in 1:n
        X[i, i] = large
    end
    index_i = 1
    tree = Int[]
    for _ in 1:n-1
        push!(tree, index_i)
        m = [minimum(view(X, :, t)) for t in tree]
        a = [argmin(view(X, :, t)) for t in tree]
        b = argmin(m)
        index_j = tree[b]
        index_i = a[b]
        N[index_i, index_j] = 1
        N[index_j, index_i] = 1
        for j in tree
            X[index_i, j] = large
            X[j, index_i] = large
        end
    end
    return N
end
