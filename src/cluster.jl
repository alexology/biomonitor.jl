# Hierarchical clustering, silhouette and indicator value analysis
# (the computations behind R `plot_indicator_taxa`).

"""
    hclust(d; method = "ward.D2") -> NamedTuple

Agglomerative hierarchical clustering of a distance matrix, following the
algorithm of R `stats::hclust` (Murtagh's nearest-neighbour Fortran code).
Methods: `"ward.D"`, `"ward.D2"`, `"single"`, `"complete"`, `"average"`,
`"mcquitty"`, `"median"`, `"centroid"`. Returns `merges` (pairs of
cluster representatives, the second absorbed into the first), `height` and
`n`; use [`cutree`](@ref) to obtain groups.
"""
function hclust(d; method = "ward.D2")
    m = string(method)
    m in ("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid") ||
        throw(ArgumentError("invalid clustering method $(method)"))
    D = copy(_D(d))
    n = size(D, 1)
    n < 2 && throw(ArgumentError("must have n >= 2 objects to cluster"))
    m == "ward.D2" && (D .= D .^ 2)
    membr = ones(n)
    flag = trues(n)
    nn = zeros(Int, n)
    disnn = fill(Inf, n)
    for i in 1:n-1
        dmin = Inf; jj = 0
        for j in i+1:n
            D[i, j] < dmin && ((dmin, jj) = (D[i, j], j))
        end
        nn[i] = jj
        disnn[i] = dmin
    end
    merges = Tuple{Int,Int}[]
    height = Float64[]
    ncl = n
    while ncl > 1
        dmin = Inf; im = 0; jm = 0
        for i in 1:n-1
            if flag[i] && disnn[i] < dmin
                dmin = disnn[i]; im = i; jm = nn[i]
            end
        end
        ncl -= 1
        i2 = min(im, jm); j2 = max(im, jm)
        push!(merges, (i2, j2))
        push!(height, m == "ward.D2" ? sqrt(dmin) : dmin)
        flag[j2] = false
        dmin = Inf; jj = 0
        for k in 1:n
            (flag[k] && k != i2) || continue
            dik, djk, dij = D[i2, k], D[j2, k], D[i2, j2]
            ni, nj, nk = membr[i2], membr[j2], membr[k]
            newd = if m == "ward.D" || m == "ward.D2"
                ((ni + nk) * dik + (nj + nk) * djk - nk * dij) / (ni + nj + nk)
            elseif m == "single"
                min(dik, djk)
            elseif m == "complete"
                max(dik, djk)
            elseif m == "average"
                (ni * dik + nj * djk) / (ni + nj)
            elseif m == "mcquitty"
                0.5 * dik + 0.5 * djk
            elseif m == "median"
                0.5 * dik + 0.5 * djk - 0.25 * dij
            else # centroid
                (ni * dik + nj * djk - ni * nj * dij / (ni + nj)) / (ni + nj)
            end
            D[i2, k] = newd
            D[k, i2] = newd
            if i2 < k
                newd < dmin && ((dmin, jj) = (newd, k))
            elseif newd < disnn[k]
                disnn[k] = newd
                nn[k] = i2
            end
        end
        membr[i2] += membr[j2]
        disnn[i2] = dmin
        nn[i2] = jj
        for i in 1:n-1
            if flag[i] && (nn[i] == i2 || nn[i] == j2)
                dmin = Inf; jj = 0
                for j in i+1:n
                    flag[j] && D[i, j] < dmin && ((dmin, jj) = (D[i, j], j))
                end
                nn[i] = jj
                disnn[i] = dmin
            end
        end
    end
    return (merges = merges, height = height, n = n, method = m)
end

"""
    cutree(tree, k) -> Vector{Int}

Groups obtained by cutting a [`hclust`](@ref) tree into `k` clusters;
clusters are numbered in order of first appearance of the objects (as R
`cutree`).
"""
function cutree(tree, k::Integer)
    n = tree.n
    1 <= k <= n || throw(ArgumentError("k must be between 1 and $(n)"))
    parent = collect(1:n)
    find(i) = (while parent[i] != i; i = parent[i]; end; i)
    for (a, b) in tree.merges[1:n-k]
        parent[find(b)] = find(a)
    end
    roots = [find(i) for i in 1:n]
    ids = Dict{Int,Int}()
    return [get!(ids, r, length(ids) + 1) for r in roots]
end

# Silhouette widths (Rousseeuw 1987, as R `cluster::silhouette`). A sample
# alone in its cluster has width 0; the R helper `silhouette_biomonitor` drops
# it instead, which inflates the average (see docs/src/r_to_julia.md).
function _silhouette_rows(D::Matrix{Float64}, cl::Vector{Int})
    rows = NamedTuple{(:cluster, :res),Tuple{Int,Float64}}[]
    labs = sort(unique(cl))
    for j in labs
        mem = findall(==(j), cl)
        others = sort(unique(cl[cl .!= j]))
        for a in mem
            if length(mem) == 1
                push!(rows, (cluster = j, res = 0.0))
                continue
            end
            md = mean(D[b, a] for b in mem if b != a)
            mn = minimum(mean(D[findall(==(o), cl), a]) for o in others)
            r = md != mn ? (mn - md) / max(md, mn) : 0.0
            push!(rows, (cluster = j, res = r))
        end
    end
    return rows
end

"""
    indval(A::AbstractMatrix, clusters; presence_absence = false) -> Matrix

Indicator value (Dufrêne & Legendre 1997) of each taxon (rows of the taxa ×
samples matrix `A`) for each cluster of samples (`clusters`, one label per
column; output columns follow the sorted labels): the square root of
specificity × fidelity. With `presence_absence = true` specificity is
computed on frequencies.
"""
function indval(A::AbstractMatrix, clusters::AbstractVector; presence_absence::Bool = false)
    labs = sort(unique(clusters))
    P = A .> 0
    np = [mean(P[i, clusters .== l]) for i in axes(A, 1), l in labs]
    npk = vec(sum(np; dims = 2))
    if presence_absence
        return sqrt.(np ./ npk .* np)
    end
    na = [mean(A[i, clusters .== l]) for i in axes(A, 1), l in labs]
    nak = vec(sum(na; dims = 2))
    return sqrt.(na ./ nak .* np)
end

"""
    indicator_taxa(x::AggregatedData; clusters = nothing, tax_lev = :Taxa, nmax = 10,
                   thresh = 25, method = "ward.D2", parent = :Order) -> NamedTuple

Indicator taxa analysis behind R `plot_indicator_taxa`. Samples are
clustered on Bray-Curtis (or Jaccard for presence-absence data)
dissimilarities with [`hclust`](@ref); the number of clusters (2 to `nmax`)
with the highest average silhouette is kept, unless `clusters` (a `Dict` or
`Pair`s sample => cluster label) is given. Indicator values
([`indval`](@ref)) of at least `thresh`% are reported, with the name of the
`parent` level of each taxon. The Sankey plot of R is not drawn; the
returned data can be plotted with any Julia plotting package.

Returns `clusters` (sample => cluster), `k`, `avg_silhouette`,
`silhouette_by_cluster`, `indval` (taxa × clusters [`LabeledMatrix`](@ref))
and `indicators` (taxon, cluster, value, parent), sorted by value.
"""
function indicator_taxa(x::AggregatedData; clusters = nothing, tax_lev = :Taxa, nmax::Int = 10, thresh::Real = 25,
                        method = "ward.D2", parent = :Order)
    _check_aggregated(x)
    smp = samples(x)
    n = length(smp)
    t = _drop_unassigned(_level(x, tax_lev))
    rank_index(parent) >= rank_index(tax_lev) && throw(ArgumentError("parent cannot have a finer taxonomic resolution than tax_lev."))
    M = permutedims(t.abundance)
    D = is_bin(x) ? jaccard_dist(M) : bray_curtis(M)
    tree = hclust(D; method)
    if clusters === nothing
        nmax >= n && throw(ArgumentError("nmax cannot be greater or equal to the number of samples."))
        best = (k = 0, s = -Inf, cl = Int[], rows = NamedTuple[])
        for k in 2:nmax
            cl = cutree(tree, k)
            rows = _silhouette_rows(D, cl)
            s = isempty(rows) ? NaN : mean(r.res for r in rows)
            s > best.s && (best = (k = k, s = s, cl = cl, rows = rows))
        end
        cl = best.cl
        k = best.k
        avg = best.s
        rows = best.rows
        labels_ = string.(cl)
    else
        cd = Dict(String(string(a)) => string(b) for (a, b) in clusters)
        all(s -> haskey(cd, s), smp) || throw(ArgumentError("clusters must give a cluster to every sample"))
        labels_ = [cd[s] for s in smp]
        u = sort(unique(labels_))
        length(u) >= n && throw(ArgumentError("the number of clusters cannot be greater or equal to the number of samples."))
        cl = [findfirst(==(l), u) for l in labels_]
        k = length(u)
        rows = _silhouette_rows(D, cl)
        avg = isempty(rows) ? NaN : mean(r.res for r in rows)
    end
    bycl = [(cluster = c, silhouette = mean(r.res for r in rows if r.cluster == c))
            for c in sort(unique(r.cluster for r in rows))]
    iv = indval(t.abundance, cl; presence_absence = is_bin(x))
    ul = sort(unique(cl))
    clnames = [labels_[findfirst(==(c), cl)] for c in ul]
    ivm = LabeledMatrix(t.taxa, clnames, iv)
    pidx = rank_index(parent)
    lidx = rank_index(tax_lev)
    d = x.data
    parentof = Dict{String,String}()
    for i in eachindex(d.taxa)
        nm = lidx == length(TAXONOMIC_LEVELS) ? d.taxa[i] : d.ranks[i, lidx]
        p = d.ranks[i, pidx]
        haskey(parentof, nm) || (parentof[nm] = isempty(p) ? "missing" : p)
    end
    inds = [(taxon = t.taxa[i], cluster = clnames[j], value = round(iv[i, j]; digits = 2),
             parent = get(parentof, t.taxa[i], "missing"))
            for i in axes(iv, 1) for j in axes(iv, 2) if iv[i, j] >= thresh / 100]
    isempty(inds) && throw(ArgumentError("no indicator taxa were found at the selected threshold value."))
    sort!(inds; by = r -> r.value)
    return (clusters = LabeledVector(smp, labels_), k = k, avg_silhouette = avg, silhouette_by_cluster = bycl,
            indval = ivm, indicators = inds, method = string(method))
end
