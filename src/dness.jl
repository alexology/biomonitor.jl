# Taxonomic distinctness (Clarke & Warwick 1998, 1999).

"""
    dness(x::AggregatedData; complete = false, tax_lev = :Genus,
          method = "delta", taxa_tree = false)

Taxonomic diversity and distinctness indices:

- `"delta"`: taxonomic diversity Δ;
- `"delta.st"`: taxonomic distinctness Δ*;
- `"delta.bin"`: Δ computed on presence-absence data (Δ+).

Path lengths between taxa are counted on the levels Phylum, Class, Order,
Family, Genus and Species down to `tax_lev`; levels with missing names and
levels that do not discriminate are dropped, and distances are scaled so
that the largest is 100 (as vegan's `taxa2dist(varstep = FALSE)`). With
`complete = false` an error is raised if some taxa are unassigned at
`tax_lev`; with `complete = true` they are excluded. `taxa_tree = true` also
returns the taxonomic table used.
"""
function dness(x::AggregatedData; complete::Bool = false, tax_lev = :Genus, method = "delta", taxa_tree::Bool = false)
    _check_aggregated(x)
    lev = _sym(tax_lev)
    lev === :Taxa && throw(ArgumentError("tax_lev cannot be equal to Taxa"))
    df = x[lev]
    (length(df.taxa) == 1 && df.taxa[1] == "unassigned") && throw(ArgumentError("The selected taxonomic level is empty"))
    used = (:Phylum, :Class, :Order, :Family, :Genus, :Species)
    lev in used || throw(ArgumentError("tax_lev must be one of $(join(used, ", "))"))
    taxpos = collect(used[1:findfirst(==(lev), used)])
    cols = [rank_index(l) for l in taxpos]
    hasun = "unassigned" in df.taxa
    (!complete && hasun) && throw(ArgumentError("Unassigned taxon at the desired taxonomic level"))
    hasun && (df = _drop_unassigned(df))

    d = x.data
    refc = unique([d.ranks[i, cols] for i in eachindex(d.taxa)])
    li = length(cols)
    # merge by tax_lev: sorted by the tax_lev name
    rows = Vector{Vector{String}}()
    abus = Vector{Vector{Float64}}()
    pos = Dict(n => i for (i, n) in enumerate(df.taxa))
    for r in sort(refc; by = r -> r[li])
        haskey(pos, r[li]) || continue
        push!(rows, r)
        push!(abus, df.abundance[pos[r[li]], :])
    end
    T = permutedims(reduce(hcat, rows))
    S = permutedims(reduce(hcat, abus))
    names = [string.(taxpos)...]
    # drop levels with empty names
    full = [all(!isempty, T[:, k]) for k in axes(T, 2)]
    if !all(full)
        T = T[:, full]
        names = names[full]
    end
    # drop non-discriminating levels (same name for all rows or all distinct)
    k = findfirst(==(string(lev)), names)
    k === nothing && throw(ArgumentError("Reference database has not enough taxonomic levels to perform the analysis"))
    keepc = trues(length(names))
    for c in 1:k-1
        u = length(unique(T[:, c]))
        (u == 1 || u == size(T, 1)) && (keepc[c] = false)
    end
    T = T[:, keepc]
    names = names[keepc]
    size(T, 2) == 1 && throw(ArgumentError("Reference database has not enough taxonomic levels to perform the analysis"))

    n = size(T, 1)
    W = [Float64(count(c -> T[i, c] != T[j, c], axes(T, 2))) for i in 1:n, j in 1:n]
    W ./= maximum(W)
    W .*= 100
    m = string(method)
    res = map(axes(S, 2)) do j
        v = S[:, j]
        m == "delta.bin" && (v = Float64.(v .> 0))
        num = 0.0
        den = 0.0
        for jj in 1:n-1, ii in jj+1:n
            a = v[ii] * v[jj]
            num += a * W[ii, jj]
            den += a
        end
        if m == "delta" || m == "delta.bin"
            num * 2 / (sum(v) * (sum(v) - 1))
        elseif m == "delta.st"
            num / den
        else
            throw(ArgumentError("method must be delta, delta.st or delta.bin"))
        end
    end
    out = LabeledVector(df.samples, res)
    taxa_tree || return out
    return (results = out, taxa_tree = (levels = names, taxa = T, abundance = LabeledMatrix(T[:, end], df.samples, S)))
end
