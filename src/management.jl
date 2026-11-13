# Abundance, richness and data management.

"""
    abundance(x::AggregatedData; tax_lev = :Taxa, unassigned = false) -> LabeledVector

Total abundance per sample at `tax_lev` (number of taxa for
presence-absence data). With `unassigned = false` the taxa without a name
at `tax_lev` are excluded.
"""
function abundance(x::AggregatedData; tax_lev = :Taxa, unassigned::Bool = false)
    t = _level(_check_aggregated(x), tax_lev)
    unassigned || (t = _drop_unassigned(t))
    return LabeledVector(t.samples, vec(sum(t.abundance; dims = 1)))
end

# total abundance at Taxa level (R helper `abu`)
_abu(x::AggregatedData) = vec(sum(x[:Taxa].abundance; dims = 1))

"""
    richness(x::AggregatedData; tax_lev = :Taxa) -> LabeledVector

Number of taxa per sample at `tax_lev`. A warning is issued when some taxa
are unassigned at that level, because richness can then be biased.
"""
function richness(x::AggregatedData; tax_lev = :Taxa, warn::Bool = true)
    t = _check_aggregated(x)[tax_lev]
    if "unassigned" in t.taxa
        warn && @warn "Unassigned taxa at $(tax_lev) level, richness calculation could be biased"
        t = _drop_unassigned(t)
    end
    return LabeledVector(t.samples, vec(sum(t.abundance .> 0; dims = 1)) .* 1.0)
end

"""
    allrich(x::AggregatedData) -> LabeledMatrix

Richness at every taxonomic level (samples × levels).
"""
function allrich(x::AggregatedData)
    _check_aggregated(x)
    res = reduce(hcat, [richness(x; tax_lev = lev, warn = false).values for lev in TAXONOMIC_LEVELS])
    return LabeledMatrix(samples(x), [string.(TAXONOMIC_LEVELS)...], reshape(res, length(samples(x)), :))
end

"""
    general_info(x::AggregatedData; taxalist = false)

Number of taxa at each taxonomic level in the whole dataset (only taxa
present in at least one sample), plus the total abundance when the data are
abundances. With `taxalist = true` returns a `Dict` level => taxon names.
"""
function general_info(x::AggregatedData; taxalist::Bool = false)
    d = _check_aggregated(x).data
    present = vec(sum(d.abundance; dims = 2)) .> 0
    M = hcat(d.ranks, d.taxa)[present, :]
    lists = Dict{Symbol,Vector{String}}()
    for (k, lev) in enumerate(TAXONOMIC_LEVELS)
        lists[lev] = filter(!isempty, unique(M[:, k]))
    end
    taxalist && return lists
    labs = [string.(TAXONOMIC_LEVELS)...]
    vals = Float64[length(lists[lev]) for lev in TAXONOMIC_LEVELS]
    if any(>(1), d.abundance)
        push!(labs, "Abundance")
        push!(vals, sum(d.abundance))
    end
    return LabeledVector(labs, vals)
end

function _tree_names(x::AggregatedData)
    v = unique(vec(x.data.ranks))
    return filter(!isempty, v)
end

"""
    get_taxa_abundance(x::AggregatedData, taxa; rel = false) -> LabeledVector

Abundance of the given taxa (any taxonomic level) per sample; with
`rel = true` relative to the total abundance. Not available for
presence-absence data.
"""
function get_taxa_abundance(x::AggregatedData, taxa = nothing; rel::Bool = false)
    _check_aggregated(x)
    is_bin(x) && throw(ArgumentError("This function cannot be applied to presence-absence data."))
    (taxa === nothing || (taxa isa AbstractString && isempty(taxa)) ||
     (taxa isa AbstractVector && length(taxa) == 1 && isempty(taxa[1]))) && throw(ArgumentError("Please provide a taxon name"))
    tx = taxa isa AbstractString ? [String(taxa)] : String.(collect(taxa))
    allnames = Set(_tree_names(x))
    missing_taxa = [t for t in tx if !(t in allnames)]
    if !isempty(missing_taxa)
        tx = [t for t in tx if t in allnames]
        isempty(tx) && throw(ArgumentError("None of the taxa provided were found in the database"))
        @info "The following taxa were not found in the database and have been excluded: $(join(missing_taxa, ", "))"
    end
    cols = Int[]
    for t in tx
        for (_, c) in positions(x.data.ranks, t)
            c in cols || push!(cols, c)
        end
    end
    tot = zeros(length(samples(x)))
    want = Set(tx)
    for c in cols
        t = x[RANKS[c]]
        for (i, n) in enumerate(t.taxa)
            n in want && (tot .+= t.abundance[i, :])
        end
    end
    rel && (tot ./= _abu(x))
    return LabeledVector(samples(x), tot)
end

"""
    get_taxa_richness(x::AggregatedData, taxa, tax_lev) -> LabeledVector

Number of taxa at `tax_lev` belonging to the given taxa (e.g. the number of
Ephemeroptera families). `tax_lev` is a single level or one level per
taxon; it cannot be coarser than the level of the taxa.
"""
function get_taxa_richness(x::AggregatedData, taxa = nothing, tax_lev = nothing)
    _check_aggregated(x)
    (taxa === nothing || tax_lev === nothing) && throw(ArgumentError("Please provide a taxon name and/or a taxonomic level"))
    tx = taxa isa AbstractString ? [String(taxa)] : String.(collect(taxa))
    tl = tax_lev isa Union{AbstractString,Symbol} ? [string(tax_lev)] : string.(collect(tax_lev))
    length(tx) > 1 && length(tl) == 1 && (tl = fill(tl[1], length(tx)))
    length(tx) != length(tl) && throw(ArgumentError("tax_lev must be of the same length of taxa"))
    tx = capwords.(tx)
    tl = capwords.(tl)
    DF = hcat(x.data.ranks, x.data.taxa)
    for t in tx
        isempty(positions(DF, t)) && throw(ArgumentError("Please provide a valid taxon name. Names provided can also be absent in your database."))
    end
    out = zeros(length(samples(x)))
    for (t, lev) in zip(tx, tl)
        li = rank_index(lev)
        pos = positions(DF, t)
        rcols = unique(c for (_, c) in positions(x.data.ranks, t))
        if !isempty(rcols) && any(>(li), rcols)
            throw(ArgumentError("Taxonomic level of taxa cannot be lower than taxonomic level of tax_lev"))
        end
        sel = Set(DF[r, li] for (r, _) in pos)
        tab = x[Symbol(lev)]
        for (i, n) in enumerate(tab.taxa)
            (n in sel && n != "unassigned") || continue
            out .+= tab.abundance[i, :] .> 0
        end
    end
    return LabeledVector(samples(x), out)
end

"""
    remove_taxa(x::AggregatedData, taxa) -> CommunityTable

The original data without the given taxa (and their children), ready to be
imported again with [`as_biomonitor`](@ref).
"""
function remove_taxa(x::AggregatedData, taxa = nothing)
    _check_aggregated(x)
    (taxa === nothing || (taxa isa AbstractString && isempty(taxa))) && throw(ArgumentError("Please provide at least taxon name"))
    tx = taxa isa AbstractString ? [String(taxa)] : String.(collect(taxa))
    allnames = Set(_tree_names(x))
    missing_taxa = [t for t in tx if !(t in allnames)]
    if !isempty(missing_taxa)
        tx = [t for t in tx if t in allnames]
        isempty(tx) && throw(ArgumentError("None of the taxa provided were found"))
        @info "The following taxa were not found and have been excluded: $(join(missing_taxa, ", "))"
    end
    drop = Set{Int}()
    for t in tx
        for (r, _) in positions(x.data.ranks, t)
            push!(drop, r)
        end
    end
    keep = [i for i in eachindex(x.data.taxa) if !(i in drop)]
    return CommunityTable(x.data.taxa[keep], samples(x), x.data.abundance[keep, :])
end

"""
    subset_taxa(x::BiomonitorData; taxa = nothing, exclude = nothing) -> BiomonitorData

Keep only the rows containing `taxa` at any level, and/or drop those
containing `exclude` (R method `subset.asb`).
"""
function subset_taxa(x::BiomonitorData; taxa = nothing, exclude = nothing)
    (taxa === nothing && exclude === nothing) && throw(ArgumentError("taxa and exclude cannot be both null"))
    tolist(v) = v === nothing ? String[] : v isa AbstractString ? [String(strip(v))] : String.(strip.(collect(v)))
    tk = tolist(taxa)
    ex = tolist(exclude)
    DF = hcat(x.ranks, x.taxa)
    allnames = Set(filter(!isempty, unique(vec(DF))))
    notfound = [t for t in [tk; ex] if !(t in allnames)]
    if !isempty(notfound)
        tk = [t for t in tk if t in allnames]
        ex = [t for t in ex if t in allnames]
        (isempty(tk) && isempty(ex)) && throw(ArgumentError("None of the taxa provided were found in the database"))
        @info "The following taxa were not found in the database and have been excluded: $(join(notfound, ", "))"
    end
    rowsof(names) = unique(r for t in names for (r, _) in positions(DF, t))
    keep = taxa === nothing ? collect(eachindex(x.taxa)) : sort(rowsof(tk))
    exclude === nothing || (keep = setdiff(keep, rowsof(ex)))
    return BiomonitorData(x.ranks[keep, :], x.taxa[keep], x.samples, x.abundance[keep, :], x.group,
                          x.presence_absence, x.corrected_names, x.suggested_names)
end

"""
    combine_taxa(x::AggregatedData; ntaxa = 2, tax_lev = :Taxa) -> CommunityTable

Abundance of every combination of `ntaxa` taxa, labelled by joining the
names with `_` (e.g. `"Baetis_Caenis"`). Useful to find indicator
combinations of taxa.
"""
function combine_taxa(x::AggregatedData; ntaxa::Int = 2, tax_lev = :Taxa)
    _check_aggregated(x)
    t = x[tax_lev]
    if is_bin(x)
        @warn "combine_taxa with presence-absence data can be meaningless"
        t = to_bin(t)
    end
    t = _drop_unassigned(t)
    combs = combinations(length(t.taxa), ntaxa)
    labs = [join(t.taxa[c], "_") for c in combs]
    A = zeros(length(combs), length(t.samples))
    for (k, c) in enumerate(combs)
        A[k, :] = vec(sum(t.abundance[c, :]; dims = 1))
    end
    return CommunityTable(labs, t.samples, A)
end

"""
    convert_to_vegan(x::AggregatedData; tax_lev = :Family) -> LabeledMatrix

Samples × taxa matrix at `tax_lev` (the layout used by community ecology
packages), without unassigned taxa.
"""
function convert_to_vegan(x::AggregatedData; tax_lev = :Family)
    _check_aggregated(x)
    _sym(tax_lev) in TAXONOMIC_LEVELS || throw(ArgumentError("Please provide a valid tax_lev"))
    t = _drop_unassigned(_level(x, tax_lev))
    return LabeledMatrix(t.samples, t.taxa, permutedims(t.abundance))
end

"""
    convert_to_biotic(x::AggregatedData) -> CommunityTable

Family-level table for biotic indices, where Oligochaeta are kept as a
single taxon (as required by many scoring systems).
"""
function convert_to_biotic(x::AggregatedData)
    _check_aggregated(x)
    tree = x.data.ranks
    pos = positions(tree, "Oligochaeta")
    fam = x[:Family]
    DF = if isempty(pos)
        fam
    else
        rows = unique(r for (r, _) in pos)
        lev = RANKS[unique(c for (_, c) in pos)[1]]
        famcol = tree[rows, rank_index(:Family)]
        if all(!isempty, famcol)
            fam
        else
            orow = x[lev]
            # aggregate_taxa groups the same tree column, so the name is present
            oi = findfirst(==("Oligochaeta"), orow.taxa)
            oi === nothing && throw(ArgumentError("inconsistent AggregatedData: Oligochaeta missing from the $(lev) table"))
            rem = Set(filter(!isempty, famcol))
            keep = [i for (i, n) in enumerate(fam.taxa) if !(n in rem) && n != "Oligochaeta"]
            CommunityTable([fam.taxa[keep]; "Oligochaeta"], fam.samples,
                           vcat(fam.abundance[keep, :], orow.abundance[oi:oi, :]))
        end
    end
    is_bin(x) && (DF = to_bin(DF))
    DF = _drop_unassigned(DF)
    o = sortperm(DF.taxa)
    return CommunityTable(DF.taxa[o], DF.samples, DF.abundance[o, :])
end

"""
    solve_ambiguous(x::AggregatedData; method = "MCWP") -> CommunityTable

Solve ambiguous (parent-child) taxa, e.g. *Baetis* and *Baetis rhodani*
in the same dataset:

- `"RPKC"`: remove the parents and keep the children;
- `"MCWP"`: merge the children with the parent.

Returns a table at `Taxa` level that can be imported again with
[`as_biomonitor`](@ref).
"""
function solve_ambiguous(x::AggregatedData; method = "MCWP")
    _check_aggregated(x)
    d = x.data
    finest = [something(findlast(!isempty, d.ranks[i, :]), 0) for i in eachindex(d.taxa)]
    o = sortperm(finest; rev = true, alg = MergeSort)
    tree = d.ranks[o, :]
    taxa = d.taxa[o]
    A = d.abundance[o, :]
    acc = String[]
    m = uppercase(string(method))
    if m == "RPKC"
        taxaset = Set(taxa)
        for i in eachindex(taxa)
            temp = filter(!isempty, tree[i, :])
            temp = temp[1:end-1]
            for t in temp
                t in taxaset && push!(acc, t)
            end
        end
        accset = Set(acc)
        tt = x[:Taxa]
        keep = [i for (i, n) in enumerate(tt.taxa) if !(n in accset)]
        names, G = group_rows(tt.taxa[keep], tt.abundance[keep, :])
    elseif m == "MCWP"
        cur = copy(taxa)
        for i in eachindex(cur)
            cur[i] in acc && continue
            ri = findfirst(==(cur[i]), d.taxa[o])
            ri === nothing && continue
            nz = findall(!isempty, tree[ri, :])
            temp = tree[ri, nz]
            levs = nz[1:end-1]
            temp = temp[1:end-1]
            present = Set(cur)
            k = findfirst(in(present), temp)
            if k === nothing
                push!(acc, cur[i])
            else
                up = temp[k]
                lev = levs[k]
                for r in eachindex(cur)
                    tree[r, lev] == up && (cur[r] = up)
                end
                push!(acc, up)
            end
        end
        names, G = group_rows(cur, A)
    else
        throw(ArgumentError("method must be RPKC or MCWP"))
    end
    res = CommunityTable(names, samples(x), G)
    return is_bin(x) ? to_bin(res) : res
end

"""
    abundance_classes(x::CommunityTable; abundance_classes = [0.1, 1, 10, 50],
                      class_labels = nothing) -> CommunityTable

Convert abundances (e.g. percent cover) into classes. With the default
breaks: (0, 0.1] → 1, (0.1, 1] → 2, (1, 10] → 3, (10, 50] → 4, > 50 → 5;
0 stays 0.
"""
function abundance_classes(x::CommunityTable; abundance_classes = [0.1, 1, 10, 50], class_labels = nothing)
    labs = class_labels === nothing ? collect(1.0:length(abundance_classes)+1) : Float64.(class_labels)
    length(labs) - 1 != length(abundance_classes) &&
        throw(ArgumentError("length of class_labels must be equal to the length of abundance_classes + 1"))
    br = [0.0; Float64.(abundance_classes); Inf]
    B = map(x.abundance) do a
        k = cut_class(a, br)
        k == 0 ? 0.0 : labs[k]
    end
    return CommunityTable(x.taxa, x.samples, B)
end
