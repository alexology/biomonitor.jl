# Biotic indices based on scores (BMWP, ASPT, WHPT, LIFE, PSI, E-PSI,
# DEHLI, IBMR) and other indices (IGOLD, EPT, EPTD, BIOCO).

# ---------------------------------------------------------------------------
# Built-in score tables
# ---------------------------------------------------------------------------

const _ASPT_FILES = Dict("ita" => "aspt_scores_fam_ita", "spa" => "aspt_scores_fam_spa",
                         "a" => "aspt_scores_fam_armitage", "armitage" => "aspt_scores_fam_armitage",
                         "uk" => "aspt_scores_fam_uk")
const _ASPT_RULES = Dict("ita" => "aspt_acc_fam_ita", "uk" => "aspt_acc_fam_uk")

_whpt_scores() = _cached(:whpt_scores) do
    t = _raw_dataset("whpt_scores_fam_uk")
    Dict((String(a), round(Int, parse_num(b))) => parse_num(c)
         for (a, b, c) in zip(column(t, "Taxon"), column(t, "ABUCLASS"), column(t, "Scores")))
end

function _class_scores(name, groupcol)
    _cached(Symbol("classscores_", name)) do
        t = _raw_dataset(name)
        Dict((round(Int, parse_num(a)), round(Int, parse_num(b))) => parse_num(c)
             for (a, b, c) in zip(column(t, groupcol), column(t, "ABUCLASS"), column(t, "SCORE")))
    end
end

_ibmr_scores() = _cached(:ibmr_scores) do
    t = _raw_dataset("ibmr_scores")
    Dict(String(a) => (csi = parse_num(b), ei = parse_num(c))
         for (a, b, c) in zip(column(t, "Taxon"), column(t, "Csi"), column(t, "Ei")))
end

"""
    show_scores(index, method) -> NamedTuple

Score tables used by a biotic index, as Julia dictionaries that can be
modified and passed back as custom `method`/`agg` arguments.

| index | methods | content |
|---|---|---|
| `"aspt"`, `"bmwp"` | `ita`, `spa`, `armitage`, `uk` | `scores` (taxon => score), `aggregation_rule` (taxon => taxon) |
| `"whpt"` | `uk` | `scores` ((taxon, abundance class) => score; class -1 is presence-only) |
| `"life"` | `extence`, `life_2017` | `fs_groups` (taxon => flow group), `scores` ((group, class) => score) |
| `"psi"` | `extence` | `fssr_groups`, `scores` ((group, class) => score) |
| `"epsi"` | `uk` | `scores` (taxon => sensitivity weight) |
| `"dehli"` | `chadd_2017` | `scores`, `aggregation_rule` |
| `"ibmr"` | `ita` | `scores` (taxon => (csi, ei)) |
"""
function show_scores(index, method)
    idx = lowercase(string(index))
    m = lowercase(string(method))
    if idx in ("aspt", "bmwp") && haskey(_ASPT_FILES, m)
        sc = copy(_score_dict(_ASPT_FILES[m]))
        haskey(_ASPT_RULES, m) && return (scores = sc, aggregation_rule = copy(_rule_dict(_ASPT_RULES[m])))
        return (scores = sc,)
    elseif idx == "whpt" && m == "uk"
        return (scores = copy(_whpt_scores()), aggregation_rule = copy(_rule_dict("whpt_acc_fam_uk")))
    elseif idx == "life" && m == "extence"
        return (fs_groups = _int_dict("life_scores_fam_extence", "FG_Score"), scores = copy(_class_scores("life_fs", "FS")),
                aggregation_rule = copy(_rule_dict("life_acc_fam_extence")))
    elseif idx == "life" && m == "life_2017"
        return (fs_groups = _int_dict("life_scores_fam_2017", "FG_Score"), scores = copy(_class_scores("life_fs", "FS")),
                aggregation_rule = copy(_rule_dict("life_acc_fam_2017")))
    elseif idx == "psi" && m == "extence"
        return (fssr_groups = _int_dict("psi_scores_fam_extence", "FSSR_Score"), scores = copy(_class_scores("psi_fssr", "FSSR")),
                aggregation_rule = copy(_rule_dict("epsi_acc_fam_uk")))
    elseif idx == "epsi" && m == "uk"
        return (scores = copy(_score_dict("epsi_scores_fam_uk")), aggregation_rule = copy(_rule_dict("epsi_acc_fam_uk")))
    elseif idx == "dehli" && m == "chadd_2017"
        return (scores = copy(_score_dict("dehli_scores_chadd_2017")), aggregation_rule = copy(_rule_dict("dehli_acc_chadd_2017")))
    elseif idx == "ibmr" && m == "ita"
        return (scores = copy(_ibmr_scores()),)
    end
    throw(ArgumentError("Method $(method) is not implemented into the index $(index) or index and methods names are wrong.\n" *
                        "Index must be one of: aspt, bmwp, whpt, psi, epsi, life, dehli, ibmr\n" *
                        "Method must be one of: ita, spa, armitage, uk, extence, life_2017, chadd_2017"))
end

_int_dict(name, col) = Dict(k => round(Int, v) for (k, v) in _score_dict(name, col))

# ---------------------------------------------------------------------------
# Shared steps
# ---------------------------------------------------------------------------

# all the taxonomic levels (Phylum ... Subspecies) stacked and summed by name
function _stack_levels(x::AggregatedData)
    names = String[]
    rows = Vector{Vector{Float64}}()
    for lev in RANKS
        t = x[lev]
        for i in eachindex(t.taxa)
            push!(names, t.taxa[i])
            push!(rows, t.abundance[i, :])
        end
    end
    A = isempty(rows) ? zeros(0, length(samples(x))) : permutedims(reduce(hcat, rows))
    return group_rows(names, A)
end

_keep(names, A, keepset) = (idx = [i for (i, n) in enumerate(names) if n in keepset]; (names[idx], A[idx, :]))

# R helper `checkBmwpFam`: rename composite taxa and aggregate again
function _apply_rules(names, A, rules::AbstractDict)
    any(n -> haskey(rules, n), names) || return names, A
    newnames = [String(get(rules, n, n)) for n in names]
    return group_rows(newnames, A)
end

# R helper `manage_inconsistencies`: when a taxon and one of its children are
# both scored, the abundance of the child is removed from the parent.
function _manage_inconsistencies(names::Vector{String}, A::Matrix{Float64}, tree::Matrix{String})
    entries = Tuple{Int,Int}[]
    for t in capwords.(names)
        entries = vcat(positions(tree, t), entries)
    end
    sort!(entries; by = e -> e[2], alg = MergeSort)
    counts = Dict{Int,Int}()
    for (r, _) in entries
        counts[r] = get(counts, r, 0) + 1
    end
    incon = [e for e in entries if counts[e[1]] > 1]
    pairs = NamedTuple{(:child, :parent),Tuple{String,String}}[]
    isempty(incon) && return names, A, pairs
    B = copy(A)
    pos = Dict(n => i for (i, n) in enumerate(names))
    stored = Set{String}()
    for (row, _) in incon
        temp = [e for e in incon if e[1] == row]
        col1 = temp[end-1][2]
        col2 = temp[end][2]
        parent = tree[row, col1]
        child = tree[row, col2]
        (child in stored && parent in stored) && continue
        if haskey(pos, parent) && haskey(pos, child)
            B[pos[parent], :] .-= B[pos[child], :]
        end
        push!(pairs, (child = child, parent = parent))
        push!(stored, child, parent)
    end
    @info "Parent-child pairs found, please see traceB for further information"
    return names, B, pairs
end

# R helper `manage_exceptions`: abundances of the exception taxa are removed
# from their scored parents and their scored children are set to zero.
function _manage_exceptions(names::Vector{String}, A::Matrix{Float64}, tree::Matrix{String}, scored, exceptions)
    exc = exceptions isa AbstractString ? [String(exceptions)] : String.(collect(exceptions))
    taxind = Tuple{Int,Int,String}[]
    for t in capwords.(exc)
        taxind = vcat([(r, c, t) for (r, c) in positions(tree, t)], taxind)
    end
    changes = NamedTuple{(:taxon, :operation, :affected),Tuple{String,String,String}}[]
    isempty(taxind) && return names, A, changes
    minrow = Dict{Tuple{Int,String},Int}()
    for (r, c, t) in taxind
        minrow[(c, t)] = min(get(minrow, (c, t), typemax(Int)), r)
    end
    ks = sort!(collect(keys(minrow)); by = k -> (k[2], k[1]))
    sort!(ks; by = k -> k[1], alg = MergeSort)
    seen = Set{String}()
    sel = Tuple{Int,Int,String}[]
    for (c, t) in ks
        t in seen && continue
        push!(seen, t)
        push!(sel, (minrow[(c, t)], c, t))
    end
    higher = [[h for h in tree[r, 1:c] if h in scored] for (r, c, _) in sel]
    lower = [[l for l in vec(tree[tree[:, c] .== t, c+1:NRANKS]) if l in scored] for (_, c, t) in sel]
    B = copy(A)
    pos = Dict(n => i for (i, n) in enumerate(names))
    for (i, (_, _, t)) in enumerate(sel)
        if !isempty(higher[i])
            rem = haskey(pos, t) ? B[pos[t], :] : zeros(size(B, 2))
            for h in unique(higher[i])
                haskey(pos, h) && (B[pos[h], :] .-= rem)
            end
        end
        for l in unique(lower[i])
            haskey(pos, l) && (B[pos[l], :] .= 0.0)
        end
    end
    for (i, (_, _, t)) in enumerate(sel)
        for h in higher[i]
            push!(changes, (taxon = t, operation = "removed from", affected = h))
        end
        for l in lower[i]
            push!(changes, (taxon = t, operation = "set to zero", affected = l))
        end
    end
    unique!(changes)
    filter!(c -> c.taxon != c.affected, changes)
    @info "Exceptions found, please see traceB for further information"
    return names, B, changes
end

_bin!(A) = (A .= Float64.(A .> 0); A)

_resolve_agg(agg, rulefile) = agg === true ? (rulefile === nothing ? nothing : _rule_dict(rulefile)) :
                              agg === false ? nothing :
                              agg isa AbstractDict ? Dict(String(capwords(k)) => String(capwords(v)) for (k, v) in agg) :
                              throw(ArgumentError("agg is not one of true, false or a Dict of aggregation rules"))

_custom_scores(d::AbstractDict) = Dict(String(capwords(k)) => Float64(v) for (k, v) in d)

# ---------------------------------------------------------------------------
# BMWP and ASPT
# ---------------------------------------------------------------------------

function _bmwp_aspt(x::AggregatedData, method, agg, exceptions, traceB, aspt::Bool)
    _check_aggregated(x)
    BIN = is_bin(x)
    if method isa AbstractDict
        (agg === false || agg isa AbstractDict) ||
            throw(ArgumentError("When method is a Dict agg needs to be false or a Dict containing the aggregation rules"))
        y = _custom_scores(method)
        z = agg === false ? nothing : _resolve_agg(agg, nothing)
    else
        m = lowercase(string(method))
        haskey(_ASPT_FILES, m) || throw(ArgumentError("Please provide a valid method"))
        (agg === true || agg === false) || throw(ArgumentError("When using the default method agg can only be true or false"))
        y = _score_dict(_ASPT_FILES[m])
        z = _resolve_agg(agg, get(_ASPT_RULES, m, nothing))
    end
    tree = ranks(x)
    names, A = _stack_levels(x)
    exce = nothing
    if exceptions !== nothing
        names, A, exce = _manage_exceptions(names, A, tree, Set(keys(y)), exceptions)
    end
    !aspt && BIN && _bin!(A)
    names, A = _keep(names, A, Set(keys(y)))
    composite = nothing
    if z !== nothing
        composite = [n for n in names if haskey(z, n)]
        names, A = _apply_rules(names, A, z)
    end
    names, A, incon = _manage_inconsistencies(names, A, tree)
    aspt && BIN && _bin!(A)
    names, A = _keep(names, A, Set(keys(y)))
    sc = [y[n] for n in names]
    P = A .> 0
    tot = vec(sum(sc .* P; dims = 1))
    res = aspt ? tot ./ vec(sum(P; dims = 1)) : tot
    out = LabeledVector(samples(x), res)
    traceB || return out
    return (results = out, taxa_df = (taxa = names, scores = sc, abundance = LabeledMatrix(names, samples(x), A)),
            composite_taxa = composite, exceptions = exce, parent_child_pairs = incon)
end

"""
    bmwp(x::AggregatedData; method = "ita", agg = false, exceptions = nothing,
         traceB = false)

Biological Monitoring Working Party score: sum of the scores of the scored
families present in each sample.

- `method`: `"ita"` (Buffagni et al. 2014), `"spa"` (Alba-Tercedor et al.
  2002), `"a"` (Armitage et al. 1983), `"uk"` (Davy-Bowker et al. 2008), or
  a custom `Dict` taxon => score (see [`show_scores`](@ref)).
- `agg`: `true` to merge composite families as in the chosen method
  (`ita`, `uk`), or a `Dict` taxon => taxon with custom rules.
- `exceptions`: taxa whose abundance is removed from their scored parents
  and whose scored children are excluded.
- `traceB = true` returns a `NamedTuple` with `results`, the scored taxa
  (`taxa_df`), the merged composite taxa, the exceptions and the
  parent-child pairs whose abundances were corrected.
"""
bmwp(x::AggregatedData; method = "ita", agg = false, exceptions = nothing, traceB::Bool = false) =
    _bmwp_aspt(x, method, agg, exceptions, traceB, false)

"""
    aspt(x::AggregatedData; method = "ita", agg = false, exceptions = nothing,
         traceB = false)

Average Score Per Taxon: [`bmwp`](@ref) divided by the number of scored
taxa. Same arguments as [`bmwp`](@ref).
"""
aspt(x::AggregatedData; method = "ita", agg = false, exceptions = nothing, traceB::Bool = false) =
    _bmwp_aspt(x, method, agg, exceptions, traceB, true)

"""
    dehli(x::AggregatedData; method = "chadd_2017", agg = true, exceptions = nothing,
          traceB = false)

Drought Effect of Habitat Loss on Invertebrates index (Chadd et al. 2017):
the average score of the scored taxa, computed like [`aspt`](@ref).
`method` can also be a custom `Dict` of scores.
"""
function dehli(x::AggregatedData; method = "chadd_2017", agg = true, exceptions = nothing, traceB::Bool = false)
    if !(method isa AbstractDict)
        string(method) == "chadd_2017" || throw(ArgumentError("Please provide a valid method"))
        method = _score_dict("dehli_scores_chadd_2017")
    end
    agg === true && (agg = _rule_dict("dehli_acc_chadd_2017"))
    return _bmwp_aspt(x, method, agg, exceptions, traceB, true)
end

# ---------------------------------------------------------------------------
# WHPT
# ---------------------------------------------------------------------------

"""
    whpt(x::AggregatedData; method = "uk", type = "ab", metric = "aspt",
         agg = false, abucl = [1, 9, 99, 999], exceptions = nothing, traceB = false)

Walley Hawkes Paisley Trigg index (UK). `type = "ab"` uses abundance-weighted
scores, `"po"` presence-only scores; `metric` is `"aspt"`, `"ntaxa"` or
`"bmwp"`. `abucl` are the lower bounds of the abundance classes. A custom
`method` is a `Dict` (taxon, abundance class) => score, with class -1 for
presence-only scores.
"""
function whpt(x::AggregatedData; method = "uk", type = "ab", metric = "aspt", agg = false,
              abucl = [1, 9, 99, 999], exceptions = nothing, traceB::Bool = false)
    _check_aggregated(x)
    BIN = is_bin(x)
    type = string(type)
    metric = string(metric)
    type in ("po", "ab") || throw(ArgumentError("Please provide a valid type: po or ab"))
    metric in ("aspt", "ntaxa", "bmwp") || throw(ArgumentError("Please provide a valid metric: aspt, ntaxa or bmwp"))
    if method isa AbstractDict
        (agg === false || agg isa AbstractDict) ||
            throw(ArgumentError("When method is a Dict agg needs to be false or a Dict containing the aggregation rules"))
        y = Dict((String(capwords(k[1])), Int(k[2])) => Float64(v) for (k, v) in method)
        z = agg === false ? nothing : _resolve_agg(agg, nothing)
    else
        string(method) == "uk" || throw(ArgumentError("Method need to be set to uk or a custom Dict"))
        (agg === true || agg === false) || throw(ArgumentError("When using the default method agg can only be true or false"))
        y = _whpt_scores()
        z = _resolve_agg(agg, "whpt_acc_fam_uk")
    end
    type == "ab" && BIN && throw(ArgumentError("Type cannot be set to abundance if presence-absence data are used"))
    ytaxa = Set(k[1] for k in keys(y))
    tree = ranks(x)
    names, A = _stack_levels(x)
    exce = nothing
    exceptions === nothing || ((names, A, exce) = _manage_exceptions(names, A, tree, ytaxa, exceptions))
    names, A = _keep(names, A, ytaxa)
    composite = nothing
    if z !== nothing
        composite = [n for n in names if haskey(z, n)]
        names, A = _apply_rules(names, A, z)
    end
    names, A, incon = _manage_inconsistencies(names, A, tree)
    BIN && _bin!(A)
    C = abundance_class_matrix(A, abucl)
    smp = samples(x)
    rows = NamedTuple{(:taxon, :sample, :abuclass, :score),Tuple{String,String,Int,Float64}}[]
    for j in eachindex(smp), i in eachindex(names)
        cl = Int(C[i, j])
        key = type == "po" ? (cl > 0 ? (names[i], -1) : nothing) : (names[i], cl)
        (key === nothing || !haskey(y, key)) && continue
        push!(rows, (taxon = names[i], sample = smp[j], abuclass = cl, score = y[key]))
    end
    res = map(smp) do s
        sc = [r.score for r in rows if r.sample == s]
        isempty(sc) && return NaN
        metric == "aspt" ? mean(sc) : metric == "ntaxa" ? Float64(length(sc)) : sum(sc)
    end
    out = LabeledVector(smp, res)
    traceB || return out
    return (results = out, taxa_df = LabeledMatrix(names, smp, A), abu_cl = LabeledMatrix(names, smp, C),
            whpt_df = rows, composite_taxa = composite, exceptions = exce, parent_child_pairs = incon)
end

# ---------------------------------------------------------------------------
# LIFE, PSI and E-PSI
# ---------------------------------------------------------------------------

function _class_index_setup(x, method, agg, groups_custom, default_methods)
    if method isa AbstractDict
        (agg === false || agg isa AbstractDict) ||
            throw(ArgumentError("When method is a Dict agg needs to be false or a Dict containing the aggregation rules"))
        groups_custom isa AbstractDict || throw(ArgumentError("class scores are needed when method is a custom Dict"))
        y = Dict(String(capwords(k)) => Int(v) for (k, v) in method)
        w = Dict((Int(k[1]), Int(k[2])) => Float64(v) for (k, v) in groups_custom)
        z = agg === false ? nothing : _resolve_agg(agg, nothing)
    else
        m = string(method)
        haskey(default_methods, m) || throw(ArgumentError("method is not one of $(join(keys(default_methods), ", ")) or a custom Dict"))
        (agg === true || agg === false) || throw(ArgumentError("When using the default method agg can only be true or false"))
        groups_custom === nothing || throw(ArgumentError("When using default methods the class scores can only be nothing"))
        gfile, gcol, wfile, wcol, rfile = default_methods[m]
        y = _int_dict(gfile, gcol)
        w = _class_scores(wfile, wcol)
        z = _resolve_agg(agg, rfile)
    end
    return y, w, z
end

function _prepare_scored(x::AggregatedData, y, z, exceptions)
    BIN = is_bin(x)
    tree = ranks(x)
    names, A = _stack_levels(x)
    exce = nothing
    exceptions === nothing || ((names, A, exce) = _manage_exceptions(names, A, tree, Set(keys(y)), exceptions))
    names, A = _keep(names, A, Set(keys(y)))
    composite = nothing
    if z !== nothing
        composite = [n for n in names if haskey(z, n)]
        names, A = _apply_rules(names, A, z)
    end
    before = (names, copy(A))
    names, A, incon = _manage_inconsistencies(names, A, tree)
    BIN && _bin!(A)
    names, A = _keep(names, A, Set(keys(y)))
    return names, A, composite, exce, incon, before
end

"""
    life(x::AggregatedData; method = "extence", abucl = [1, 9, 99, 999, 9999],
         agg = false, fs_scores = nothing, exceptions = nothing, traceB = false)

Lotic-invertebrate Index for Flow Evaluation (Extence et al. 1999).
`method` is `"extence"` (family scores of Extence et al. 1999),
`"life_2017"` (2017 family scores) or a custom `Dict` taxon => flow group;
in that case `fs_scores` must be a `Dict` (flow group, abundance class) =>
score. Abundance classes have lower bounds `abucl`.
"""
function life(x::AggregatedData; method = "extence", abucl = [1, 9, 99, 999, 9999], agg = false,
              fs_scores = nothing, exceptions = nothing, traceB::Bool = false)
    _check_aggregated(x)
    defaults = Dict("extence" => ("life_scores_fam_extence", "FG_Score", "life_fs", "FS", "life_acc_fam_extence"),
                    "life_2017" => ("life_scores_fam_2017", "FG_Score", "life_fs", "FS", "life_acc_fam_2017"))
    y, w, z = _class_index_setup(x, method, agg, fs_scores, defaults)
    names, A, composite, exce, incon, _ = _prepare_scored(x, y, z, exceptions)
    C = abundance_class_matrix(A, abucl)
    smp = samples(x)
    S = fill(NaN, length(names), length(smp))
    for j in eachindex(smp), i in eachindex(names)
        key = (y[names[i]], Int(C[i, j]))
        haskey(w, key) && (S[i, j] = w[key])
    end
    res = [(v = filter(!isnan, S[:, j]); isempty(v) ? NaN : mean(v)) for j in eachindex(smp)]
    out = LabeledVector(smp, res)
    traceB || return out
    return (results = out, taxa_df = LabeledMatrix(names, smp, A), abu_cl = LabeledMatrix(names, smp, C),
            life_df = LabeledMatrix(names, smp, replace(S, NaN => 0.0)),
            composite_taxa = composite, exceptions = exce, parent_child_pairs = incon)
end

"""
    psi(x::AggregatedData; method = "extence", abucl = [1, 9, 99, 999], agg = false,
        fssr_scores = nothing, exceptions = nothing, traceB = false)

Proportion of Sediment-sensitive Invertebrates (Extence et al. 2013): the
percentage of the scores of sensitive taxa (FSSR groups 1 and 2) over the
scores of all scored taxa. A custom `method` is a `Dict` taxon => FSSR
group with `fssr_scores` a `Dict` (group, abundance class) => score.

Samples without sensitive taxa get 0 (R biomonitoR returns `NA` when some,
but not all, samples lack sensitive taxa); samples without scored taxa get
`NaN`.
"""
function psi(x::AggregatedData; method = "extence", abucl = [1, 9, 99, 999], agg = false,
             fssr_scores = nothing, exceptions = nothing, traceB::Bool = false)
    _check_aggregated(x)
    defaults = Dict("extence" => ("psi_scores_fam_extence", "FSSR_Score", "psi_fssr", "FSSR", "epsi_acc_fam_uk"))
    y, w, z = _class_index_setup(x, method, agg, fssr_scores, defaults)
    names, A, composite, exce, incon, before = _prepare_scored(x, y, z, exceptions)
    C = abundance_class_matrix(A, abucl)
    smp = samples(x)
    S = zeros(Float64, length(names), length(smp))
    M = falses(length(names), length(smp))
    for j in eachindex(smp), i in eachindex(names)
        key = (y[names[i]], Int(C[i, j]))
        if haskey(w, key)
            S[i, j] = w[key]
            M[i, j] = true
        end
    end
    sens = [y[n] in (1, 2) for n in names]
    res = map(eachindex(smp)) do j
        any(M[:, j]) || return NaN
        100 * sum(S[sens, j]) / sum(S[:, j])
    end
    out = LabeledVector(smp, res)
    traceB || return out
    return (results = out, taxa_df = LabeledMatrix(before[1], smp, before[2]), abu_cl = LabeledMatrix(names, smp, C),
            psi_df = LabeledMatrix(names, smp, S), composite_taxa = composite, exceptions = exce,
            parent_child_pairs = incon)
end

"""
    epsi(x::AggregatedData; method = "uk", agg = false, abucl = [1, 9, 99, 999],
         exceptions = nothing, traceB = false)

Empirically-weighted Proportion of Sediment-sensitive Invertebrates (Turley
et al. 2016): percentage of the abundance-class-weighted scores of sensitive
taxa (score ≥ 0.5) over all scored taxa. A custom `method` is a `Dict`
taxon => score.
"""
function epsi(x::AggregatedData; method = "uk", agg = false, abucl = [1, 9, 99, 999],
              exceptions = nothing, traceB::Bool = false)
    _check_aggregated(x)
    if method isa AbstractDict
        y = _custom_scores(method)
        z = agg === false ? nothing : _resolve_agg(agg, nothing)
    else
        string(method) == "uk" || throw(ArgumentError("Please provide a valid method"))
        (agg === true || agg === false) || throw(ArgumentError("When using the default method agg can only be true or false"))
        y = _score_dict("epsi_scores_fam_uk")
        z = _resolve_agg(agg, "epsi_acc_fam_uk")
    end
    names, A, composite, exce, incon, _ = _prepare_scored(x, y, z, exceptions)
    C = abundance_class_matrix(A, abucl)
    sc = [y[n] for n in names]
    E = sc .* C
    sensrows = sc .>= 0.5
    res = 100 .* vec(sum(E[sensrows, :]; dims = 1)) ./ vec(sum(E; dims = 1))
    out = LabeledVector(samples(x), res)
    traceB || return out
    return (results = out, taxa_df = LabeledMatrix(names, samples(x), A), epsi_df = (taxa = names, scores = sc,
            classes = LabeledMatrix(names, samples(x), C)), composite_taxa = composite, exceptions = exce,
            parent_child_pairs = incon)
end

# ---------------------------------------------------------------------------
# IBMR
# ---------------------------------------------------------------------------

"""
    ibmr(x::AggregatedData; method = "ita", exceptions = nothing,
         coverage_coef = false, coverage_classes = [0.1, 1, 10, 50], traceB = false)

Indice Biologique Macrophytique en Rivière (AFNOR 2003), Italian scores.
Data are percent covers converted to coverage coefficients 1-5 with
`coverage_classes`, unless `coverage_coef = true` (data already
coefficients). A custom `method` is a `Dict` taxon => (csi, ei), i.e.
specific score and stenoecy coefficient.
"""
function ibmr(x::AggregatedData; method = "ita", exceptions = nothing, coverage_coef::Bool = false,
              coverage_classes = [0.1, 1, 10, 50], traceB::Bool = false)
    _check_aggregated(x)
    y = method isa AbstractDict ? Dict(String(capwords(k)) => (csi = Float64(v[1]), ei = Float64(v[2])) for (k, v) in method) :
        _ibmr_scores()
    tree = ranks(x)
    names, A = _stack_levels(x)
    exce = nothing
    exceptions === nothing || ((names, A, exce) = _manage_exceptions(names, A, tree, Set(keys(y)), exceptions))
    names, A = _keep(names, A, Set(keys(y)))
    names, A, incon = _manage_inconsistencies(names, A, tree)
    all(a -> a in (0, 1, 2, 3, 4), A) ||
        @info "Data different from 0, 1, 2, 3 and 4 detected. Probably they are not coverage coefficients."
    if !coverage_coef
        A = abundance_classes(CommunityTable(names, samples(x), A); abundance_classes = coverage_classes).abundance
    end
    csi = [y[n].csi for n in names]
    ei = [y[n].ei for n in names]
    num = vec(sum(A .* (csi .* ei); dims = 1))
    den = vec(sum(A .* ei; dims = 1))
    out = LabeledVector(samples(x), num ./ den)
    traceB || return out
    return (results = out, taxa_df = (taxa = names, csi = csi, ei = ei, coverage = LabeledMatrix(names, samples(x), A)),
            exceptions = exce, parent_child_pairs = incon)
end

# ---------------------------------------------------------------------------
# Other indices
# ---------------------------------------------------------------------------

"""
    igold(x::AggregatedData; traceB = false)

1 - relative abundance of Gastropoda, Oligochaeta and Diptera (Pinto et al.
2004).
"""
function igold(x::AggregatedData; traceB::Bool = false)
    _check_aggregated(x)
    d = x.data
    ci, si, oi = rank_index(:Class), rank_index(:Subclass), rank_index(:Order)
    rows = [i for i in eachindex(d.taxa) if d.ranks[i, ci] == "Gastropoda" || d.ranks[i, si] == "Oligochaeta" || d.ranks[i, oi] == "Diptera"]
    G = d.abundance[rows, :]
    is_bin(x) && (G = Float64.(G .> 0))
    res = 1 .- vec(sum(G; dims = 1)) ./ abundance(x; tax_lev = :Taxa, unassigned = true).values
    out = LabeledVector(samples(x), res)
    traceB || return out
    return (results = out, taxa_df = CommunityTable(d.taxa[rows], samples(x), d.abundance[rows, :]))
end

"""
    ept(x::AggregatedData; tax_lev = :Taxa) -> LabeledVector

Number of Ephemeroptera, Plecoptera and Trichoptera taxa at `tax_lev`
(finer than Order).
"""
function ept(x::AggregatedData; tax_lev = :Taxa)
    _check_aggregated(x)
    li = rank_index(tax_lev)
    li <= 4 && throw(ArgumentError("Taxonomic level cannot be equal or higher than Order"))
    d = x.data
    oi = rank_index(:Order)
    rows = [i for i in eachindex(d.taxa) if d.ranks[i, oi] in ("Plecoptera", "Ephemeroptera", "Trichoptera")]
    isempty(rows) && return LabeledVector(samples(x), zeros(length(samples(x))))
    keys = li == length(TAXONOMIC_LEVELS) ? d.taxa[rows] : d.ranks[rows, li]
    keys = [isempty(k) ? "unassigned" : k for k in keys]
    names, G = group_rows(keys, d.abundance[rows, :])
    keep = names .!= "unassigned"
    return LabeledVector(samples(x), vec(sum(G[keep, :] .> 0; dims = 1)) .* 1.0)
end

const EPTD_FAMILIES = ["Heptageniidae", "Ephemeridae", "Leptophlebiidae", "Brachycentridae", "Goeridae",
                       "Polycentropodidae", "Limnephilidae", "Odontoceridae", "Dolichopodidae", "Stratiomyidae",
                       "Dixidae", "Empididae", "Athericidae", "Nemouridae"]

"""
    eptd(x::AggregatedData; base = 10, eptd_families = nothing, traceB = false)

``\\log_{base}(1 + n)``, where n is the abundance of the EPTD families
(Buffagni et al. 2004): Heptageniidae, Ephemeridae, Leptophlebiidae,
Brachycentridae, Goeridae, Polycentropodidae, Limnephilidae, Odontoceridae,
Dolichopodidae, Stratiomyidae, Dixidae, Empididae, Athericidae, Nemouridae.
"""
function eptd(x::AggregatedData; base = 10, eptd_families = nothing, traceB::Bool = false)
    _check_aggregated(x)
    fams = eptd_families === nothing ? EPTD_FAMILIES : capwords.(string.(collect(eptd_families)))
    t = _level(x, :Family)
    rows = [i for (i, n) in enumerate(t.taxa) if n in fams]
    res = log.(base, vec(sum(t.abundance[rows, :]; dims = 1)) .+ 1)
    out = LabeledVector(t.samples, res)
    traceB || return out
    return (results = out, taxa_df = CommunityTable(t.taxa[rows], t.samples, t.abundance[rows, :]))
end

"""
    bioco(x::AggregatedData; alien, dfref, digits = 2) -> LabeledMatrix

Biocontamination indices (Arbačiauskas et al. 2008): abundance
contamination (`aci`), richness contamination (`rci`, at family level) and
site-specific biocontamination (`sbci`). `alien` are the alien taxa (any
level); `dfref` is the reference used to import the data (`:mi`, `:mf`,
`:fi` or a [`ReferenceDB`](@ref)). For presence-absence data `aci` and
`sbci` are `NaN`.
"""
function bioco(x::AggregatedData; alien = nothing, dfref = nothing, digits::Int = 2)
    _check_aggregated(x)
    alien === nothing && throw(ArgumentError("Please provide a vector containing the names of alien taxa"))
    dfref === nothing && throw(ArgumentError("Please set dfref"))
    if is_custom(x) && !(dfref isa ReferenceDB)
        throw(ArgumentError("Please provide the dfref you used for as_biomonitor"))
    end
    ref = dfref isa ReferenceDB ? dfref : reference_db(dfref)
    al = capwords.(alien isa AbstractString ? [String(alien)] : String.(collect(alien)))
    smp = samples(x)
    vec_ = Set(vec(ref.ranks))
    absent = [a for a in al if !(a in vec_)]
    if !isempty(absent)
        if length(absent) == length(al)
            @warn "$(join(absent, ", ")) are not part of the taxonomic tree of the reference database"
            return LabeledMatrix(smp, ["aci", "rci", "sbci"], zeros(length(smp), 3))
        end
        al = [a for a in al if a in vec_]
        @warn "The following taxa are absent from the reference database: $(join(absent, ", "))"
    end
    taxind = Tuple{Int,Int}[]
    for a in al
        taxind = vcat(positions(ref.ranks, a), taxind)
    end
    alienall = String[]
    for (r, c) in taxind
        append!(alienall, ref.ranks[r, c:NRANKS])
    end
    alienall = Set(filter(!isempty, unique(alienall)))
    d = x.data
    rows = [i for i in eachindex(d.taxa) if d.taxa[i] in alienall]
    if isempty(rows)
        v = is_bin(x) ? fill(NaN, length(smp)) : zeros(length(smp))
        return LabeledMatrix(smp, ["aci", "rci", "sbci"], hcat(v, v, v))
    end
    names, G = group_rows(d.ranks[rows, rank_index(:Family)], d.abundance[rows, :])
    abu_alien = vec(sum(G; dims = 1))
    tax_alien = vec(sum(G .> 0; dims = 1))
    aci = round.(abu_alien ./ abundance(x; tax_lev = :Taxa, unassigned = true).values; digits)
    rci = round.(tax_alien ./ richness(x; tax_lev = :Family, warn = false).values; digits)
    br = [0, 0.01, 0.1, 0.2, 0.5, 1]
    cls(v) = (k = cut_class(v, br; include_lowest = true); k == 0 ? NaN : Float64(k - 1))
    sbci = [max(cls(a), cls(r)) for (a, r) in zip(aci, rci)]
    if is_bin(x)
        aci = fill(NaN, length(smp))
        sbci = fill(NaN, length(smp))
    end
    return LabeledMatrix(smp, ["aci", "rci", "sbci"], hcat(aci, rci, sbci))
end
