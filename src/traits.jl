# Trait assignment and trait-based indices.

# level of each name in the reference (finest level first, R `ref_long`)
function _name_levels(ref::ReferenceDB)
    lev = Dict{String,Int}()
    for c in NRANKS:-1:1
        for v in view(ref.ranks, :, c)
            isempty(v) || haskey(lev, v) || (lev[v] = c)
        end
    end
    return lev
end

"""
    assign_traits(x::AggregatedData; trait_db = nothing, group = :mi, tax_lev = :Taxa,
                  dfref = nothing, filter_by_distance = nothing) -> AssignedTraits

Assign traits to the taxa of `x`, also borrowing them across taxonomic
levels (from Subspecies to Family). For instance a genus present in the
data receives the traits of its species in the trait database, and a
species missing from the trait database receives the traits of its genus.
When traits exist at the matched level, traits of other taxa at that level
are discarded (e.g. *Anabolia lombarda* gets the traits of *Anabolia*, not
those of *Anabolia nervosa*).

- `trait_db`: a [`TraitTable`](@ref); the default is
  [`traits_tachet`](@ref) (macroinvertebrates).
- `dfref`: reference database to place the trait taxa in the taxonomy
  (default: the built-in one for `group`).
- `tax_lev`: work at a coarser level (`:Family` at most).
- `filter_by_distance`: `:pos` (traits from equal or finer taxa), `:neg`
  (equal or coarser) or an integer (maximum absolute distance).

Use [`manage_traits`](@ref) to select among multiple assignments and
[`average_traits`](@ref) or [`sample_traits`](@ref) to obtain one row per
taxon.
"""
function assign_traits(x::AggregatedData; trait_db::Union{Nothing,TraitTable} = nothing, group = :mi,
                       tax_lev = :Taxa, dfref::Union{Nothing,ReferenceDB} = nothing, filter_by_distance = nothing)
    _check_aggregated(x)
    if trait_db === nothing
        is_custom(x) && @warn "It seems that you used your own reference database. Please check the consistency of the taxonomy used for calculating the index with those of your reference database to have reliable results."
        tdb = traits_tachet()
    else
        tdb = TraitTable(String.(strip.(trait_db.taxa)), trait_db.traits, trait_db.values)
    end
    ref = dfref === nothing ? reference_db(group) : dfref

    # trait taxa placed in the reference taxonomy
    tidx = Dict{String,Vector{Int}}()
    for (i, t) in enumerate(tdb.taxa)
        push!(get!(tidx, t, Int[]), i)
    end
    refrows = Tuple{Int,Int}[]   # (reference row, trait row)
    for (r, t) in enumerate(ref.taxa)
        haskey(tidx, t) || continue
        for k in tidx[t]
            push!(refrows, (r, k))
        end
    end
    # expanded trait table: one row per (level name, trait taxon)
    exp_name = String[]
    exp_real = String[]
    exp_row = Int[]
    fam = rank_index(:Family)
    for lev in NRANKS:-1:fam
        for (r, k) in refrows
            nm = ref.ranks[r, lev]
            isempty(nm) && continue
            fr = findlast(!isempty, view(ref.ranks, r, fam:NRANKS))
            push!(exp_name, nm)
            push!(exp_real, fr === nothing ? "" : ref.ranks[r, fam+fr-1])
            push!(exp_row, k)
        end
    end
    byname = Dict{String,Vector{Int}}()
    for (e, nm) in enumerate(exp_name)
        push!(get!(byname, nm, Int[]), e)
    end

    lev = _sym(tax_lev)
    lev in (:Family, :Subfamily, :Tribus, :Genus, :Species, :Subspecies, :Taxa) ||
        throw(ArgumentError("Maximum taxonomic level is family."))
    d = x.data
    DFranks = copy(d.ranks)
    DFtaxa = copy(d.taxa)
    if lev !== :Taxa
        li = rank_index(lev)
        DFranks[:, li+1:end] .= ""
        DFtaxa = DFranks[:, li]
        seen = Set{Vector{String}}()
        keep = Int[]
        for i in eachindex(DFtaxa)
            key = [DFranks[i, :]; DFtaxa[i]]
            key in seen || (push!(seen, key); push!(keep, i))
        end
        DFranks = DFranks[keep, :]
        DFtaxa = DFtaxa[keep]
    end

    rows = Tuple{String,String,String,Int}[]   # (user taxon, trait taxon, matched name, trait row)
    for i in eachindex(DFtaxa)
        t = DFtaxa[i]
        isempty(t) && continue
        same = findall(==(t), DFtaxa)
        tree = filter(!isempty, vec(DFranks[same, :]))
        for nm in reverse(tree)
            haskey(byname, nm) || continue
            for e in byname[nm]
                push!(rows, (t, exp_real[e], exp_name[e], exp_row[e]))
            end
            break
        end
    end
    # when traits exist at the matched level, drop the other taxa of that level
    keep = trues(length(rows))
    for u in unique(r[1] for r in rows)
        idx = findall(r -> r[1] == u, rows)
        s = [(r[1] == r[2]) + (r[1] == r[3]) + (r[2] == r[3]) for r in rows[idx]]
        if any(==(0), s) && sum(s) != 0
            keep[idx[s .== 0]] .= false
        end
    end
    rows = rows[keep]

    levels = _name_levels(ref)
    rows = [r for r in rows if haskey(levels, r[1]) && haskey(levels, r[2])]
    ul = [levels[r[1]] for r in rows]
    tl = [levels[r[2]] for r in rows]
    dist = tl .- ul
    V = isempty(rows) ? zeros(0, length(tdb.traits)) : permutedims(reduce(hcat, [tdb.values[r[4], :] for r in rows]))
    res = AssignedTraits([r[1] for r in rows], [RANKS[l] for l in ul], [r[2] for r in rows], [RANKS[l] for l in tl],
                         dist, copy(tdb.traits), V)
    filter_by_distance === nothing && return res
    if filter_by_distance isa Union{Symbol,AbstractString}
        f = string(filter_by_distance)
        f == "pos" && return _subset(res, findall(>=(0), res.distance))
        f == "neg" && return _subset(res, findall(<=(0), res.distance))
        throw(ArgumentError("pos, neg or an integer are needed when filter_by_distance is not nothing"))
    end
    return _subset(res, findall(v -> abs(v) <= abs(filter_by_distance), res.distance))
end

"""
    manage_traits(x::AssignedTraits; method = "nearest+-", traceB = false) -> AssignedTraits

Select, for each taxon, the traits assigned from the taxonomically nearest
taxa:

- `"nearest"`: smallest distance, irrespective of the direction;
- `"nearest+"`: nearest among taxa of equal or finer level;
- `"nearest-"`: nearest among taxa of equal or coarser level;
- `"nearest+-"`: nearest giving priority to equal or finer taxa;
- `"nearest-+"`: nearest giving priority to equal or coarser taxa.

With `traceB = true` also returns the taxa excluded by the selection.

For `"nearest"` the R implementation keeps the minimum signed distance
(i.e. coarser taxa first); this behaviour is kept for compatibility.
"""
function manage_traits(x::AssignedTraits; method = "nearest+-", traceB::Bool = false)
    m = string(method)
    groups = sort(unique(x.taxa))
    keep = Int[]
    for g in groups
        idx = findall(==(g), x.taxa)
        dv = x.distance[idx]
        sel = if m == "nearest"
            idx[dv .== minimum(dv)]
        elseif m == "nearest+"
            p = idx[dv .>= 0]
            isempty(p) ? Int[] : p[x.distance[p] .== minimum(x.distance[p])]
        elseif m == "nearest-"
            q = idx[dv .<= 0]
            isempty(q) ? Int[] : q[x.distance[q] .== maximum(x.distance[q])]
        elseif m == "nearest+-"
            p = idx[dv .>= 0]
            if !isempty(p)
                p[x.distance[p] .== minimum(x.distance[p])]
            else
                idx[dv .== maximum(dv)]
            end
        elseif m == "nearest-+"
            q = idx[dv .<= 0]
            if !isempty(q)
                q[x.distance[q] .== maximum(x.distance[q])]
            else
                idx[dv .== minimum(dv)]
            end
        else
            throw(ArgumentError("method must be one of nearest, nearest+, nearest-, nearest+- or nearest-+"))
        end
        append!(keep, sel)
    end
    isempty(keep) && m in ("nearest+", "nearest-") &&
        throw(ArgumentError(m == "nearest+" ? "No traits assigned at lower level" : "No traits assigned at higher level"))
    res = _subset(x, keep)
    traceB || return res
    excl = setdiff(unique(x.taxa), unique(res.taxa))
    return (results = res, taxa_excluded = excl)
end

"""
    average_traits(x::AssignedTraits; col_blocks = nothing, type = "F") -> TraitTable

Average the traits assigned to the same taxon. For fuzzy traits
(`type = "F"`) the averages are rescaled with [`prep_fuzzy`](@ref) using
`col_blocks` (default: the blocks of [`traits_tachet`](@ref)); for
continuous traits (`type = "C"`) the plain means are returned. Missing
values count as 0.
"""
function average_traits(x::AssignedTraits; col_blocks = nothing, type = "F")
    string(type) in ("F", "C") || throw(ArgumentError("type must be one of F or C"))
    V = replace(x.values, NaN => 0.0)
    names, M = group_rows(x.taxa, V; f = :mean)
    if string(type) == "F"
        blocks = col_blocks === nothing ? DEFAULT_TACHET_BLOCKS : col_blocks
        P, _, keep = _prep_fuzzy(M, blocks)
        return TraitTable(names, x.traits[keep], P)
    end
    return TraitTable(names, x.traits, M)
end

"""
    sample_traits(x::AssignedTraits; rng = Random.default_rng()) -> TraitTable

For each taxon, pick at random one of the assigned trait rows.
"""
function sample_traits(x::AssignedTraits; rng::AbstractRNG = Random.default_rng())
    u = unique(x.taxa)
    rows = [rand(rng, findall(==(t), x.taxa)) for t in u]
    return TraitTable(u, x.traits, x.values[rows, :])
end

"""
    add_bias_to_traits(x; fuzzy = true, col_blocks = nothing, SD = 0.001,
                       rng = Random.default_rng())

Add a small positive noise (|N(0, SD)|) to the trait values of a
[`TraitTable`](@ref) or [`AssignedTraits`](@ref), so that taxa with
identical traits become distinguishable. With `fuzzy = true` the values are
rescaled with [`prep_fuzzy`](@ref) (`col_blocks` required).
"""
function add_bias_to_traits(x::Union{TraitTable,AssignedTraits}; fuzzy::Bool = true, col_blocks = nothing,
                            SD::Real = 0.001, rng::AbstractRNG = Random.default_rng())
    fuzzy && col_blocks === nothing && throw(ArgumentError("Please set col_blocks"))
    !fuzzy && col_blocks !== nothing && @warn "col_blocks will be ignored because fuzzy = false"
    V = x.values .+ abs.(SD .* randn(rng, size(x.values)))
    traits = x.traits
    if fuzzy
        V, _, keep = _prep_fuzzy(V, col_blocks)
        traits = traits[keep]
    end
    x isa TraitTable && return TraitTable(x.taxa, traits, V)
    return AssignedTraits(x.taxa, x.taxa_level, x.trait_taxa, x.trait_level, x.distance, traits, V)
end

function _trait_cols(t::TraitTable, names)
    idx = Int[]
    for n in names
        k = findfirst(==(String(n)), t.traits)
        k === nothing && throw(ArgumentError("trait $(n) not found"))
        push!(idx, k)
    end
    return idx
end

_apply_trans(trans, A) = trans === nothing ? A : map(trans, A)

"""
    fuzzy_trait_ratio(x::AggregatedData; tax_lev = :Taxa, trait_db, type = "ab",
                      numerator, denominator, trans = nothing) -> LabeledVector

Community-level ratio between the summed affinities of the `numerator`
modalities and of the `denominator` modalities, weighted by abundance
(`type = "ab"`, optionally transformed elementwise by `trans`, e.g. `log1p`)
or presence (`type = "po"`).
"""
function fuzzy_trait_ratio(x::AggregatedData; tax_lev = :Taxa, trait_db::TraitTable, type = "ab",
                           numerator, denominator, trans = nothing)
    taxo = x[tax_lev]
    common = sort(intersect(taxo.taxa, trait_db.taxa))
    ti = [findfirst(==(c), taxo.taxa) for c in common]
    ri = [findfirst(==(c), trait_db.taxa) for c in common]
    num = vec(sum(trait_db.values[ri, _trait_cols(trait_db, numerator)]; dims = 2))
    den = vec(sum(trait_db.values[ri, _trait_cols(trait_db, denominator)]; dims = 2))
    A = taxo.abundance[ti, :]
    string(type) == "ab" && (A = _apply_trans(trans, A))
    (string(type) == "po" || is_bin(x)) && (A = Float64.(A .> 0))
    return LabeledVector(taxo.samples, vec(sum(A .* num; dims = 1)) ./ vec(sum(A .* den; dims = 1)))
end

"""
    flow_t(x::AggregatedData; tax_lev = :Taxa, type = "ab", trans = nothing) -> LabeledVector

Flow index based on the current-velocity preferences of Tachet et al.
(2010): affinity for medium and fast currents (CURRENT_3, CURRENT_4) over
all current modalities.
"""
function flow_t(x::AggregatedData; tax_lev = :Taxa, type = "ab", trans = nothing)
    ts = average_traits(assign_traits(x; tax_lev))
    return fuzzy_trait_ratio(x; tax_lev, trait_db = ts, type, trans,
                             numerator = ["CURRENT_3", "CURRENT_4"],
                             denominator = ["CURRENT_1", "CURRENT_2", "CURRENT_3", "CURRENT_4"])
end

# a trait name or a collection of names as Vector{String}
_names_vector(x::Union{AbstractString,Symbol}) = [String(x)]
_names_vector(x) = String.(collect(x))

const _BIO_INDICES = Dict{String,Tuple{Vector{String},Vector{String},Vector{Int},Vector{String}}}(
    "trophic_t" => (["TROPHIC_1"], ["TROPHIC_1", "TROPHIC_2", "TROPHIC_3"], [3], ["TROPHIC_1", "TROPHIC_2", "TROPHIC_3"]),
    "saprobity_t" => (["SAPROBITY_1", "SAPROBITY_2"], ["SAPROBITY_$(i)" for i in 1:5], [5], ["SAPROBITY_$(i)" for i in 1:5]),
    "temperature_t" => (["TEMPERATURE_1"], ["TEMPERATURE_1", "TEMPERATURE_2", "TEMPERATURE_3"], [3],
                        ["TEMPERATURE_1", "TEMPERATURE_2", "TEMPERATURE_3"]))

"""
    bio_index(x::AggregatedData; tax_lev = :Taxa, trait_db = nothing, type = "ab",
              trans = nothing, index = ["trophic_t"], rm_taxa = true, group = :mi,
              num = nothing, den = nothing, col_blocks = nothing, col_order = nothing) -> NamedTuple

Trait-based indices computed with [`fuzzy_trait_ratio`](@ref):

- `"trophic_t"`: oligotrophic affinity (TROPHIC_1) over all trophic modalities;
- `"saprobity_t"`: xeno- and oligosaprobic (SAPROBITY_1, _2) over all;
- `"temperature_t"`: psychrophilic (TEMPERATURE_1) over all;
- `"custom"`: `num` over `den`, using the traits in `col_order` grouped in
  blocks of sizes `col_blocks`.

Traits are assigned with [`assign_traits`](@ref) and averaged. With
`rm_taxa = true` taxa with an empty block are dropped. Returns a
`NamedTuple` index => [`LabeledVector`](@ref).
"""
function bio_index(x::AggregatedData; tax_lev = :Taxa, trait_db = nothing, type = "ab", trans = nothing,
                   index = ["trophic_t"], rm_taxa::Bool = true, group = :mi, num = nothing, den = nothing,
                   col_blocks = nothing, col_order = nothing)
    idx = index isa Union{AbstractString,Symbol} ? [string(index)] : string.(collect(index))
    valid = ["trophic_t", "saprobity_t", "temperature_t", "custom"]
    all(in(valid), idx) || throw(ArgumentError("Invalid 'index' value(s). Must be one or more of: $(join(valid, ", "))"))
    if any(!=("custom"), idx) && (num !== nothing || den !== nothing)
        throw(ArgumentError("You can set numerator and denominator ONLY if all indices are 'custom'"))
    end
    custom = nothing
    if "custom" in idx
        (num === nothing || den === nothing) && throw(ArgumentError("For custom index, both numerator and denominator must be provided"))
        col_order === nothing && throw(ArgumentError("For custom index, please provide a col_order"))
        col_blocks === nothing && throw(ArgumentError("For custom index, please provide a col_blocks"))
        custom_blocks = col_blocks isa Integer ? [Int(col_blocks)] : Int.(collect(col_blocks))
        custom_order = _names_vector(col_order)
        _check_blocks(custom_blocks, length(custom_order); what = "col_blocks", of = "col_order")
        custom = (_names_vector(num), _names_vector(den), custom_blocks, custom_order)
    end
    data_ts = trait_db === nothing ? assign_traits(x; tax_lev, group) : assign_traits(x; tax_lev, trait_db)
    results = map(idx) do id
        n, d, blocks, order = id == "custom" ? something(custom) : _BIO_INDICES[id]
        ci = [findfirst(==(c), data_ts.traits) for c in order]
        any(isnothing, ci) && throw(ArgumentError("some traits of col_order are not in the trait database"))
        V = data_ts.values[:, ci]
        starts = cumsum([1; blocks[1:end-1]])
        ends = cumsum(blocks)
        nz = [vec(sum(V[:, s:e]; dims = 2)) .> 0 for (s, e) in zip(starts, ends)]
        rowsok = rm_taxa ? reduce(.&, nz) : reduce(.|, nz)
        sub = AssignedTraits(data_ts.taxa[rowsok], data_ts.taxa_level[rowsok], data_ts.trait_taxa[rowsok],
                             data_ts.trait_level[rowsok], data_ts.distance[rowsok], String.(order), V[rowsok, :])
        av = average_traits(sub; col_blocks = blocks)
        av = TraitTable(av.taxa, av.traits, replace(av.values, NaN => 0.0))
        fuzzy_trait_ratio(x; tax_lev, trait_db = av, type, numerator = n, denominator = d, trans)
    end
    return NamedTuple{Tuple(Symbol.(idx))}(Tuple(results))
end

function _abundances_for_traits(x, tax_lev, trait_db)
    t = _drop_unassigned(_level(x, tax_lev))
    common = sort(intersect(t.taxa, trait_db.taxa))
    ti = [findfirst(==(c), t.taxa) for c in common]
    ri = [findfirst(==(c), trait_db.taxa) for c in common]
    return t, common, ti, ri
end

"""
    cwm(x::AggregatedData; trait_db = traits_tachet(), tax_lev = :Taxa,
        trans = log1p, traceB = false) -> LabeledMatrix

Community-weighted mean of each trait modality (samples × modalities), with
weights `trans(abundance)`. Only taxa present in `trait_db` are used.
"""
function cwm(x::AggregatedData; trait_db::Union{Nothing,TraitTable} = nothing, tax_lev = :Taxa,
             trans = log1p, traceB::Bool = false)
    _check_aggregated(x)
    tdb = trait_db === nothing ? traits_tachet() : trait_db
    _sym(tax_lev) in (:Family, :Genus, :Species, :Taxa) ||
        throw(ArgumentError("tax_lev should be one of the following: Family, Genus, Species or Taxa"))
    t, common, ti, ri = _abundances_for_traits(x, tax_lev, tdb)
    W = map(trans, t.abundance[ti, :])
    tot = vec(sum(W; dims = 1))
    T = tdb.values[ri, :]
    res = zeros(length(t.samples), length(tdb.traits))
    for j in eachindex(t.samples), k in eachindex(tdb.traits)
        s = 0.0
        for i in eachindex(common)
            v = W[i, j] * T[i, k] / tot[j]
            isnan(v) || (s += v)
        end
        res[j, k] = s
    end
    out = LabeledMatrix(t.samples, tdb.traits, res)
    traceB || return out
    return (results = out, traits = TraitTable(common, tdb.traits, T),
            taxa = CommunityTable(common, t.samples, t.abundance[ti, :]),
            taxa_excluded = setdiff(t.taxa, common))
end

function _tsi_matrix(tdb::TraitTable, rows)
    groups = String[]
    for tr in tdb.traits
        g = split(tr, "_")[1]
        g in groups || push!(groups, g)
    end
    gid = [findfirst(==(split(tr, "_")[1]), groups) for tr in tdb.traits]
    M = zeros(Float64, length(rows), length(groups))
    for (k, g) in enumerate(groups)
        cols = findall(==(k), gid)
        kk = length(cols)
        for (a, r) in enumerate(rows)
            M[a, k] = (sum(abs2, tdb.values[r, cols]) - 1 / kk) / (1 - 1 / kk)
        end
    end
    return groups, M
end

"""
    tsi(x::AggregatedData; trait_db = traits_tachet(), tax_lev = :Taxa) -> LabeledMatrix

Taxon Specialisation Index of each taxon and trait (taxa × traits)
(Mondy & Usseglio-Polatera 2014): ``(\\sum a_k^2 - 1/K) / (1 - 1/K)`` for
the K modalities of a trait. Trait names are the modality prefixes
before `_`.
"""
function tsi(x::AggregatedData; trait_db::Union{Nothing,TraitTable} = nothing, tax_lev = :Taxa)
    _check_aggregated(x)
    tdb = trait_db === nothing ? traits_tachet() : trait_db
    t = _level(x, tax_lev)
    common = sort(intersect(t.taxa, tdb.taxa))
    ri = [findfirst(==(c), tdb.taxa) for c in common]
    groups, M = _tsi_matrix(tdb, ri)
    return LabeledMatrix(common, groups, M)
end

"""
    csi(x::AggregatedData; trait_db = traits_tachet(), tax_lev = :Taxa,
        trans = log1p, traceB = false) -> LabeledMatrix

Community Specialisation Index (samples × traits): the mean
[`tsi`](@ref) of the taxa, weighted by `trans(abundance)`.
"""
function csi(x::AggregatedData; trait_db::Union{Nothing,TraitTable} = nothing, tax_lev = :Taxa,
             trans = log1p, traceB::Bool = false)
    _check_aggregated(x)
    tdb = trait_db === nothing ? traits_tachet() : trait_db
    t, common, ti, ri = _abundances_for_traits(x, tax_lev, tdb)
    groups, M = _tsi_matrix(tdb, ri)
    W = map(trans, t.abundance[ti, :])
    tot = vec(sum(W; dims = 1))
    res = zeros(length(t.samples), length(groups))
    for j in eachindex(t.samples), k in eachindex(groups)
        s = 0.0
        for i in eachindex(common)
            v = W[i, j] * M[i, k] / tot[j]
            isnan(v) || (s += v)
        end
        res[j, k] = s
    end
    out = LabeledMatrix(t.samples, groups, res)
    traceB || return out
    return (results = out, traits = TraitTable(common, tdb.traits, tdb.values[ri, :]),
            taxa = CommunityTable(common, t.samples, t.abundance[ti, :]),
            taxa_excluded = setdiff(t.taxa, common))
end
