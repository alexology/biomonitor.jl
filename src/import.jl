# Import of community data, aggregation at taxonomic levels and reference
# databases built from user trees.

_to_change_pairs(tc::AbstractDict) = [String(k) => String(v) for (k, v) in tc]
_to_change_pairs(tc::AbstractVector{<:Pair}) = [String(first(p)) => String(last(p)) for p in tc]
function _to_change_pairs(tc)
    throw(ArgumentError("to_change needs to be nothing, :default or a collection of Pairs taxon => correct taxon"))
end

function _default_to_change()
    t = _raw_dataset("to_change_mi")
    return [String(a) => String(b) for (a, b) in zip(column(t, "Taxon"), column(t, "Correct_Taxon"))]
end

"""
    as_biomonitor(x::CommunityTable; group = :mi, dfref = nothing,
                  to_change = :default, FUN = sum, correct_names = false) -> BiomonitorData

Merge the community table `x` with a taxonomic reference database.

- `group`: `:mi` macroinvertebrates, `:mf` macrophytes, `:fi` fish,
  `:di` diatoms (ignored when `dfref` is given).
- `dfref`: a custom [`ReferenceDB`](@ref).
- `to_change`: names to replace before the check, as `Pair`s
  `"Taxon" => "Correct_Taxon"` (a `Dict` or a vector of pairs). With
  `:default`, Hydracarina, Hydracnidia and Acariformes become
  Trombidiformes. `nothing` disables the replacement.
- `FUN`: `sum` for abundances, [`bin`](@ref) for presence-absence data.
- `correct_names`: `false` removes the taxa not found in the reference
  (spelling suggestions are stored, see [`excluded_taxa`](@ref)); a `Dict`
  `wrong => correct` renames them; `true` asks interactively on the
  terminal for each unknown name.

Taxon names are capitalised (first letter upper case), rows with the same
name are aggregated with `FUN` and missing values are set to 0.
"""
function as_biomonitor(x::CommunityTable; group = :mi, dfref::Union{Nothing,ReferenceDB} = nothing,
                       to_change = :default, FUN = sum, correct_names = false)
    fun = _funsym(FUN)
    A = copy(x.abundance)
    check_na = any(isnan, A)
    A[isnan.(A)] .= 0.0

    tc = if to_change === nothing
        nothing
    elseif to_change === :default || to_change == "default"
        _default_to_change()
    else
        p = _to_change_pairs(to_change)
        isempty(p) && throw(ArgumentError("to_change must have at least one entry"))
        length(unique(first.(p))) == length(p) ||
            throw(ArgumentError("the same name cannot be present twice in the Taxon column of to_change"))
        p
    end

    check_pa = false
    isint = all(a -> a % 1 == 0, A)
    if !any(>(1), A) && isint && fun !== :bin
        @warn "Presence-absence data detected but FUN is not set to bin. Is it this what you want?"
        check_pa = true
    end
    if any(>(1), A) && isint && fun !== :sum
        @warn "Abundance data detected but FUN is not set to sum. Is it this what you want?"
    end
    isint || @warn "Decimal numbers detected. Please check carefully which FUN to use."

    ref = dfref === nothing ? reference_db(group) : dfref
    grp = dfref === nothing ? _sym(group) : :custom

    names = [capwords(t) for t in x.taxa]

    corrected = Pair{String,String}[]
    if tc !== nothing
        tcc = [capwords(a) => capwords(b) for (a, b) in tc]
        present = Set(names)
        for (a, b) in tcc
            a in present && push!(corrected, a => b)
        end
        lookup = Dict(tcc)
        names = [get(lookup, t, t) for t in names]
    end

    refset = Set(ref.taxa)
    wrong = String[]
    for t in names
        (t in refset || t in wrong) || push!(wrong, t)
    end
    suggestions = [w => (s = suggest_names(w, ref.taxa); isempty(s) ? ["none"] : s) for w in wrong]

    keep = trues(length(names))
    if correct_names === false
        keep .= .!in.(names, Ref(Set(wrong)))
    else
        mapping = correct_names === true ? _ask_corrections(suggestions) : Dict(String(k) => String(v) for (k, v) in correct_names)
        for (i, t) in enumerate(names)
            if t in wrong
                if haskey(mapping, t)
                    newname = capwords(mapping[t])
                    push!(corrected, t => newname)
                    names[i] = newname
                else
                    keep[i] = false
                end
            end
        end
        unique!(corrected)
    end

    gnames, G = group_rows(names[keep], A[keep, :]; f = fun)
    if check_pa
        G = Float64.(G .> 0)
        @info "data imported as presence absence"
    end
    check_na && @info "NA detected, transformed to 0"

    refidx = Dict{String,Int}()
    dups = Set{String}()
    for (i, t) in enumerate(ref.taxa)
        haskey(refidx, t) ? push!(dups, t) : (refidx[t] = i)
    end
    rows = [i for (i, t) in enumerate(gnames) if haskey(refidx, t)]
    dupused = [gnames[i] for i in rows if gnames[i] in dups]
    isempty(dupused) || @warn "Taxa present more than once in the reference database (the first occurrence is used): $(join(dupused, ", "))"
    R = Matrix{String}(undef, length(rows), NRANKS)
    for (k, i) in enumerate(rows)
        R[k, :] = ref.ranks[refidx[gnames[i]], :]
    end

    if !isempty(wrong) && correct_names === false
        @info "Some taxa were excluded, check them with excluded_taxa()"
    end

    return BiomonitorData(R, gnames[rows], copy(x.samples), G[rows, :], grp, fun === :bin,
                          corrected, correct_names === false ? suggestions : Pair{String,Vector{String}}[])
end

function _ask_corrections(suggestions)
    mapping = Dict{String,String}()
    for (wrongname, choices) in suggestions
        opts = [filter(!=("none"), choices); "exit"]
        println("\n", wrongname)
        for (k, o) in enumerate(opts)
            println("  ", k, ": ", o)
        end
        print("Selection: ")
        ans = readline()
        k = tryparse(Int, strip(ans))
        (k === nothing || k < 1 || k > length(opts)) && throw(ArgumentError("invalid selection"))
        opts[k] == "exit" && throw(InterruptException())
        mapping[wrongname] = opts[k]
    end
    return mapping
end

"""
    excluded_taxa(x::BiomonitorData) -> Vector{Pair{String,Vector{String}}}

Taxa removed by [`as_biomonitor`](@ref) because they are not in the
reference database, with spelling suggestions (`"none"` when there are
none). This is the R `suggested_taxa_names` element of `traceB = TRUE`.
"""
excluded_taxa(x::BiomonitorData) = x.suggested_names
excluded_taxa(x::AggregatedData) = excluded_taxa(x.data)

"""
    corrected_names(x::BiomonitorData) -> Vector{Pair{String,String}}

Names changed by [`as_biomonitor`](@ref) (`wrong => correct`), through
`to_change` or `correct_names`.
"""
corrected_names(x::BiomonitorData) = x.corrected_names
corrected_names(x::AggregatedData) = corrected_names(x.data)

"""
    to_bin(x)

Presence-absence version of a [`CommunityTable`](@ref) (or of a
[`BiomonitorData`](@ref)): positive values become 1.
"""
to_bin(x::CommunityTable) = CommunityTable(x.taxa, x.samples, Float64.(x.abundance .> 0))
to_bin(x::BiomonitorData) = BiomonitorData(x.ranks, x.taxa, x.samples, Float64.(x.abundance .> 0),
                                           x.group, x.presence_absence, x.corrected_names, x.suggested_names)

"""
    aggregate_taxa(x::BiomonitorData; FUN = sum) -> AggregatedData

Aggregate the data at every taxonomic level ([`TAXONOMIC_LEVELS`](@ref)).
Taxa without a name at a level are pooled in the row `"unassigned"`. This
step is needed by all the index functions.
"""
function aggregate_taxa(x::BiomonitorData; FUN = sum)
    fun = _funsym(FUN)
    levels = Dict{Symbol,CommunityTable}()
    for (k, lev) in enumerate(TAXONOMIC_LEVELS)
        keys = lev === :Taxa ? x.taxa : x.ranks[:, k]
        names, G = group_rows(keys, x.abundance; f = fun)
        names = [isempty(n) ? "unassigned" : n for n in names]
        levels[lev] = CommunityTable(names, x.samples, G)
    end
    return AggregatedData(x, levels)
end

function _check_aggregated(x)
    x isa AggregatedData || throw(ArgumentError("The object you provided is not an object created with aggregate_taxa"))
    return x
end

# level table converted to presence-absence for bin objects
function _level(x::AggregatedData, tax_lev)
    t = x[tax_lev]
    return is_bin(x) ? to_bin(t) : t
end

function _drop_unassigned(t::CommunityTable)
    keep = t.taxa .!= "unassigned"
    all(keep) && return t
    return CommunityTable(t.taxa[keep], t.samples, t.abundance[keep, :])
end

# ---------------------------------------------------------------------------
# quick_rename
# ---------------------------------------------------------------------------

"""
    quick_rename(x::CommunityTable; group = :mi, corrections = nothing,
                 interactive = false) -> CommunityTable

Fix misspelled taxon names without importing the data. Names not found in
the reference database of `group` are replaced according to `corrections`
(a `Dict` `wrong => correct`; map a name to `""` to delete the row), or by
asking on the terminal when `interactive = true`. Unknown names without a
correction are kept unchanged. Underscores are converted to spaces.
"""
function quick_rename(x::CommunityTable; group = :mi, corrections = nothing, interactive::Bool = false)
    ref = reference_db(group)
    refset = Set(ref.taxa)
    wrongidx = [i for (i, t) in enumerate(x.taxa) if !(capwords(t) in refset)]
    if isempty(wrongidx)
        @info "All names are correct"
        return x
    end
    mapping = Dict{String,String}()
    if corrections !== nothing
        for (k, v) in corrections
            mapping[String(k)] = String(v)
        end
    elseif interactive
        for i in wrongidx
            w = x.taxa[i]
            opts = [suggest_names(capwords(w), ref.taxa); w; "User taxon input"; "Delete taxon"; "exit"]
            println("\n", w)
            for (k, o) in enumerate(opts)
                println("  ", k, ": ", o)
            end
            print("Selection: ")
            k = tryparse(Int, strip(readline()))
            (k === nothing || k < 1 || k > length(opts)) && throw(ArgumentError("invalid selection"))
            choice = opts[k]
            choice == "exit" && throw(InterruptException())
            if choice == "User taxon input"
                print("Enter taxon name: ")
                choice = strip(readline())
            elseif choice == "Delete taxon"
                choice = ""
            end
            mapping[w] = choice
        end
    else
        return x
    end
    names = String[]
    rows = Int[]
    for (i, t) in enumerate(x.taxa)
        newname = get(mapping, t, t)
        isempty(newname) && continue
        push!(names, replace(newname, "_" => " "))
        push!(rows, i)
    end
    return CommunityTable(names, x.samples, x.abundance[rows, :])
end

# ---------------------------------------------------------------------------
# Reference database from a user taxonomic tree
# ---------------------------------------------------------------------------

"""
    ref_from_tree(levels, values; group = :none) -> ReferenceDB

Build a reference database from a taxonomic table. `levels` are the column
names (a subset of `Phylum`, `Class`, `Subclass`, `Order`, `Family`,
`Subfamily`, `Tribus`, `Genus`, `Species`, `Subspecies`, `Taxa`, in any
case) and `values` the matrix of names (empty string for missing). Every
name becomes a row of the database, with its parent ranks. With
`group = :mi` or `:mf` the built-in reference is appended (user rows take
precedence).

```julia
ref = ref_from_tree(example_tree()...)
```
"""
function ref_from_tree(levels::AbstractVector, values::AbstractMatrix; group = :none)
    size(values, 2) == 1 && throw(ArgumentError("tables with 1 column are not allowed!"))
    size(values, 2) == length(levels) || throw(DimensionMismatch("levels and columns of values differ"))
    names = [capwords(string(l)) for l in levels]
    canon = [string.(TAXONOMIC_LEVELS)...]
    all(in(canon), names) || throw(ArgumentError("Provide valid column names"))
    V = [v === missing || v === nothing ? "" : String(strip(string(v))) for v in values]
    # remove duplicated rows
    seen = Set{Vector{String}}()
    rows = Int[]
    for i in axes(V, 1)
        r = V[i, :]
        r in seen || (push!(seen, r); push!(rows, i))
    end
    V = V[rows, :]
    order = [findfirst(==(c), names) for c in canon if c in names]
    V = V[:, order]
    names = names[order]

    out = Vector{Vector{String}}()
    for c in eachindex(names)
        sub = V[:, 1:c]
        seenc = Set{Vector{String}}()
        for i in axes(sub, 1)
            isempty(sub[i, c]) && continue
            r = sub[i, :]
            r in seenc && continue
            push!(seenc, r)
            row = fill("", length(TAXONOMIC_LEVELS))
            for (k, n) in enumerate(names[1:c])
                row[rank_index(n)] = r[k]
            end
            row[end] = r[c]
            push!(out, row)
        end
    end
    # remove duplicated rows again
    uniq = Vector{Vector{String}}()
    seen2 = Set{Vector{String}}()
    for r in out
        r in seen2 || (push!(seen2, r); push!(uniq, r))
    end
    M = isempty(uniq) ? fill("", 0, length(TAXONOMIC_LEVELS)) : permutedims(reduce(hcat, uniq))
    ref = ReferenceDB(M[:, 1:NRANKS], M[:, end])
    msg = check_tree(ref)
    msg === nothing || error(msg)
    g = _sym(group)
    if g in (:mi, :mf)
        base = reference_db(g)
        have = Set(ref.taxa)
        extra = [i for i in eachindex(base.taxa) if !(base.taxa[i] in have)]
        # drop duplicates within the built-in reference as R `duplicated` does
        seen3 = Set{String}()
        extra = [i for i in extra if !(base.taxa[i] in seen3) && (push!(seen3, base.taxa[i]); true)]
        ref = ReferenceDB(vcat(ref.ranks, base.ranks[extra, :]), vcat(ref.taxa, base.taxa[extra]))
    end
    return ref
end

ref_from_tree(t::Tuple; kw...) = ref_from_tree(t...; kw...)

"""
    check_tree(ref::ReferenceDB) -> Union{Nothing,String}

Consistency checks of a reference database (R `checkTree`): returns an
error message, or `nothing` when no problem is found.
"""
function check_tree(ref::ReferenceDB)
    M = hcat(ref.ranks, ref.taxa)
    cols = [string.(TAXONOMIC_LEVELS)...]
    mes = nothing
    tx = ref.taxa
    if length(unique(tx)) != length(tx)
        seen = Set{String}()
        dups = String[]
        for t in tx
            t in seen ? push!(dups, t) : push!(seen, t)
        end
        mes = "Duplicate names are present in the column Taxa: " * join(dups, ",")
    end
    # replicate the R checks column by column, from Subspecies to Phylum
    rowcount = Int[]
    for c in NRANKS:-1:1
        keep = isempty(rowcount) ? collect(axes(M, 1)) : setdiff(collect(axes(M, 1)), rowcount)
        sub = M[keep, c]
        nz = findall(!isempty, sub)
        names = sub[nz]
        if length(unique(names)) < length(names)
            mes = "Duplicate names are present in $(cols[c]) column"
        end
        append!(rowcount, nz)
    end
    return mes
end
