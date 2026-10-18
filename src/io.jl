# Reading and writing delimited text files, and access to the bundled data.

const DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))

# A parsed delimited file: header plus columns of strings.
struct RawTable
    names::Vector{String}
    columns::Vector{Vector{String}}
end

nrows(t::RawTable) = isempty(t.columns) ? 0 : length(t.columns[1])

function column(t::RawTable, name::AbstractString)
    i = findfirst(==(name), t.names)
    i === nothing && throw(ArgumentError("column $(name) not found"))
    return t.columns[i]
end
hascolumn(t::RawTable, name::AbstractString) = name in t.names

parse_num(s::AbstractString) = isempty(s) || s == "NA" ? NaN : parse(Float64, s)

# split one line, honouring double quotes ("" is an escaped quote)
function split_line(line::AbstractString, delim::Char)
    fields = String[]
    buf = IOBuffer()
    inq = false
    chars = collect(line)
    i = 1
    while i <= length(chars)
        c = chars[i]
        if inq
            if c == '"'
                if i < length(chars) && chars[i+1] == '"'
                    write(buf, '"')
                    i += 1
                else
                    inq = false
                end
            else
                write(buf, c)
            end
        elseif c == '"' && position(buf) == 0
            inq = true
        elseif c == delim
            push!(fields, String(take!(buf)))
        else
            write(buf, c)
        end
        i += 1
    end
    push!(fields, String(take!(buf)))
    return fields
end

function read_raw(io::IO; delim::Char = '\t')
    lines = [chomp(l) for l in eachline(io)]
    filter!(!isempty, lines)
    isempty(lines) && return RawTable(String[], Vector{String}[])
    header = split_line(replace(lines[1], '﻿' => ""), delim)
    cols = [String[] for _ in header]
    for l in lines[2:end]
        f = split_line(l, delim)
        length(f) < length(header) && append!(f, fill("", length(header) - length(f)))
        length(f) > length(header) && throw(ArgumentError("row with $(length(f)) fields, header has $(length(header))"))
        for (j, v) in enumerate(f)
            push!(cols[j], v)
        end
    end
    return RawTable(header, cols)
end
read_raw(path::AbstractString; delim::Char = '\t') = open(io -> read_raw(io; delim), path)

guess_delim(path::AbstractString) = endswith(lowercase(path), ".csv") ? ',' : '\t'

# ---------------------------------------------------------------------------
# Lazy cache of the bundled data (read at run time, so editing a data file
# never leaves stale precompiled constants behind)
# ---------------------------------------------------------------------------

const _CACHE = Dict{Symbol,Any}()
const _CACHE_LOCK = ReentrantLock()

function _cached(f, key::Symbol)
    lock(_CACHE_LOCK) do
        get!(f, _CACHE, key)
    end
end

_raw_dataset(name::AbstractString) = _cached(Symbol("raw_", name)) do
    path = joinpath(DATA_DIR, name * ".tsv")
    isfile(path) || throw(ArgumentError("dataset $(name) not found in $(DATA_DIR)"))
    read_raw(path)
end

# ---------------------------------------------------------------------------
# Conversions from RawTable
# ---------------------------------------------------------------------------

function community_from_raw(t::RawTable; taxa_col::AbstractString = "Taxa")
    j = findfirst(==(taxa_col), t.names)
    j === nothing && throw(ArgumentError("A column called $(taxa_col) is needed"))
    others = [k for k in eachindex(t.names) if k != j]
    A = Matrix{Float64}(undef, nrows(t), length(others))
    for (c, k) in enumerate(others)
        col = t.columns[k]
        for i in eachindex(col)
            v = strip(col[i])
            x = tryparse(Float64, v)
            if x === nothing
                (isempty(v) || v == "NA") || throw(ArgumentError("Non-numeric columns other than Taxa are not allowed (column $(t.names[k]))"))
                x = NaN
            end
            A[i, c] = x
        end
    end
    return CommunityTable(t.columns[j], t.names[others], A)
end

function traits_from_raw(t::RawTable; taxa_col::AbstractString = "Taxa")
    j = findfirst(==(taxa_col), t.names)
    j === nothing && throw(ArgumentError("A column called $(taxa_col) is needed"))
    others = [k for k in eachindex(t.names) if k != j]
    V = Matrix{Float64}(undef, nrows(t), length(others))
    for (c, k) in enumerate(others)
        V[:, c] = parse_num.(strip.(t.columns[k]))
    end
    return TraitTable(strip.(t.columns[j]), t.names[others], V)
end

function reference_from_raw(t::RawTable)
    n = nrows(t)
    R = fill("", n, NRANKS)
    for (k, r) in enumerate(RANKS)
        hascolumn(t, string(r)) || throw(ArgumentError("reference database without column $(r)"))
        col = column(t, string(r))
        for i in 1:n
            v = strip(col[i])
            R[i, k] = v == "NA" ? "" : v
        end
    end
    tx = [v == "NA" ? "" : String(strip(v)) for v in column(t, "Taxa")]
    return ReferenceDB(R, tx)
end

# ---------------------------------------------------------------------------
# Public readers and writers
# ---------------------------------------------------------------------------

"""
    read_community(path; delim = auto, taxa_col = "Taxa") -> CommunityTable

Read a taxa × samples table: one column with taxon names (`Taxa`) and one
numeric column per sample. Empty cells and `NA` become `NaN` (converted to 0
by [`as_biomonitor`](@ref)). The delimiter is `,` for `.csv` files and tab
otherwise.
"""
read_community(path::AbstractString; delim::Char = guess_delim(path), taxa_col::AbstractString = "Taxa") =
    community_from_raw(read_raw(path; delim); taxa_col)

"""
    read_traits(path; delim = auto, taxa_col = "Taxa") -> TraitTable

Read a trait table with a `Taxa` column and one numeric column per trait
modality.
"""
read_traits(path::AbstractString; delim::Char = guess_delim(path), taxa_col::AbstractString = "Taxa") =
    traits_from_raw(read_raw(path; delim); taxa_col)

"""
    read_reference(path; delim = auto) -> ReferenceDB

Read a custom reference database with the columns `Phylum`, `Class`,
`Subclass`, `Order`, `Family`, `Subfamily`, `Tribus`, `Genus`, `Species`,
`Subspecies` and `Taxa` (see [`ref_from_tree`](@ref) to build one from a
simpler taxonomic table).
"""
read_reference(path::AbstractString; delim::Char = guess_delim(path)) =
    reference_from_raw(read_raw(path; delim))

"""
    read_dist(path; delim = auto) -> DistMatrix

Read a square distance matrix whose first column and header hold the labels.
"""
function read_dist(path::AbstractString; delim::Char = guess_delim(path))
    t = read_raw(path; delim)
    labs = t.columns[1]
    D = reduce(hcat, [parse_num.(c) for c in t.columns[2:end]])
    return DistMatrix(labs, D)
end

"""
    write_table(path, x; delim = '\\t')

Write a [`CommunityTable`](@ref), [`TraitTable`](@ref),
[`ReferenceDB`](@ref), [`LabeledMatrix`](@ref), [`LabeledVector`](@ref) or
[`DistMatrix`](@ref) as delimited text.
"""
function write_table(path::AbstractString, x; delim::Char = '\t')
    header, rows = _table_rows(x)
    open(path, "w") do io
        println(io, join(header, delim))
        for r in rows
            println(io, join(r, delim))
        end
    end
    return path
end

_cell(v::AbstractFloat) = isnan(v) ? "NA" : string(v)
_cell(v) = string(v)

_table_rows(x::CommunityTable) = (["Taxa"; x.samples], [[t; _cell.(x.abundance[i, :])] for (i, t) in enumerate(x.taxa)])
_table_rows(x::TraitTable) = (["Taxa"; x.traits], [[t; _cell.(x.values[i, :])] for (i, t) in enumerate(x.taxa)])
_table_rows(x::ReferenceDB) = ([string.(RANKS)...; "Taxa"], [[x.ranks[i, :]; x.taxa[i]] for i in eachindex(x.taxa)])
_table_rows(x::LabeledMatrix) = (["";  x.colnames], [[r; _cell.(x.values[i, :])] for (i, r) in enumerate(x.rownames)])
_table_rows(x::LabeledVector) = (["name", "value"], [[l, _cell(v)] for (l, v) in zip(x.labels, x.values)])
_table_rows(x::DistMatrix) = (["";  x.labels], [[l; _cell.(x.d[i, :])] for (i, l) in enumerate(x.labels)])
_table_rows(x::BiomonitorData) = ([string.(RANKS)...; "Taxa"; x.samples],
                                  [[x.ranks[i, :]; x.taxa[i]; _cell.(x.abundance[i, :])] for i in eachindex(x.taxa)])

# ---------------------------------------------------------------------------
# Bundled datasets
# ---------------------------------------------------------------------------

"""
    macro_ex() -> CommunityTable

Small macroinvertebrate dataset with 34 taxa and 2 samples (R `macro_ex`).
"""
macro_ex() = _cached(:macro_ex) do
    community_from_raw(_raw_dataset("macro_ex"))
end

"""
    mi_prin() -> CommunityTable

Macroinvertebrates from 8 Italian rivers, 5 replicates each (R `mi_prin`).
"""
mi_prin() = _cached(:mi_prin) do
    community_from_raw(_raw_dataset("mi_prin"))
end

"""
    oglio() -> CommunityTable

Presence-absence of macrophytes in 30 sites of the Oglio river (R `oglio`).
"""
oglio() = _cached(:oglio) do
    community_from_raw(_raw_dataset("oglio"))
end

"""
    traits_mf() -> TraitTable

Traits of the macrophytes of the `oglio` dataset (R `traits_mf`).
"""
traits_mf() = _cached(:traits_mf) do
    traits_from_raw(_raw_dataset("traits_mf"))
end

"""
    traits_tachet() -> TraitTable

Fuzzy-coded macroinvertebrate traits of Tachet et al. (2010), 484 taxa and
113 modalities in 21 traits (R internal `traitsTachet`). The trait blocks
are `[8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8]`.
"""
traits_tachet() = _cached(:traitsTachet) do
    traits_from_raw(_raw_dataset("traitsTachet"))
end

"""
    diatom_traits() -> TraitTable

Diatom traits and ecological values (R internal `diatom_traits`).
"""
diatom_traits() = _cached(:diatom_traits) do
    traits_from_raw(_raw_dataset("diatom_traits"))
end

"""
    example_tree() -> (levels, values)

Small taxonomic table (R dataset `Tree`) for [`ref_from_tree`](@ref):
`ref_from_tree(example_tree()...)`.
"""
function example_tree()
    t = _raw_dataset("Tree")
    return (t.names, reduce(hcat, t.columns))
end

"""
    reference_db(group = :mi) -> ReferenceDB

Built-in reference database: `:mi` macroinvertebrates, `:mf` macrophytes,
`:fi` fish, `:di` diatoms. They rely on www.freshwaterecology.info
(Schmidt-Kloiber & Hering, 2015).
"""
function reference_db(group = :mi)
    g = _sym(group)
    g in (:mi, :mf, :fi, :di) || throw(ArgumentError("group must be one of mi, mf, fi or di"))
    return _cached(Symbol(g, "_ref")) do
        reference_from_raw(_raw_dataset(string(g, "_ref")))
    end
end

"""
    reference_datasets() -> Vector{NamedTuple}

Datasets available online through [`get_reference_db`](@ref): name,
download link, citation and version.
"""
function reference_datasets()
    t = _raw_dataset("dfMain")
    return [(dataset = column(t, "dataset")[i], link = column(t, "link")[i],
             citation = column(t, "citation_APA")[i], version = column(t, "version")[i]) for i in 1:nrows(t)]
end

"""
    get_reference_db(name = "mi_cuba_2023") -> ReferenceDB

Download one of the reference databases listed by
[`reference_datasets`](@ref) (OSF repository of biomonitoR). Requires an
internet connection. Datasets that are not taxonomic trees (score tables)
are returned as a [`CommunityTable`]-like raw table of strings, i.e. a
`Vector{NamedTuple}`.
"""
function get_reference_db(name::AbstractString = "mi_cuba_2023")
    ds = reference_datasets()
    i = findfirst(d -> d.dataset == name, ds)
    i === nothing && throw(ArgumentError("Please provide a valid db_name"))
    buf = IOBuffer()
    Downloads.download(ds[i].link, buf)
    seekstart(buf)
    t = read_raw(buf; delim = ',')
    if all(r -> hascolumn(t, string(r)), RANKS) && hascolumn(t, "Taxa")
        return reference_from_raw(t)
    end
    return [NamedTuple{Tuple(Symbol.(t.names))}(Tuple(c[k] for c in t.columns)) for k in 1:nrows(t)]
end

# score tables -----------------------------------------------------------------

_score_dict(name, col = "Scores") = _cached(Symbol("scores_", name, col)) do
    t = _raw_dataset(name)
    Dict(String(strip(a)) => parse_num(b) for (a, b) in zip(column(t, "Taxon"), column(t, col)))
end

_rule_dict(name) = _cached(Symbol("rules_", name)) do
    t = _raw_dataset(name)
    Dict(String(strip(a)) => String(strip(b)) for (a, b) in zip(column(t, "Taxon"), column(t, "Correct_Taxon")))
end
