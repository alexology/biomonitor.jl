# Core data types.

"""
Taxonomic ranks used by the reference databases, from the coarsest to the
finest. Every reference row has a value (possibly empty) for each rank.
"""
const RANKS = (:Phylum, :Class, :Subclass, :Order, :Family, :Subfamily,
               :Tribus, :Genus, :Species, :Subspecies)

"""
All taxonomic levels at which [`aggregate_taxa`](@ref) aggregates data:
the ten [`RANKS`](@ref) plus `:Taxa` (the name used in the data).
"""
const TAXONOMIC_LEVELS = (RANKS..., :Taxa)

const NRANKS = length(RANKS)

function rank_index(level)
    lev = Symbol(level)
    i = findfirst(==(lev), TAXONOMIC_LEVELS)
    i === nothing && throw(ArgumentError("unknown taxonomic level $(level); valid levels are $(join(TAXONOMIC_LEVELS, ", "))"))
    return i
end

# ---------------------------------------------------------------------------
# Labeled vectors and matrices (the Julia counterpart of R named vectors and
# data.frames of results)
# ---------------------------------------------------------------------------

"""
    LabeledVector{T}(labels, values)

A vector whose elements also have names, used for per-sample results
(e.g. `shannon(x)` returns one value per sample). Index it by position or by
label: `v[1]`, `v["Sample_1"]`. `labels(v)` returns the names.
"""
struct LabeledVector{T} <: AbstractVector{T}
    labels::Vector{String}
    values::Vector{T}
    function LabeledVector{T}(labels::AbstractVector, values::AbstractVector) where {T}
        length(labels) == length(values) || throw(DimensionMismatch("labels and values must have the same length"))
        new{T}(String.(labels), convert(Vector{T}, collect(values)))
    end
end
LabeledVector(labels::AbstractVector, values::AbstractVector{T}) where {T} = LabeledVector{T}(labels, values)

labels(v::LabeledVector) = v.labels
Base.size(v::LabeledVector) = size(v.values)
Base.getindex(v::LabeledVector, i::Int) = v.values[i]
Base.getindex(v::LabeledVector, name::AbstractString) = v.values[_findlabel(v.labels, name)]
Base.setindex!(v::LabeledVector, x, i::Int) = (v.values[i] = x)
Base.setindex!(v::LabeledVector, x, name::AbstractString) = (v.values[_findlabel(v.labels, name)] = x)
Base.IndexStyle(::Type{<:LabeledVector}) = IndexLinear()
Base.haskey(v::LabeledVector, name::AbstractString) = name in v.labels
Base.pairs(v::LabeledVector) = (l => x for (l, x) in zip(v.labels, v.values))
Base.similar(v::LabeledVector, ::Type{S}, dims::Dims{1}) where {S} =
    dims == size(v) ? LabeledVector{S}(copy(v.labels), Vector{S}(undef, dims)) : Vector{S}(undef, dims)

function _findlabel(labs, name)
    i = findfirst(==(name), labs)
    i === nothing && throw(KeyError(name))
    return i
end

function Base.show(io::IO, ::MIME"text/plain", v::LabeledVector{T}) where {T}
    println(io, length(v), "-element LabeledVector{", T, "}:")
    w = isempty(v.labels) ? 0 : maximum(textwidth, v.labels)
    for (k, (l, x)) in enumerate(zip(v.labels, v.values))
        print(io, " ", rpad(l, w), " => ", x)
        k < length(v) && println(io)
    end
end

"""
    LabeledMatrix{T}(rownames, colnames, values)

A matrix with row and column names, used for tabular results (e.g.
`allindices(x)` returns samples × indices). Index it by position or by
names: `m["Sample_1", "shannon"]`; `m[:, "shannon"]` returns a
[`LabeledVector`](@ref) indexed by row names.
"""
struct LabeledMatrix{T} <: AbstractMatrix{T}
    rownames::Vector{String}
    colnames::Vector{String}
    values::Matrix{T}
    function LabeledMatrix{T}(rn::AbstractVector, cn::AbstractVector, values::AbstractMatrix) where {T}
        size(values) == (length(rn), length(cn)) || throw(DimensionMismatch("matrix size $(size(values)) does not match the names ($(length(rn)), $(length(cn)))"))
        new{T}(String.(rn), String.(cn), convert(Matrix{T}, values))
    end
end
LabeledMatrix(rn::AbstractVector, cn::AbstractVector, values::AbstractMatrix{T}) where {T} = LabeledMatrix{T}(rn, cn, values)

rownames(m::LabeledMatrix) = m.rownames
colnames(m::LabeledMatrix) = m.colnames
Base.size(m::LabeledMatrix) = size(m.values)
Base.getindex(m::LabeledMatrix, i::Int, j::Int) = m.values[i, j]
Base.getindex(m::LabeledMatrix, r::AbstractString, c::AbstractString) =
    m.values[_findlabel(m.rownames, r), _findlabel(m.colnames, c)]
Base.getindex(m::LabeledMatrix, ::Colon, c::AbstractString) =
    LabeledVector(m.rownames, m.values[:, _findlabel(m.colnames, c)])
Base.getindex(m::LabeledMatrix, r::AbstractString, ::Colon) =
    LabeledVector(m.colnames, m.values[_findlabel(m.rownames, r), :])
Base.setindex!(m::LabeledMatrix, x, i::Int, j::Int) = (m.values[i, j] = x)
Base.IndexStyle(::Type{<:LabeledMatrix}) = IndexCartesian()

function Base.show(io::IO, ::MIME"text/plain", m::LabeledMatrix{T}) where {T}
    println(io, size(m, 1), "×", size(m, 2), " LabeledMatrix{", T, "}:")
    _print_table(io, m.rownames, m.colnames, m.values)
end

function _fmt(x)
    x isa AbstractFloat && return isnan(x) ? "NaN" : string(round(x; sigdigits = 6))
    return string(x)
end

function _print_table(io, rn, cn, vals; maxrows = 30, maxcols = 10)
    nr, nc = size(vals)
    showr = min(nr, maxrows)
    showc = min(nc, maxcols)
    cells = [_fmt(vals[i, j]) for i in 1:showr, j in 1:showc]
    rw = isempty(rn) ? 0 : maximum(textwidth, rn[1:showr]; init = 0)
    cw = [max(textwidth(cn[j]), maximum(textwidth, cells[:, j]; init = 0)) for j in 1:showc]
    print(io, " ", " "^rw)
    for j in 1:showc
        print(io, "  ", lpad(cn[j], cw[j]))
    end
    showc < nc && print(io, "  … (", nc - showc, " more columns)")
    for i in 1:showr
        println(io)
        print(io, " ", rpad(rn[i], rw))
        for j in 1:showc
            print(io, "  ", lpad(cells[i, j], cw[j]))
        end
    end
    showr < nr && print(io, "\n ⋮ (", nr - showr, " more rows)")
end

# ---------------------------------------------------------------------------
# Community data
# ---------------------------------------------------------------------------

"""
    CommunityTable(taxa, samples, abundance)

Taxa × samples table of abundances (or presence-absence values). This is the
input of [`as_biomonitor`](@ref) (the R `data.frame` with a `Taxa` column)
and also the per-level tables stored in an [`AggregatedData`](@ref).
Missing values are stored as `NaN`.

Construct it from vectors and a matrix, or from a `NamedTuple`/`Dict` of
columns where the `Taxa` column holds the names and every other column is a
sample, e.g. `CommunityTable((Taxa = ["Baetis"], S1 = [3], S2 = [0]))`.
"""
struct CommunityTable
    taxa::Vector{String}
    samples::Vector{String}
    abundance::Matrix{Float64}
    function CommunityTable(taxa::AbstractVector, samples::AbstractVector, abundance::AbstractMatrix)
        size(abundance) == (length(taxa), length(samples)) ||
            throw(DimensionMismatch("abundance is $(size(abundance)) but there are $(length(taxa)) taxa and $(length(samples)) samples"))
        A = Matrix{Float64}(undef, size(abundance))
        for i in eachindex(abundance)
            v = abundance[i]
            A[i] = v === missing || v === nothing ? NaN : Float64(v)
        end
        new(String.(string.(taxa)), String.(string.(samples)), A)
    end
end

function CommunityTable(cols::Union{NamedTuple,AbstractDict})
    ks = collect(keys(cols))
    tk = findfirst(k -> string(k) == "Taxa", ks)
    tk === nothing && throw(ArgumentError("A column called Taxa is needed"))
    taxacol = [string(t) for t in cols[ks[tk]]]
    others = [k for k in ks if string(k) != "Taxa"]
    for k in others
        all(v -> v === missing || v isa Real, cols[k]) ||
            throw(ArgumentError("Non-numeric columns other than Taxa are not allowed"))
    end
    A = Matrix{Union{Missing,Float64}}(undef, length(taxacol), length(others))
    for (j, k) in enumerate(others)
        A[:, j] = [v === missing ? missing : Float64(v) for v in cols[k]]
    end
    return CommunityTable(taxacol, string.(others), A)
end

taxa(x::CommunityTable) = x.taxa
samples(x::CommunityTable) = x.samples
abundances(x::CommunityTable) = x.abundance
Base.size(x::CommunityTable) = (length(x.taxa), length(x.samples))
Base.copy(x::CommunityTable) = CommunityTable(copy(x.taxa), copy(x.samples), copy(x.abundance))
Base.:(==)(a::CommunityTable, b::CommunityTable) =
    a.taxa == b.taxa && a.samples == b.samples && isequal(a.abundance, b.abundance)

"""
    getindex(x::CommunityTable, taxon::AbstractString) -> LabeledVector

Abundances of one taxon across samples.
"""
Base.getindex(x::CommunityTable, taxon::AbstractString) =
    LabeledVector(x.samples, x.abundance[_findlabel(x.taxa, taxon), :])

function Base.show(io::IO, ::MIME"text/plain", x::CommunityTable)
    println(io, "CommunityTable with ", length(x.taxa), " taxa and ", length(x.samples), " samples:")
    _print_table(io, x.taxa, x.samples, x.abundance)
end
Base.show(io::IO, x::CommunityTable) = print(io, "CommunityTable(", length(x.taxa), " taxa × ", length(x.samples), " samples)")

# ---------------------------------------------------------------------------
# Reference database
# ---------------------------------------------------------------------------

"""
    ReferenceDB(ranks, taxa)

Taxonomic reference database: `ranks` is an n × 10 matrix of names for the
[`RANKS`](@ref) (empty string when a rank is not defined) and `taxa` the n
taxon names that can appear in the data. Built-in databases are returned by
[`reference_db`](@ref); custom ones can be read with [`read_reference`](@ref)
or built with [`ref_from_tree`](@ref).
"""
struct ReferenceDB
    ranks::Matrix{String}
    taxa::Vector{String}
    function ReferenceDB(ranks::AbstractMatrix, taxa::AbstractVector)
        size(ranks, 2) == NRANKS || throw(DimensionMismatch("ranks must have $(NRANKS) columns ($(join(RANKS, ", ")))"))
        size(ranks, 1) == length(taxa) || throw(DimensionMismatch("ranks and taxa must have the same number of rows"))
        R = Matrix{String}(undef, size(ranks))
        for i in eachindex(ranks)
            v = ranks[i]
            R[i] = v === missing || v === nothing ? "" : String(strip(string(v)))
        end
        new(R, [String(strip(string(t))) for t in taxa])
    end
end

ranks(r::ReferenceDB) = r.ranks
taxa(r::ReferenceDB) = r.taxa
Base.length(r::ReferenceDB) = length(r.taxa)
Base.:(==)(a::ReferenceDB, b::ReferenceDB) = a.ranks == b.ranks && a.taxa == b.taxa
Base.show(io::IO, r::ReferenceDB) = print(io, "ReferenceDB(", length(r.taxa), " taxa)")
function Base.show(io::IO, ::MIME"text/plain", r::ReferenceDB)
    println(io, "ReferenceDB with ", length(r.taxa), " taxa:")
    _print_table(io, r.taxa, [string.(RANKS)...], r.ranks; maxrows = 10)
end

# ---------------------------------------------------------------------------
# Imported data (R class "asb") and aggregated data (R class "biomonitoR")
# ---------------------------------------------------------------------------

"""
    BiomonitorData

Result of [`as_biomonitor`](@ref) (R class `asb`): the user data merged with
the reference database. Fields:

- `ranks`: n × 10 matrix with the taxonomic tree of each taxon;
- `taxa`: taxon names; `samples`: sample names;
- `abundance`: n × samples matrix;
- `group`: `:mi`, `:mf`, `:fi`, `:di` or `:custom`;
- `presence_absence`: `true` when imported with `FUN = bin`;
- `corrected_names`: names changed through `to_change` (`wrong => correct`);
- `suggested_names`: names not found in the reference with spelling
  suggestions (`excluded => suggestions`).
"""
struct BiomonitorData
    ranks::Matrix{String}
    taxa::Vector{String}
    samples::Vector{String}
    abundance::Matrix{Float64}
    group::Symbol
    presence_absence::Bool
    corrected_names::Vector{Pair{String,String}}
    suggested_names::Vector{Pair{String,Vector{String}}}
end

taxa(x::BiomonitorData) = x.taxa
samples(x::BiomonitorData) = x.samples
ranks(x::BiomonitorData) = x.ranks
abundances(x::BiomonitorData) = x.abundance
is_custom(x::BiomonitorData) = x.group === :custom

function Base.show(io::IO, ::MIME"text/plain", x::BiomonitorData)
    kind = x.presence_absence ? "presence-absence" : "abundance"
    println(io, "BiomonitorData (", kind, ", group ", x.group, ") with ", length(x.taxa),
            " taxa and ", length(x.samples), " samples:")
    _print_table(io, x.taxa, x.samples, x.abundance)
end
Base.show(io::IO, x::BiomonitorData) = print(io, "BiomonitorData(", length(x.taxa), " taxa × ", length(x.samples), " samples)")

"""
    AggregatedData

Result of [`aggregate_taxa`](@ref) (R class `biomonitoR`): the imported data
plus one [`CommunityTable`](@ref) per taxonomic level. Index it with a level:
`x[:Family]`, `x["Genus"]`. Taxa without a name at a level are pooled in the
row `"unassigned"`. `x.data` is the underlying [`BiomonitorData`](@ref) (the
R `Tree` element).
"""
struct AggregatedData
    data::BiomonitorData
    levels::Dict{Symbol,CommunityTable}
end

Base.getindex(x::AggregatedData, level::Symbol) = x.levels[TAXONOMIC_LEVELS[rank_index(level)]]
Base.getindex(x::AggregatedData, level::AbstractString) = x[Symbol(level)]
samples(x::AggregatedData) = x.data.samples
taxa(x::AggregatedData) = x.data.taxa
ranks(x::AggregatedData) = x.data.ranks
is_bin(x::AggregatedData) = x.data.presence_absence
is_bin(x::BiomonitorData) = x.presence_absence
is_custom(x::AggregatedData) = is_custom(x.data)

function Base.show(io::IO, ::MIME"text/plain", x::AggregatedData)
    kind = x.data.presence_absence ? "presence-absence" : "abundance"
    println(io, "AggregatedData (", kind, ", group ", x.data.group, "), ", length(x.data.samples), " samples")
    for lev in TAXONOMIC_LEVELS
        t = x.levels[lev]
        n = count(!=("unassigned"), t.taxa)
        print(io, "  ", rpad(string(lev), 11), n, " taxa")
        lev === :Taxa || println(io)
    end
end
Base.show(io::IO, x::AggregatedData) = print(io, "AggregatedData(", length(x.data.taxa), " taxa × ", length(x.data.samples), " samples)")

# ---------------------------------------------------------------------------
# Traits
# ---------------------------------------------------------------------------

"""
    TraitTable(taxa, traits, values)

Taxa × traits table (fuzzy-coded affinities or continuous values). `NaN`
marks missing values. Built-in tables: [`traits_tachet`](@ref),
[`traits_mf`](@ref), [`diatom_traits`](@ref). Read your own with
[`read_traits`](@ref).
"""
struct TraitTable
    taxa::Vector{String}
    traits::Vector{String}
    values::Matrix{Float64}
    function TraitTable(taxa::AbstractVector, traits::AbstractVector, values::AbstractMatrix)
        size(values) == (length(taxa), length(traits)) ||
            throw(DimensionMismatch("values must be $(length(taxa)) × $(length(traits))"))
        V = Matrix{Float64}(undef, size(values))
        for i in eachindex(values)
            v = values[i]
            V[i] = v === missing || v === nothing ? NaN : Float64(v)
        end
        new(String.(string.(taxa)), String.(string.(traits)), V)
    end
end

function TraitTable(cols::Union{NamedTuple,AbstractDict})
    ks = collect(keys(cols))
    tk = findfirst(k -> string(k) == "Taxa", ks)
    tk === nothing && throw(ArgumentError("A column called Taxa is needed"))
    others = [k for k in ks if string(k) != "Taxa"]
    V = Matrix{Union{Missing,Float64}}(undef, length(cols[ks[tk]]), length(others))
    for (j, k) in enumerate(others)
        V[:, j] = [v === missing ? missing : Float64(v) for v in cols[k]]
    end
    return TraitTable([string(t) for t in cols[ks[tk]]], string.(others), V)
end

taxa(t::TraitTable) = t.taxa
Base.size(t::TraitTable) = (length(t.taxa), length(t.traits))
Base.:(==)(a::TraitTable, b::TraitTable) = a.taxa == b.taxa && a.traits == b.traits && isequal(a.values, b.values)
Base.getindex(t::TraitTable, taxon::AbstractString) = LabeledVector(t.traits, t.values[_findlabel(t.taxa, taxon), :])
function Base.show(io::IO, ::MIME"text/plain", t::TraitTable)
    println(io, "TraitTable with ", length(t.taxa), " taxa and ", length(t.traits), " traits:")
    _print_table(io, t.taxa, t.traits, t.values)
end
Base.show(io::IO, t::TraitTable) = print(io, "TraitTable(", length(t.taxa), " taxa × ", length(t.traits), " traits)")

"""
    AssignedTraits

Result of [`assign_traits`](@ref): one row for each pair (taxon in the data,
taxon of the trait database whose traits were assigned to it). Fields:

- `taxa`: taxon of the user data;
- `taxa_level`: its taxonomic level;
- `trait_taxa`: taxon of the trait database providing the traits;
- `trait_level`: its taxonomic level;
- `distance`: taxonomic distance (positive when the trait taxon is finer,
  e.g. traits of a species assigned to a genus);
- `traits`, `values`: trait names and values.

In the R package these columns are called `Taxa`, `Taxa_db`,
`Traits_taxlev`, `Traits_real` and `Taxonomic_distance`.
"""
struct AssignedTraits
    taxa::Vector{String}
    taxa_level::Vector{Symbol}
    trait_taxa::Vector{String}
    trait_level::Vector{Symbol}
    distance::Vector{Int}
    traits::Vector{String}
    values::Matrix{Float64}
end

taxa(a::AssignedTraits) = a.taxa
Base.length(a::AssignedTraits) = length(a.taxa)

function _subset(a::AssignedTraits, idx)
    AssignedTraits(a.taxa[idx], a.taxa_level[idx], a.trait_taxa[idx], a.trait_level[idx],
                   a.distance[idx], copy(a.traits), a.values[idx, :])
end

function Base.show(io::IO, ::MIME"text/plain", a::AssignedTraits)
    println(io, "AssignedTraits with ", length(a.taxa), " rows (", length(unique(a.taxa)),
            " taxa) and ", length(a.traits), " traits:")
    n = min(length(a.taxa), 20)
    w1 = maximum(textwidth, a.taxa[1:n]; init = 4)
    w2 = maximum(textwidth, a.trait_taxa[1:n]; init = 10)
    print(io, " ", rpad("taxon", w1), "  ", rpad("level", 10), "  ", rpad("trait taxon", w2), "  ", rpad("level", 10), "  distance")
    for i in 1:n
        print(io, "\n ", rpad(a.taxa[i], w1), "  ", rpad(string(a.taxa_level[i]), 10), "  ",
              rpad(a.trait_taxa[i], w2), "  ", rpad(string(a.trait_level[i]), 10), "  ", a.distance[i])
    end
    n < length(a.taxa) && print(io, "\n ⋮ (", length(a.taxa) - n, " more rows)")
end
Base.show(io::IO, a::AssignedTraits) = print(io, "AssignedTraits(", length(a.taxa), " rows)")

# ---------------------------------------------------------------------------
# Distance matrices
# ---------------------------------------------------------------------------

"""
    DistMatrix(labels, d)

Symmetric dissimilarity matrix with labels (the R `dist` class). `d` is the
full n × n matrix with zeros on the diagonal.
"""
struct DistMatrix
    labels::Vector{String}
    d::Matrix{Float64}
    function DistMatrix(labels::AbstractVector, d::AbstractMatrix)
        n = length(labels)
        size(d) == (n, n) || throw(DimensionMismatch("distance matrix must be $(n) × $(n)"))
        D = Matrix{Float64}(d)
        for i in 1:n
            D[i, i] = 0.0
        end
        new(String.(string.(labels)), D)
    end
end

labels(d::DistMatrix) = d.labels
Base.size(d::DistMatrix) = size(d.d)
Base.length(d::DistMatrix) = length(d.labels)
Base.Matrix(d::DistMatrix) = copy(d.d)
Base.getindex(d::DistMatrix, i::Int, j::Int) = d.d[i, j]
Base.getindex(d::DistMatrix, a::AbstractString, b::AbstractString) = d.d[_findlabel(d.labels, a), _findlabel(d.labels, b)]

"lower-triangle values of a distance matrix, in R `dist` order"
function lowertri(D::AbstractMatrix)
    n = size(D, 1)
    out = Vector{Float64}(undef, n * (n - 1) ÷ 2)
    k = 0
    for j in 1:n-1, i in j+1:n
        out[k += 1] = D[i, j]
    end
    return out
end
lowertri(d::DistMatrix) = lowertri(d.d)

function Base.show(io::IO, ::MIME"text/plain", d::DistMatrix)
    println(io, length(d.labels), "×", length(d.labels), " DistMatrix:")
    _print_table(io, d.labels, d.labels, d.d)
end
Base.show(io::IO, d::DistMatrix) = print(io, "DistMatrix(", length(d.labels), " objects)")

# ---------------------------------------------------------------------------
# Accessor docstrings
# ---------------------------------------------------------------------------

"""
    labels(v::LabeledVector)
    labels(d::DistMatrix)

Names of the elements of a [`LabeledVector`](@ref) or of the objects of a
[`DistMatrix`](@ref).
"""
labels

"""
    rownames(m::LabeledMatrix)

Row names of a [`LabeledMatrix`](@ref).
"""
rownames

"""
    colnames(m::LabeledMatrix)

Column names of a [`LabeledMatrix`](@ref).
"""
colnames

"""
    taxa(x)

Taxon names of a [`CommunityTable`](@ref), [`ReferenceDB`](@ref),
[`BiomonitorData`](@ref), [`AggregatedData`](@ref), [`TraitTable`](@ref)
or [`AssignedTraits`](@ref).
"""
taxa

"""
    samples(x)

Sample names of a [`CommunityTable`](@ref), [`BiomonitorData`](@ref) or
[`AggregatedData`](@ref).
"""
samples

"""
    abundances(x)

Taxa × samples abundance matrix of a [`CommunityTable`](@ref) or
[`BiomonitorData`](@ref).
"""
abundances

"""
    ranks(x)

n × 10 matrix with the names of the [`RANKS`](@ref) of each taxon of a
[`ReferenceDB`](@ref), [`BiomonitorData`](@ref) or
[`AggregatedData`](@ref) (empty string when undefined).
"""
ranks
