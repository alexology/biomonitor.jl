# Functional diversity indices (FRic, FDis, FEve, Rao's Q, functional
# redundancy) and the quality of functional spaces.

const _DEFAULT_SET_PARAM = (max_nbdim = 15, prec = "Qt", tol = 1e-7, cor_zero = true)

_set_param(sp) = sp === nothing ? _DEFAULT_SET_PARAM : merge(_DEFAULT_SET_PARAM, sp)

# Common preparation of f_rich, f_disp, f_eve, f_divs and f_red: aligns the
# community and trait data, computes the trait dissimilarity, removes zero
# distances and applies the correction.
function _functional_setup(x::AggregatedData, trait_db, tax_lev, type, col_blocks, distance,
                           zerodist_rm, correction, sp)
    _check_aggregated(x)
    BIN = is_bin(x)
    trait_db === nothing && throw(ArgumentError("Please provide trait_db"))
    trait_db isa Union{TraitTable,DistMatrix} || throw(ArgumentError("trait_db must be a TraitTable or a DistMatrix"))
    distance = string(distance)
    tp = type === nothing ? nothing : string(type)
    if trait_db isa TraitTable
        tp === nothing && throw(ArgumentError("Please specify a type when trait_db is a TraitTable"))
        tp in ("F", "C") || throw(ArgumentError("type must be C or F when trait_db is a TraitTable"))
        tp == "C" && distance == "gower" && throw(ArgumentError("Using gower distance when type is C is currently not allowed"))
        tp == "F" && distance == "euclidean" && @warn "Are you sure to use euclidean distance when type is F?"
        if tp == "F"
            col_blocks === nothing && throw(ArgumentError("Please provide col_blocks"))
            sum(col_blocks) == length(trait_db.traits) ||
                throw(ArgumentError("The number of traits in trait_db is not equal to the sum of col_blocks"))
        end
    end
    t = x[tax_lev]
    labels_ = trait_db isa TraitTable ? [capwords(s) for s in trait_db.taxa] : trait_db.labels
    common = sort(intersect(t.taxa, labels_))
    ti = [findfirst(==(c), t.taxa) for c in common]
    li = [findfirst(==(c), labels_) for c in common]
    DF = CommunityTable(common, t.samples, t.abundance[ti, :])
    BIN && (DF = to_bin(DF))
    tr_prep = nothing
    traits_used = nothing
    if trait_db isa TraitTable
        V = trait_db.values[li, :]
        traits_used = TraitTable(common, trait_db.traits, V)
        blocks = [size(V, 2)]
        if tp == "F"
            P, blocks, keep = _prep_fuzzy(V, col_blocks)
            tr_prep = TraitTable(common, trait_db.traits[keep], P)
        else
            tr_prep = traits_used
        end
        D = distance == "gower" ? gower_fuzzy(tr_prep.values, blocks) :
            distance == "euclidean" ? euclidean_dist(tr_prep.values) :
            throw(ArgumentError("distance must be gower or euclidean"))
    else
        D = trait_db.d[li, li]
    end
    off = [D[i, j] for j in axes(D, 2) for i in axes(D, 1) if i > j]
    dup_pairs = nothing
    if zerodist_rm && any(<(sp.tol), off)
        DF, D, dup_pairs = _zero_dist_traits(DF, D, BIN)
    end
    dm = DistMatrix(DF.taxa, D)
    distance == "gower" && (dm = _correct(dm, correction; tol = sp.tol, cor_zero = sp.cor_zero))
    euclid = is_euclid(dm; tol = sp.tol)
    off = lowertri(dm)
    if any(<(sp.tol), off)
        mes = "At least a pair of species has the same traits. Depending on your needs, this could be an issue."
        @info mes
    else
        mes = "no taxa with the same traits"
    end
    na = if tr_prep === nothing
        "NAs cannot be detected when trait_db is a DistMatrix"
    else
        ids = [(taxon = tr_prep.taxa[i], trait = tr_prep.traits[j]) for i in axes(tr_prep.values, 1)
               for j in axes(tr_prep.values, 2) if isnan(tr_prep.values[i, j])]
        isempty(ids) ? "No NAs detected" : sort(ids; by = r -> r.taxon)
    end
    return (DF = DF, dist = dm, euclid = euclid, message = mes, dup_pairs = dup_pairs, na = na,
            traits = trait_db isa TraitTable ? traits_used : dm)
end

# R helper `zero_dist_traits`: taxa with identical traits are merged. With
# missing traits a zero distance is not transitive (A-B and B-C zero, A-C
# positive), so a taxon joins the first group whose members are all at zero
# distance from it, named after its first taxon (R merges along the pairs and
# stops with an error in that case).
function _zero_dist_traits(DF::CommunityTable, D::Matrix{Float64}, BIN::Bool)
    n = length(DF.taxa)
    pairs = Tuple{String,String}[]
    for i in 1:n, j in i+1:n
        D[i, j] == 0 && push!(pairs, (DF.taxa[i], DF.taxa[j]))
    end
    groups = Vector{Int}[]
    for i in 1:n
        g = findfirst(m -> all(j -> D[j, i] == 0, m), groups)
        g === nothing ? push!(groups, [i]) : push!(groups[g], i)
    end
    names = copy(DF.taxa)
    removed = Set{String}()
    for m in groups, j in m[2:end]
        names[j] = DF.taxa[m[1]]
        push!(removed, DF.taxa[j])
    end
    gn, G = group_rows(names, DF.abundance)
    keep = [i for i in 1:n if !(DF.taxa[i] in removed)]
    DF2 = CommunityTable(gn, DF.samples, G)
    BIN && (DF2 = to_bin(DF2))
    DF.taxa[keep] == gn || error("Something went wrong with the zero distance removal.")
    return DF2, D[keep, keep], [(taxon = a, same_traits = b) for (a, b) in pairs]
end

function _nbdim(dm, nbdim, sp)
    nbdim === :auto || nbdim == "auto" || return Int(nbdim)
    q = qfs(dm; nbdim = sp.max_nbdim)
    k = findfirst(<(0.01), q.meanSD.values)
    k === nothing && throw(ArgumentError("there is no optimal number of dimension, please increase the number of dimensions"))
    return k + 1
end

function _functional_trace(res, s, m, correction)
    base = (results = res, traits = s.traits, taxa = s.DF)
    extra = m === nothing ? (;) : (nbdim = m,)
    return merge(base, extra, (correction = string(correction), NA_detection = s.na,
                              duplicated_traits = s.dup_pairs === nothing ? s.message : s.dup_pairs))
end

_noneuclid_msg() = @info "Negative eigenvalues found, please consider to add a correction to the pcoa"

"""
    f_rich(x::AggregatedData; trait_db, tax_lev = :Taxa, type = nothing,
           col_blocks = nothing, nbdim = 2, distance = "gower", zerodist_rm = false,
           correction = "none", traceB = false, set_param = nothing) -> LabeledVector

Functional richness (Villéger et al. 2008): volume of the convex hull of the
taxa of each sample in a PCoA functional space with `nbdim` axes, relative
to the volume of all taxa. Samples with fewer than `nbdim + 1` taxa get 0.

- `trait_db`: a [`TraitTable`](@ref) (`type = "F"` fuzzy with
  `col_blocks`, or `"C"` continuous) or a [`DistMatrix`](@ref).
- `distance`: `"gower"` (fuzzy traits) or `"euclidean"` (scaled traits).
- `nbdim`: number of axes, or `:auto` to choose it with [`qfs`](@ref)
  (first space with mean squared deviation < 0.01).
- `zerodist_rm`: merge taxa with identical traits.
- `correction`: `"none"`, `"cailliez"`, `"lingoes"`, `"sqrt"`, `"quasi"`
  to make the distance Euclidean.
- `set_param`: `NamedTuple` with `max_nbdim` (15), `tol` (1e-7) and
  `cor_zero` (true).
"""
function f_rich(x::AggregatedData; trait_db = nothing, tax_lev = :Taxa, type = nothing, col_blocks = nothing,
                nbdim = 2, distance = "gower", zerodist_rm::Bool = false, correction = "none",
                traceB::Bool = false, set_param = nothing)
    sp = _set_param(set_param)
    s = _functional_setup(x, trait_db, tax_lev, type, col_blocks, distance, zerodist_rm, correction, sp)
    s.euclid || _noneuclid_msg()
    m = _nbdim(s.dist, nbdim, sp)
    pc = pcoa(s.dist; nf = m)
    size(pc.coords, 2) < m && throw(ArgumentError("the functional space has only $(size(pc.coords, 2)) axes, fewer than nbdim = $(m)"))
    C = pc.coords[:, 1:m]
    P = s.DF.abundance .> 0
    vmax = convex_hull_volume(C)
    res = map(axes(P, 2)) do j
        rows = findall(view(P, :, j))
        length(rows) >= m + 1 || return 0.0
        v = convex_hull_volume(C[rows, :]) / vmax
        isnan(v) ? 0.0 : v
    end
    out = LabeledVector(samples(x), res)
    traceB || return out
    return _functional_trace(out, s, m, correction)
end

"""
    f_disp(x::AggregatedData; trait_db, ...) -> LabeledVector

Functional dispersion (Laliberté & Legendre 2010): abundance-weighted mean
distance of the taxa to the community centroid in the PCoA space. Same
arguments as [`f_rich`](@ref); with `nbdim` axes treated as positive and
the others as negative (as in the R implementation).
"""
function f_disp(x::AggregatedData; trait_db = nothing, tax_lev = :Taxa, type = nothing, col_blocks = nothing,
                nbdim = 2, distance = "gower", zerodist_rm::Bool = false, correction = "none",
                traceB::Bool = false, set_param = nothing)
    sp = _set_param(set_param)
    s = _functional_setup(x, trait_db, tax_lev, type, col_blocks, distance, zerodist_rm, correction, sp)
    s.euclid || _noneuclid_msg()
    m = _nbdim(s.dist, nbdim, sp)
    out = LabeledVector(samples(x), _fdisp_k(s.dist.d, s.DF.abundance, m))
    traceB || return out
    return _functional_trace(out, s, m, correction)
end

# translation of the R helper `fdisp_k` (modified from FD::fdisp)
function _fdisp_k(D::Matrix{Float64}, A::Matrix{Float64}, m::Int)
    tol = 1e-7
    n = size(D, 1)
    Am = -0.5 .* D .^ 2
    for i in 1:n
        Am[i, i] = 0.0
    end
    G = bicenter(Am)
    eig, vectors = _eigen_desc(G)
    w0 = eig[n] / eig[1]
    r = w0 > -tol ? count(>(eig[1] * tol), eig) : length(eig)
    eig = eig[1:r]
    vectors = vectors[:, 1:r] .* sqrt.(abs.(eig))'
    pos = eig .> 0
    m > 0 && (pos = [pos[1:min(m, r)]; falses(max(r - m, 0))])
    any(vec(sum(A; dims = 1)) .== 0) && @warn "At least one community has zero-sum abundances (no species)."
    any(vec(sum(A; dims = 2)) .== 0) && @warn "At least one species does not occur in any community (zero total abundance across all communities)."
    res = zeros(size(A, 2))
    for j in axes(A, 2)
        pres = findall(>(0), view(A, :, j))
        vec_ = vectors[pres, :]
        nsp = length(unique(eachrow(vec_)))
        if nsp >= 2
            w = A[pres, j]
            centroid = vec(sum(vec_ .* w; dims = 1)) ./ sum(w)
            dpos = vec(sum((vec_[:, pos] .- centroid[pos]') .^ 2; dims = 2))
            dneg = any(.!pos) ? vec(sum((vec_[:, .!pos] .- centroid[.!pos]') .^ 2; dims = 2)) : zeros(length(pres))
            zij = sqrt.(abs.(dpos .- dneg))
            res[j] = sum(zij .* w) / sum(w)
        end
    end
    return res
end

"""
    f_eve(x::AggregatedData; trait_db, ...) -> LabeledVector

Functional evenness (Villéger et al. 2008): regularity of the
abundance-weighted minimum spanning tree linking the taxa in the PCoA
space. Samples with fewer than 3 taxa get `NaN`. Same arguments as
[`f_rich`](@ref).
"""
function f_eve(x::AggregatedData; trait_db = nothing, tax_lev = :Taxa, type = nothing, col_blocks = nothing,
               nbdim = 2, distance = "gower", zerodist_rm::Bool = false, correction = "none",
               traceB::Bool = false, set_param = nothing)
    sp = _set_param(set_param)
    s = _functional_setup(x, trait_db, tax_lev, type, col_blocks, distance, zerodist_rm, correction, sp)
    s.euclid || _noneuclid_msg()
    m = _nbdim(s.dist, nbdim, sp)
    pc = pcoa(s.dist; nf = m)
    C = pc.coords[:, 1:min(m, size(pc.coords, 2))]
    out = LabeledVector(samples(x), _feve_k(C, s.DF.abundance))
    traceB || return out
    return _functional_trace(out, s, m, correction)
end

# translation of the R helper `feve_k` (modified from FD::dbFD)
function _feve_k(C::Matrix{Float64}, A::Matrix{Float64})
    res = fill(NaN, size(A, 2))
    for j in axes(A, 2)
        pres = findall(>(0), view(A, :, j))
        S = length(pres)
        S > 2 || continue
        ab = A[pres, j]
        rel = ab ./ sum(ab)
        tr = C[pres, :]
        D = [sqrt(sum(abs2, view(tr, a, :) .- view(tr, b, :))) for a in 1:S, b in 1:S]
        N = _mst(D)
        EW = Float64[]
        for b in 1:S-1, a in b+1:S
            N[a, b] != 0 && push!(EW, D[a, b] / (rel[a] + rel[b]))
        end
        odsmo = 1 / (S - 1)
        minpew = [min(e / sum(EW), odsmo) for e in EW]
        res[j] = (sum(minpew) - odsmo) / (1 - odsmo)
    end
    return res
end

"""
    f_divs(x::AggregatedData; trait_db, tax_lev = :Taxa, type = nothing, col_blocks = nothing,
           distance = "gower", zerodist_rm = false, traceB = false,
           correction = "none", set_param = nothing) -> LabeledVector

Rao's quadratic entropy of each sample, scaled by its maximum
([`divc`](@ref) with `scale = true`). The trait distance must be Euclidean
(use `correction` otherwise).
"""
function f_divs(x::AggregatedData; trait_db = nothing, tax_lev = :Taxa, type = nothing, col_blocks = nothing,
                distance = "gower", zerodist_rm::Bool = false, traceB::Bool = false, correction = "none",
                set_param = nothing)
    sp = _set_param(set_param)
    s = _functional_setup(x, trait_db, tax_lev, type, col_blocks, distance, zerodist_rm, correction, sp)
    s.euclid || throw(ArgumentError("Non euclidean trait distance. Euclidean property is needed. Please use the correction options otherwise consider to remove taxa with the same traits."))
    out = LabeledVector(samples(x), divc(s.DF.abundance, s.dist; scale = true))
    traceB || return out
    return _functional_trace(out, s, nothing, correction)
end

"""
    f_red(x::AggregatedData; trait_db, ...) -> LabeledMatrix

Functional redundancy (de Bello et al. 2007): Gini-Simpson diversity
(`GS_rich`), scaled Rao's Q (`raoQ`) and their difference (`fred`,
negative values set to 0), samples × 3. Same arguments as
[`f_divs`](@ref).
"""
function f_red(x::AggregatedData; trait_db = nothing, tax_lev = :Taxa, type = nothing, col_blocks = nothing,
               distance = "gower", zerodist_rm::Bool = false, traceB::Bool = false, correction = "none",
               set_param = nothing)
    sp = _set_param(set_param)
    s = _functional_setup(x, trait_db, tax_lev, type, col_blocks, distance, zerodist_rm, correction, sp)
    s.euclid || throw(ArgumentError("Non euclidean trait distance. Euclidean property is needed. Please use the correction options otherwise consider to remove taxa with the same traits."))
    gs = divc(s.DF.abundance)
    rao = divc(s.DF.abundance, s.dist; scale = true)
    fred = max.(gs .- rao, 0.0)
    out = LabeledMatrix(samples(x), ["GS_rich", "raoQ", "fred"], hcat(gs, rao, fred))
    traceB || return out
    return _functional_trace(out, s, nothing, correction)
end

# ---------------------------------------------------------------------------
# Quality of functional spaces
# ---------------------------------------------------------------------------

"""
    qfs(d; nbdim = 15) -> (meanSD, coords, eig)

Quality of functional spaces (Maire et al. 2015): for 2 to `nbdim` PCoA
axes, the mean squared deviation between the initial distances and the
scaled Euclidean distances in the space (`meanSD`, labelled `m_2D`, ...).
"""
function qfs(d; nbdim::Integer = 15)
    D = _D(d)
    S = size(D, 1)
    S < 3 && throw(ArgumentError("there must be at least 3 species"))
    any(isnan, D) && throw(ArgumentError("NA are not allowed"))
    nbdim < 2 && throw(ArgumentError("'nbdim' must be higher than 1"))
    pc = pcoa(D; nf = nbdim)
    nb = min(nbdim, size(pc.coords, 2))
    C = pc.coords[:, 1:nb]
    x = lowertri(D)
    mx = maximum(x)
    vals = Float64[]
    for k in 2:nb
        y = lowertri([sqrt(sum(abs2, view(C, a, 1:k) .- view(C, b, 1:k))) for a in 1:S, b in 1:S])
        yst = y ./ maximum(y) .* mx
        push!(vals, round(sum((x .- yst) .^ 2) / (S * (S - 1) / 2); digits = 6))
    end
    return (meanSD = LabeledVector(["m_$(k)D" for k in 2:nb], vals), coords = C, eig = pc.eig)
end

function _r2_eig(eig, euclid)
    n = length(eig)
    euclid && return cumsum(eig) ./ sum(eig)
    neg = maximum(abs.(eig[eig .< 0]); init = 0.0)
    return (cumsum(eig) .+ (1:n) .* neg) ./ (sum(eig) + (n - 1) * neg)
end

function _pcoa_quality(D::Matrix{Float64}, type::String, method::String, tresh::Real, nbdim::Integer, nrepet, rng)
    euclid0 = is_euclid(D)
    Dt = euclid0 ? D : _D(_correct(D, type))
    trans = euclid0 ? "none" : type
    euclid = is_euclid(Dt)
    eu = euclid ? "yes" : "no"
    first_above(v) = (k = findfirst(>(tresh), v); k === nothing ? missing : k)
    if method == "cor"
        cm = _cmdscale(Dt)
        r2 = _r2_eig(cm.eig, euclid)
        rr = Float64[]
        pp = Float64[]
        for i in axes(cm.points, 2)
            pd = [sqrt(sum(abs2, view(cm.points, a, 1:i) .- view(cm.points, b, 1:i))) for a in axes(D, 1), b in axes(D, 1)]
            mt = mantel_test(pd, D; nrepet, rng)
            push!(rr, mt.obs)
            push!(pp, mt.pvalue)
        end
        k = if trans == "cailliez"
            tm = max(rr[1], tresh)
            count(<=(tm), rr)
        else
            first_above(rr)
        end
        k === missing && return (Transformation = trans, Euclidean = eu, nbdim = missing, cor = missing, p_value = missing, r2 = missing)
        return (Transformation = trans, Euclidean = eu, nbdim = k, cor = rr[k], p_value = pp[k], r2 = r2[k])
    elseif method == "legendre"
        cm = _cmdscale(Dt)
        r2 = _r2_eig(cm.eig, euclid)
        k = first_above(r2)
        return (Transformation = trans, Euclidean = eu, nbdim = k, r2 = k === missing ? missing : r2[k])
    elseif method == "maire"
        q = qfs(Dt; nbdim)
        k0 = findfirst(<(tresh), q.meanSD.values)
        k0 === nothing && throw(ArgumentError("there is no optimal number of dimension, please increase the number of dimensions"))
        k = k0 + 1
        eig = euclid0 ? q.eig : _cmdscale(Dt).eig
        r2 = euclid ? cumsum(eig) ./ sum(eig) : _r2_eig(eig, false)
        return (Transformation = trans, Euclidean = eu, nbdim = k, SD = q.meanSD.values[k0], r2 = r2[k])
    end
    throw(ArgumentError("method must be one of cor, legendre or maire"))
end

"""
    select_pcoa_axes(d; method = "legendre", tresh = nothing, nbdim = 15,
                     nrepet = 999, rng = Random.default_rng()) -> Vector{NamedTuple}

Choose the number of PCoA axes for a trait distance. For a non-Euclidean
distance one row is returned for each transformation (none, cailliez,
lingoes, sqrt, quasi).

- `"legendre"`: first axis where the cumulative eigenvalues exceed `tresh`
  (default 0.7; Legendre & Legendre 2012, with the correction for negative
  eigenvalues);
- `"cor"`: first axis where the Mantel correlation between the original
  and the PCoA distances exceeds `tresh` (0.7); p-values come from
  `nrepet` permutations, so they vary between runs;
- `"maire"`: first space with mean squared deviation below `tresh` (0.01;
  Maire et al. 2015).
"""
function select_pcoa_axes(d; method = "legendre", tresh = nothing, nbdim::Integer = 15, nrepet::Integer = 999,
                          rng::AbstractRNG = Random.default_rng())
    m = string(method)
    t = tresh === nothing ? (m == "maire" ? 0.01 : 0.7) : tresh
    D = _D(d)
    types = is_euclid(D) ? ["none"] : ["none", "cailliez", "lingoes", "sqrt", "quasi"]
    return [_pcoa_quality(D, ty, m, t, nbdim, nrepet, rng) for ty in types]
end
