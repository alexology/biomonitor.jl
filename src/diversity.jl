# Diversity indices (R helper `Pi` and its wrappers).

function _pi_index(x::AbstractVector{<:Real}, index::Symbol; base::Real = ℯ)
    ntaxa = count(>(0), x)
    abu = sum(x)
    rel = x ./ abu
    if index === :shannon
        s = 0.0
        for p in rel
            v = p * log(base, p)
            isfinite(v) && (s += v)
        end
        return -s
    elseif index === :simpson
        return 1 - sum(abs2, rel)
    elseif index === :esimpson
        return (1 / sum(abs2, rel)) / ntaxa
    elseif index === :invsimpson
        return 1 / sum(abs2, rel)
    elseif index === :margalef
        return (ntaxa - 1) / log(abu)
    elseif index === :menhinick
        return ntaxa / sqrt(abu)
    elseif index === :berpar
        return maximum(x) / abu
    elseif index === :invberpar
        return 1 / (maximum(x) / abu)
    elseif index === :brill
        return (lfactorial(abu) - sum(lfactorial, x)) / abu
    elseif index === :mcintosh
        return (abu - sqrt(sum(abs2, x))) / (abu - sqrt(abu))
    elseif index === :fisher
        n = abu
        r = ntaxa
        fn(t) = abs((1 - t) * n / t * (-log(1 - t)) - r)
        t = brent_fmin(fn, 0.0, 1.0, 1e-9)
        return n * (1 - t) / t
    end
    throw(ArgumentError("unknown index $(index)"))
end

function _diversity(x::AggregatedData, index::Symbol, tax_lev; base = ℯ)
    t = _drop_unassigned(_level(_check_aggregated(x), tax_lev))
    vals = [_pi_index(view(t.abundance, :, j), index; base) for j in eachindex(t.samples)]
    return LabeledVector(t.samples, vals)
end

"""
    shannon(x::AggregatedData; base = ℯ, tax_lev = :Taxa) -> LabeledVector

Shannon index ``H = -\\sum p_i \\log p_i`` per sample.
"""
shannon(x::AggregatedData; base = ℯ, tax_lev = :Taxa) = _diversity(x, :shannon, tax_lev; base)

"""
    simpson(x::AggregatedData; tax_lev = :Taxa) -> LabeledVector

Simpson's index of diversity ``1 - \\sum p_i^2``.
"""
simpson(x::AggregatedData; tax_lev = :Taxa) = _diversity(x, :simpson, tax_lev)

"""
    invsimpson(x::AggregatedData; tax_lev = :Taxa) -> LabeledVector

Inverse Simpson index ``1 / \\sum p_i^2``.
"""
invsimpson(x::AggregatedData; tax_lev = :Taxa) = _diversity(x, :invsimpson, tax_lev)

"""
    esimpson(x::AggregatedData; tax_lev = :Taxa) -> LabeledVector

Simpson's evenness: inverse Simpson divided by richness.
"""
esimpson(x::AggregatedData; tax_lev = :Taxa) = _diversity(x, :esimpson, tax_lev)

"""
    berpar(x::AggregatedData; tax_lev = :Taxa) -> LabeledVector

Berger-Parker dominance: abundance of the most abundant taxon over the total.
"""
berpar(x::AggregatedData; tax_lev = :Taxa) = _diversity(x, :berpar, tax_lev)

"""
    invberpar(x::AggregatedData; tax_lev = :Taxa) -> LabeledVector

Inverse of the Berger-Parker index.
"""
invberpar(x::AggregatedData; tax_lev = :Taxa) = _diversity(x, :invberpar, tax_lev)

"""
    brill(x::AggregatedData; tax_lev = :Taxa) -> LabeledVector

Brillouin index ``(\\ln N! - \\sum \\ln n_i!) / N``.
"""
brill(x::AggregatedData; tax_lev = :Taxa) = _diversity(x, :brill, tax_lev)

"""
    fisher(x::AggregatedData; tax_lev = :Taxa) -> LabeledVector

Fisher's alpha, found by minimising ``|(1-x) N/x (-\\ln(1-x)) - S|`` on
(0, 1) with Brent's method, as R `optimize()` does.
"""
fisher(x::AggregatedData; tax_lev = :Taxa) = _diversity(x, :fisher, tax_lev)

"""
    margalef(x::AggregatedData; tax_lev = :Taxa) -> LabeledVector

Margalef index ``(S - 1) / \\ln N``.
"""
margalef(x::AggregatedData; tax_lev = :Taxa) = _diversity(x, :margalef, tax_lev)

"""
    mcintosh(x::AggregatedData; tax_lev = :Taxa) -> LabeledVector

McIntosh index ``(N - \\sqrt{\\sum n_i^2}) / (N - \\sqrt{N})``.
"""
mcintosh(x::AggregatedData; tax_lev = :Taxa) = _diversity(x, :mcintosh, tax_lev)

"""
    menhinick(x::AggregatedData; tax_lev = :Taxa) -> LabeledVector

Menhinick index ``S / \\sqrt{N}``.
"""
menhinick(x::AggregatedData; tax_lev = :Taxa) = _diversity(x, :menhinick, tax_lev)

"""
    pielou(x::AggregatedData; base = ℯ, tax_lev = :Taxa) -> LabeledVector

Pielou's evenness: Shannon index divided by ``\\log S``.
"""
function pielou(x::AggregatedData; base = ℯ, tax_lev = :Taxa)
    t = _drop_unassigned(_level(_check_aggregated(x), tax_lev))
    vals = map(eachindex(t.samples)) do j
        col = view(t.abundance, :, j)
        _pi_index(col, :shannon; base) / log(base, count(>(0), col))
    end
    return LabeledVector(t.samples, vals)
end

"""
    allindices(x::AggregatedData; tax_lev = :Taxa, base = ℯ) -> LabeledMatrix

All the diversity indices (samples × indices): shannon, berpar, brill,
invberpar, invsimpson, margalef, mcintosh, menhinick, pielou, simpson,
esimpson, fisher.
"""
function allindices(x::AggregatedData; tax_lev = :Taxa, base = ℯ)
    _check_aggregated(x)
    cols = [shannon(x; tax_lev, base), berpar(x; tax_lev), brill(x; tax_lev), invberpar(x; tax_lev),
            invsimpson(x; tax_lev), margalef(x; tax_lev), mcintosh(x; tax_lev), menhinick(x; tax_lev),
            pielou(x; tax_lev, base), simpson(x; tax_lev), esimpson(x; tax_lev), fisher(x; tax_lev)]
    M = reduce(hcat, [c.values for c in cols])
    return LabeledMatrix(samples(x), ["shannon", "berpar", "brill", "invberpar", "invsimpson", "margalef",
                                      "mcintosh", "menhinick", "pielou", "simpson", "esimpson", "fisher"],
                         reshape(M, length(samples(x)), :))
end
