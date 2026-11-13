# Small helpers shared by the package.

"""
    capwords(s)

First character upper case, the rest lower case, surrounding spaces removed
(the R helper `capWords` followed by `trimws`).
"""
function capwords(s::AbstractString)
    t = strip(s)
    isempty(t) && return ""
    c = first(t)
    return string(uppercase(c), lowercase(t[nextind(t, 1):end]))
end

# Group rows of `A` by `keys` and reduce each group column-wise with `f`
# (`:sum` or `:bin`). Groups are returned sorted by key.
function group_rows(keys::AbstractVector{<:AbstractString}, A::AbstractMatrix; f::Symbol = :sum)
    idx = Dict{String,Vector{Int}}()
    for (i, k) in enumerate(keys)
        push!(get!(idx, String(k), Int[]), i)
    end
    names = sort!(collect(Base.keys(idx)))
    out = zeros(Float64, length(names), size(A, 2))
    for (g, k) in enumerate(names)
        rows = idx[k]
        for j in axes(A, 2)
            if f === :sum
                s = 0.0
                for r in rows
                    s += A[r, j]
                end
                out[g, j] = s
            elseif f === :bin
                out[g, j] = any(r -> A[r, j] > 0, rows) ? 1.0 : 0.0
            elseif f === :mean
                s = 0.0
                for r in rows
                    s += A[r, j]
                end
                out[g, j] = s / length(rows)
            else
                throw(ArgumentError("unknown aggregation $(f)"))
            end
        end
    end
    return names, out
end

"""
    bin(x)

Aggregation function for presence-absence data: `1.0` if any value of `x`
is positive, `0.0` otherwise. Pass it as `FUN = bin` to
[`as_biomonitor`](@ref).
"""
bin(x) = any(>(0), x) ? 1.0 : 0.0

_funsym(f) = f === sum ? :sum : f === bin ? :bin : f isa Symbol ? f :
    throw(ArgumentError("FUN must be sum or bin"))

# R `cut(x, breaks, labels = 1:k, right = TRUE)`: returns the label (1-based
# interval index) or 0 when `x` falls outside all intervals (R gives NA, which
# biomonitoR then replaces with 0).
function cut_class(x::Real, breaks::AbstractVector{<:Real}; include_lowest::Bool = false)
    isnan(x) && return 0
    b = sort(Float64.(breaks))
    for k in 1:length(b)-1
        lo, hi = b[k], b[k+1]
        if (x > lo || (include_lowest && k == 1 && x == lo)) && x <= hi
            return k
        end
    end
    return 0
end

"abundance classes as used by WHPT, LIFE, PSI and E-PSI: breaks `[abucl; 1e18]`, lowest included"
abundance_class_matrix(A::AbstractMatrix, abucl::AbstractVector) =
    Float64[cut_class(a, [abucl; 1e18]; include_lowest = true) for a in A]

# ---------------------------------------------------------------------------
# log-gamma (Stirling series with upward recurrence), accurate to ~1e-15
# ---------------------------------------------------------------------------

function _lgamma(x::Float64)
    x <= 0 && return x == 0 ? Inf : throw(DomainError(x, "lgamma is defined here only for x > 0"))
    isinf(x) && return Inf
    shift = 0.0
    while x < 15.0
        shift += log(x)
        x += 1.0
    end
    z = 1 / (x * x)
    series = (1 / 12 + z * (-1 / 360 + z * (1 / 1260 + z * (-1 / 1680 + z * (1 / 1188 + z * (-691 / 360360 + z * (1 / 156))))))) / x
    return (x - 0.5) * log(x) - x + 0.5 * log(2π) + series - shift
end

"`log(x!)` for real `x ≥ 0` (R `lfactorial`)"
lfactorial(x::Real) = _lgamma(Float64(x) + 1.0)

# ---------------------------------------------------------------------------
# Brent minimisation, a line-by-line translation of R's `Brent_fmin`
# (src/library/stats/src/optimize.c) used by `optimize()`.
# ---------------------------------------------------------------------------

function brent_fmin(f, ax::Float64, bx::Float64, tol::Float64)
    c = (3.0 - sqrt(5.0)) * 0.5
    eps = sqrt(eps_double())
    a = ax; b = bx
    v = a + c * (b - a)
    w = v; x = v
    d = 0.0; e = 0.0
    fx = f(x)
    fv = fx; fw = fx
    tol3 = tol / 3.0
    while true
        xm = (a + b) * 0.5
        tol1 = eps * abs(x) + tol3
        t2 = tol1 * 2.0
        abs(x - xm) <= t2 - (b - a) * 0.5 && break
        p = 0.0; q = 0.0; r = 0.0
        if abs(e) > tol1
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = (q - r) * 2.0
            if q > 0.0
                p = -p
            else
                q = -q
            end
            r = e
            e = d
        end
        if abs(p) >= abs(q * 0.5 * r) || p <= q * (a - x) || p >= q * (b - x)
            e = x < xm ? b - x : a - x
            d = c * e
        else
            d = p / q
            u = x + d
            if u - a < t2 || b - u < t2
                d = tol1
                x >= xm && (d = -d)
            end
        end
        u = abs(d) >= tol1 ? x + d : (d > 0.0 ? x + tol1 : x - tol1)
        fu = f(u)
        if fu <= fx
            if u < x
                b = x
            else
                a = x
            end
            v = w; w = x; x = u
            fv = fw; fw = fx; fx = fu
        else
            if u < x
                a = u
            else
                b = u
            end
            if fu <= fw || w == x
                v = w; fv = fw
                w = u; fw = fu
            elseif fu <= fv || v == x || v == w
                v = u; fv = fu
            end
        end
    end
    return x
end
eps_double() = eps(Float64)

# all (row, col) positions of `name` in `M` in column-major order (R `which(M == name, arr.ind = TRUE)`)
function positions(M::AbstractMatrix{<:AbstractString}, name::AbstractString)
    out = Tuple{Int,Int}[]
    for j in axes(M, 2), i in axes(M, 1)
        M[i, j] == name && push!(out, (i, j))
    end
    return out
end

# k-combinations of 1:n in lexicographic order (R `combn`)
function combinations(n::Int, k::Int)
    out = Vector{Vector{Int}}()
    (k < 1 || k > n) && return out
    c = collect(1:k)
    while true
        push!(out, copy(c))
        i = k
        while i >= 1 && c[i] == n - k + i
            i -= 1
        end
        i == 0 && break
        c[i] += 1
        for j in i+1:k
            c[j] = c[j-1] + 1
        end
    end
    return out
end

_sym(x) = x isa Symbol ? x : Symbol(x)
_str(x) = x isa AbstractString ? String(x) : string(x)

# block sizes must be positive and cover the `n` columns exactly
function _check_blocks(blocks::AbstractVector{<:Integer}, n::Integer; what = "blocks", of = "columns")
    remaining = n
    for b in blocks
        0 < b <= remaining || throw(ArgumentError("$what must be positive and add up to the length of $of ($n)"))
        remaining -= b
    end
    remaining == 0 || throw(ArgumentError("$what must add up to the length of $of ($n)"))
    return nothing
end
