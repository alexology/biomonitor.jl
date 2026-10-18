# Spelling suggestions for taxon names not found in the reference database.
#
# The R package uses hunspell with a dictionary of the reference taxa. Here
# the reference names are ranked by optimal string alignment distance
# (Damerau-Levenshtein with adjacent transpositions), compared case
# insensitively. Suggestions can therefore differ from hunspell in order and
# number, but they are always names of the reference database.

function osa_distance(a::AbstractString, b::AbstractString; maxdist::Int = typemax(Int))
    s = collect(a)
    t = collect(b)
    m, n = length(s), length(t)
    abs(m - n) > maxdist && return maxdist + 1
    prev2 = zeros(Int, n + 1)
    prev = collect(0:n)
    cur = zeros(Int, n + 1)
    for i in 1:m
        cur[1] = i
        rowmin = cur[1]
        for j in 1:n
            cost = s[i] == t[j] ? 0 : 1
            v = min(prev[j+1] + 1, cur[j] + 1, prev[j] + cost)
            if i > 1 && j > 1 && s[i] == t[j-1] && s[i-1] == t[j]
                v = min(v, prev2[j-1] + 1)
            end
            cur[j+1] = v
            rowmin = min(rowmin, v)
        end
        rowmin > maxdist && return maxdist + 1
        prev2, prev, cur = prev, cur, prev2
    end
    return prev[n+1]
end

"""
    suggest_names(name, candidates; maxdist = auto, n = 10) -> Vector{String}

Spelling suggestions for `name` among `candidates` (by default the taxa of
a reference database), ordered by edit distance and then alphabetically.
The default maximum distance is 1 for names up to 4 characters, 2 up to 8
and 3 for longer names.
"""
function suggest_names(name::AbstractString, candidates::AbstractVector{<:AbstractString};
                       maxdist::Union{Nothing,Int} = nothing, n::Int = 10)
    q = lowercase(strip(name))
    md = maxdist === nothing ? (length(q) <= 4 ? 1 : length(q) <= 8 ? 2 : 3) : maxdist
    scored = Tuple{Int,String}[]
    for c in unique(candidates)
        isempty(c) && continue
        d = osa_distance(q, lowercase(c); maxdist = md)
        d <= md && push!(scored, (d, String(c)))
    end
    sort!(scored)
    return [s[2] for s in scored[1:min(n, length(scored))]]
end
suggest_names(name::AbstractString, ref::ReferenceDB; kw...) = suggest_names(name, ref.taxa; kw...)
suggest_names(name::AbstractString, group::Union{Symbol,AbstractString} = :mi; kw...) =
    suggest_names(name, reference_db(group); kw...)
