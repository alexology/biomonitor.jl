# Taxonomic trees from online services (GBIF, WoRMS, NBN Atlas).
#
# The R package also queries the IUCN Red List API v3, which has been
# retired by IUCN; that function is not ported.

function _get_json(url::AbstractString)
    buf = IOBuffer()
    Downloads.download(url, buf)
    return parse_json(String(take!(buf)))
end

_urlname(s) = replace(strip(String(s)), " " => "%20")

# R helpers `om`, `only_char`, `countWords`, `specieName`
_om(x) = x === nothing ? missing : x
_only_char(x) = x === missing ? missing : (occursin(r"[^A-Za-z]", replace(string(x), " " => "")) ? missing : string(x))
_nwords(x) = x === missing ? 1 : length(split(string(x), ' '))
function _species_name(x)
    s = _only_char(_om(x))
    n = _nwords(s)
    n == 1 && return missing
    n == 2 && return string(x)
    w = split(string(s), ' ')
    return string(w[1], " ", w[2])
end

const _TreeRow = NamedTuple{(:Phylum, :Class, :Subclass, :Order, :Family, :Subfamily, :Tribus, :Genus, :Species, :Subspecies, :Taxa),
                            NTuple{11,Union{Missing,String}}}

_s(x) = x === missing || x === nothing ? missing : string(x)

function _tree_result(rows::Vector{_TreeRow}, ref_from_tree_flag::Bool)
    ref_from_tree_flag || return rows
    levels = [string.(RANKS)...]
    V = [something(r[k], "") for r in rows, k in 1:NRANKS]
    return ref_from_tree(levels, V)
end

"""
    get_gbif_taxa_tree(taxa; ref_from_tree = false) -> NamedTuple

Retrieve the taxonomy of `taxa` from the GBIF backbone
(`api.gbif.org/v1/species/match`). Synonyms are replaced by the accepted
names. Returns `taxonomy` (rows with Phylum ... Taxa, or a
[`ReferenceDB`](@ref) when `ref_from_tree = true`), `not_found`,
`synonyms` (synonym => accepted) and `status`. Requires internet access.
"""
function get_gbif_taxa_tree(taxa; ref_from_tree::Bool = false, verbose::Bool = true)
    tx = taxa isa AbstractString ? [String(taxa)] : String.(collect(taxa))
    rows = _TreeRow[]
    notfound = String[]
    synonyms = Pair{String,String}[]
    status = Pair{String,String}[]
    for (j, t) in enumerate(tx)
        verbose && println("---- ", j, " of ", length(tx), " ---- ", t, " ----")
        c1 = _get_json("https://api.gbif.org/v1/species/match?&name=" * _urlname(t))
        if get(c1, "matchType", "NONE") != "NONE"
            if get(c1, "synonym", false) === true
                c2 = _get_json("https://api.gbif.org/v1/species/" * string(c1["acceptedUsageKey"]))
                pushfirst!(synonyms, t => string(get(c2, "canonicalName", "")))
            else
                c2 = c1
            end
            g(k) = _only_char(_om(get(c2, k, nothing)))
            sub = get(c2, "rank", "") == "SUBSPECIES" ? replace(string(get(c1, "scientificName", "")), "subsp. " => "") : missing
            push!(rows, _TreeRow((g("phylum"), g("class"), missing, g("order"), g("family"), missing, missing,
                                  g("genus"), g("species"), sub, _s(get(c2, "canonicalName", nothing)))))
            push!(status, t => string(haskey(c2, "status") ? c2["status"] : get(c2, "taxonomicStatus", "")))
        else
            push!(notfound, t)
        end
    end
    verbose && !isempty(notfound) && println("Taxa not found: ", length(notfound))
    verbose && !isempty(synonyms) && println("Synonym detected: ", length(synonyms))
    return (taxonomy = _tree_result(rows, ref_from_tree), not_found = notfound, synonyms = synonyms, status = status)
end

"""
    get_worms_taxa_tree(taxa; ref_from_tree = false) -> NamedTuple

Retrieve the taxonomy of `taxa` from the World Register of Marine Species
(`marinespecies.org/rest`). Unaccepted names are replaced by the valid ones.
Returns `taxonomy`, `not_found`, `synonyms` and `multiple_matches` (names
with several AphiaIDs). Requires internet access.
"""
function get_worms_taxa_tree(taxa; ref_from_tree::Bool = false, verbose::Bool = true)
    tx = taxa isa AbstractString ? [String(taxa)] : String.(collect(taxa))
    rows = _TreeRow[]
    notfound = String[]
    synonyms = Pair{String,String}[]
    multiple = String[]
    for (j, t) in enumerate(tx)
        verbose && println("---- ", j, " of ", length(tx), " ---- ", t, " ----")
        aphia = try
            _get_json("https://www.marinespecies.org/rest/AphiaIDByName/" * _urlname(t) * "?marine_only=false")
        catch
            nothing
        end
        if aphia === nothing
            push!(notfound, t)
        elseif aphia == -999
            push!(multiple, t)
        else
            c1 = _get_json("https://www.marinespecies.org/rest/AphiaRecordByAphiaID/" * string(aphia))
            if get(c1, "status", "") != "accepted"
                c2 = _get_json("https://www.marinespecies.org/rest/AphiaRecordByAphiaID/" * string(c1["valid_AphiaID"]))
                pushfirst!(synonyms, t => string(get(c2, "valid_name", "")))
            else
                c2 = c1
            end
            g(k) = _only_char(_om(get(c2, k, nothing)))
            vn = get(c2, "valid_name", nothing)
            sp = _species_name(_only_char(_om(vn)))
            sub = _nwords(_only_char(_om(vn))) == 3 ? _s(vn) : missing
            push!(rows, _TreeRow((g("phylum"), g("class"), missing, g("order"), g("family"), missing, missing,
                                  g("genus"), _s(sp), sub, _s(vn))))
        end
    end
    verbose && !isempty(notfound) && println("Taxa not found: ", length(notfound))
    verbose && !isempty(synonyms) && println("Synonym detected: ", length(synonyms))
    verbose && !isempty(multiple) && println("Taxa with multiple match detected: ", length(multiple))
    return (taxonomy = _tree_result(rows, ref_from_tree), not_found = notfound, synonyms = synonyms,
            multiple_matches = multiple)
end

function _nbn_row(classification)
    byrank = Dict{String,String}()
    for c in classification
        r = get(c, "rank", nothing)
        r === nothing && continue
        haskey(byrank, r) || (byrank[r] = string(get(c, "scientificName", "")))
    end
    g(k) = haskey(byrank, k) ? _only_char(byrank[k]) : missing
    last_ = isempty(classification) ? missing : _s(get(classification[end], "scientificName", nothing))
    return _TreeRow((g("phylum"), g("class"), missing, g("order"), g("family"), missing, missing,
                     g("genus"), g("species"), g("subspecies"), last_))
end

"""
    get_nbn_taxa_tree(taxa; ref_from_tree = false, interactive = false) -> NamedTuple

Retrieve the taxonomy of `taxa` from the NBN Atlas
(`species-ws.nbnatlas.org`). Synonyms are replaced by the accepted names.
Names with several matches are listed in `multiple_matches`; with
`interactive = true` the accepted name is chosen on the terminal (as in R).
Requires internet access.
"""
function get_nbn_taxa_tree(taxa; ref_from_tree::Bool = false, interactive::Bool = false, verbose::Bool = true)
    tx = taxa isa AbstractString ? [String(taxa)] : String.(collect(taxa))
    rows = _TreeRow[]
    notfound = String[]
    synonyms = Pair{String,String}[]
    multiple = String[]
    for (j, t) in enumerate(tx)
        verbose && println("---- ", j, " of ", length(tx), " ---- ", t, " ----")
        c1 = try
            _get_json("https://species-ws.nbnatlas.org/guid/batch?q=" * _urlname(t))
        catch
            nothing
        end
        matches = c1 isa AbstractDict && !isempty(c1) ? first(values(c1)) : nothing
        if matches isa AbstractVector && length(matches) == 1
            m = matches[1]
            key = get(m, "identifier", nothing) == get(m, "acceptedIdentifier", nothing) ? m["identifier"] : m["acceptedIdentifier"]
            cl = _get_json("https://species-ws.nbnatlas.org/classification/" * string(key))
            row = _nbn_row(cl)
            if get(m, "identifier", nothing) != get(m, "acceptedIdentifier", nothing)
                pushfirst!(synonyms, t => something(row.Taxa, ""))
            end
            push!(rows, row)
        elseif matches isa AbstractVector && length(matches) > 1
            push!(multiple, t)
            if interactive
                names = [string(get(m, "acceptedName", get(m, "name", ""))) for m in matches]
                println("Choose options for ", t)
                for (k, nm) in enumerate(names)
                    println("  ", k, ": ", nm)
                end
                print("Selection: ")
                k = tryparse(Int, strip(readline()))
                (k === nothing || k < 1 || k > length(names)) && throw(ArgumentError("invalid selection"))
                cl = _get_json("https://species-ws.nbnatlas.org/classification/" * string(matches[k]["acceptedIdentifier"]))
                push!(rows, _nbn_row(cl))
            end
        else
            push!(notfound, t)
        end
    end
    verbose && !isempty(notfound) && println("Taxa not found: ", length(notfound))
    verbose && !isempty(synonyms) && println("Synonym detected: ", length(synonyms))
    verbose && !isempty(multiple) && println("Taxa with multiple match detected: ", length(multiple))
    return (taxonomy = _tree_result(rows, ref_from_tree), not_found = notfound, synonyms = synonyms,
            multiple_matches = multiple)
end
