# Taxonomic trees from online services (GBIF, WoRMS, NBN Atlas).
#
# The R package also queries the IUCN Red List API v3, which has been
# retired by IUCN; that function is not ported.

# Download and parse a JSON response. An empty body or an HTTP 4xx status
# (e.g. WoRMS answers 204 and NBN 404 for unknown names) gives `nothing`;
# network failures and server errors are rethrown, so that they are not
# reported as taxa not found.
function _download_json(url::AbstractString)
    buf = IOBuffer()
    try
        Downloads.download(url, buf)
    catch err
        err isa Downloads.RequestError && 400 <= err.response.status < 500 && return nothing
        rethrow()
    end
    body = String(take!(buf))
    isempty(strip(body)) && return nothing
    return parse_json(body)
end

# Source of the parsed responses (url -> value). Tests replace it with
# simulated responses, so that they make no network requests.
const _JSON_SOURCE = Ref{Function}(_download_json)
_get_json(url::AbstractString) = _JSON_SOURCE[](url)

# Percent-encode a taxon name for a URL path or query (RFC 3986 unreserved
# characters are kept, the UTF-8 bytes of the others are encoded).
function _urlname(s)
    io = IOBuffer()
    for byte in codeunits(strip(String(s)))
        c = Char(byte)
        if byte < 0x80 && (isletter(c) || isdigit(c) || c in ('-', '.', '_', '~'))
            write(io, c)
        else
            print(io, '%', uppercase(string(byte; base = 16, pad = 2)))
        end
    end
    return String(take!(io))
end

# Shape checks on the parsed responses, with the service and the taxon in the message
_bad_response(service, taxon, what) =
    ArgumentError("$service: unexpected response for \"$taxon\": $what")
_describe(x) = x === nothing ? "an empty response" : string(typeof(x))
function _jobject(x, service, taxon)
    x isa AbstractDict && return x
    throw(_bad_response(service, taxon, "expected a JSON object, got $(_describe(x))"))
end
function _jid(d::AbstractDict, key, service, taxon)
    v = get(d, key, nothing)
    v isa Integer && !(v isa Bool) && return string(v)
    v isa AbstractString && !isempty(strip(v)) && return String(strip(v))
    throw(_bad_response(service, taxon, "missing or invalid \"$key\""))
end
# optional identifier: nothing when absent or null
_jid_or_nothing(d::AbstractDict, key, service, taxon) =
    get(d, key, nothing) === nothing ? nothing : _jid(d, key, service, taxon)

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
        c1 = _jobject(_get_json("https://api.gbif.org/v1/species/match?&name=" * _urlname(t)), "GBIF", t)
        if get(c1, "matchType", "NONE") != "NONE"
            if get(c1, "synonym", false) === true
                key = _jid(c1, "acceptedUsageKey", "GBIF", t)
                c2 = _jobject(_get_json("https://api.gbif.org/v1/species/" * key), "GBIF", t)
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
        aphia = _get_json("https://www.marinespecies.org/rest/AphiaIDByName/" * _urlname(t) * "?marine_only=false")
        aphia === nothing || (aphia isa Integer && !(aphia isa Bool)) ||
            throw(_bad_response("WoRMS", t, "expected an AphiaID, got $(_describe(aphia))"))
        if aphia === nothing
            push!(notfound, t)
        elseif aphia == -999
            push!(multiple, t)
        else
            c1 = _jobject(_get_json("https://www.marinespecies.org/rest/AphiaRecordByAphiaID/" * string(aphia)), "WoRMS", t)
            if get(c1, "status", "") != "accepted"
                valid = _jid_or_nothing(c1, "valid_AphiaID", "WoRMS", t)
                if valid === nothing
                    # unaccepted record without a valid name (e.g. nomen dubium)
                    @warn "WoRMS: \"$t\" is not accepted and has no valid AphiaID; reported as not found"
                    push!(notfound, t)
                    continue
                end
                c2 = _jobject(_get_json("https://www.marinespecies.org/rest/AphiaRecordByAphiaID/" * valid), "WoRMS", t)
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

function _nbn_row(classification, taxon)
    classification isa AbstractVector ||
        throw(_bad_response("NBN Atlas", taxon, "expected a classification array, got $(_describe(classification))"))
    all(c -> c isa AbstractDict, classification) ||
        throw(_bad_response("NBN Atlas", taxon, "classification entries must be JSON objects"))
    byrank = Dict{String,String}()
    for c in classification
        r = get(c, "rank", nothing)
        r === nothing && continue
        haskey(byrank, string(r)) || (byrank[string(r)] = string(get(c, "scientificName", "")))
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
        c1 = _get_json("https://species-ws.nbnatlas.org/guid/batch?q=" * _urlname(t))
        c1 === nothing || c1 isa AbstractDict ||
            throw(_bad_response("NBN Atlas", t, "expected a JSON object, got $(_describe(c1))"))
        matches = c1 isa AbstractDict && !isempty(c1) ? first(values(c1)) : nothing
        matches === nothing || matches isa AbstractVector ||
            throw(_bad_response("NBN Atlas", t, "expected an array of matches, got $(_describe(matches))"))
        matches isa AbstractVector && !all(m -> m isa AbstractDict, matches) &&
            throw(_bad_response("NBN Atlas", t, "matches must be JSON objects"))
        if matches isa AbstractVector && length(matches) == 1
            m = matches[1]
            id = _jid(m, "identifier", "NBN Atlas", t)
            accepted = something(_jid_or_nothing(m, "acceptedIdentifier", "NBN Atlas", t), id)
            cl = _get_json("https://species-ws.nbnatlas.org/classification/" * accepted)
            row = _nbn_row(cl, t)
            if id != accepted
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
                mk = matches[k]
                accepted = something(_jid_or_nothing(mk, "acceptedIdentifier", "NBN Atlas", t),
                                     _jid(mk, "identifier", "NBN Atlas", t))
                cl = _get_json("https://species-ws.nbnatlas.org/classification/" * accepted)
                push!(rows, _nbn_row(cl, t))
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
