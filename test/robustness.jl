# Regressions for BM6 (unchecked invariants and input validation) and BM7
# (shape of web-service responses). No test makes network requests: the
# responses are simulated through Biomonitor._JSON_SOURCE.

# Run `f` with simulated responses: `table` maps a URL prefix to a JSON
# string (parsed as the services would) or to `nothing` (empty response or
# HTTP 4xx). The longest matching prefix wins; unknown URLs fail the test.
function with_responses(f, table)
    source = function (url)
        keys_ = sort!([k for k in keys(table) if startswith(url, k)]; by = length, rev = true)
        isempty(keys_) && error("unexpected request: $url")
        v = table[first(keys_)]
        v isa Exception && throw(v)
        return v === nothing ? nothing : Biomonitor.parse_json(v)
    end
    old = Biomonitor._JSON_SOURCE[]
    Biomonitor._JSON_SOURCE[] = source
    try
        return f()
    finally
        Biomonitor._JSON_SOURCE[] = old
    end
end

const GBIF = "https://api.gbif.org/v1/species/"
const WORMS = "https://www.marinespecies.org/rest/"
const NBN = "https://species-ws.nbnatlas.org/"

const GBIF_RHODANI = """{"matchType": "EXACT", "synonym": false, "rank": "SPECIES", "status": "ACCEPTED",
    "phylum": "Arthropoda", "class": "Insecta", "order": "Ephemeroptera", "family": "Baetidae",
    "genus": "Baetis", "species": "Baetis rhodani", "canonicalName": "Baetis rhodani"}"""
const WORMS_RHODANI = """{"status": "accepted", "valid_AphiaID": 101, "valid_name": "Baetis rhodani",
    "phylum": "Arthropoda", "class": "Insecta", "order": "Ephemeroptera", "family": "Baetidae",
    "genus": "Baetis"}"""
const NBN_CLASSIFICATION = """[{"rank": "phylum", "scientificName": "Arthropoda"},
    {"rank": "class", "scientificName": "Insecta"}, {"rank": "order", "scientificName": "Ephemeroptera"},
    {"rank": "family", "scientificName": "Baetidae"}, {"rank": "genus", "scientificName": "Baetis"},
    {"rank": "species", "scientificName": "Baetis rhodani"}]"""

@testset "Web services: names in URLs are percent-encoded" begin
    @test Biomonitor._urlname(" Baetis rhodani ") == "Baetis%20rhodani"
    @test Biomonitor._urlname("A&b=c?d/e") == "A%26b%3Dc%3Fd%2Fe"
    @test Biomonitor._urlname("Ä-x_y.z~") == "%C3%84-x_y.z~"
end

@testset "Web services: empty and missing responses" begin
    mktempdir() do dir
        url(p) = "file:///" * replace(p, '\\' => '/')
        empty = joinpath(dir, "empty.json"); write(empty, "  \n")
        full = joinpath(dir, "full.json"); write(full, """{"a": [1, 2]}""")
        @test Biomonitor._download_json(url(empty)) === nothing
        @test Biomonitor._download_json(url(full)) == Dict{String,Any}("a" => Any[1, 2])
        # a failed transfer is not an HTTP 4xx answer: it must not become "not found"
        @test_throws Exception Biomonitor._download_json(url(joinpath(dir, "missing.json")))
    end
end

@testset "GBIF: simulated responses" begin
    r = with_responses(Dict(GBIF * "match" => GBIF_RHODANI)) do
        get_gbif_taxa_tree("Baetis rhodani"; verbose = false)
    end
    row = only(r.taxonomy)
    @test (row.Phylum, row.Family, row.Genus, row.Species, row.Taxa) ==
          ("Arthropoda", "Baetidae", "Baetis", "Baetis rhodani", "Baetis rhodani")
    @test r.status == ["Baetis rhodani" => "ACCEPTED"] && isempty(r.not_found) && isempty(r.synonyms)

    syn = """{"matchType": "EXACT", "synonym": true, "acceptedUsageKey": 42, "rank": "SPECIES"}"""
    r = with_responses(Dict(GBIF * "match" => syn, GBIF * "42" => GBIF_RHODANI)) do
        get_gbif_taxa_tree("Baetis vernus old"; verbose = false)
    end
    @test r.synonyms == ["Baetis vernus old" => "Baetis rhodani"]
    @test only(r.taxonomy).Taxa == "Baetis rhodani"

    r = with_responses(Dict(GBIF * "match" => """{"matchType": "NONE"}""")) do
        get_gbif_taxa_tree(["Nonexistent"]; verbose = false)
    end
    @test r.not_found == ["Nonexistent"] && isempty(r.taxonomy)

    for bad in ("""{"matchType": "EXACT", "synonym": true}""",            # no accepted key
                """{"matchType": "EXACT", "synonym": true, "acceptedUsageKey": null}""",
                """{"matchType": "EXACT", "synonym": true, "acceptedUsageKey": [1]}""",
                """[1, 2]""", """ "text" """, nothing)
        err = try
            with_responses(() -> get_gbif_taxa_tree("Baetis"; verbose = false), Dict(GBIF * "match" => bad))
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("GBIF", err.msg) && occursin("Baetis", err.msg)
    end
    # accepted record that is not an object
    @test_throws ArgumentError with_responses(Dict(GBIF * "match" => """{"matchType": "EXACT", "synonym": true,
                                                       "acceptedUsageKey": "7"}""", GBIF * "7" => "[]")) do
        get_gbif_taxa_tree("Baetis"; verbose = false)
    end
end

@testset "WoRMS: simulated responses" begin
    byname = WORMS * "AphiaIDByName/"
    record = WORMS * "AphiaRecordByAphiaID/"
    r = with_responses(Dict(byname => "101", record * "101" => WORMS_RHODANI)) do
        get_worms_taxa_tree("Baetis rhodani"; verbose = false)
    end
    row = only(r.taxonomy)
    @test (row.Order, row.Genus, row.Species, row.Taxa) == ("Ephemeroptera", "Baetis", "Baetis rhodani", "Baetis rhodani")

    r = with_responses(Dict(byname => "-999")) do
        get_worms_taxa_tree("Ambiguous"; verbose = false)
    end
    @test r.multiple_matches == ["Ambiguous"] && isempty(r.taxonomy)

    r = with_responses(Dict(byname => nothing)) do
        get_worms_taxa_tree("Nonexistent"; verbose = false)
    end
    @test r.not_found == ["Nonexistent"]

    unaccepted = """{"status": "unaccepted", "valid_AphiaID": 101, "valid_name": "Baetis rhodani"}"""
    r = with_responses(Dict(byname => "100", record * "100" => unaccepted, record * "101" => WORMS_RHODANI)) do
        get_worms_taxa_tree("Old name"; verbose = false)
    end
    @test r.synonyms == ["Old name" => "Baetis rhodani"] && only(r.taxonomy).Taxa == "Baetis rhodani"

    # unaccepted without a valid AphiaID: warning and not found, the batch goes on
    dubium = """{"status": "nomen dubium", "valid_AphiaID": null}"""
    table = Dict(byname * "Dubium" => "100", record * "100" => dubium,
                 byname * "Baetis" => "101", record * "101" => WORMS_RHODANI)
    run_dubium() = with_responses(() -> get_worms_taxa_tree(["Dubium", "Baetis rhodani"]; verbose = false), table)
    r = @test_logs (:warn, r"no valid AphiaID") run_dubium()
    @test r.not_found == ["Dubium"] && length(r.taxonomy) == 1

    for (table, what) in ((Dict(byname => "\"101\""), "AphiaID"),
                          (Dict(byname => "{\"id\": 1}"), "AphiaID"),
                          (Dict(byname => "true"), "AphiaID"),
                          (Dict(byname => "101", record * "101" => "[1]"), "JSON object"),
                          (Dict(byname => "101", record * "101" => nothing), "empty response"),
                          (Dict(byname => "100", record * "100" => """{"status": "unaccepted", "valid_AphiaID": true}"""),
                           "valid_AphiaID"))
        err = try
            with_responses(() -> get_worms_taxa_tree("Baetis"; verbose = false), table)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("WoRMS", err.msg) && occursin(what, err.msg)
    end
    # connection failures are no longer reported as taxa not found
    @test_throws ErrorException with_responses(Dict(byname => ErrorException("connection refused"))) do
        get_worms_taxa_tree("Baetis"; verbose = false)
    end
end

@testset "NBN Atlas: simulated responses" begin
    batch = NBN * "guid/batch"
    cl = NBN * "classification/"
    one = """{"Baetis rhodani": [{"identifier": "NBNSYS1", "acceptedIdentifier": "NBNSYS1"}]}"""
    r = with_responses(Dict(batch => one, cl * "NBNSYS1" => NBN_CLASSIFICATION)) do
        get_nbn_taxa_tree("Baetis rhodani"; verbose = false)
    end
    row = only(r.taxonomy)
    @test (row.Class, row.Family, row.Species, row.Taxa) == ("Insecta", "Baetidae", "Baetis rhodani", "Baetis rhodani")
    @test isempty(r.synonyms)

    syn = """{"Old": [{"identifier": "NBNSYS0", "acceptedIdentifier": "NBNSYS1"}]}"""
    r = with_responses(Dict(batch => syn, cl * "NBNSYS1" => NBN_CLASSIFICATION)) do
        get_nbn_taxa_tree("Old"; verbose = false)
    end
    @test r.synonyms == ["Old" => "Baetis rhodani"]

    # no acceptedIdentifier: the identifier itself is used
    noacc = """{"Baetis rhodani": [{"identifier": "NBNSYS1"}]}"""
    r = with_responses(Dict(batch => noacc, cl * "NBNSYS1" => NBN_CLASSIFICATION)) do
        get_nbn_taxa_tree("Baetis rhodani"; verbose = false)
    end
    @test only(r.taxonomy).Taxa == "Baetis rhodani" && isempty(r.synonyms)

    two = """{"Baetis": [{"identifier": "A", "acceptedIdentifier": "A"}, {"identifier": "B", "acceptedIdentifier": "B"}]}"""
    r = with_responses(Dict(batch => two)) do
        get_nbn_taxa_tree("Baetis"; verbose = false)
    end
    @test r.multiple_matches == ["Baetis"] && isempty(r.taxonomy)

    for none in (nothing, "{}", """{"X": null}""", """{"X": []}""")
        r = with_responses(Dict(batch => none)) do
            get_nbn_taxa_tree("X"; verbose = false)
        end
        @test r.not_found == ["X"]
    end

    for table in (Dict(batch => "[1]"),
                  Dict(batch => """{"X": {"identifier": "A"}}"""),
                  Dict(batch => """{"X": ["A"]}"""),
                  Dict(batch => """{"X": [{"acceptedIdentifier": "A"}]}"""),
                  Dict(batch => one, cl * "NBNSYS1" => """{"rank": "phylum"}"""),
                  Dict(batch => one, cl * "NBNSYS1" => """["phylum"]"""),
                  Dict(batch => one, cl * "NBNSYS1" => nothing))
        err = try
            with_responses(() -> get_nbn_taxa_tree("X"; verbose = false), table)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("NBN Atlas", err.msg)
    end
end

@testset "BM6: invariants and custom trait blocks" begin
    q(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
    agg = q(() -> aggregate_taxa(as_biomonitor(macro_ex())))

    # bio_index: a custom index equal to trophic_t gives the same values
    order = ["TROPHIC_1", "TROPHIC_2", "TROPHIC_3"]
    ref = bio_index(agg; index = "trophic_t").trophic_t
    cus = bio_index(agg; index = "custom", num = ["TROPHIC_1"], den = order, col_blocks = [3], col_order = order).custom
    @test cus.values == ref.values
    # a single name is a one-element vector, a single block an integer
    one = bio_index(agg; index = "custom", num = "TROPHIC_1", den = order, col_blocks = 3, col_order = order).custom
    @test one.values == ref.values
    for blocks in ([2], [4], [0, 3], [1, 1, 2], [3, -1, 1], Int[])
        err = try
            bio_index(agg; index = "custom", num = ["TROPHIC_1"], den = order, col_blocks = blocks, col_order = order)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("col_blocks", err.msg)
    end

    # convert_to_biotic: Oligochaeta without a family are kept as one taxon;
    # an AggregatedData built by hand without the Oligochaeta row is rejected
    d = agg.data
    pos = Biomonitor.positions(d.ranks, "Oligochaeta")
    @test !isempty(pos)
    R = copy(d.ranks)
    R[unique(r for (r, _) in pos), Biomonitor.rank_index(:Family)] .= ""
    agg2 = aggregate_taxa(Biomonitor.BiomonitorData(R, d.taxa, d.samples, d.abundance, d.group,
                                                    d.presence_absence, d.corrected_names, d.suggested_names))
    bt = convert_to_biotic(agg2)
    @test "Oligochaeta" in bt.taxa
    lev = Biomonitor.RANKS[first(pos)[2]]
    levels = copy(agg2.levels)
    t = levels[lev]
    keep = t.taxa .!= "Oligochaeta"
    levels[lev] = CommunityTable(t.taxa[keep], t.samples, t.abundance[keep, :])
    @test_throws ArgumentError convert_to_biotic(Biomonitor.AggregatedData(agg2.data, levels))
    @test convert_to_biotic(agg) isa CommunityTable

    # dness: unsupported levels are rejected before any computation
    @test_throws ArgumentError dness(agg; tax_lev = :Subfamily)
    @test_throws ArgumentError dness(agg; tax_lev = :Taxa)
end
