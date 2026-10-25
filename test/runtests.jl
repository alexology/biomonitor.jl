using biomonitor
using Test
using Random

const quietlog = Base.CoreLogging.NullLogger()
q(f) = Base.CoreLogging.with_logger(f, quietlog)

@testset "Presence-absence conversion preserves analysis semantics" begin
    original = as_biomonitor(macro_ex())
    converted = to_bin(original)
    imported = q(() -> as_biomonitor(macro_ex(); FUN = bin))
    @test converted.presence_absence
    @test !original.presence_absence
    @test converted.abundance == imported.abundance
    @test to_bin(converted).abundance == converted.abundance
    @test to_bin(converted).presence_absence
    actual = aggregate_taxa(converted)
    expected = aggregate_taxa(imported)
    @test shannon(actual; tax_lev = :Family) ≈ shannon(expected; tax_lev = :Family)
    @test_throws ArgumentError get_taxa_abundance(actual, "Ephemeroptera")
end

@testset "Silhouette counts singleton clusters as 0" begin
    D = [0.0 0.1 0.9; 0.1 0.0 0.8; 0.9 0.8 0.0]
    rows = biomonitor._silhouette_rows(D, [1, 1, 2])
    @test [r.cluster for r in rows] == [1, 1, 2]
    @test sum(r.res for r in rows) / 3 ≈ 0.587963 atol = 1e-6
    # Standard silhouette (R cluster 2.1.8.1 on the same Ward.D2 tree): averages
    # 0.281627, 0.334973, 0.309103 for k = 2, 3, 4. The R helper chose k = 4.
    A = Float64[0 7 14 8 12; 14 9 6 0 3; 12 1 19 5 17]
    x = aggregate_taxa(as_biomonitor(CommunityTable(["Baetis", "Caenis", "Hydracarina"],
                                                   ["S$i" for i in 1:5], A)))
    res = indicator_taxa(x; nmax = 4, thresh = 0)
    @test res.k == 3
    @test res.avg_silhouette ≈ 0.334973 atol = 1e-6
    @test length(res.silhouette_by_cluster) == 3
end

@testset "biomonitor.jl" begin
    @testset "types and helpers" begin
        v = LabeledVector(["a", "b"], [1.0, 2.0])
        @test v["b"] == 2.0
        @test labels(v) == ["a", "b"]
        @test sum(v) == 3.0
        m = LabeledMatrix(["r1", "r2"], ["c1"], [1.0; 2.0;;])
        @test m["r2", "c1"] == 2.0
        @test m[:, "c1"]["r1"] == 1.0
        @test biomonitor.capwords("  bAETIS rhodani ") == "Baetis rhodani"
        @test biomonitor.lfactorial(10) ≈ log(factorial(10))
        @test biomonitor.lfactorial(0.5) ≈ 0.5 * log(π) - log(2) atol = 1e-14  # log(Γ(1.5))
        @test biomonitor.combinations(4, 2) == [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]]
        ct = CommunityTable((Taxa = ["Baetis", "Caenis"], S1 = [1, missing], S2 = [0, 3]))
        @test isnan(ct.abundance[2, 1])
        @test ct.samples == ["S1", "S2"]
        @test biomonitor.parse_json("""{"a": [1, 2.5, "x\\u00e9"], "b": null, "c": true}""") ==
              Dict("a" => Any[1, 2.5, "xé"], "b" => nothing, "c" => true)
    end

    @testset "import" begin
        x = macro_ex()
        asb = as_biomonitor(x)
        @test size(asb.abundance) == (34, 2)
        @test asb.group === :mi
        wrong = CommunityTable([x.taxa[1:end-1]; "Acentrela"], x.samples, x.abundance)
        wrong.taxa[findfirst(==("Acentrella"), wrong.taxa)] = "Acentrela"
        asbw = @test_logs (:info, r"Some taxa were excluded") match_mode = :any as_biomonitor(wrong)
        @test first(excluded_taxa(asbw)).first == "Acentrela"
        @test first(excluded_taxa(asbw)).second[1] == "Acentrella"
        fixed = q(() -> as_biomonitor(wrong; correct_names = Dict("Acentrela" => "Acentrella")))
        @test "Acentrella" in fixed.taxa
        @test_throws ArgumentError CommunityTable((A = ["x"], S = [1]))
        pa = to_bin(x)
        @test_logs (:warn, r"Presence-absence data detected") match_mode = :any as_biomonitor(pa)
        asbb = as_biomonitor(pa; FUN = bin)
        @test asbb.presence_absence
        # to_change
        tc = as_biomonitor(CommunityTable(["Hydracarina", "Baetis"], ["S"], [2.0; 3.0;;]))
        @test "Trombidiformes" in tc.taxa
        @test suggest_names("Baetiss") == ["Baetis"] || "Baetis" in suggest_names("Baetiss")
    end

    agg = q(() -> aggregate_taxa(as_biomonitor(macro_ex())))
    @testset "richness and abundance (R tests)" begin
        @test q(() -> richness(agg; tax_lev = :Genus)).values == [13, 20]
        @test richness(agg; tax_lev = :Taxa).values == vec(sum(macro_ex().abundance .> 0; dims = 1))
        @test abundance(agg).values == vec(sum(macro_ex().abundance; dims = 1))
        @test get_taxa_abundance(agg, "Ancylus").values == macro_ex()["Ancylus"].values
        @test_throws ArgumentError get_taxa_abundance(agg)
        @test_throws ArgumentError get_taxa_abundance(agg, "Ergo")
        @test_throws ArgumentError get_taxa_richness(agg, ["Ephemeroptera", "Trichoptera", "Diptera"], ["Family", "Genus"])
        @test_throws ArgumentError get_taxa_richness(agg, "Baetis", "Family")
        @test_throws ArgumentError get_taxa_richness(agg, "Epemeroptera", "Family")
        @test ept(agg).values == get_taxa_richness(agg, ["Ephemeroptera", "Plecoptera", "Trichoptera"], "Taxa").values
        @test ept(agg; tax_lev = :Family).values ==
              get_taxa_richness(agg, ["Ephemeroptera", "Plecoptera", "Trichoptera"], "Family").values
        @test igold(agg).values ≈ 1 .- get_taxa_abundance(agg, ["Oligochaeta", "Diptera", "Gastropoda"]; rel = true).values
        fams = ["Heptageniidae", "Leptophlebiidae", "Polycentropodidae", "Empididae"]
        @test eptd(agg).values ≈ log10.(get_taxa_abundance(agg, fams).values .+ 1)
        @test eptd(agg; eptd_families = fams, base = 2).values ≈ log2.(get_taxa_abundance(agg, fams).values .+ 1)
    end

    @testset "data management (R tests)" begin
        asb = q(() -> as_biomonitor(macro_ex()))
        h = subset_taxa(asb; taxa = "Hydropsychidae")
        @test all(==("Hydropsychidae"), h.ranks[:, 5])
        both = subset_taxa(asb; taxa = "Hydropsychidae", exclude = "Cheumatopsyche")
        @test all(==("Hydropsyche"), both.ranks[:, 8])
        @test_throws ArgumentError subset_taxa(asb)
        r = remove_taxa(agg, "Hydroptilidae")
        @test !("Hydroptila" in r.taxa) && length(r.taxa) == 32
        @test_throws ArgumentError remove_taxa(agg, "ignita")
        sub = CommunityTable(macro_ex().taxa[2:5], macro_ex().samples, macro_ex().abundance[2:5, :])
        ag2 = q(() -> aggregate_taxa(as_biomonitor(sub)))
        ct = combine_taxa(ag2; tax_lev = :Genus)
        @test length(ct.taxa) == 3
    end

    @testset "biotic indices (R tests)" begin
        tdir = joinpath(@__DIR__, "data")
        sp = q(() -> aggregate_taxa(as_biomonitor(biomonitor.community_from_raw(biomonitor.read_raw(joinpath(tdir, "aspt_bmwp_spain.tsv"))))))
        @test q(() -> bmwp(sp; method = "spa")).values == [35, 18]
        @test q(() -> aspt(sp; method = "spa")).values ≈ [35, 18] ./ [6, 3]
        @test q(() -> bmwp(sp; method = "spa", exceptions = "Ferrissia")).values == [29, 12]
        tb = q(() -> bmwp(sp; method = "spa", traceB = true))
        @test tb.parent_child_pairs == [(child = "Ferrissia", parent = "Planorbidae")]
        @test_logs (:info, r"Parent-child pairs found") match_mode = :any bmwp(sp; method = "spa")
        @test q(() -> bmwp(sp; method = "ita")).values == [35, 21]
        @test q(() -> bmwp(sp; method = "ita", agg = true)).values == [21, 21]
        @test q(() -> aspt(sp; method = "uk", agg = true)).values ≈ [21, 21] ./ [4, 4]
        @test q(() -> bmwp(sp; method = "uk", agg = true, traceB = true)).composite_taxa == ["Apataniidae", "Glossosomatidae"]
        sc = show_scores("aspt", "uk")
        @test q(() -> bmwp(agg; method = sc.scores, agg = sc.aggregation_rule)).values == q(() -> bmwp(agg; method = "uk", agg = true)).values
        @test epsi(agg).values ≈ [96.848485, 94.339623] atol = 1e-6
        @test_throws ArgumentError bmwp(agg; method = "xx")
        @test_throws ArgumentError whpt(agg; type = "zz")
        aggb = q(() -> aggregate_taxa(as_biomonitor(to_bin(macro_ex()); FUN = bin)))
        @test_throws ArgumentError whpt(aggb; type = "ab")
    end

    @testset "traits" begin
        ts = q(() -> assign_traits(agg))
        @test length(ts) == 71
        @test count(==(0), ts.distance) == 21
        @test_throws ArgumentError assign_traits(agg; tax_lev = :Order)
        mt = manage_traits(ts)
        av = average_traits(mt)
        @test length(av.taxa) == 34
        cb = biomonitor.DEFAULT_TACHET_BLOCKS
        @test all(isapprox.(sum(av.values[:, 1:8]; dims = 2), 1) .| isnan.(sum(av.values[:, 1:8]; dims = 2)))
        st = sample_traits(ts; rng = MersenneTwister(1))
        @test length(st.taxa) == length(unique(ts.taxa))
        b = add_bias_to_traits(av; col_blocks = cb, rng = MersenneTwister(2))
        @test size(b.values) == size(av.values)
        @test_throws ArgumentError add_bias_to_traits(av)
    end

    @testset "functional diversity errors" begin
        av = q(() -> average_traits(assign_traits(agg)))
        @test_throws ArgumentError f_rich(agg)
        @test_throws ArgumentError f_rich(agg; trait_db = av)
        @test_throws ArgumentError f_rich(agg; trait_db = av, type = "X")
        @test_throws ArgumentError f_rich(agg; trait_db = av, type = "C", distance = "gower")
        @test_throws ArgumentError f_rich(agg; trait_db = av, type = "F")
        @test_throws ArgumentError f_rich(agg; trait_db = av, type = "F", col_blocks = [1, 1])
        tb = q(() -> f_rich(agg; trait_db = av, type = "F", col_blocks = biomonitor.DEFAULT_TACHET_BLOCKS, traceB = true))
        @test length(tb) == 7
        # averaged Tachet traits contain empty blocks (NaN), as in R
        @test tb.NA_detection isa Vector
        @test first(tb.NA_detection).taxon == "Beraeamyia"
        cont = TraitTable(av.taxa, av.traits[1:3], replace(av.values[:, 1:3], NaN => 0.0))
        tc = q(() -> f_rich(agg; trait_db = cont, type = "C", distance = "euclidean", traceB = true))
        @test tc.NA_detection == "No NAs detected"
        dc = DistMatrix(cont.taxa, euclidean_dist(cont.values))
        td = q(() -> f_rich(agg; trait_db = dc, traceB = true))
        @test td.results ≈ tc.results
        @test td.NA_detection == "NAs cannot be detected when trait_db is a DistMatrix"
    end

    @testset "zero distances" begin
        av = q(() -> average_traits(assign_traits(agg)))
        cb = biomonitor.DEFAULT_TACHET_BLOCKS
        G = gower_fuzzy(prep_fuzzy(av, cb), cb)
        D = copy(G.d)
        D[4, 1] = D[1, 4] = 0
        G0 = DistMatrix(G.labels, D)
        r0 = q(() -> f_rich(agg; trait_db = G0, nbdim = 3, zerodist_rm = true, traceB = true))
        @test r0.duplicated_traits == [(taxon = G.labels[1], same_traits = G.labels[4])]
        keep = [i for i in eachindex(G.labels) if i != 4]
        mx = macro_ex()
        A = copy(mx.abundance)
        i1 = findfirst(==(G.labels[1]), mx.taxa)
        i4 = findfirst(==(G.labels[4]), mx.taxa)
        A[i1, :] .+= A[i4, :]
        rows = [i for i in eachindex(mx.taxa) if i != i4]
        agg0 = q(() -> aggregate_taxa(as_biomonitor(CommunityTable(mx.taxa[rows], mx.samples, A[rows, :]))))
        r1 = q(() -> f_rich(agg0; trait_db = DistMatrix(G.labels[keep], G.d[keep, keep]), nbdim = 3, traceB = true))
        @test r0.results.values ≈ r1.results.values
        @test r0.taxa == r1.taxa
    end

    @testset "Euclidean corrections, degenerate cases" begin
        # Lingoes corrects distances above tol and leaves the others unchanged
        # (ade4 compares d^2 with tol and returns d^2 for those it skips)
        D = [0.0 0.00005 1.0; 0.00005 0.0 2.0; 1.0 2.0 0.0]
        L = q(() -> lingoes(D))
        c = (L[1, 3]^2 - 1) / 2
        @test c > 0
        @test L[1, 2] ≈ sqrt(0.00005^2 + 2c)
        @test L[2, 3] ≈ sqrt(4 + 2c)
        D[1, 2] = D[2, 1] = 5e-8
        @test q(() -> lingoes(D))[1, 2] == 5e-8
        # Coincident points are Euclidean; the maximal Rao entropy is 0
        @test is_euclid(zeros(3, 3))
        m = divcmax(zeros(3, 3))
        @test m.value == 0.0
        @test m.weights ≈ fill(1 / 3, 3)
        # An explicit tol also decides whether a correction is needed
        D = [0.0 1.0 2.00000001; 1.0 0.0 1.0; 2.00000001 1.0 0.0]
        @test is_euclid(D) && !is_euclid(D; tol = 1e-10)
        @test is_euclid(lingoes(D; tol = 1e-10); tol = 1e-10)
        @test is_euclid(cailliez(D; tol = 1e-10); tol = 1e-10)
    end

    @testset "zero distances with missing traits" begin
        # A-B and B-C are at zero distance, A-C is not: only A and B are merged
        names = ["A", "B", "C"]
        ref = ReferenceDB(fill("", 3, 10), names)
        x = q(() -> aggregate_taxa(as_biomonitor(CommunityTable(names, ["S"], reshape([1.0, 2.0, 3.0], 3, 1));
                                                 dfref = ref, to_change = nothing)))
        tr = TraitTable(names, ["f1a", "f1b", "f2a", "f2b"], [1.0 0.0 1.0 0.0; 1.0 0.0 NaN NaN; 1.0 0.0 0.0 1.0])
        r = q(() -> f_disp(x; trait_db = tr, type = "F", col_blocks = [2, 2], zerodist_rm = true,
                           nbdim = 1, traceB = true))
        @test r.taxa.taxa == ["A", "C"]
        @test vec(r.taxa.abundance) == [3.0, 3.0]
        @test r.results.values ≈ [0.35355339059327373]
    end

    @testset "convex hull" begin
        cube = [x for x in Iterators.product(0:1, 0:1, 0:1)]
        P = reduce(vcat, [[c[1] c[2] c[3]] for c in vec(cube)])
        P = vcat(P, [0.5 0.5 0.5; 0.2 0.3 0.4])
        @test convex_hull_volume(P) ≈ 1.0
        @test convex_hull_volume([0 0; 1 0; 0 1; 1 1; 0.5 0.5]) ≈ 1.0
        @test convex_hull_volume([0 0 0; 1 0 0; 0 1 0; 1 1 0]) == 0.0
        rng = MersenneTwister(3)
        S = randn(rng, 300, 3)
        S ./= sqrt.(sum(abs2, S; dims = 2))
        @test isapprox(convex_hull_volume(S), 4π / 3; rtol = 0.05)
    end

    include("parity.jl")
end
