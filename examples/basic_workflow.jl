# Basic biomonitor.jl workflow on the bundled Italian macroinvertebrate data.
#
#   julia --project=. examples/basic_workflow.jl [output_dir]
#
# Writes diversity, biotic and functional indices as TSV tables.

using biomonitor

outdir = length(ARGS) >= 1 ? ARGS[1] : "biomonitoR_results"
mkpath(outdir)

# 1. import and aggregate
comm = mi_prin()
asb = as_biomonitor(comm; group = :mi)
agg = aggregate_taxa(asb)
println(agg)

# 2. diversity
write_table(joinpath(outdir, "richness_by_level.tsv"), allrich(agg))
write_table(joinpath(outdir, "diversity_indices.tsv"), allindices(agg))

# 3. biotic indices
samples_ = samples(agg)
biotic = LabeledMatrix(samples_, ["BMWP_ita", "ASPT_ita", "WHPT_ASPT", "LIFE", "PSI", "EPSI", "EPT_families"],
                       hcat(bmwp(agg; method = "ita", agg = true).values,
                            aspt(agg; method = "ita", agg = true).values,
                            whpt(agg; agg = true).values,
                            life(agg).values,
                            psi(agg).values,
                            epsi(agg).values,
                            ept(agg; tax_lev = :Family).values))
write_table(joinpath(outdir, "biotic_indices.tsv"), biotic)

# 4. traits and functional diversity
traits = average_traits(manage_traits(assign_traits(agg)))
blocks = [8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8]
write_table(joinpath(outdir, "cwm.tsv"), cwm(agg; trait_db = traits))
fd = LabeledMatrix(samples_, ["FRic", "FDis", "FEve"],
                   hcat(f_rich(agg; trait_db = traits, type = "F", col_blocks = blocks, nbdim = 3, correction = "cailliez").values,
                        f_disp(agg; trait_db = traits, type = "F", col_blocks = blocks, nbdim = 3, correction = "cailliez").values,
                        f_eve(agg; trait_db = traits, type = "F", col_blocks = blocks, nbdim = 3, correction = "cailliez").values))
write_table(joinpath(outdir, "functional_diversity.tsv"), fd)
println("Results written to ", abspath(outdir))
