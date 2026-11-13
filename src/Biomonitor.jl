"""
    Biomonitor

Native Julia port of the R package biomonitoR (version 0.9.5): import and
validation of community data against taxonomic reference databases,
aggregation at taxonomic levels, diversity indices, biotic indices
(BMWP, ASPT, WHPT, LIFE, PSI, E-PSI, DEHLI, IBMR, ...), trait assignment and
functional diversity (FRic, FDis, FEve, Rao's Q, functional redundancy).

All data are stored in Julia types ([`CommunityTable`](@ref),
[`ReferenceDB`](@ref), [`BiomonitorData`](@ref), [`AggregatedData`](@ref),
[`TraitTable`](@ref), [`AssignedTraits`](@ref), [`DistMatrix`](@ref),
[`LabeledVector`](@ref), [`LabeledMatrix`](@ref)); no R runtime is needed.
"""
module Biomonitor

using LinearAlgebra: LinearAlgebra, I, Symmetric, det, dot, eigen, eigvals, norm, nullspace, tr
using Random: Random, AbstractRNG, randperm
using Statistics: Statistics, cor, mean, std
import Downloads

# core types and helpers
export LabeledVector, LabeledMatrix, labels, rownames, colnames
export CommunityTable, ReferenceDB, BiomonitorData, AggregatedData
export TraitTable, AssignedTraits, DistMatrix
export taxa, samples, abundances, ranks, TAXONOMIC_LEVELS, RANKS

# data and I/O
export read_community, write_table, read_traits, read_reference, read_dist
export macro_ex, mi_prin, oglio, traits_mf, traits_tachet, diatom_traits, example_tree
export reference_db, get_reference_db, reference_datasets

# import and data management
export as_biomonitor, aggregate_taxa, bin, to_bin, quick_rename, suggest_names
export excluded_taxa, corrected_names, ref_from_tree, check_tree
export subset_taxa, remove_taxa, combine_taxa, convert_to_vegan, convert_to_biotic
export solve_ambiguous, abundance_classes, general_info
export get_taxa_abundance, get_taxa_richness

# diversity
export abundance, richness, allrich, allindices
export shannon, simpson, invsimpson, esimpson, berpar, invberpar, brill, fisher
export margalef, mcintosh, menhinick, pielou

# biotic indices
export bmwp, aspt, whpt, life, psi, epsi, dehli, ibmr, igold, ept, eptd, bioco
export show_scores, dness

# traits
export assign_traits, manage_traits, average_traits, sample_traits, add_bias_to_traits
export prep_fuzzy, fuzzy_trait_ratio, flow_t, bio_index, cwm, tsi, csi

# functional diversity and distances
export f_rich, f_disp, f_eve, f_divs, f_red, select_pcoa_axes, qfs
export gower_fuzzy, euclidean_dist, is_euclid, cailliez, lingoes, quasieuclid
export pcoa, divc, divcmax, mantel_test, bray_curtis, jaccard_dist, convex_hull_volume
export hclust, cutree, indicator_taxa, indval

# web services
export get_gbif_taxa_tree, get_worms_taxa_tree, get_nbn_taxa_tree

include("types.jl")
include("utils.jl")
include("io.jl")
include("spelling.jl")
include("import.jl")
include("management.jl")
include("diversity.jl")
include("biotic.jl")
include("dness.jl")
include("distances.jl")
include("hull.jl")
include("traits.jl")
include("functional.jl")
include("cluster.jl")
include("json.jl")
include("webservices.jl")

end # module
