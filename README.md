# Biomonitor.jl

[![CI](https://github.com/alexology/Biomonitor.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/alexology/Biomonitor.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://alexology.github.io/Biomonitor.jl/dev/)

Native Julia port of the R package
[biomonitoR](https://github.com/alexology/biomonitoR) 0.9.5 (Laini et al.
2022): management of taxonomic and trait data and calculation of
biomonitoring indices and metrics for running waters.

This package is part of an experimental effort to use AI to port
bioinformatics and quantitative ecology packages to Julia, with the goal of
creating a unified environment for high performance scientific computing.

The package does not call R. Data are held in Julia types, and the
algorithms from ade4, ape, FD, geometry (Qhull) and stats that biomonitoR
relies on are reimplemented in Julia. It depends only on Julia standard
libraries.

## Installation

The package is not registered. Install it from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/alexology/Biomonitor.jl")
```

or, from a local clone, `Pkg.develop(path = "path/to/Biomonitor.jl")`.
Julia ≥ 1.10.

Up to the first registration the package was called `biomonitor`; the
module is now `Biomonitor` (same UUID), so replace `using biomonitor` with
`using Biomonitor`. The function names are unchanged.

## Quick start

```julia
using Biomonitor

x   = macro_ex()                 # CommunityTable: 34 taxa × 2 samples
asb = as_biomonitor(x)           # check names against the reference database
agg = aggregate_taxa(asb)        # tables at every taxonomic level

richness(agg; tax_lev = :Family) # LabeledVector, one value per sample
allindices(agg)                  # LabeledMatrix, samples × 12 diversity indices
aspt(agg; method = "ita", agg = true)
whpt(agg; type = "ab", metric = "aspt")
life(agg; method = "life_2017")

traits = average_traits(manage_traits(assign_traits(agg)))
cwm(agg; trait_db = traits)
f_rich(agg; trait_db = traits, type = "F",
       col_blocks = [8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8],
       nbdim = 3)
```

Your own data: `read_community("my_data.csv")`, with a `Taxa` column and one
numeric column per sample, or `CommunityTable((Taxa = [...], S1 = [...], ...))`
built from any `NamedTuple`/`Dict` of columns (for example
`Tables.columntable(df)` of a DataFrame).

## Data structures

| R object | Julia type |
|---|---|
| `data.frame` with `Taxa` and sample columns | `CommunityTable` (taxa, samples, abundance matrix; `NaN` = missing) |
| reference database (`mi_ref`, `dfref`) | `ReferenceDB` (n × 10 rank matrix + taxa) |
| `as_biomonitor()` result (class `asb`) | `BiomonitorData` |
| `aggregate_taxa()` result (class `biomonitoR`) | `AggregatedData`; `agg[:Genus]` is a `CommunityTable` |
| trait `data.frame` | `TraitTable` |
| `assign_traits()` result | `AssignedTraits` |
| `dist` object | `DistMatrix` |
| named vector of results | `LabeledVector` (`v["Sample_1"]`, `labels(v)`) |
| `data.frame` of results | `LabeledMatrix` (`m["Sample_1", "shannon"]`) |
| list returned with `traceB = TRUE` | `NamedTuple` |

Taxonomic levels are symbols (`:Family`) or strings (`"Family"`); see
`TAXONOMIC_LEVELS`.

## Functions

- **Import and data management**: `as_biomonitor`, `aggregate_taxa`,
  `quick_rename`, `suggest_names`, `excluded_taxa`, `corrected_names`,
  `ref_from_tree`, `check_tree`, `subset_taxa`, `remove_taxa`,
  `combine_taxa`, `solve_ambiguous`, `convert_to_vegan`,
  `convert_to_biotic`, `abundance_classes`, `general_info`, `to_bin`, `bin`.
- **Abundance and richness**: `abundance`, `richness`, `allrich`,
  `get_taxa_abundance`, `get_taxa_richness`.
- **Diversity**: `shannon`, `simpson`, `invsimpson`, `esimpson`, `berpar`,
  `invberpar`, `brill`, `fisher`, `margalef`, `mcintosh`, `menhinick`,
  `pielou`, `allindices`, `dness` (taxonomic distinctness).
- **Biotic indices**: `bmwp`, `aspt`, `whpt`, `life`, `psi`, `epsi`,
  `dehli`, `ibmr`, `igold`, `ept`, `eptd`, `bioco`; score tables with
  `show_scores` (editable `Dict`s accepted as custom methods).
- **Traits**: `assign_traits`, `manage_traits`, `average_traits`,
  `sample_traits`, `add_bias_to_traits`, `prep_fuzzy`, `cwm`, `tsi`,
  `csi`, `fuzzy_trait_ratio`, `flow_t`, `bio_index`.
- **Functional diversity**: `f_rich`, `f_disp`, `f_eve`, `f_divs`,
  `f_red`, `qfs`, `select_pcoa_axes`; tools `gower_fuzzy`,
  `euclidean_dist`, `is_euclid`, `cailliez`, `lingoes`, `quasieuclid`,
  `pcoa`, `divc`, `divcmax`, `mantel_test`, `convex_hull_volume`.
- **Clustering and indicator taxa**: `hclust`, `cutree`, `indval`,
  `indicator_taxa` (the analysis behind `plot_indicator_taxa`),
  `bray_curtis`, `jaccard_dist`.
- **Online taxonomy**: `get_gbif_taxa_tree`, `get_worms_taxa_tree`,
  `get_nbn_taxa_tree`, `get_reference_db`, `reference_datasets`.
- **Data**: `macro_ex`, `mi_prin`, `oglio`, `traits_mf`, `traits_tachet`,
  `diatom_traits`, `example_tree`, `reference_db(:mi | :mf | :fi | :di)`;
  I/O with `read_community`, `read_traits`, `read_reference`, `read_dist`,
  `write_table`.

Every exported function has a docstring (`?f_rich` in the REPL).
[docs/src/r_to_julia.md](docs/src/r_to_julia.md) maps each R function to its Julia
counterpart and lists the intentional differences.

## Validation

`test/parity.jl` compares the Julia results with those of the R package on
the bundled datasets (diversity, biotic indices with aggregation rules,
exceptions and parent-child corrections, custom reference databases,
presence-absence data, trait assignment, CWM/TSI/CSI, Gower distances,
Euclidean corrections, FRic/FDis/FEve/Rao/FRed, PCoA axis selection, convex
hull volumes against Qhull, hierarchical clustering). The reference values
are produced by `tools/make_reference.R` and stored in `test/reference/`.

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

The bundled data (`data/`, `test/data/`) were exported from the R package
with `Rscript tools/export_data.R ../biomonitoR`.

## Citation

Laini A., Guareschi S., Bolpagni R., Burgazzi G., Bruno D., Gutiérrez-Cánovas
C., Miranda R., Mondy C., Várbíró G., Cancellario T. (2022). biomonitoR: an R
package for managing ecological data and calculating biomonitoring indices.
PeerJ 10:e14183. https://doi.org/10.7717/peerj.14183

## License

GPL (≥ 2), as the R package it is derived from. See [LICENSE.md](LICENSE.md)
and [NOTICE.md](NOTICE.md).
