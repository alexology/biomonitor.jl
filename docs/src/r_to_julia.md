# From biomonitoR (R) to Biomonitor.jl

This page maps the exported functions of biomonitoR 0.9.5 to Biomonitor.jl
and lists where the Julia version behaves differently on purpose.

## General conventions

- R `data.frame`s become Julia types: `CommunityTable`, `ReferenceDB`,
  `BiomonitorData` (class `asb`), `AggregatedData` (class `biomonitoR`),
  `TraitTable`, `AssignedTraits`, `DistMatrix`.
- Per-sample results are `LabeledVector`s; tabular results are
  `LabeledMatrix`es; `traceB = TRUE` lists are `NamedTuple`s.
- Arguments keep the R names and defaults (`tax_lev`, `traceB`,
  `col_blocks`, ...), passed as keywords: `richness(agg; tax_lev = :Family)`.
- Custom score tables (R `data.frame`s passed as `method`, `agg`,
  `fs_scores`, `fssr_scores`) are `Dict`s; `show_scores` returns the
  built-in ones in that form.
- `NA` results of R are `NaN` (numbers) in Julia.
- R messages and warnings are Julia `@info` and `@warn` log records.
- Names are sorted by Unicode code point; R sorts with the locale
  collation, so the row order of some tables can differ (values do not).

## Function map

| R | Julia | Notes |
|---|---|---|
| `as_biomonitor(x, group, dfref, to_change, FUN, correct_names, traceB)` | `as_biomonitor(x; group, dfref, to_change, FUN, correct_names)` | `to_change`: `Dict`/`Pair`s; `correct_names`: `false`, `true` (terminal) or a `Dict` wrong => correct. Excluded and corrected names are always stored: `excluded_taxa(asb)`, `corrected_names(asb)` |
| `aggregate_taxa` | `aggregate_taxa` | |
| `as.data.frame.asb`, `print.asb` | fields of `BiomonitorData`, `show` | `write_table` writes it to a file |
| `subset.asb` | `subset_taxa(x; taxa, exclude)` | |
| `quick_rename(x, group, write_table)` | `quick_rename(x; group, corrections, interactive)` | write the result with `write_table` |
| `ref_from_tree` | `ref_from_tree(levels, values; group)` | `ref_from_tree(example_tree()...)` |
| `get_reference_db` | `get_reference_db`, `reference_datasets` | |
| `get_gbif_taxa_tree`, `get_worms_taxa_tree`, `get_nbn_taxa_tree` | same names | rows in a `TaxonomyTable` (a vector of `NamedTuple`s shown as a table); GBIF also returns `fuzzy_matches`; NBN multiple matches are interactive only with `interactive = true` |
| `get_iucn_taxa_tree` | — | not ported: it targets the IUCN Red List API v3; IUCN now serves API v4 with different endpoints |
| `abundance`, `richness`, `allrich`, `general_info` | same | |
| `get_taxa_abundance(x, taxa, rel)`, `get_taxa_richness(x, taxa, tax_lev)` | `get_taxa_abundance(x, taxa; rel)`, `get_taxa_richness(x, taxa, tax_lev)` | |
| `remove_taxa`, `combine_taxa`, `convert_to_vegan`, `convert_to_biotic`, `solve_ambiguous`, `abundance_classes`, `to_bin`, `bin` | same | `remove_taxa` returns a `CommunityTable` |
| `shannon`, `simpson`, `invsimpson`, `esimpson`, `berpar`, `invberpar`, `brill`, `fisher`, `margalef`, `mcintosh`, `menhinick`, `pielou`, `allindices` | same | Fisher's alpha uses a translation of R's `optimize` (Brent) |
| `bmwp`, `aspt`, `whpt`, `life`, `psi`, `epsi`, `dehli`, `ibmr`, `igold`, `ept`, `eptd`, `bioco`, `show_scores` | same | see differences below |
| `dness` | `dness` | |
| `assign_traits`, `manage_traits`, `average_traits`, `sample_traits`, `add_bias_to_traits` | same | `AssignedTraits` fields: `taxa`, `taxa_level`, `trait_taxa`, `trait_level`, `distance` (R: `Taxa`, `Taxa_db`, `Traits_taxlev`, `Traits_real`, `Taxonomic_distance`) |
| `cwm`, `tsi`, `csi`, `fuzzy_trait_ratio`, `flow_t`, `bio_index` | same | `bio_index` returns a `NamedTuple` index => values |
| `f_rich`, `f_disp`, `f_eve`, `f_divs`, `f_red`, `select_pcoa_axes` | same | `traitSel` (interactive trait selection) is not ported: subset the `TraitTable` instead |
| `plot.asb`, `plot.biomonitoR` | — | plotting with plotly is not ported; the data are in `BiomonitorData` and `AggregatedData` |
| `plot_indicator_taxa` | `indicator_taxa` | returns clusters, silhouette and indicator values instead of a Sankey plot |
| `macro_ex`, `mi_prin`, `oglio`, `traits_mf`, `Tree` | `macro_ex()`, `mi_prin()`, `oglio()`, `traits_mf()`, `example_tree()` | |
| internal `traitsTachet`, `diatom_traits`, `mi_ref`, `mf_ref`, `fi_ref`, `di_ref` | `traits_tachet()`, `diatom_traits()`, `reference_db(:mi)`, ... | |
| ade4 helpers used internally | `prep_fuzzy`, `gower_fuzzy`, `is_euclid`, `cailliez`, `lingoes`, `quasieuclid`, `pcoa`, `divc`, `divcmax`, `mantel_test` | exported for direct use |
| defunct names (`asBiomonitor`, `aggregatoR`, ...) | — | not ported |

## Intentional differences

1. **Spelling suggestions.** R uses hunspell with a dictionary of the
   reference taxa; Julia ranks the reference taxa by edit distance
   (`suggest_names`). Suggestions may differ in number and order. No
   dictionary files are written in the working directory.
2. **PSI.** When a sample has scored taxa but no sensitive ones, R returns
   `NA` (when the other samples do have sensitive taxa); Julia returns 0, the
   value of the index. Samples without scored taxa give `NaN`.
3. **Convex hulls.** FRic uses a Julia incremental convex hull instead of
   Qhull. Volumes agree with Qhull to floating-point precision. When the
   points do not span the space, the volume is 0 (Qhull stops with an
   error); the `prec` element of `set_param` is ignored.
4. **`manage_traits`, methods `"nearest+-"` and `"nearest-+"`.** The
   documented priority (finer, respectively coarser, taxa first) is applied.
   The R code takes the minimum over all distances, which differs only for
   taxa assigned traits from both finer and coarser taxa (this does not
   happen with `assign_traits` output from the built-in data).
5. **Indicator taxa.** R computes Bray-Curtis dissimilarities dropping the
   first taxon column (a side effect of `lapply` on a matrix, which removes
   the `unassigned` row when present); Julia excludes `unassigned` and keeps
   every other taxon.
6. **`cwm`, `tsi`, `csi`** keep the trait order of the trait database
   (R sorts the columns alphabetically).
7. **`select_pcoa_axes(method = "cor")`**: p-values come from random
   permutations, as in R, so they differ between runs and from R.
8. **Duplicated taxa in a reference database**: the first occurrence is
   used (R duplicates the rows). The built-in `mi_ref` lists *Telmaturgus*
   twice, `fi_ref` *Squaliobarbinae*.
9. **`cwm` with an invalid `tax_lev`** throws an `ArgumentError` (R returns
   the message as a string).
10. **WHPT with `traceB`** returns the exceptions element correctly (the R
    code fails when exceptions are used together with `traceB = TRUE`).
11. **`lingoes` with `cor_zero = true`** corrects the distances larger than
    `tol` and leaves the others unchanged, as `cailliez`. ade4 compares the
    squared distance with `tol`, so distances up to about 3e-4 are skipped
    and returned squared (5e-5 becomes 2.5e-9).
12. **`is_euclid` on a matrix of zero distances** returns `true` (ade4
    returns `NA`, from 0/0); `divcmax` then returns 0 with uniform weights.
13. **`cailliez` and `lingoes` use `tol` also to decide whether a correction
    is needed** (ade4 checks `is.euclid` with its default tolerance, so a
    smaller `tol` could return an uncorrected non-Euclidean matrix).
14. **Zero-distance removal (`zerodist_rm = true`) with missing traits.** A
    zero distance is then not transitive (A–B and B–C zero, A–C positive):
    a taxon is merged into the first group whose members are all at zero
    distance from it. R merges along each zero pair and stops with "Something
    went wrong with the zero distance removal"; when zero distances are
    transitive the result is the same.
15. **Silhouette in `indicator_taxa`.** A sample alone in its cluster has
    silhouette width 0 (Rousseeuw 1987; R `cluster::silhouette`). The R
    helper `silhouette_biomonitor` assigns it `NaN` and drops it, so the
    average silhouette is inflated and the number of clusters can change
    (an example with five samples gives k = 4 in R, k = 3 in Julia).
    Singleton clusters also appear in `silhouette_by_cluster`.
16. **Web services (`get_gbif_taxa_tree`, `get_worms_taxa_tree`,
    `get_nbn_taxa_tree`).** R catches every error of the WoRMS and NBN
    requests (`tryCatch`), so a lost connection reports all the taxa as not
    found; Julia reports as not found only empty responses and HTTP 4xx
    answers, and rethrows network and server errors, timeouts (408) and rate
    limiting (429). Responses with an
    unexpected shape raise an `ArgumentError` that names the service and the
    taxon. A WoRMS record that is not accepted and has no `valid_AphiaID` is
    reported as not found with a warning (R stops with an error). Names are
    fully percent-encoded in the URLs (R encodes only spaces).
17. **`bio_index` with `index = "custom"`.** `col_blocks` must be positive
    and add up to the length of `col_order`; otherwise an `ArgumentError`
    is thrown before any calculation.
18. **GBIF matches.** A name that GBIF matches only at a higher rank
    (`matchType = "HIGHERRANK"`, e.g. the kingdom *Animalia* for an unknown
    name) is reported as not found, with a warning; R returns the tree of that
    higher taxon under the name. Names corrected for spelling
    (`matchType = "FUZZY"`) are kept and listed in `fuzzy_matches`.
    Synonyms are recognised by `acceptedUsageKey`: GBIF no longer returns the
    `synonym` flag that R checks (verified on 2026-11-13), so R keeps the
    synonym as the taxon name and reports no synonyms.
