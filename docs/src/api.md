# API

```@meta
CurrentModule = Biomonitor
```

```@docs
Biomonitor
```

## Types and accessors

```@docs
LabeledVector
LabeledMatrix
labels
rownames
colnames
CommunityTable
ReferenceDB
BiomonitorData
AggregatedData
TraitTable
AssignedTraits
DistMatrix
taxa
samples
abundances
ranks
TAXONOMIC_LEVELS
RANKS
```

## Data and I/O

```@docs
read_community
write_table
read_traits
read_reference
read_dist
macro_ex
mi_prin
oglio
traits_mf
traits_tachet
diatom_traits
example_tree
reference_db
get_reference_db
reference_datasets
```

## Import and data management

```@docs
as_biomonitor
aggregate_taxa
bin
to_bin
quick_rename
suggest_names
excluded_taxa
corrected_names
ref_from_tree
check_tree
subset_taxa
remove_taxa
combine_taxa
convert_to_vegan
convert_to_biotic
solve_ambiguous
abundance_classes
general_info
get_taxa_abundance
get_taxa_richness
```

## Abundance, richness and diversity

```@docs
abundance
richness
allrich
allindices
shannon
simpson
invsimpson
esimpson
berpar
invberpar
brill
fisher
margalef
mcintosh
menhinick
pielou
dness
```

## Biotic indices

```@docs
bmwp
aspt
whpt
life
psi
epsi
dehli
ibmr
igold
ept
eptd
bioco
show_scores
```

## Traits

```@docs
assign_traits
manage_traits
average_traits
sample_traits
add_bias_to_traits
prep_fuzzy
fuzzy_trait_ratio
flow_t
bio_index
cwm
tsi
csi
```

## Functional diversity and distances

```@docs
f_rich
f_disp
f_eve
f_divs
f_red
select_pcoa_axes
qfs
gower_fuzzy
euclidean_dist
is_euclid
cailliez
lingoes
quasieuclid
pcoa
divc
divcmax
mantel_test
bray_curtis
jaccard_dist
convex_hull_volume
```

## Clustering and indicator taxa

```@docs
hclust
cutree
indicator_taxa
indval
```

## Online taxonomy

```@docs
get_gbif_taxa_tree
get_worms_taxa_tree
get_nbn_taxa_tree
TaxonomyTable
```
