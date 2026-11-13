# Validation

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

