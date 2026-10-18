# Notice

biomonitor.jl is a Julia translation of the R package **biomonitoR**
(version 0.9.5, GPL ≥ 2) by Alex Laini, Tommaso Cancellario, Simone
Guareschi, Daniel Bruno, Cayetano Gutiérrez-Cánovas and Cédric Mondy, with
contributions by Rossano Bolpagni, Gemma Burgazzi, Rafael Miranda and Gábor
Várbíró.

The files in `data/` and `test/data/` are the datasets distributed with
biomonitoR (reference taxonomies, index scores, traits and examples),
converted to TSV by `tools/export_data.R`. Their sources are cited in the
biomonitoR documentation; the reference taxonomies and the macroinvertebrate
traits rely on www.freshwaterecology.info (Schmidt-Kloiber & Hering 2015)
and Tachet et al. (2010).

Parts of the code translate algorithms from other GPL-licensed R packages
used by biomonitoR:

- **ade4** (Dray, Dufour, Chessel et al.): `prep.fuzzy`, `dist.ktab`
  (fuzzy variables), `dist.prop`, `bicenter.wt`, `is.euclid`, `cailliez`,
  `lingoes`, `quasieuclid`, `dudi.pco`, `divc`, `divcmax`,
  `mantel.randtest` (`src/distances.jl`);
- **ape** (Paradis et al.): `mst` (`src/distances.jl`);
- **FD** (Laliberté, Legendre, Shipley): functional dispersion and evenness,
  through the biomonitoR helpers `fdisp_k` and `feve_k`
  (`src/functional.jl`);
- the quality of functional spaces of Maire et al. (2015), adapted by
  C. Gutiérrez-Cánovas (`qfs`);
- **stats** (R Core Team): Brent minimisation of `optimize` and the
  hierarchical clustering algorithm of `hclust` (F. Murtagh)
  (`src/utils.jl`, `src/cluster.jl`).

The convex hull (`src/hull.jl`), the spelling suggestions and the JSON
parser are original code written for this package.

The test reference values in `test/reference/` were computed with the R
packages biomonitoR 0.9.5, ade4, FD and geometry by `tools/make_reference.R`.
