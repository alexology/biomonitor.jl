# Compute reference results with the R package biomonitoR, used by
# test/parity.jl to check that biomonitor.jl gives the same numbers.
#
# Usage (from the biomonitor.jl directory, with biomonitoR, ade4, FD and
# vegan installed):
#   Rscript tools/make_reference.R
#
# Output: test/reference/<case>.tsv. Vectors are written as name/value
# pairs, tables with their row names in the first column.

suppressPackageStartupMessages({
  library(biomonitoR)
})
out <- file.path("test", "reference")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
td <- function(f) system.file("testdata", f, package = "biomonitoR")

num <- function(v) ifelse(is.na(v), "NaN", sprintf("%.17g", v))

w <- function(name, x) {
  path <- file.path(out, paste0(name, ".tsv"))
  if (is.null(dim(x))) {
    nm <- names(x)
    if (is.null(nm)) nm <- seq_along(x)
    df <- data.frame(name = nm, value = if (is.numeric(x)) num(x) else as.character(x), stringsAsFactors = FALSE)
  } else {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    df <- data.frame(row = rownames(x), stringsAsFactors = FALSE)
    for (j in seq_along(x)) {
      v <- x[[j]]
      df[[names(x)[j]]] <- if (is.numeric(v)) num(v) else as.character(v)
    }
  }
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
}
q <- function(expr) suppressMessages(suppressWarnings(expr))

# --- macro_ex --------------------------------------------------------------
data(macro_ex)
mex <- aggregate_taxa(as_biomonitor(macro_ex))
w("macro_ex_allrich", q(allrich(mex)))
w("macro_ex_allindices", q(allindices(mex)))
w("macro_ex_allindices_family", q(allindices(mex, tax_lev = "Family")))
w("macro_ex_allindices_base2", q(allindices(mex, base = 2)))
w("macro_ex_general_info", q(general_info(mex)))
w("macro_ex_abundance_genus", q(abundance(mex, tax_lev = "Genus")))
for (m in c("ita", "spa", "a", "uk")) {
  w(paste0("macro_ex_bmwp_", m), q(bmwp(mex, method = m)))
  w(paste0("macro_ex_aspt_", m), q(aspt(mex, method = m)))
}
w("macro_ex_bmwp_ita_agg", q(bmwp(mex, method = "ita", agg = TRUE)))
w("macro_ex_aspt_uk_agg", q(aspt(mex, method = "uk", agg = TRUE)))
w("macro_ex_whpt_aspt", q(whpt(mex)))
w("macro_ex_whpt_po_bmwp", q(whpt(mex, type = "po", metric = "bmwp")))
w("macro_ex_whpt_ntaxa_agg", q(whpt(mex, metric = "ntaxa", agg = TRUE)))
w("macro_ex_life", q(life(mex)))
w("macro_ex_life_2017", q(life(mex, method = "life_2017", agg = TRUE)))
w("macro_ex_psi", q(psi(mex)))
w("macro_ex_epsi", q(epsi(mex)))
w("macro_ex_epsi_agg", q(epsi(mex, agg = TRUE)))
w("macro_ex_dehli", q(dehli(mex)))
w("macro_ex_igold", q(igold(mex)))
w("macro_ex_ept", q(ept(mex)))
w("macro_ex_ept_family", q(ept(mex, tax_lev = "Family")))
w("macro_ex_eptd", q(eptd(mex)))
w("macro_ex_taxa_abundance", q(get_taxa_abundance(mex, taxa = c("Ephemeroptera", "Baetis"), rel = TRUE)))
w("macro_ex_taxa_richness", q(get_taxa_richness(mex, taxa = "Ephemeroptera", tax_lev = "Genus")))
w("macro_ex_convert_to_biotic", q(convert_to_biotic(mex)))
w("macro_ex_combine_taxa", q(combine_taxa(mex, ntaxa = 2, tax_lev = "Family")))
w("macro_ex_remove_taxa", q(remove_taxa(mex, taxa = c("Hydroptilidae", "Acentrella"))))
w("macro_ex_vegan_family", q(convert_to_vegan(mex, tax_lev = "Family")))

# traits
ts <- q(assign_traits(mex))
w("macro_ex_assign_traits", ts[, 1:5])
w("macro_ex_assign_traits_sum", data.frame(s = rowSums(ts[, -(1:5)])))
mt <- q(manage_traits(ts, method = "nearest+-"))
w("macro_ex_manage_traits", mt[, 1:5])
av <- q(average_traits(ts))
rownames(av) <- av$Taxa
w("macro_ex_average_traits", av[, -1])
w("macro_ex_cwm", q(cwm(mex, trait_db = av)))
w("macro_ex_csi", q(csi(mex)))
tsi_res <- q(tsi(mex))
rownames(tsi_res) <- tsi_res$Taxa
w("macro_ex_tsi", tsi_res[, -1])
w("macro_ex_flow_t", q(flow_t(mex)))
bi <- q(bio_index(mex, index = c("trophic_t", "saprobity_t", "temperature_t")))
w("macro_ex_bio_index", do.call(cbind, bi))
ts_fam <- q(assign_traits(mex, tax_lev = "Family"))
w("macro_ex_assign_traits_family", ts_fam[, 1:5])

cb <- c(8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8)
w("macro_ex_f_rich_2", q(f_rich(mex, trait_db = av, type = "F", col_blocks = cb)))
w("macro_ex_f_rich_3", q(f_rich(mex, trait_db = av, type = "F", col_blocks = cb, nbdim = 3)))
w("macro_ex_f_rich_4", q(f_rich(mex, trait_db = av, type = "F", col_blocks = cb, nbdim = 4)))
w("macro_ex_f_disp_3", q(f_disp(mex, trait_db = av, type = "F", col_blocks = cb, nbdim = 3)))
w("macro_ex_f_eve_3", q(f_eve(mex, trait_db = av, type = "F", col_blocks = cb, nbdim = 3)))
w("macro_ex_f_divs", q(f_divs(mex, trait_db = av, type = "F", col_blocks = cb)))
w("macro_ex_f_red", q(f_red(mex, trait_db = av, type = "F", col_blocks = cb)))
w("macro_ex_qfs", q(biomonitoR:::qfs(ade4::dist.ktab(ade4::ktab.list.df(list(ade4::prep.fuzzy(data.frame(av[, -1], row.names = av$Taxa), col.blocks = cb))), type = "F"), nbdim = 10))$meanSD)
avd <- data.frame(av[, -1], row.names = av$Taxa)
gd <- ade4::dist.ktab(ade4::ktab.list.df(list(ade4::prep.fuzzy(avd, col.blocks = cb))), type = "F")
w("macro_ex_gower", as.matrix(gd))
w("macro_ex_select_legendre", q(select_pcoa_axes(gd, method = "legendre")))
w("macro_ex_select_maire", q(select_pcoa_axes(gd, method = "maire")))

# --- mi_prin (non-Euclidean trait distance) ---------------------------------
data(mi_prin)
mip <- aggregate_taxa(as_biomonitor(mi_prin))
w("mi_prin_allindices", q(allindices(mip)))
w("mi_prin_allrich", q(allrich(mip)))
w("mi_prin_bmwp_ita_agg", q(bmwp(mip, method = "ita", agg = TRUE)))
w("mi_prin_aspt_ita", q(aspt(mip, method = "ita")))
w("mi_prin_aspt_spa", q(aspt(mip, method = "spa")))
w("mi_prin_whpt", q(whpt(mip, agg = TRUE)))
w("mi_prin_life", q(life(mip)))
w("mi_prin_psi", q(psi(mip)))
w("mi_prin_epsi", q(epsi(mip)))
w("mi_prin_dehli", q(dehli(mip)))
w("mi_prin_ept", q(ept(mip, tax_lev = "Family")))
w("mi_prin_igold", q(igold(mip)))
w("mi_prin_eptd", q(eptd(mip)))
w("mi_prin_dness_family", q(dness(mip, tax_lev = "Family", complete = TRUE)))
w("mi_prin_dness_family_st", q(dness(mip, tax_lev = "Family", complete = TRUE, method = "delta.st")))
w("mi_prin_dness_family_bin", q(dness(mip, tax_lev = "Family", complete = TRUE, method = "delta.bin")))
avp <- q(average_traits(assign_traits(mip)))
w("mi_prin_cwm", q(cwm(mip, trait_db = avp)))
w("mi_prin_csi", q(csi(mip)))
for (cr in c("cailliez", "lingoes", "sqrt", "quasi")) {
  w(paste0("mi_prin_f_rich_3_", cr), q(f_rich(mip, trait_db = avp, type = "F", col_blocks = cb, nbdim = 3, correction = cr)))
  w(paste0("mi_prin_f_divs_", cr), q(f_divs(mip, trait_db = avp, type = "F", col_blocks = cb, correction = cr)))
  w(paste0("mi_prin_f_red_", cr), q(f_red(mip, trait_db = avp, type = "F", col_blocks = cb, correction = cr)))
}
w("mi_prin_f_rich_3", q(f_rich(mip, trait_db = avp, type = "F", col_blocks = cb, nbdim = 3)))
w("mi_prin_f_disp_3", q(f_disp(mip, trait_db = avp, type = "F", col_blocks = cb, nbdim = 3)))
w("mi_prin_f_eve_3", q(f_eve(mip, trait_db = avp, type = "F", col_blocks = cb, nbdim = 3)))
w("mi_prin_f_rich_auto", q(f_rich(mip, trait_db = avp, type = "F", col_blocks = cb, nbdim = "auto", correction = "cailliez")))
w("mi_prin_f_rich_zero", q(f_rich(mip, trait_db = avp, type = "F", col_blocks = cb, nbdim = 3, zerodist_rm = TRUE, correction = "lingoes")))
avpd <- data.frame(avp[, -1], row.names = avp$Taxa)
gdp <- ade4::dist.ktab(ade4::ktab.list.df(list(ade4::prep.fuzzy(avpd, col.blocks = cb))), type = "F")
w("mi_prin_gower", as.matrix(gdp))
w("mi_prin_cailliez", as.matrix(q(ade4::cailliez(gdp))))
w("mi_prin_lingoes", as.matrix(q(ade4::lingoes(gdp))))
w("mi_prin_quasi", as.matrix(q(ade4::quasieuclid(gdp))))
w("mi_prin_select_legendre", q(select_pcoa_axes(gdp, method = "legendre")))
# select_pcoa_axes(method = "maire") fails in R for mi_prin (sqrt transformation)

# --- presence-absence -------------------------------------------------------
mexb <- aggregate_taxa(q(as_biomonitor(macro_ex, FUN = bin)))
w("macro_ex_bin_allindices", q(allindices(mexb)))
w("macro_ex_bin_bmwp_ita", q(bmwp(mexb, method = "ita")))
w("macro_ex_bin_aspt_ita", q(aspt(mexb, method = "ita")))
w("macro_ex_bin_whpt_po", q(whpt(mexb, type = "po")))
w("macro_ex_bin_life", q(life(mexb)))
w("macro_ex_bin_igold", q(igold(mexb)))

# --- UK data with custom reference -----------------------------------------
load(td("ref_def.rda")); load(td("taxa_uk_data.rda")); load(td("to_change.rda"))
uk <- aggregate_taxa(q(as_biomonitor(taxa_uk_data, dfref = ref_def, to_change = to_change)))
w("uk_life_2017", q(life(uk, method = "life_2017")))
w("uk_whpt_agg", q(whpt(uk, method = "uk", agg = TRUE)))
w("uk_whpt_bmwp_agg", q(whpt(uk, method = "uk", metric = "bmwp", agg = TRUE)))
w("uk_psi", q(psi(uk)))
w("uk_bmwp_uk_agg", q(bmwp(uk, method = "uk", agg = TRUE)))
w("uk_epsi", q(epsi(uk)))

load(td("aspt_bmwp_spain.rda"))
spa <- aggregate_taxa(as_biomonitor(aspt_bmwp_spain))
w("spain_bmwp_exc", q(bmwp(spa, method = "spa", exceptions = "Ferrissia")))
w("spain_aspt_spa", q(aspt(spa, method = "spa")))

load(td("data_uk_bioco.rda"))
alien <- c("Potamopyrgus antipodarum", "Crangonyx pseudogracilis", "Crangonyx", "Dikerogammarus",
           "Dikerogammarus haemobaphes", "Chelicorophium curvispinum", "Corbicula fluminea", "Hypania invalida",
           "Dreissena polymorpha", "Dugesia tigrina")
bio <- aggregate_taxa(q(as_biomonitor(data_uk_bioco, dfref = ref_def)))
w("bioco", q(bioco(bio, alien = alien, dfref = ref_def, digits = 3)))

# --- macrophytes ------------------------------------------------------------
load(td("ibmr_taxa.rda"))
ib <- aggregate_taxa(q(as_biomonitor(ibmr_taxa, group = "mf")))
to_remove <- c("Spirogyra dubia", "Spirogyra gracilis", "Spirogyra majuscola", "Spirogyra parva")
w("ibmr_coef_exc", q(ibmr(ib, coverage_coef = TRUE, exceptions = to_remove)))
w("ibmr_coef", q(ibmr(ib, coverage_coef = TRUE)))
data(oglio)
og <- aggregate_taxa(q(as_biomonitor(oglio, group = "mf", FUN = bin)))
w("oglio_dness", q(dness(og, method = "delta", complete = TRUE)))
w("oglio_dness_st", q(dness(og, method = "delta.st", complete = TRUE)))
w("oglio_allrich", q(allrich(og)))
data(traits_mf)
w("oglio_assign_traits", q(assign_traits(og, trait_db = traits_mf, group = "mf"))[, 1:5])

# --- ambiguous taxa ---------------------------------------------------------
load(td("ambiguous_df.rda"))
amb <- aggregate_taxa(q(as_biomonitor(ambiguous_df, FUN = sum)))
w("ambiguous_mcwp", q(solve_ambiguous(amb, method = "MCWP")))
w("ambiguous_rpkc", q(solve_ambiguous(amb, method = "RPKC")))

# --- reference from tree ------------------------------------------------------
data(Tree)
w("ref_from_tree", ref_from_tree(Tree))

# --- distances and hull checks ------------------------------------------------
set.seed(1)
pts <- matrix(rnorm(40 * 4), 40, 4)
w("hull_points", pts)
w("hull_volumes", c(d2 = geometry::convhulln(pts[, 1:2], "FA")$vol,
                    d3 = geometry::convhulln(pts[, 1:3], "FA")$vol,
                    d4 = geometry::convhulln(pts, "FA")$vol))
hc <- hclust(dist(pts), method = "ward.D2")
w("hclust_cutree", sapply(2:6, function(k) cutree(hc, k)))
hc2 <- hclust(dist(pts), method = "average")
w("hclust_cutree_average", sapply(2:6, function(k) cutree(hc2, k)))
# as_biomonitor(dfref = ...) writes a hunspell dictionary in the working directory
unlink(c("custom_dictio.dic", "custom_dictio.aff"))
cat("reference results written to", out, "\n")
