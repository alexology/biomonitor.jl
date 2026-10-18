# Export the datasets shipped with the R package biomonitoR to plain TSV files
# used by biomonitor.jl.
#
# Usage (from the biomonitor.jl directory):
#   Rscript tools/export_data.R ../biomonitoR
#
# Output:
#   data/<name>.tsv       internal reference data (R/sysdata.rda) and the
#                         example datasets (data/*.rda)
#   test/data/<name>.tsv  test fixtures (inst/testdata/*.rda)
#
# Numbers are written with 17 significant digits, missing values as empty
# cells, text as UTF-8. Factors are written as their labels.

args <- commandArgs(trailingOnly = TRUE)
src <- if (length(args) >= 1) args[1] else "../biomonitoR"

write_tsv <- function(x, path) {
  if (is.matrix(x)) {
    x <- data.frame(Row = rownames(x), x, check.names = FALSE, stringsAsFactors = FALSE)
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  for (j in seq_along(x)) {
    v <- x[[j]]
    if (is.factor(v)) v <- as.character(v)
    if (is.numeric(v)) {
      out <- ifelse(is.na(v), "", sprintf("%.17g", v))
      # keep integers readable
      int <- !is.na(v) & v == round(v) & abs(v) < 1e15
      out[int] <- sprintf("%.0f", v[int])
      v <- out
    } else {
      v <- as.character(v)
      v[is.na(v)] <- ""
      # some citations contain Mac Roman bytes (e.g. "\x87" for "á")
      bad <- !validUTF8(v)
      if (any(bad)) v[bad] <- iconv(v[bad], from = "macintosh", to = "UTF-8")
    }
    if (any(grepl("[\t\n\r]", v))) stop("tab or newline inside a field of ", path)
    x[[j]] <- v
  }
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con))
  writeLines(paste(names(x), collapse = "\t"), con)
  if (nrow(x) > 0) {
    writeLines(do.call(paste, c(unname(as.list(x)), sep = "\t")), con)
  }
}

export_rda <- function(file, outdir) {
  e <- new.env()
  load(file, envir = e)
  for (n in ls(e)) {
    write_tsv(e[[n]], file.path(outdir, paste0(n, ".tsv")))
    cat("wrote", file.path(outdir, paste0(n, ".tsv")), "\n")
  }
}

dir.create("data", showWarnings = FALSE)
dir.create(file.path("test", "data"), recursive = TRUE, showWarnings = FALSE)

export_rda(file.path(src, "R", "sysdata.rda"), "data")
for (f in list.files(file.path(src, "data"), pattern = "\\.rda$", full.names = TRUE)) {
  export_rda(f, "data")
}
for (f in list.files(file.path(src, "inst", "testdata"), pattern = "\\.rda$", full.names = TRUE)) {
  export_rda(f, file.path("test", "data"))
}
