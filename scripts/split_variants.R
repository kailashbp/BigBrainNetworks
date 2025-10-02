#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(optparse); library(readr); library(dplyr) })

opt_list <- list(
  make_option("--input",    type="character"),
  make_option("--n-splits", type="integer", default = 1000),
  make_option("--outdir",   type="character"),
  make_option("--prefix",   type="character", default = "variants.part")
)
opt <- parse_args(OptionParser(option_list = opt_list))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# If your BED has NO header (typical): col_names = FALSE
bed <- read_tsv(opt$input, col_names = FALSE, show_col_types = FALSE)
# Expect X1=chr, X2=start(0-based), X3=end(1-based), X4=variant_id
if (ncol(bed) < 3) stop("BED must have at least 3 columns (chr, start, end)")

# Round-robin split by variant (row)
bed$.__bucket__ <- (seq_len(nrow(bed)) - 1L) %% opt$`n-splits` + 1L

for (i in seq_len(opt$`n-splits`)) {
  outp <- file.path(opt$outdir, sprintf("%s_%05d.bed", opt$prefix, i))
  bed %>%
    filter(.__bucket__ == i) %>%
    select(-.__bucket__) %>%
    write_tsv(outp, col_names = FALSE)  # <-- no header
}
