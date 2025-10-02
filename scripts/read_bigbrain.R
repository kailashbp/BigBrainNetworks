#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr); library(doParallel)
})

# === adjust defaults to your columns ===
opt_list <- list(
  make_option("--sumstats", type="character", help="bgzip+tabix subset per split"),
  make_option("--out",      type="character", help="output RDS"),
  make_option("--beta_col", default="fixed_beta"),
  make_option("--p_col",    default="Random_P"),
  make_option("--id_col",   default="variant_id"),
  make_option("--se_col",   default="fixed_sd"),
  make_option("--p_thresh", type="double", default=5e-8),
  make_option("--chr_col",  default="chr"),
  make_option("--pos_col",  default="pos"),
  make_option("--ref_col",  default="ref"),
  make_option("--alt_col",  default="alt"),
  make_option("--feature_col", default="#feature"),
  make_option("--delim", default="\t"),
  make_option("--filter_col", default="crossmap"),
  make_option("--filter_min", type="double", default=0),
  make_option("--filter_max", type="double", default=0.1),
  make_option("--cores", type="integer", default=1)
)
opt <- parse_args(OptionParser(option_list = opt_list))

# --- setup parallel backend ---
# --- setup parallel backend (use CLI cores) ---
if (opt$cores > 1) {
  doParallel::registerDoParallel(cores = opt$cores)
  } else {
    foreach::registerDoSEQ()
  }

# Loading the function
source("scripts/process_bigbrain_functions.R")

# The function call

res <- read_bigbrain(
  filename    = opt$sumstats,
  beta_col    = opt$beta_col,
  p_col       = opt$p_col,
  id_col      = opt$id_col,
  se_col      = opt$se_col,
  p_thresh    = opt$p_thresh,
  chr_col     = opt$chr_col,
  pos_col     = opt$pos_col,
  ref_col     = opt$ref_col,
  alt_col     = opt$alt_col,
  feature_col = opt$feature_col,
  delim       = opt$delim,
  filter_col  = opt$filter_col,
  filter_min  = opt$filter_min,
  filter_max  = opt$filter_max
)

saveRDS(res, file = opt$out)
doParallel::stopImplicitCluster()
