#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(purrr)
})

opt_list <- list(
  make_option("--inputs",  type="character", help="space-separated RDS files", action="store"),
  make_option("--inst",    type="character", help="output inst_df .tsv.gz"),
  make_option("--result",  type="character", help="output result_df .rds")
)
opt <- parse_args(OptionParser(option_list = opt_list))

rds_paths <- strsplit(opt$inputs, "\\s+")[[1]]
objs <- map(rds_paths, readRDS)

# Loading the function
source("scripts/process_bigbrain_functions.R")

# The function call

joined <- do.call(join_bigbrain, objs)

# write outputs
joined$inst_df %>%
  distinct() %>%
  write_tsv(opt$inst)   # readr auto-detects .gz and bgzips

saveRDS(joined$result_df, file = opt$result)
