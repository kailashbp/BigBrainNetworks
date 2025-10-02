#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(purrr)
  library(MendelianRandomization)
  library(mr.raps)
  library(MRPRESSO)
  library(MRMix)
  library(inspre)
})

# ---- CLI ----
opt_list <- list(
  make_option("--input_rds",  type="character", help="join_bigbrain RDS"),
  make_option("--method",     type="character", default="ivw",
              help="MR method: ivw, ps, aps, raps, egger_p, median, mbe, mr_presso, mr_mix"),
  make_option("--min_inst",   type="integer", default=3,
              help="minimum instruments per row"),
  make_option("--out_tsv_gz", type="character", help="output TSV.GZ")
)
opt <- parse_args(OptionParser(option_list = opt_list))

# ---- load helpers ----

obj <- readRDS(opt$input_rds)
res_df <- obj$result_df

# keep rows with enough instruments
res_df <- res_df %>%
  mutate(n_inst = map_int(inst, length)) %>%
  filter(n_inst >= opt$min_inst)

# run MR per row
mr_tbl <- fit_R_MR(res_df, mr_method = opt$method) %>%
  mutate(method = tolower(opt$method)) %>%
  relocate(method, .before = 1)

# write gzipped TSV
con <- gzfile(opt$out_tsv_gz, open = "wb")
readr::write_tsv(mr_tbl, con)
close(con)
