#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(purrr); library(doParallel); library(foreach)
})

opt_list <- list(
  make_option("--inputs",  type="character", help="Text file containing RDS paths (one per line)"),
  make_option("--inst",    type="character", help="output inst_df .tsv.gz"),
  make_option("--result",  type="character", help="output result_df .rds"),
  make_option("--group_size", type="integer",   help="Number of files per worker per round", default = 4),
  make_option("--cores",      type="integer",   help="Number of parallel workers", default = parallel::detectCores())
)
opt <- parse_args(OptionParser(option_list = opt_list))

rds_files <- readLines(opt$inputs)
rds_files <- rds_files[nzchar(rds_files)] # Remove empty lines
stopifnot(length(opt$inputs) > 0)

message("Merging ", length(rds_files), " input RDS files.")
message("group_size = ", opt$group_size, ", cores = ", opt$cores)

# Loading the function
source("scripts/process_bigbrain_functions.R")

# The function call

## Parallel tree reduction
current <- as.list(rds_files)

# spin up cluster
cores <- max(1L, opt$cores)
cl <- parallel::makeCluster(cores)
doParallel::registerDoParallel(cl)
on.exit(parallel::stopCluster(cl), add = TRUE)

while (length(current) > 1) {
  # make groups of size group_size (last group may be smaller)
  groups <- split(current, ceiling(seq_along(current) / opt$group_size))

  # each group is a character vector of paths in the first round,
  # and a list of *either* paths or already-merged objects in later rounds.
  # Normalize each group into a list of objects for join_bigbrain.
  results <- foreach(g = groups, .packages = c("tibble","dplyr","purrr")) %dopar% {
    # if group elements are already merged objects, pass them through;
    # otherwise readRDS paths. Then merge all with join_bigbrain.
    objs <- lapply(g, function(x) {
      if (is.list(x) && !is.null(x$inst_df) && !is.null(x$result_df)) {
        x
      } else if (is.character(x)) {
        readRDS(x)
      } else {
        x
      }
    })
    do.call(join_bigbrain, objs)
  }

  current <- results
  message("Rounds remaining: ~", ceiling(log2(length(current))), " (chunks now: ", length(current), ")")
}

joined <- current[[1L]]

# write outputs
joined$inst_df %>%
  distinct() %>%
  write_tsv(opt$inst) 

saveRDS(joined$result_df, file = opt$result)
message("Done. Wrote:\n  inst_df -> ", opt$inst, "\n  result_df -> ", opt$result)
