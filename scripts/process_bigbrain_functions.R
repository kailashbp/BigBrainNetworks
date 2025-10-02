# The function

read_bigbrain <- function(filename, beta_col, p_col, id_col = "variant_id",
                          se_col = NULL, p_thresh = 5e-08, chr_col = "chr",
                          pos_col = "pos", ref_col = "ref", alt_col = "alt",
                          feature_col = "#feature", delim = "\t",
                          filter_col = NULL, filter_max = NULL,
                          filter_min = NULL){

  if (is.null(se_col)) stop("se_col must be provided (used to enforce se > 0).")

  # ---- helpers -------------------------------------------------------------
  as_chr_or_na <- function(x) if (is.na(x)) NA_character_ else as.character(x)
  as_num_or_na <- function(x) if (is.na(x)) NA_real_      else as.numeric(x)

  # initialize
  last_snp  <- NULL
  effect_df <- tibble::tibble()
  inst_df   <- tibble::tibble(inst=character(), chr=character(), pos=numeric(),
                              ref=character(), alt=character())
  result_df <- tibble::tibble(exp_out=character(), Exposure=character(), Outcome=character(),
                              inst=list(), beta_exp=list(), se_exp=list(),
                              beta_out=list(), se_out=list())

  # Build exp-outcome rows for a single SNP's effects
  build_exp_data <- function(effect_df, inst_exposures) {
    if (length(inst_exposures) == 0) return(NULL)

    use_parallel <- requireNamespace("foreach", quietly = TRUE) &&
      foreach::getDoParWorkers() > 1

    cross_it <- function(exp) {
      dplyr::select(effect_df, -p) %>%
        dplyr::filter(.data[[feature_col]] == exp) %>%
        dplyr::rename(Exposure = !!feature_col, beta_exp = beta, se_exp = se) %>%
        dplyr::cross_join(
          dplyr::select(effect_df, -inst, -p) %>%
            dplyr::rename(Outcome = !!feature_col, beta_out = beta, se_out = se)
        ) %>%
        dplyr::filter(Exposure != Outcome) %>%                 # ← drop exposure=outcome
        dplyr::mutate(
          exp_out  = paste(Exposure, Outcome, sep = "_"),
          inst     = list(as.character(inst)),
          beta_exp = list(as.numeric(beta_exp)),
          se_exp   = list(as.numeric(se_exp)),
          beta_out = list(as.numeric(beta_out)),
          se_out   = list(as.numeric(se_out))
        )
    }

    if (use_parallel) {
      foreach::foreach(exp = inst_exposures, .combine = dplyr::bind_rows) %dopar% cross_it(exp)
    } else {
      if (!requireNamespace("purrr", quietly = TRUE)) stop("Needs purrr or foreach")
      purrr::map_dfr(inst_exposures, cross_it)
    }
  }

  # per-row handler
  f <- function(x, pos){
    if (is.null(last_snp)) {
      message("Initializing")
      inst_df  <<- dplyr::bind_rows(inst_df, tibble::tibble(
        inst = as_chr_or_na(x[[id_col]]),
        chr  = as_chr_or_na(x[[chr_col]]),
        pos  = as_num_or_na(x[[pos_col]]),
        ref  = as_chr_or_na(x[[ref_col]]),
        alt  = as_chr_or_na(x[[alt_col]])
      ))
      effect_df <<- tibble::tibble(
        inst    = as_chr_or_na(x[[id_col]]),
        feature = as_chr_or_na(x[[feature_col]]),
        beta    = as_num_or_na(x[[beta_col]]),
        se      = as_num_or_na(x[[se_col]]),
        p       = as_num_or_na(x[[p_col]])
      )
      last_snp <<- x[[id_col]]

    } else if (!identical(x[[id_col]], last_snp)) {
      cat(sprintf("\nProcessing SNP %s.\n", as_chr_or_na(last_snp)))
      effect_df <- dplyr::distinct(effect_df)
      inst_exposures <- dplyr::filter(effect_df, p < p_thresh)[[feature_col]]
      message("It is an instrument for ", length(inst_exposures), " exposures.")
      exp_data <- build_exp_data(effect_df, inst_exposures)

      if (!is.null(exp_data) && nrow(exp_data) > 0) {
        match_idx <- match(exp_data$exp_out, result_df$exp_out)
        has_match <- which(!is.na(match_idx))
        if (length(has_match) > 0) {
          m <- match_idx[has_match]; i <- has_match
          result_df$inst[m]     <<- purrr::map2(result_df$inst[m],     exp_data$inst[i],     ~c(.x, .y))
          result_df$beta_exp[m] <<- purrr::map2(result_df$beta_exp[m], exp_data$beta_exp[i], ~c(.x, .y))
          result_df$se_exp[m]   <<- purrr::map2(result_df$se_exp[m],   exp_data$se_exp[i],   ~c(.x, .y))
          result_df$beta_out[m] <<- purrr::map2(result_df$beta_out[m], exp_data$beta_out[i], ~c(.x, .y))
          result_df$se_out[m]   <<- purrr::map2(result_df$se_out[m],   exp_data$se_out[i],   ~c(.x, .y))
        }
        if (any(is.na(match_idx))) {
          result_df <<- dplyr::bind_rows(result_df, exp_data[is.na(match_idx), ])
        }
      }

      # reset for new SNP
      effect_df <<- tibble::tibble(
        inst    = as_chr_or_na(x[[id_col]]),
        feature = as_chr_or_na(x[[feature_col]]),
        beta    = as_num_or_na(x[[beta_col]]),
        se      = as_num_or_na(x[[se_col]]),
        p       = as_num_or_na(x[[p_col]])
      )
      inst_df   <<- dplyr::bind_rows(inst_df, tibble::tibble(
        inst = as_chr_or_na(x[[id_col]]),
        chr  = as_chr_or_na(x[[chr_col]]),
        pos  = as_num_or_na(x[[pos_col]]),
        ref  = as_chr_or_na(x[[ref_col]]),
        alt  = as_chr_or_na(x[[alt_col]])
      ))
      last_snp  <<- x[[id_col]]

    } else {
      # same SNP; optional filter
      skip <- FALSE
      if (!is.null(filter_col) && !is.na(x[[filter_col]])) {
        if (is.null(filter_min) && is.null(filter_max)) stop("Filter col specified but no values given.")
        if (!is.null(filter_min) && x[[filter_col]] < filter_min) skip <- TRUE
        if (!is.null(filter_max) && x[[filter_col]] > filter_max) skip <- TRUE
      }
      if (!skip && !is.na(x[[se_col]]) && x[[se_col]] > 0) {
        effect_df <<- dplyr::bind_rows(effect_df, tibble::tibble(
          inst    = as_chr_or_na(x[[id_col]]),
          feature = as_chr_or_na(x[[feature_col]]),
          beta    = as_num_or_na(x[[beta_col]]),
          se      = as_num_or_na(x[[se_col]]),
          p       = as_num_or_na(x[[p_col]])
        ))
      }
    }
  }

  finalize_current_snp <- function() {
    if (is.null(last_snp) || nrow(effect_df) == 0) return(invisible(NULL))
    effect_df <<- dplyr::distinct(effect_df)
    inst_exposures <- dplyr::filter(effect_df, p < p_thresh)[[feature_col]]
    message("It is an instrument for ", length(inst_exposures), " exposures.")
    exp_data <- build_exp_data(effect_df, inst_exposures)
    if (!is.null(exp_data) && nrow(exp_data) > 0) {
      match_idx <- match(exp_data$exp_out, result_df$exp_out)
      has_match <- which(!is.na(match_idx))
      if (length(has_match) > 0) {
        m <- match_idx[has_match]; i <- has_match
        result_df$inst[m]     <<- purrr::map2(result_df$inst[m],     exp_data$inst[i],     ~c(.x, .y))
        result_df$beta_exp[m] <<- purrr::map2(result_df$beta_exp[m], exp_data$beta_exp[i], ~c(.x, .y))
        result_df$se_exp[m]   <<- purrr::map2(result_df$se_exp[m],   exp_data$se_exp[i],   ~c(.x, .y))
        result_df$beta_out[m] <<- purrr::map2(result_df$beta_out[m], exp_data$beta_out[i], ~c(.x, .y))
        result_df$se_out[m]   <<- purrr::map2(result_df$se_out[m],   exp_data$se_out[i],   ~c(.x, .y))
      }
      if (any(is.na(match_idx))) {
        result_df <<- dplyr::bind_rows(result_df, exp_data[is.na(match_idx), ])
      }
    }
  }

  # ---- read with fixed column types ----------------------------------------
  # Build a col_types spec using the provided column names
  ct <- readr::cols()
  force_char <- c(id_col, chr_col, ref_col, alt_col, feature_col)
  for (nm in unique(force_char)) ct$cols[[nm]] <- readr::col_character()
  force_num  <- c(pos_col, beta_col, se_col, p_col, filter_col)
  for (nm in unique(force_num[!is.na(force_num)])) ct$cols[[nm]] <- readr::col_double()

  cb <- readr::SideEffectChunkCallback$new(f)
  readr::read_delim_chunked(
    file        = filename,
    callback    = cb,
    chunk_size  = 1,                 # process 1 row at a time
    delim       = delim,
    col_types   = ct,
    progress    = FALSE
  )

  finalize_current_snp()

  # standardize list-column types
  result_df <- result_df %>%
    dplyr::mutate(
      beta_exp = purrr::map(beta_exp, as.numeric),
      se_exp   = purrr::map(se_exp,   as.numeric),
      beta_out = purrr::map(beta_out, as.numeric),
      se_out   = purrr::map(se_out,   as.numeric),
      inst     = purrr::map(inst,     as.character)
    )

  return(list(result_df = result_df, inst_df = inst_df))
}

join_bigbrain <- function(...) {
  input_res <- list(...)
  inst_out   <- tibble::tibble(inst=character(), chr=character(), pos=numeric(),
                               ref=character(), alt=character())
  result_out <- tibble::tibble(exp_out=character(), Exposure=character(), Outcome=character(),
                               inst=list(), beta_exp=list(), se_exp=list(),
                               beta_out=list(), se_out=list())
  
  aggregate_res <- function(res){
    inst_out <<- dplyr::bind_rows(inst_out, res$inst_df)
    
    idx <- seq_len(nrow(res$result_df))
    match_idx <- match(res$result_df$exp_out, result_out$exp_out)
    has_match <- which(!is.na(match_idx))
    if (length(has_match) > 0) {
      i <- has_match; m <- match_idx[has_match]
      result_out$inst[m]     <<- purrr::map2(result_out$inst[m],     res$result_df$inst[i],     ~append(.x, .y))
      result_out$beta_exp[m] <<- purrr::map2(result_out$beta_exp[m], res$result_df$beta_exp[i], ~append(.x, .y))
      result_out$se_exp[m]   <<- purrr::map2(result_out$se_exp[m],   res$result_df$se_exp[i],   ~append(.x, .y))
      result_out$beta_out[m] <<- purrr::map2(result_out$beta_out[m], res$result_df$beta_out[i], ~append(.x, .y))
      result_out$se_out[m]   <<- purrr::map2(result_out$se_out[m],   res$result_df$se_out[i],   ~append(.x, .y))
    }
    if (any(is.na(match_idx))) {
      result_out <<- dplyr::bind_rows(result_out, res$result_df[is.na(match_idx), ])
    }
  }
  
  purrr::walk(input_res, aggregate_res)
  list(inst_df = inst_out, result_df = result_out)
}
