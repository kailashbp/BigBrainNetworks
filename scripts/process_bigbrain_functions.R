# The function

read_bigbrain <- function(filename, beta_col, p_col, id_col = "variant_id",
                          se_col, p_thresh = 5e-08, chr_col = "chr",
                          pos_col = "pos", ref_col = "ref", alt_col = "alt",
                          feature_col = "feature", delim = "\t",
                          filter_col = NULL, filter_max = NULL,
                          filter_min = NULL){
  
  if (missing(se_col) || is.null(se_col))
    stop("`se_col` must be provided (cannot be NULL).")
  
  last_snp <- NULL
  inst_df <- NULL
  effect_df <- NULL
  result_df <- NULL
  
  # helper: merge exp_data into result_df by exp_out, appending list-cols
  merge_exp_data <- function(exp_data) {
    if (is.null(exp_data) || nrow(exp_data) == 0) return(invisible(NULL))
    if (nrow(result_df) == 0) {
      result_df <<- exp_data
      return(invisible(NULL))
    }
    idx <- seq_len(nrow(exp_data))
    match_idx <- match(exp_data$exp_out, result_df$exp_out)
    # update matches
    if (any(!is.na(match_idx))) {
      i <- idx[!is.na(match_idx)]
      m <- match_idx[!is.na(match_idx)]
      result_df[m, ] <<- result_df[m, ] %>% dplyr::mutate(
        inst     = purrr::map2(inst,     exp_data$inst[i],     append),
        beta_exp = purrr::map2(beta_exp, exp_data$beta_exp[i], append),
        se_exp   = purrr::map2(se_exp,   exp_data$se_exp[i],   append),
        beta_out = purrr::map2(beta_out, exp_data$beta_out[i], append),
        se_out   = purrr::map2(se_out,   exp_data$se_out[i],   append)
      )
    }
    # bind new pairs
    result_df <<- dplyr::bind_rows(result_df, exp_data[is.na(match_idx), ])
  }
  
  process_effect_df <- function(effect_df) {
    effect_df <- dplyr::distinct(effect_df)
    inst_exposures <- effect_df %>%
      dplyr::filter(p < p_thresh) %>%
      dplyr::pull(feature)
    if (length(inst_exposures) == 0) return(NULL)
    
    exp_data <- purrr::map_dfr(inst_exposures, function(exp) {
      left  <- effect_df %>%
        dplyr::filter(feature == exp) %>%
        dplyr::select(-p) %>%
        dplyr::rename(Exposure = feature, beta_exp = beta, se_exp = se)
      
      right <- effect_df %>%
        dplyr::select(-inst, -p) %>%
        dplyr::rename(Outcome = feature, beta_out = beta, se_out = se)
      
      dplyr::cross_join(left, right) %>%
        # dplyr::filter(Exposure != Outcome) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          exp_out  = paste(Exposure, Outcome, sep = "_"),
          inst     = list(inst),
          beta_exp = list(beta_exp),
          se_exp   = list(se_exp),
          beta_out = list(beta_out),
          se_out   = list(se_out)
        ) %>%
        dplyr::ungroup()
    })
    
    if (nrow(exp_data) == 0) return(NULL)
    exp_data
  }
  
  f <- function(x, pos){
    # print(colnames(x))
    # print(pos)
    if (is.null(last_snp)) {
      message("Initializing")
      inst_df <<- tibble::tibble(
        inst = x[[id_col]], chr = x[[chr_col]], pos = x[[pos_col]],
        ref  = x[[ref_col]], alt = x[[alt_col]]
      )
      effect_df <<- tibble::tibble(
        inst = x[[id_col]], feature = x[[feature_col]],
        beta = x[[beta_col]], se = x[[se_col]], p = x[[p_col]]
      )
      last_snp <<- x[[id_col]]
      result_df <<- tibble::tibble(
        exp_out = character(), Exposure = character(), Outcome = character(),
        inst = list(), beta_exp = list(), se_exp = list(),
        beta_out = list(), se_out = list()
      )
      
    } else if (x[[id_col]] != last_snp) {
      cat(paste0("\nProcessing SNP ", last_snp, ".\n"))
      exp_data <- process_effect_df(effect_df)
      merge_exp_data(exp_data)
      
      # Reset effect_df for new SNP
      effect_df <<- tibble::tibble(
        inst = x[[id_col]], feature = x[[feature_col]],
        beta = x[[beta_col]], se = x[[se_col]], p = x[[p_col]]
      )
      # Add new SNP info to inst_df
      inst_df <<- rbind(inst_df, tibble::tibble(
        inst = x[[id_col]], chr = x[[chr_col]],
        pos = x[[pos_col]], ref = x[[ref_col]],
        alt = x[[alt_col]]
      ))
      last_snp <<- x[[id_col]]
      
    } else {
      # Continue adding SNP effects to effect_df.
      # Optionally skip SNPs based on filer_col.
      skip <- FALSE
      if (!is.null(filter_col)) {
        val <- x[[filter_col]]
        if (is.na(val)) {
          skip <- TRUE
        } else {
          if (is.null(filter_min) && is.null(filter_max))
            stop("Filter col specified but no values (min/max) provided.")
          if (!is.null(filter_min) && val < filter_min) skip <- TRUE
          if (!is.null(filter_max) && val > filter_max) skip <- TRUE
        }
      }
      
      if (!skip && is.finite(x[[se_col]]) && x[[se_col]] > 0) {
        effect_df <<- rbind(effect_df, tibble::tibble(
          inst = x[[id_col]], feature = x[[feature_col]],
          beta = x[[beta_col]], se = x[[se_col]], p = x[[p_col]]
        ))
      }
    }
  }
  
  callback_fn <- readr::SideEffectChunkCallback$new(f) # Call f for every chunk, and not return anything, but mutate global variable effect_df & result_df outside callback
  readr::read_delim_chunked(filename, callback_fn, chunk_size = 1, delim = delim, show_col_types = FALSE)
  
  # final flush for the last SNP
  exp_data <- process_effect_df(effect_df)
  merge_exp_data(exp_data)
  
  # de-duplicate inst_df
  inst_df <- dplyr::distinct(inst_df)
  
  list(result_df = result_df, inst_df = inst_df)
}

join_bigbrain <- function(...) {
  input_res <- list(...)
  # filter out NULLs or empty results
  input_res <- purrr::compact(input_res)
  if (length(input_res) == 0) {
    return(list(inst_df = tibble::tibble(), result_df = tibble::tibble()))
  }
  
  inst_out   <- tibble::tibble(inst=character(), chr=character(), pos=numeric(),
                               ref=character(), alt=character())
  result_out <- tibble::tibble(exp_out=character(), Exposure=character(), Outcome=character(),
                               inst=list(), beta_exp=list(), se_exp=list(),
                               beta_out=list(), se_out=list())
  
  aggregate_res <- function(res){
    # skip if malformed
    if (is.null(res$result_df) || is.null(res$inst_df) || nrow(res$result_df) == 0) return(invisible())
 
    # enforce column types
    res$inst_df <- res$inst_df %>%
    dplyr::mutate(
      inst = as.character(inst),
      chr  = as.character(chr),
      ref  = as.character(ref),
      alt  = as.character(alt),
      pos  = as.numeric(pos)
    )
   
    # add instruments   
    inst_out <<- dplyr::bind_rows(inst_out, res$inst_df)

    # collapse current chunk to one row per key with list-cols
    res_c <- res$result_df %>%
      dplyr::group_by(exp_out, Exposure, Outcome) %>%
      dplyr::summarise(
        inst     = list(inst),
        beta_exp = list(beta_exp),
        se_exp   = list(se_exp),
        beta_out = list(beta_out),
        se_out   = list(se_out),
        .groups = "drop"
      )
    
    # Use first chunk as is
    if (nrow(result_out) == 0) {
      result_out <<- res_c
      return(invisible())
    }

    # match on all keys (safer than exp_out alone)
    key_join <- dplyr::left_join(
      res_c %>% dplyr::mutate(.row_id = dplyr::row_number()),
      result_out %>% dplyr::mutate(.idx = dplyr::row_number()),
      by = c("exp_out","Exposure","Outcome")
    )

    has_match <- stats::na.omit(key_join$.idx)
    if (length(has_match)) {
      i <- key_join$.row_id[!is.na(key_join$.idx)]
      m <- key_join$.idx[!is.na(key_join$.idx)]
      result_out$inst[m]     <<- purrr::map2(result_out$inst[m],     res_c$inst[i],     ~c(.x, .y))
      result_out$beta_exp[m] <<- purrr::map2(result_out$beta_exp[m], res_c$beta_exp[i], ~c(.x, .y))
      result_out$se_exp[m]   <<- purrr::map2(result_out$se_exp[m],   res_c$se_exp[i],   ~c(.x, .y))
      result_out$beta_out[m] <<- purrr::map2(result_out$beta_out[m], res_c$beta_out[i], ~c(.x, .y))
      result_out$se_out[m]   <<- purrr::map2(result_out$se_out[m],   res_c$se_out[i],   ~c(.x, .y))
    }

    if (any(is.na(key_join$.idx))) {
      result_out <<- dplyr::bind_rows(
        result_out,
        res_c[key_join$.row_id[is.na(key_join$.idx)], ]
      )
    }
  }

  purrr::walk(input_res, aggregate_res)
  inst_out <- dplyr::distinct(inst_out)

  list(inst_df = inst_out, result_df = result_out) 
}
