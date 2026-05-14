

choose_contrast_mode <- function(meta,
                                 treatment_col = "treatment",
                                 treatment_type_col = "treatment_type") {
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  
  if (!treatment_col %in% colnames(meta)) {
    stop("Missing treatment column: ", treatment_col)
  }
  
  if (treatment_type_col %in% colnames(meta)) {
    tt <- trimws(tolower(as.character(meta[[treatment_type_col]])))
    has_control <- any(tt == "control", na.rm = TRUE)
    has_treatment <- any(tt == "treatment", na.rm = TRUE)
    
    if (has_control && has_treatment) {
      return(list(
        mode = "treatment_vs_control",
        control_treatment = unique(as.character(meta[[treatment_col]][tt == "control"]))[1]
      ))
    }
  }
  

  choices <- c(
    "treatment_vs_control",
    "pairwise",
    "all_samples_as_is",
    "average_treatment_vs_control",
    "manual_control"
  )
  
  mode <- utils::select.list(
    choices = choices,
    title = "Choose contrast mode",
    multiple = FALSE, graphics = T
  )
  
  if (!nzchar(mode)) {
    stop("No contrast mode selected.")
  }
  
  if (mode == "manual_control") {
    control_treatment <- utils::select.list(
      choices = unique(as.character(meta[[treatment_col]])),
      title = "Choose the control treatment",
      multiple = FALSE, graphics = T
    )
    if (!nzchar(control_treatment)) {
      stop("No control treatment selected.")
    }
    return(list(
      mode = "treatment_vs_control",
      control_treatment = control_treatment
    ))
  }
  
  return(list(mode = mode, control_treatment = NULL))
  
  
  # list(mode = "pairwise", control_treatment = NULL)
}


build_group_contrast_matrix <- function(score_matrix,
                                        group_meta,
                                        mode = c("treatment_vs_control",
                                                 "pairwise",
                                                 "all_samples_as_is",
                                                 "average_treatment_vs_control"),
                                        treatment_col = "treatment",
                                        treatment_type_col = "treatment_type",
                                        sample_col = "sample",
                                        control_treatment = NULL
                                        ) {
  mode <- match.arg(mode)
  group_meta <- as.data.frame(group_meta, stringsAsFactors = FALSE)
  
  required_cols <- c(sample_col, treatment_col)
  missing_cols <- setdiff(required_cols, colnames(group_meta))
  if (length(missing_cols) > 0) {
    stop("group_meta is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  sample_ids <- as.character(group_meta[[sample_col]])
  sample_ids <- sample_ids[sample_ids %in% colnames(score_matrix)]
  
  if (length(sample_ids) == 0) {
    stop("No group samples found in score_matrix.")
  }
  
  submat <- score_matrix[, sample_ids, drop = FALSE]
  treatments <- as.character(group_meta[[treatment_col]])
  names(treatments) <- as.character(group_meta[[sample_col]])
  
  if (mode == "all_samples_as_is") {
    return(list(
      contrast_matrix = submat,
      contrast_info = data.frame(
        contrast = colnames(submat),
        mode = "all_samples_as_is",
        stringsAsFactors = FALSE
      )
    ))
  }
  
  if (mode == "pairwise") {
    trt_levels <- unique(treatments)
    trt_levels <- trt_levels[!is.na(trt_levels) & nzchar(trimws(trt_levels))]
    
    if (length(trt_levels) < 2) {
      stop("Need at least 2 treatment values for pairwise contrasts.")
    }
    
    pairs <- combn(trt_levels, 2, simplify = FALSE)
    
    contrast_list <- lapply(pairs, function(p) {
      a <- p[1]
      b <- p[2]
      
      a_samples <- names(treatments)[treatments == a]
      b_samples <- names(treatments)[treatments == b]
      
      a_mean <- rowMeans(submat[, a_samples, drop = FALSE], na.rm = TRUE)
      b_mean <- rowMeans(submat[, b_samples, drop = FALSE], na.rm = TRUE)
      
      out <- a_mean - b_mean
      out
    })
    
    contrast_names <- vapply(pairs, function(p) {
      paste0(p[1], "_vs_", p[2])
    }, character(1))
    
    contrast_matrix <- do.call(cbind, contrast_list)
    colnames(contrast_matrix) <- contrast_names
    rownames(contrast_matrix) <- rownames(submat)
    
    return(list(
      contrast_matrix = contrast_matrix,
      contrast_info = data.frame(
        contrast = contrast_names,
        mode = "pairwise",
        stringsAsFactors = FALSE
      )
    ))
  }
  
  if (mode %in% c("treatment_vs_control", "average_treatment_vs_control")) {
    if (is.null(control_treatment) || !nzchar(control_treatment)) {
      if (treatment_type_col %in% colnames(group_meta)) {
        tt <- trimws(tolower(as.character(group_meta[[treatment_type_col]])))
        controls <- unique(as.character(group_meta[[treatment_col]][tt == "control"]))
        controls <- controls[!is.na(controls) & nzchar(trimws(controls))]
        if (length(controls) >= 1) {
          control_treatment <- controls[1]
        }
      }
    }
    
    if (is.null(control_treatment) || !nzchar(control_treatment)) {
      stop("Could not determine control treatment.")
    }
    
    ctrl_samples <- names(treatments)[treatments == control_treatment]
    if (length(ctrl_samples) == 0) {
      stop("No samples found for control treatment: ", control_treatment)
    }
    
    other_treatments <- setdiff(unique(treatments), control_treatment)
    other_treatments <- other_treatments[!is.na(other_treatments) & nzchar(trimws(other_treatments))]
    
    if (length(other_treatments) == 0) {
      stop("No non-control treatments found in this group.")
    }
    
    
    if (mode == "average_treatment_vs_control") {
      trt_samples <- names(treatments)[treatments %in% other_treatments]
      
      ctrl_mean <- rowMeans(submat[, ctrl_samples, drop = FALSE], na.rm = TRUE)
      trt_mean <- rowMeans(submat[, trt_samples, drop = FALSE], na.rm = TRUE)
      
      contrast_matrix <- cbind(treatment_vs_control = trt_mean - ctrl_mean)
      rownames(contrast_matrix) <- rownames(submat)
      
      return(list(
        contrast_matrix = contrast_matrix,
        contrast_info = data.frame(
          contrast = "treatment_vs_control",
          treatment = paste(other_treatments, collapse = ";"),
          control = control_treatment,
          mode = "average_treatment_vs_control",
          stringsAsFactors = FALSE
        )
      ))
    }
    

    trt_samples <- names(treatments)[treatments %in% other_treatments]
    contrast_matrix = submat[, trt_samples, drop = FALSE] - rowMeans(submat[, ctrl_samples, drop = FALSE])
    contrast_names <- paste0(trt_samples, "_vs_", ctrl_samples)
    colnames(contrast_matrix) <- contrast_names
    rownames(contrast_matrix) <- rownames(submat)
    
    return(list(
      contrast_matrix = contrast_matrix,
      contrast_info = data.frame(
        contrast = contrast_names,
        treatment = other_treatments,
        control = control_treatment,
        mode = "treatment_vs_control",
        stringsAsFactors = FALSE
      )
    ))
  }
  
  stop("Unsupported contrast mode: ", mode)
}

run_single_mcd <- function(analysis_matrix, B = 25) {
  X <- as.matrix(analysis_matrix)
  
  if (!is.numeric(X)) {
    stop("analysis_matrix must be numeric.")
  }
  
  if (is.null(rownames(X))) {
    stop("analysis_matrix must have rownames corresponding to genes/features.")
  }
  
  if (nrow(X) < 2) {
    stop("analysis_matrix must have at least 2 rows.")
  }
  
  if (ncol(X) < 1) {
    stop("analysis_matrix must have at least 1 column.")
  }
  
  genes <- rownames(X)
  m <- nrow(X)
  k <- ncol(X)
  
  if (k == 1) {
    # one-dimensional fallback
    z <- as.numeric(scale(X[, 1]))
    p.two.sided <- 2 * pnorm(abs(z), lower.tail = FALSE)
    
    res <- data.frame(
      gene = genes,
      mahalanobis = z^2,
      p.chisq = pchisq(z^2, df = 1, lower.tail = FALSE),
      p.mdn = pchisq(z^2, df = 1, lower.tail = FALSE),
      p.sim = p.two.sided,
      rank_mahalanobis = rank(-(z^2), ties.method = "min"),
      stringsAsFactors = FALSE
    )
    rownames(res) <- genes
    return(res)
  }
  
  mcd.res <- robustbase::covMcd(X)
  mhd <- mahalanobis(X, center = mcd.res$center, cov = mcd.res$cov)
  
  p.chisq <- pchisq(mhd, df = k, lower.tail = FALSE)
  
  qqr <- sort(mhd) / qchisq(seq_along(mhd) / (length(mhd) + 1), df = k)
  mdn.ratio <- median(qqr, na.rm = TRUE)
  mr.mhd <- mhd / mdn.ratio
  p.mdn.scl <- pchisq(mr.mhd, df = k, lower.tail = FALSE)
  
  gren.ebp <- grenander.ebp(p.mdn.scl)
  
  boot_prob <- 1 - gren.ebp$ebp.null
  boot_prob[!is.finite(boot_prob)] <- 0
  if (sum(boot_prob) <= 0) {
    boot_prob <- rep(1 / m, m)
  } else {
    boot_prob <- boot_prob / sum(boot_prob)
  }
  
  sim.list <- replicate(B, {
    idx <- sample.int(n = m, size = m, replace = TRUE, prob = boot_prob)
    sim.X <- X[idx, , drop = FALSE]
    sim.mcd <- robustbase::covMcd(sim.X, raw = TRUE)
    mahalanobis(sim.X, center = sim.mcd$center, cov = sim.mcd$cov)
  }, simplify = FALSE)
  
  sim.dist.all <- unlist(sim.list, use.names = FALSE)
  p.sim <- vapply(mhd, function(d) mean(sim.dist.all >= d), numeric(1))
  
  res <- data.frame(
    gene = genes,
    mahalanobis = mhd,
    p.chisq = p.chisq,
    p.mdn = p.mdn.scl,
    p.sim = p.sim,
    rank_mahalanobis = rank(-mhd, ties.method = "min"),
    stringsAsFactors = FALSE
  )
  
  rownames(res) <- genes
  res
}



grenander =  function(F, type=c("decreasing", "increasing"))  # from fdrtool package
{
  if( !any(class(F) == "ecdf") ) stop("ecdf object required as input!")
  type = match.arg(type)
  if (type == "decreasing")
  {
    # find least concave majorant of ECDF
    ll = fdrtool::gcmlcm(environment(F)$x, environment(F)$y, type="lcm")
  }
  else
  {
    # find greatest convex minorant of ECDF
    l = length(environment(F)$y)
    ll = fdrtool::gcmlcm(environment(F)$x, c(0,environment(F)$y[-l]), type="gcm")
  }
  f.knots = ll$slope.knots
  f.knots = c(f.knots, f.knots[length(f.knots)])
  g = list(F=F,
           x.knots=ll$x.knots,
           F.knots=ll$y.knots,
           f.knots=f.knots)
  class(g) = "grenander"
  return(g)
}

grenander.ebp =function(p)     # Compute the grenander.ebp from a vector of p-values
{
  na=is.na(p)
  p.edf=ecdf(p[!na])
  gren.res=grenander(p.edf)
  gren.pdf=approx(gren.res$x.knots,gren.res$f.knots,xout=p)$y
  gren.ebp=min(gren.res$f.knots)/gren.pdf
  ebp.null=pval.pdf=rep(NA,length(p))
  ebp.null[!na]=gren.ebp
  pval.pdf[!na]=gren.pdf   
  return(cbind.data.frame(pval = p,pval.pdf=pval.pdf,ebp.null=ebp.null))
}

# Langaas, M., Lindqvist, B., 2005. Estimating the proportion of true null hypotheses, 
# with application to DNA microarray data. J.R. Statist. Soc. B67, part4, 555-572. Strimmer, K. 2008. A
# unified approach to false discovery rate estimation. BMC Bioinformatics 9: 303. Strimmer, K.
# 2008. fdrtool: a versatile R package for estimating local and tail area-based false discovery rates.
# Bioinformatics 24: 1461-1462. Grenander, U. (1956). On the theory of mortality measurement, 
# Part II. Skan. Aktuarietidskr, 39, 125–153.

cwMCD=function(X,alpha=0.75,quant=0.99,
               crit=1e-4,noCits=100,lmin=1e-4,
               fixedCenter = FALSE, checkPars=list(silent = TRUE)
               )
{
  if (is.null(rownames(X))) rownames(X)=paste0("row_",1:nrow(X))
  if (is.null(colnames(X))) colnames(X)=paste0("clm_",1:ncol(X))
  
  # cellMCD
  genes = rownames(X)
  out=cellMCD(X,alpha,quant,crit,noCits,lmin,fixedCenter = F,checkPars)

  # convert a matrix to long data.frame
  mat_to_df <- function(mat, value_name) {
    df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    colnames(df) <- c("gene", "contrast", value_name)
    df
  }
  
  # build long tables
  df_X     <- mat_to_df(out$X,     "X")
  df_W     <- mat_to_df(out$W,     "W")
  df_flag  <- mat_to_df(out$W==0,   "flag")
  df_pred  <- mat_to_df(out$preds, "pred")
  df_csd   <- mat_to_df(out$csds,  "csd")
  df_Ximp  <- mat_to_df(out$Ximp,  "Ximp")
  df_Zres  <- mat_to_df(out$Zres,  "Zres")
  
  # combine
  df_long <- Reduce(
    function(x, y) merge(x, y, by = c("gene", "contrast"), sort = FALSE),
    list(df_X, df_W, df_flag, df_pred, df_csd, df_Ximp, df_Zres)
  )
  
  
  
  list(res.df = df_long,
       cw.res.obj = out[c('mu', 'S', 'raw.S', 
                             'locsca', 'nosteps', 
                             'quant')])
}
  
  
# helper function applied to each (group_id, gene) subset
summarize_cmcd_group <- function(df_long) {

  gene_summary <- df_long %>%
    mutate(abs_Z = abs(Zres)) %>%
    group_by(group_id, gene) %>%
    summarise(
      n_flag = sum(flag, na.rm = TRUE),
      prop_flag = mean(flag, na.rm = TRUE),
      ssZ = sum(Zres^2, na.rm = TRUE),
      max_abs_Z = max(abs_Z, na.rm = TRUE),
      mean_abs_Z_flag = if (any(flag, na.rm = TRUE)) {
        mean(abs_Z[flag], na.rm = TRUE)
      } else {
        0
      },
      n_nonmiss = sum(!is.na(X)),
      dominant_sign = if (n_nonmiss == 0) {
        NA_real_
      } else {
        sum(X > 0, na.rm = TRUE) / n_nonmiss -
          sum(X < 0, na.rm = TRUE) / n_nonmiss
      },
      local_dominant_sign_meaning = {
        idx <- which.max(abs(Zres))
        if (length(idx) == 0 || is.na(X[idx])) {
          NA_character_
        } else if (X[idx] > 0) {
          "Positive"
        } else if (X[idx] < 0) {
          "Negative"
        } else {
          "Mixed"
        }
      },
      .groups = "drop"
    ) %>%
    mutate(
      max_sq_Z = max_abs_Z^2,
      eff_n = if_else(max_sq_Z > 0, ssZ / max_sq_Z, 0),
      spread_score = if_else(ssZ > 0, 1 - (max_sq_Z / ssZ), 0),
      
      # keep your original scores too
      global_score_flag = prop_flag * mean_abs_Z_flag,
      local_score_flag  = (1 - prop_flag) * max_abs_Z,
      p_global_flag = if_else(
        global_score_flag + local_score_flag > 0,
        global_score_flag / (global_score_flag + local_score_flag),
        0
      ),
      
      # new Z-based scores
      global_score = ssZ * spread_score,
      local_score  = max_abs_Z^2 / ssZ,
      p_global = if_else(eff_n > 0, (eff_n - 1) / pmax(1, n_nonmiss - 1), 0),
      
      pattern_flag = case_when(
        n_flag == 0 ~ "none",
        p_global_flag >= 0.6~ "global",
        TRUE ~ "local"
      ),

      pattern = case_when(
        n_flag == 0 ~ "none",
        p_global >= 0.5 ~ "global",
        TRUE ~ "local"
      ),
      
      sign_meaning = case_when(
        is.na(dominant_sign) ~ NA_character_,
        dominant_sign > 0.5 ~ "Positive",
        dominant_sign < -0.5 ~ "Negative",
        TRUE ~ "Mixed"
      ),
      
      local_dominant_sign_meaning = case_when(
        pattern == "local" & sign_meaning == "Mixed" ~ local_dominant_sign_meaning,
        TRUE ~ sign_meaning
      )
    ) %>%
    mutate(
      local_dominant_sign = case_when(
        is.na(local_dominant_sign_meaning) ~ NA_real_,
        pattern == "local" & local_dominant_sign_meaning == "Positive" ~ 1,
        pattern == "local" & local_dominant_sign_meaning == "Negative" ~ -1,
        local_dominant_sign_meaning == "Positive" ~ abs(dominant_sign),
        local_dominant_sign_meaning == "Negative" ~ -abs(dominant_sign),
        sign_meaning == "Mixed" ~ 0,
        TRUE ~ dominant_sign
      ),
      
      gene_class = case_when(
        pattern == "none" ~ "None",
        pattern == "global" & sign_meaning == "Positive" ~ "Global Positive",
        pattern == "global" & sign_meaning == "Negative" ~ "Global Negative",
        pattern == "global" & sign_meaning == "Mixed" ~ "Global Mixed",
        pattern == "local" & local_dominant_sign_meaning == "Positive" ~ "Local Positive",
        pattern == "local" & local_dominant_sign_meaning == "Negative" ~ "Local Negative",
        TRUE ~ "None"
      )
    ) %>%
    select(
      group_id, gene,
      ssZ, eff_n, p_global,
      gene_class, pattern, local_dominant_sign_meaning,
      spread_score, n_flag, prop_flag, max_abs_Z,
      global_score, local_score,  
      global_score_flag, local_score_flag, p_global_flag, pattern_flag,
      dominant_sign, sign_meaning, local_dominant_sign
    )  
  gene_summary
}


# ============================================================
# MCD Pipeline: Metadata-driven comparison definition + export
# ============================================================




run_mcd_pipeline <- function(final_input,
                             B = 20) {
  
  if (!(exists("score_matrix", final_input)) |(is.null(final_input$score_matrix))) {
    stop("final_input$score_matrix is missing.")
  }
  
  if (!exists("metadata", final_input) | is.null(final_input$metadata)) {
    stop("final_input$metadata is missing.")
  }
  
  score_matrix <- as.matrix(final_input$score_matrix)
  meta <- as.data.frame(final_input$metadata, stringsAsFactors = FALSE)
  
  if (!is.numeric(score_matrix)) {
    stop("score_matrix must be numeric.")
  }
  
  if (is.null(rownames(score_matrix))) {
    if (!is.null(final_input$gene_data) &&
        length(final_input$gene_data) == nrow(score_matrix)) {
      rownames(score_matrix) <- as.character(final_input$gene_data)
    } else {
      stop("score_matrix must have rownames or matching gene_data.")
    }
  }
  
  group_sep = ","
  treatment_col = "treatment"
  treatment_type_col = "treatment_type"
  sample_col = "sample"
  interactive_if_needed = TRUE
  contrast_mode = NULL
  control_treatment = NULL
  
  required_meta_cols <- c(sample_col, "Group_By", treatment_col)
  missing_meta <- setdiff(required_meta_cols, colnames(meta))
  if (length(missing_meta) > 0) {
    stop("metadata is missing required columns: ",
         paste(missing_meta, collapse = ", "))
  }
  
  meta <- meta[meta[[sample_col]] %in% colnames(score_matrix), , drop = FALSE]
  meta <- meta[match(colnames(score_matrix), meta[[sample_col]]), , drop = FALSE]
  
  if (any(is.na(meta[[sample_col]]))) {
    stop("Not all score_matrix columns matched metadata.")
  }
  
  group_entries <- unique(trimws(as.character(meta$Group_By)))
  group_entries <- group_entries[!is.na(group_entries) & nzchar(group_entries)]
  
  if (length(group_entries) == 0 || all(tolower(group_entries) == "none")) {
    group_by_cols <- character(0)
  } else {
    if (length(group_entries) > 1) {
      stop("Conflicting Group_By specifications in metadata: ",
           paste(group_entries, collapse = ", "))
    }
    group_by_cols <- trimws(unlist(strsplit(group_entries[1], split = group_sep, fixed = TRUE)))
    group_by_cols <- group_by_cols[nzchar(group_by_cols)]
    
    invalid_group_cols <- setdiff(group_by_cols, colnames(meta))
    if (length(invalid_group_cols) > 0) {
      stop("Invalid Group_By columns: ",
           paste(invalid_group_cols, collapse = ", "))
    }
  }
  
  if (length(group_by_cols) == 0) {
    group_keys <- rep("all_samples", nrow(meta))
  } else {
    group_keys <- apply(meta[, group_by_cols, drop = FALSE], 1, function(x) {
      paste(paste(group_by_cols, x, sep = "="), collapse = "__")
    })
  }
  
  split_idx <- split(seq_len(nrow(meta)), group_keys)
  
  if (is.null(contrast_mode)) {
    contrast_choice <- choose_contrast_mode(
      meta = meta,
      treatment_col = treatment_col,
      treatment_type_col = treatment_type_col
    )
    contrast_mode <- contrast_choice$mode
    control_treatment <- contrast_choice$control_treatment
  }
  
  cat("Using contrast mode: ", contrast_mode, "\n")
  if (!is.null(control_treatment) && nzchar(control_treatment)) {
    cat("Control treatment: ", control_treatment, "\n\n")
  }
  cat(" Running ", length(split_idx), " group(s). \n")
  
  mcd_results_by_group <- vector("list", length(split_idx))
  names(mcd_results_by_group) <- names(split_idx)
  
  cw_mcd_results_by_group <- vector("list", length(split_idx))
  names(cw_mcd_results_by_group) <- names(split_idx)
  
  contrast_matrices <- vector("list", length(split_idx))
  names(contrast_matrices) <- names(split_idx)
  
  contrast_info_by_group <- vector("list", length(split_idx))
  names(contrast_info_by_group) <- names(split_idx)
  
  for (group_name in names(split_idx)) {
    
    idx <- split_idx[[group_name]]
    group_meta <- meta[idx, , drop = FALSE]
    group_samples <- as.character(group_meta[[sample_col]])
    
    cat("Building contrasts for group: ", group_name, "\n")
    
    built <- build_group_contrast_matrix(
      score_matrix = score_matrix,
      group_meta = group_meta,
      mode = contrast_mode,
      treatment_col = treatment_col,
      treatment_type_col = treatment_type_col,
      sample_col = sample_col,
      control_treatment = control_treatment
    )
    
    contrast_matrix <- built$contrast_matrix
    contrast_info <- built$contrast_info
    
    if (ncol(contrast_matrix) < 1) {
      stop("Contrast matrix for group '", group_name, "' has no columns.")
    }
    
    contrast_matrices[[group_name]] <- contrast_matrix
    contrast_info_by_group[[group_name]] <- contrast_info
    
    # run mcd
    mcd_res <- run_single_mcd(
      analysis_matrix = contrast_matrix,
      B = B
    )
    
    mcd_res$group_id <- group_name
    rownames(mcd_res) <- NULL
    mcd_results_by_group[[group_name]] <- mcd_res
    
    # run cellMCD
    cw_mcd_list <- cwMCD(contrast_matrix)
    cw_mcd_res <- cw_mcd_list$res.df
    cw_mcd_res = merge(cw_mcd_res, contrast_info, by = 'contrast')
    cw_mcd_res = merge(cw_mcd_res, group_meta, by = 'treatment')
    
    
    remove.cols = c('control', 'mode', 'sample', 'treatment_type', 'rep',
                    'include_sample', 'Group_By', 'original_samples')
    
    keep.cols = setdiff(colnames(cw_mcd_res), remove.cols)
    keep.cols = c(
      colnames(cw_mcd_list$res.df), 
      setdiff(keep.cols, 
              colnames(cw_mcd_list$res.df))
    )
    
    cw_mcd_res=cw_mcd_res[, keep.cols]
    cw_mcd_res$group_id <- group_name

    cw_gene_sum = summarize_cmcd_group(cw_mcd_res)

    cw_mcd_results_by_group[[group_name]] <- list(cw_mcd_res, cw_gene_sum,
                                                  cw_mcd_list$cw.res.obj)
    
  }
  
  mcd_combined_results <- do.call(rbind, mcd_results_by_group)
  rownames(mcd_combined_results) <- NULL
  
  cw_dfs <- lapply(cw_mcd_results_by_group, `[[`, 1)
  cw_mcd_combined_results <- do.call(rbind, cw_dfs)
  rownames(cw_mcd_combined_results) <- NULL
  
  list(
    group_by = group_by_cols,
    contrast_mode = contrast_mode,
    control_treatment = control_treatment,
    metadata_used = meta,
    contrast_matrices = contrast_matrices,
    contrast_info_by_group = contrast_info_by_group,
    mcd_results_by_group = mcd_results_by_group,
    cw_mcd_results_by_group = cw_mcd_results_by_group,
    mcd_combined_results = mcd_combined_results,
    cw_mcd_combined_results = cw_mcd_combined_results
    
  )
}

# ============================================================
# Make MCD Report Tables (run with main.R: )
# ============================================================



build_outlier_summary_table <- function(res, alpha = 0.05, top_n = 25) {
  out <- lapply(names(res$mcd_results_by_group), function(g) {
    tbl <- res$mcd_results_by_group[[g]]
    cm <- res$contrast_matrices[[g]]
    
    
    get_primary_p_col <- function(tbl, default = NULL) {
      
      if(is.null(default)){
        for (nm in c("p.sim", "p.mdn", "p.chisq", "MCD.pval")) {
          if (nm %in% colnames(tbl)) return(nm)
        }
      } else {
        if (default %in% colnames(tbl)) return(default)
      }
      
      stop("No supported p-value column found.")
    }
    
    p_col <- get_primary_p_col(tbl)
    x <- tbl
    x$p_value_used <- p_col
    x$p_value <- x[[p_col]]
    x$outlier_flag <- x$p_value < alpha
    x$top_gene_flag <- x$outlier_flag | (x$rank_mahalanobis <= top_n)
    
    keep_cols <- c(
      "gene",
      "mahalanobis",
      "rank_mahalanobis",
      "p_value",
      "top_gene_flag",
      intersect(c("group_id", "timepoint", "cell_line", "treatment"), 
                colnames(x))
    )
    
    x[, keep_cols, drop=FALSE]
  })
  
  out.df = do.call(rbind, out)
  cm.all = do.call(cbind, res$contrast_matrices)
  cm.all = as.data.frame(cm.all)
  cm.all$gene = rownames(cm.all)
  out.df = merge(cm.all, out.df, by = 'gene')


}

run_outlier_hdbscan <- function(outlier_summary_tbl,
                                         contrast_cols,
                                         alpha = 0.05,
                                         minPts = 10,
                                         scale_rows = TRUE) {
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required.")
  }
  
  missing_cols <- setdiff(contrast_cols, colnames(outlier_summary_tbl))
  if (length(missing_cols) > 0) {
    stop(
      "Missing contrast columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  row_scale <- function(mat) {
    t(apply(mat, 1, function(x) {
      s <- stats::sd(x, na.rm = TRUE)
      m <- mean(x, na.rm = TRUE)
      
      if (is.na(s) || s == 0) {
        rep(0, length(x))
      } else {
        as.numeric((x - m) / s)
      }
    }))
  }
  
  # start with full table
  cluster_table <- outlier_summary_tbl
  
  # initialize new columns for all rows
  cluster_table$hdbscan_cluster <- NA_integer_
  cluster_table$hdbscan_is_noise <- NA
  cluster_table$hdbscan_membership_prob <- NA_real_
  
  # rows to cluster
  keep_idx <- outlier_summary_tbl$flag | outlier_summary_tbl$p_value < alpha
  
  if (!any(keep_idx, na.rm = TRUE)) {
    return(cluster_table)
  }
  
  dat_sub <- outlier_summary_tbl[keep_idx, , drop = FALSE]
  
  # keep original row indices so we can write results back
  dat_sub$.orig_row_id <- which(keep_idx)
  
  split_tbl <- split(dat_sub, dat_sub$group_id)
  
  cluster_list <- lapply(split_tbl, function(df_group) {
    mat <- as.matrix(df_group[, contrast_cols, drop = FALSE])
    
    if (scale_rows) {
      mat <- row_scale(mat)
      colnames(mat) <- contrast_cols
    }
    
    # if too few rows, mark all as noise
    if (nrow(mat) < minPts) {
      df_group$hdbscan_cluster <- 0L
      df_group$hdbscan_is_noise <- TRUE
      df_group$hdbscan_membership_prob <- NA_real_
      return(df_group)
    }
    
    hdb <- dbscan::hdbscan(mat, minPts = minPts)
    
    df_group$hdbscan_cluster <- hdb$cluster
    df_group$hdbscan_is_noise <- hdb$cluster == 0L
    df_group$hdbscan_membership_prob <- hdb$membership_prob
    
    df_group
  })
  
  clustered_sub <- do.call(rbind, cluster_list)
  rownames(clustered_sub) <- NULL
  
  # write clustering results back into full table
  cluster_table$hdbscan_cluster[clustered_sub$.orig_row_id] <- clustered_sub$hdbscan_cluster
  cluster_table$hdbscan_is_noise[clustered_sub$.orig_row_id] <- clustered_sub$hdbscan_is_noise
  cluster_table$hdbscan_membership_prob[clustered_sub$.orig_row_id] <- clustered_sub$hdbscan_membership_prob
  
  cluster_table
}
