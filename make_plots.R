

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

get_top_genes_for_group <- function(dt, group_name, top_n = 10,
                                    rank_col = "mahalanobis",
                                    keep_mode = c("interesting_only", "rowwise_only", "all")) {
  keep_mode <- match.arg(keep_mode)
  
  sub <- dt[group_id == group_name]
  
  if (keep_mode == "rowwise_only") {
    keep_genes <- unique(sub[top_gene_flag == TRUE, gene])
  } else if (keep_mode == "interesting_only") {
    keep_genes <- unique(sub[top_gene_flag == TRUE | pattern != "none", gene])
  } else {
    keep_genes <- unique(sub$gene)
  }
  
  sub <- sub[gene %in% keep_genes]
  
  # one row per gene
  gene_dt <- unique(sub[, .(gene, rank_value = get(rank_col))])
  gene_dt <- gene_dt[!is.na(rank_value)]
  
  if (nrow(gene_dt) == 0L) return(character(0))
  
  gene_dt <- gene_dt[order(-rank_value)]
  head(gene_dt$gene, top_n)
}


get_flag_summary_text <- function(dt, group_name, mode = c("cellwise", "rowwise")) {
  mode <- match.arg(mode)
  
  sub <- dt[group_id == group_name]
  
  if (mode == "cellwise") {
    gene_flag <- unique(sub[, .(gene, flag)])[
      , .(flagged = any(flag, na.rm = TRUE)), by = gene
    ]
    n_flag <- sum(gene_flag$flagged, na.rm = TRUE)
    n_total <- nrow(gene_flag)
    txt <- paste0("Cell-wise flagged genes: ", n_flag, " / ", n_total)
  } else {
    gene_flag <- unique(sub[, .(gene, top_gene_flag)])[
      , .(flagged = any(top_gene_flag, na.rm = TRUE)), by = gene
    ]
    n_flag <- sum(gene_flag$flagged, na.rm = TRUE)
    n_total <- nrow(gene_flag)
    txt <- paste0("Row-wise flagged genes: ", n_flag, " / ", n_total)
  }
  
  txt
}


# ------------------------------------------------------------
# 7. Plot ssZ and p_global or eff_n
# ------------------------------------------------------------

plot_global_local_summary <- function(
    dt,
    x_var = c("p_global", "mahalanobis", "eff_n"),
    y_var = c("ssZ", "mahalanobis"),
    pattern_col = NULL,
    x_cutoff = NULL,
    facet_scales = "free",
    top_n_labels = 5,
    label_rank_var = NULL,
    point_alpha = 0.5,
    xlabel = NULL,
    ylabel = NULL
) {
  
  dt <- as.data.table(copy(dt))
  x_var <- match.arg(x_var)
  y_var <- match.arg(y_var)
  
  if (is.null(pattern_col)) {
    pattern_col <- if ("pattern" %in% names(dt)) {"pattern"
      } else if ("pattern_flag" %in% names(dt)) {"pattern_flag"
    } else {stop("Could not find pattern or pattern_flag")
    }
  }
  
  required_cols <- c(pattern_col, x_var, y_var, "group_id", 
                     "local_dominant_sign_meaning","gene")
  
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  if (is.null(label_rank_var)) label_rank_var <- y_var
  if (!(label_rank_var %in% names(dt))) stop("label_rank_var not found in dt")
  
  plot_dt <- dt[get(pattern_col) != "none"]
  
  if (nrow(plot_dt) == 0L) {
    stop("No rows left after filtering out 'none'")
  }
  
  plot_dt <- plot_dt[
    is.finite(get(x_var)) &
      is.finite(get(y_var)) &
      get(x_var) > 0 &
      get(y_var) > 0
  ]
  
  if (nrow(plot_dt) == 0L) {
    stop("No rows left after removing non-positive or non-finite x/y values")
  }
  
  label_dt <- plot_dt[!is.na(get(label_rank_var))]
  
  if (nrow(label_dt) > 0L && top_n_labels > 0) {
    label_dt <- label_dt[
      order(-get(label_rank_var)),
      head(.SD, top_n_labels),
      by = .(group_id, pattern, local_dominant_sign_meaning)
    ]
  }
  
  p <- ggplot(
    plot_dt,
    aes(
      x = .data[[x_var]],
      y = .data[[y_var]],
      color = .data[[pattern_col]],
      shape = .data[["local_dominant_sign_meaning"]]
    )
  ) +
    geom_point(alpha = point_alpha) +
    scale_x_log10() +
    scale_y_log10() +
    facet_wrap(
      stats::as.formula("group_id ~ local_dominant_sign_meaning"),
      scales = facet_scales
    ) +
    theme_bw() +
    labs(
      x = ifelse(!is.null(xlabel), xlabel, paste0("log10(", x_var, ")")),
      y = ifelse(!is.null(ylabel), ylabel, paste0("log10(", y_var, ")")),
      color = pattern_col,
      shape = "Dominant Sign"
    ) +
    theme(legend.position = "top")
  
  if (!is.null(x_cutoff)) {
    p <- p + geom_vline(xintercept = x_cutoff, linetype = "dashed")
  }
  
  if (nrow(label_dt) > 0L && top_n_labels > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = label_dt,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE
    )
  }
  
  p
}
# ------------------------------------------------------------
# Make Heatmaps
# ------------------------------------------------------------
make_class_heatmap <- function(
    tout,
    group_name,
    value_col = c("X", "Zres"),
    gene_classes = c(
      "Local Positive",
      "Local Negative",
      "Global Mixed",
      "Global Positive",
      "Global Negative"
    ),
    top_n_genes = 100,
    rank_var = "mahalanobis",
    max_rows = 100,
    rows_per_cluster = 10,
    row_km_repeats = 50,
    transform = c("signed_log10", "none", "log10"),
    seed = 123
) {
  
  value_col <- match.arg(value_col)
  transform <- match.arg(transform)
  
  tout <- data.table::as.data.table(data.table::copy(tout))
  
  required_cols <- c(
    "group_id",
    "gene",
    "contrast",
    "gene_class",
    value_col
  )
  
  if (!is.null(rank_var)) {
    required_cols <- c(required_cols, rank_var)
  }
  
  missing_cols <- setdiff(required_cols, names(tout))
  if (length(missing_cols) > 0L) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  signed_log10 <- function(x) {
    sign(x) * log10(abs(x) + 1)
  }
  
  transform_values <- function(x) {
    if (transform == "none") {
      return(x)
    }
    
    if (transform == "signed_log10") {
      return(signed_log10(x))
    }
    
    if (transform == "log10") {
      if (any(x <= 0, na.rm = TRUE)) {
        stop("log10 transform requires all values to be positive. Use transform = 'signed_log10' instead.")
      }
      return(log10(x))
    }
  }
  
  tout[, gene_class := factor(gene_class, levels = gene_classes)]
  
  group_dt <- tout[
    group_id == group_name &
      gene_class %in% gene_classes &
      !is.na(gene) &
      !is.na(contrast) &
      is.finite(get(value_col))
  ]
  
  if (nrow(group_dt) == 0L) {
    message("Skipping group with no heatmap data: ", group_name)
    return(NULL)
  }
  
  out <- list()
  
  for (class_name in gene_classes) {
    
    class_dt <- group_dt[gene_class == class_name]
    
    if (nrow(class_dt) == 0L) {
      message("Skipping empty class: ", group_name, " | ", class_name)
      next
    }
    
    if (!is.null(top_n_genes) && !is.null(rank_var)) {
      rank_dt <- class_dt[
        !is.na(get(rank_var)),
        .(rank_value = max(get(rank_var), na.rm = TRUE)),
        by = gene
      ]
      
      if (nrow(rank_dt) > 0L) {
        data.table::setorder(rank_dt, -rank_value)
        keep_genes <- head(rank_dt$gene, top_n_genes)
        class_dt <- class_dt[gene %in% keep_genes]
      }
    }
    
    if (nrow(class_dt) == 0L) {
      message("Skipping class after gene filtering: ", group_name, " | ", class_name)
      next
    }
    
    sub <- class_dt[
      ,
      .(value = mean(get(value_col), na.rm = TRUE)),
      by = .(gene, contrast)
    ]
    
    sub[, value := transform_values(value)]
    
    wide <- data.table::dcast(
      sub,
      gene ~ contrast,
      value.var = "value"
    )
    
    if (nrow(wide) == 0L || ncol(wide) < 3L) {
      message("Skipping class with insufficient matrix dimensions: ", group_name, " | ", class_name)
      next
    }
    
    mat <- as.matrix(wide[, -1, with = FALSE])
    rownames(mat) <- wide$gene
    
    mat <- mat[stats::complete.cases(mat), , drop = FALSE]
    
    if (nrow(mat) < 2L || ncol(mat) < 2L) {
      message("Skipping class with too few complete rows/columns: ", group_name, " | ", class_name)
      next
    }
    
    if (nrow(mat) > max_rows) {
      if (!is.null(rank_var)) {
        rank_dt <- class_dt[
          gene %in% rownames(mat) & !is.na(get(rank_var)),
          .(rank_value = max(get(rank_var), na.rm = TRUE)),
          by = gene
        ]
        
        if (nrow(rank_dt) > 0L) {
          data.table::setorder(rank_dt, -rank_value)
          keep_genes <- head(rank_dt$gene, max_rows)
        } else {
          rv <- apply(mat, 1, stats::var, na.rm = TRUE)
          keep_genes <- names(sort(rv, decreasing = TRUE))[seq_len(max_rows)]
        }
      } else {
        rv <- apply(mat, 1, stats::var, na.rm = TRUE)
        keep_genes <- names(sort(rv, decreasing = TRUE))[seq_len(max_rows)]
      }
      
      mat <- mat[intersect(keep_genes, rownames(mat)), , drop = FALSE]
    }
    
    lim <- max(abs(mat), na.rm = TRUE)
    
    if (!is.finite(lim) || lim == 0) {
      message("Skipping class with invalid heatmap range: ", group_name, " | ", class_name)
      next
    }
    
    col_fun <- circlize::colorRamp2(
      c(-lim, 0, lim),
      c("navy", "white", "firebrick")
    )
    
    
    choose_row_km <- function(mat, rows_per_cluster = 10, min_rows_for_km = 20) {
      n <- nrow(mat)
      
      if (n < min_rows_for_km) {
        return(NULL)
      }
      
      distinct_n <- nrow(unique(as.data.frame(mat)))
      
      if (distinct_n < 4L) {
        return(NULL)
      }
      
      k <- ceiling(n / rows_per_cluster)
      k <- min(k, floor(n / 5), distinct_n - 1L)
      k <- max(k, 2L)
      
      if (!is.finite(k) || k < 2L) {
        return(NULL)
      }
      
      k
    }
    
    row_km_value <- choose_row_km(
      mat = mat,
      rows_per_cluster = rows_per_cluster,
      min_rows_for_km = 20
    )
    
    set.seed(seed)
    
    ht <- ComplexHeatmap::Heatmap(
      mat,
      name = value_col,
      col = col_fun,
      row_km = row_km_value,
      row_km_repeats = row_km_repeats,
      cluster_rows = TRUE,
      cluster_columns = FALSE,      
      show_column_dend = FALSE,
      border = TRUE,
      show_row_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 6),
      column_names_gp = grid::gpar(fontsize = 6),
      column_names_rot = 45,
      column_title = paste(group_name, "|", class_name, "|", value_col, "|", transform),
      row_title = NULL,
      show_row_dend = TRUE
    )
    
    ht_drawn <- ComplexHeatmap::draw(ht)
    
    row_clusters <- ComplexHeatmap::row_order(ht_drawn)
    
    if (!is.list(row_clusters)) {
      row_clusters <- list(row_clusters)
    }
    
    genes_per_cluster <- lapply(row_clusters, function(idx) {
      rownames(mat)[idx]
    })
    
    cluster_df <- data.table::rbindlist(
      lapply(seq_along(genes_per_cluster), function(i) {
        data.table::data.table(
          group_id = group_name,
          gene_class = class_name,
          km_cluster = paste0("cluster_", i),
          gene = genes_per_cluster[[i]]
        )
      }),
      use.names = TRUE,
      fill = TRUE
    )
    
    out[[class_name]] <- list(
      # heatmap = ht,
      heatmap_drawn = ht_drawn,
      cluster_df = cluster_df
      # matrix = mat
    )
  }
  
  if (length(out) == 0L) {
    return(NULL)
  }
  
  out
}

# ------------------------------------------------------------
# Make Summary Plots
# ------------------------------------------------------------

plot_summary_bar <- function(group_class_counts) {
  library(data.table)
  library(ggplot2)
  library(gt)
  
  plot_dt <- copy(group_class_counts)
  
  plot_dt[, gene_class := factor(
    gene_class,
    levels = c(
      "None",
      "Local Negative", "Local Positive",
      "Global Negative", "Global Mixed", "Global Positive"
    )
  )]
  
  # proportions
  plot_dt[gene_class != 'None', prop := N / sum(N), by = group_id]
  
  class_cols <- c(
    "None" = "grey80",
    "Local Negative" = "#4575b4",
    "Local Positive" = "#d73027",
    "Global Negative" = "#313695",
    "Global Mixed" = "#969696",
    "Global Positive" = "#a50026"
  )
  
  # proportion plot (no scales package)
  prop_plot <- ggplot(plot_dt[gene_class != 'None'], aes(
    x = group_id, y = prop, fill = gene_class)) +
    geom_col(
      color = "black",
      width = 0.95   # increase width → less gap between bars
    ) +
    geom_text(
      aes(label = ifelse(prop > 0.02, sprintf("%.1f%%", 100 * prop), "")),
      position = position_stack(vjust = 0.5),
      size = 3
    ) +
    scale_y_continuous(
      limits = c(0, 1),   # force 0–100%
      expand = c(0, 0),   # remove padding top/bottom
      labels = function(x) paste0(round(100 * x), "%")
    ) +
    scale_x_discrete(
      expand = c(0.01, 0)  # remove side padding between groups
    ) +
    scale_fill_manual(values = class_cols, drop = FALSE) +
    labs(
      title = "Gene Class Proportions by Group",
      x = "Group",
      y = "Proportion",
      fill = "Gene Class"
    ) +
    coord_flip() +
    # theme_bw() +
    theme(
      axis.text = element_text(color = 'black', face = 'bold'),
      axis.text.x = element_text(angle = 45, hjust = 1, color = 'black'),
      panel.grid.major.x = element_blank(),  # cleaner look
      panel.grid.minor.x = element_blank()
    )
  #-----------------------------------------
  # keep only displayed classes
  show_cols <- c(
    "Local Negative", "Local Positive",
    "Global Negative", "Global Mixed", "Global Positive"
  )
  # wide counts table
  counts_dt <- dcast(
    plot_dt[gene_class != 'None'],
    group_id ~ gene_class,
    value.var = "N",
    fill = 0
  )
  
  setcolorder(counts_dt, c("group_id", show_cols))
  
  # compute totals
  counts_dt[, Total := rowSums(.SD), .SDcols = show_cols]
  grand_total <- counts_dt[, c(
    list(group_id = "Total"),
    lapply(.SD, sum)
  ), .SDcols = c(show_cols, "Total")]
  
  # bind on bottom
  counts_dt2 <- rbindlist(
    list(counts_dt[, c("group_id", show_cols, "Total"), with = FALSE], grand_total),
    use.names = TRUE,
    fill = TRUE
  )
  setnames(counts_dt2, 'group_id', 'Group')
  
  # percentages within each row total
  for (col in c(show_cols, "Total")) {
    if (col != "Total") {
      counts_dt2[, (col) := sprintf(
        "%d<br>(%.1f%%)</span>",
        get(col),
        fifelse(Total > 0, 100 * get(col) / Total, 0)
      )]
    } else {
      counts_dt2[, (col) := as.character(get(col))]
    }
  }
  
  # convert to gt table
  counts_table <- counts_dt2 |>
    as.data.frame() |>
    gt() |>
    tab_header(
      title = md("**Gene Class Distribution**")
    ) |>
    cols_label(
      `Local Negative` = md("Local<br>Negative"),
      `Local Positive` = md("Local<br>Positive"),
      `Global Negative` = md("Global<br>Negative"),
      `Global Mixed`    = md("Global<br>Mixed"),
      `Global Positive` = md("Global<br>Positive"),
      Total             = md("Row<br>Total")
    ) |>
    fmt_markdown(columns = everything()) |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = list(
        cells_stub(rows = Group %in% counts_dt2$Group)
      )
    ) |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = list(
        cells_body(rows = nrow(counts_dt2)),
        cells_stub(rows = nrow(counts_dt2))
      )
    )
  
  list(
    prop_plot = prop_plot,
    counts_table_graph = counts_table,
    counts_table = counts_dt
  )
}


# ------------------------------------------------------------
# 7. t-SNE per group
# 8. UMAP per group
# Updated: general color column, top gene labels, cluster labels
# ------------------------------------------------------------

plot_tsne_one_group <- function(dt, annot_dt, group_name, value_col = "X",
                                keep_mode = c("interesting_only", "rowwise_only", "all"),
                                seed = 1, perplexity = NULL,
                                color_col = "plot_class",
                                label_top_genes = FALSE,
                                top_n_genes = 10,
                                gene_rank_col = "mahalanobis",
                                cluster_col = "hdbscan_cluster") {
  keep_mode <- match.arg(keep_mode)
  
  mat <- prep_embedding_matrix(dt, group_name, value_col = value_col, keep_mode = keep_mode)
  if (is.null(mat)) return(NULL)
  
  n <- nrow(mat)
  
  if (is.null(perplexity)) {
    perplexity <- max(2, min(30, floor((n - 1) / 3)))
  }
  if (perplexity >= n) {
    perplexity <- max(2, floor((n - 1) / 3))
  }
  if (perplexity < 2 || n < 4) {
    return(NULL)
  }
  
  set.seed(seed)
  tsne_out <- Rtsne(
    mat,
    pca = TRUE,
    check_duplicates = FALSE,
    perplexity = perplexity
  )
  
  plot_dt <- data.table(
    gene = rownames(mat),
    dim1 = tsne_out$Y[, 1],
    dim2 = tsne_out$Y[, 2]
  )
  
  annot_sub <- annot_dt[group_id == group_name]
  plot_dt <- annot_sub[plot_dt, on = "gene"]
  
  p <- ggplot(plot_dt, aes(x = dim1, y = dim2, color = as.factor(get(color_col)))) +
    geom_point(alpha = 0.8, size = 1.6) +
    theme_bw() +
    labs(
      title = paste("t-SNE:", group_name, "|", value_col, "|", keep_mode, "| color =", color_col),
      x = "tSNE-1",
      y = "tSNE-2",
      color = color_col
    )
  
  if (label_top_genes) {
    top_genes <- get_top_genes_for_group(
      dt = dt,
      group_name = group_name,
      top_n = top_n_genes,
      rank_col = gene_rank_col,
      keep_mode = keep_mode
    )
    
    label_dt <- plot_dt[gene %in% top_genes]
    
    if (nrow(label_dt) > 0L) {
      p <- p + ggrepel::geom_text_repel(
        data = label_dt,
        aes(label = gene),
        size = 3,
        max.overlaps = Inf,
        show.legend = FALSE
      )
    }
  }
  
  p
}


plot_umap_one_group <- function(dt, annot_dt, group_name, value_col = "X",
                                keep_mode = c("interesting_only", "rowwise_only", "all"),
                                seed = 1, n_neighbors = 15,
                                color_col = "plot_class",
                                label_top_genes = FALSE,
                                top_n_genes = 10,
                                gene_rank_col = "mahalanobis",
                                cluster_col = "hdbscan_cluster") {
  keep_mode <- match.arg(keep_mode)
  
  mat <- prep_embedding_matrix(dt, group_name, value_col = value_col, keep_mode = keep_mode)
  if (is.null(mat)) return(NULL)
  
  n <- nrow(mat)
  n_neighbors <- min(n_neighbors, max(2, n - 1))
  if (n < 3) return(NULL)
  
  set.seed(seed)
  umap_out <- uwot::umap(
    mat,
    n_components = 2,
    n_neighbors = n_neighbors,
    verbose = FALSE
  )
  
  plot_dt <- data.table(
    gene = rownames(mat),
    dim1 = umap_out[, 1],
    dim2 = umap_out[, 2]
  )
  
  annot_sub <- annot_dt[group_id == group_name]
  plot_dt <- annot_sub[plot_dt, on = "gene"]
  
  p <- ggplot(plot_dt, aes(x = dim1, y = dim2, color = as.factor(get(color_col)))) +
    geom_point(alpha = 0.8, size = 1.6) +
    theme_bw() +
    labs(
      title = paste("UMAP:", group_name, "|", value_col, "|", keep_mode, "| color =", color_col),
      x = "UMAP-1",
      y = "UMAP-2",
      color = color_col
    )
  
  if (label_top_genes) {
    top_genes <- get_top_genes_for_group(
      dt = dt,
      group_name = group_name,
      top_n = top_n_genes,
      rank_col = gene_rank_col,
      keep_mode = keep_mode
    )
    
    label_dt <- plot_dt[gene %in% top_genes]
    
    if (nrow(label_dt) > 0L) {
      p <- p + ggrepel::geom_text_repel(
        data = label_dt,
        aes(label = gene),
        size = 3, 
        max.overlaps = Inf,
        show.legend = FALSE
      )
    }
  }

  p
}

# ------------------------------------------------------------
# 10. Pairs plots 2
# Updated: supports cellwise or rowwise, adds counts + top gene labels
# ------------------------------------------------------------




add_distortion_jitter <- function(
    x_raw,
    x_plot,
    affected = NULL,
    distance = NULL,
    amount = 0.02,
    direction = c("outward", "random"),
    seed = NULL
) {
  direction <- match.arg(direction)
  
  if (!is.null(seed)) {
    old_seed <- .Random.seed
    on.exit({
      if (exists("old_seed", inherits = FALSE)) .Random.seed <<- old_seed
    }, add = TRUE)
    set.seed(seed)
  }
  
  x_out <- x_plot
  
  ok <- is.finite(x_raw) & is.finite(x_plot)
  
  if (is.null(affected)) {
    affected <- ok & abs(x_raw - x_plot) > sqrt(.Machine$double.eps)
  } else {
    affected <- affected & ok
  }
  
  if (!any(affected)) {
    return(x_out)
  }
  
  if (is.null(distance)) {
    distance <- abs(x_raw - x_plot)
  }
  
  distance[!is.finite(distance)] <- NA_real_
  
  max_dist <- max(distance[affected], na.rm = TRUE)
  
  if (!is.finite(max_dist) || max_dist == 0) {
    return(x_out)
  }
  
  rel_dist <- distance / max_dist
  
  plot_range <- diff(range(x_plot[ok], na.rm = TRUE))
  
  if (!is.finite(plot_range) || plot_range == 0) {
    plot_range <- 1
  }
  
  jitter_size <- amount * plot_range * rel_dist
  
  if (direction == "random") {
    jitter <- runif(length(x_raw), min = -1, max = 1) * jitter_size
  } else {
    center <- median(x_raw[ok], na.rm = TRUE)
    jitter <- sign(x_raw - center) * jitter_size
  }
  
  x_out[affected] <- x_out[affected] + jitter[affected]
  
  x_out
}

transform_for_plot <- function(
    x,
    option = c("none", "winsorize", "rank_quantile", "mad_limit", 
               "signed_log", "asinh", "robust_z"),
    probs = c(0.005, 0.995),
    k = 6,
    base = 10,
    scale = c("mad", "sd", "none"),
    add_jitter = FALSE,
    jitter_amount = 0.02,
    jitter_direction = c("outward", "random"),
    seed = NULL
) {
  option <- match.arg(option)
  scale <- match.arg(scale)
  jitter_direction <- match.arg(jitter_direction)
  
  x_raw <- x
  x_out <- x
  
  get_scale <- function(x, scale) {
    if (scale == "none") return(1)
    
    s <- switch(
      scale,
      mad = mad(x, na.rm = TRUE),
      sd  = sd(x, na.rm = TRUE)
    )
    
    if (is.na(s) || s == 0) {
      s <- sd(x, na.rm = TRUE)
    }
    
    if (is.na(s) || s == 0) {
      s <- 1
    }
    
    s
  }
  
  if (option == "none") {
    return(x_out)
  }
  
  if (option == "winsorize") {
    q <- quantile(x, probs = probs, na.rm = TRUE)
    lower <- q[1]
    upper <- q[2]
    
    x_out <- pmin(pmax(x, lower), upper)
    
    if (add_jitter) {
      affected <- is.finite(x) & (x < lower | x > upper)
      distance <- ifelse(
        x < lower,
        lower - x,
        ifelse(x > upper, x - upper, 0)
      )
      
      x_out <- add_distortion_jitter(
        x_raw = x_raw,
        x_plot = x_out,
        affected = affected,
        distance = distance,
        amount = jitter_amount,
        direction = jitter_direction,
        seed = seed
      )
    }
    
    return(x_out)
  }
  
  if (option == "mad_limit") {
    med <- median(x, na.rm = TRUE)
    s <- mad(x, na.rm = TRUE)
    
    if (is.na(s) || s == 0) {
      s <- sd(x, na.rm = TRUE)
    }
    
    if (is.na(s) || s == 0) {
      return(x_out)
    }
    
    lower <- med - k * s
    upper <- med + k * s
    
    x_out <- pmin(pmax(x, lower), upper)
    
    if (add_jitter) {
      affected <- is.finite(x) & (x < lower | x > upper)
      distance <- ifelse(
        x < lower,
        lower - x,
        ifelse(x > upper, x - upper, 0)
      )
      
      x_out <- add_distortion_jitter(
        x_raw = x_raw,
        x_plot = x_out,
        affected = affected,
        distance = distance,
        amount = jitter_amount,
        direction = jitter_direction,
        seed = seed
      )
    }
    
    return(x_out)
  }
  
  if (option == "rank_quantile") {
    ok <- !is.na(x)
    x_out[ok] <- (rank(x[ok], ties.method = "average") - 0.5) / sum(ok)
    
    if (add_jitter) {
      affected <- is.finite(x)
      distance <- abs(x - median(x, na.rm = TRUE))
      
      x_out <- add_distortion_jitter(
        x_raw = x_raw,
        x_plot = x_out,
        affected = affected,
        distance = distance,
        amount = jitter_amount,
        direction = jitter_direction,
        seed = seed
      )
    }
    
    return(x_out)
  }
  
  if (option == "signed_log") {
    x_out <- sign(x) * log1p(abs(x)) / log(base)
    
    if (add_jitter) {
      affected <- is.finite(x) & abs(x) > 1
      distance <- abs(x)
      
      x_out <- add_distortion_jitter(
        x_raw = x_raw,
        x_plot = x_out,
        affected = affected,
        distance = distance,
        amount = jitter_amount,
        direction = jitter_direction,
        seed = seed
      )
    }
    
    return(x_out)
  }
  
  if (option == "asinh") {
    s <- get_scale(x, scale = scale)
    x_out <- asinh(x / s)
    
    if (add_jitter) {
      affected <- is.finite(x) & abs(x) > s
      distance <- abs(x)
      
      x_out <- add_distortion_jitter(
        x_raw = x_raw,
        x_plot = x_out,
        affected = affected,
        distance = distance,
        amount = jitter_amount,
        direction = jitter_direction,
        seed = seed
      )
    }
    
    return(x_out)
  }
  
  if (option == "robust_z") {
    med <- median(x, na.rm = TRUE)
    s <- get_scale(x, scale = "mad")
    x_out <- (x - med) / s
    
    if (add_jitter) {
      affected <- is.finite(x) & abs(x_out) > k
      distance <- abs(x - med)
      
      x_out <- add_distortion_jitter(
        x_raw = x_raw,
        x_plot = x_out,
        affected = affected,
        distance = distance,
        amount = jitter_amount,
        direction = jitter_direction,
        seed = seed
      )
    }
    
    return(x_out)
  }
}




plot_pairs <- function(
    dt,
    group_name,
    value_col = "Zres",
    keep_mode = c("interesting_only", "rowwise_only", "all"),
    main_title = NULL,
    pch = 16, cex = 0.5,
    label_top_genes = TRUE,
    top_n_genes = 10,
    gene_rank_col = "mahalanobis",
    transform_option = "winsorize",
    max_vars = 4
    ) {
  
  keep_mode <- match.arg(keep_mode)

  sub <- dt[group_id == group_name]
  
  if (nrow(sub) == 0L) return(invisible(NULL))
  

  # --- Full matrix: all genes ---
  wide_dt <- unique(
    sub[, .(gene, contrast, value = get(value_col))]
  )
  
  wide <- data.table::dcast(
    wide_dt,
    gene ~ contrast,
    value.var = "value"
  )
  
  if (is.null(wide) || ncol(wide) < 3L) {
    return(invisible(NULL))
  }
  
  Xmat <- as.matrix(wide[, -1, with = FALSE])
  rownames(Xmat) <- wide$gene
  
  keep_complete <- complete.cases(Xmat)
  Xmat <- Xmat[keep_complete, , drop = FALSE]
  
  if (nrow(Xmat) < 3L) {
    return(invisible(NULL))
  }
  
  # --- Subset matrix: keep_mode genes ---
  # --- used only for highlighting ---
  
  if (keep_mode == "rowwise_only") {
    keep_genes0 <- unique(sub[top_gene_flag == TRUE, gene])
  } else if (keep_mode == "interesting_only") {
    keep_genes0 <- unique(sub[top_gene_flag == TRUE | pattern != "none", gene])
  } else {
    keep_genes0 <- unique(sub$gene)
  }
  
  sub_keep <- sub[gene %in% keep_genes0]
  
  sub_wide_dt <- unique(
    sub_keep[, .(gene, contrast, value = get(value_col))]
  )
  
  sub_wide <- data.table::dcast(
    sub_wide_dt,
    gene ~ contrast,
    value.var = "value"
  )
  
  if (is.null(sub_wide) || ncol(sub_wide) < 3L) {
    keep_genes <- character(0)
  } else {
    sub_mat <- as.matrix(sub_wide[, -1, with = FALSE])
    rownames(sub_mat) <- sub_wide$gene
    
    sub_mat <- sub_mat[complete.cases(sub_mat), , drop = FALSE]
    
    keep_genes <- rownames(sub_mat)
  }
  
  flag_vec <- rownames(Xmat) %in% keep_genes
  
  Xmat <- apply(
    Xmat,
    2,
    transform_for_plot,
    option = transform_option,
    probs = c(0.0001, 0.9995),
    add_jitter = TRUE,
    jitter_amount = 0.5,
    jitter_direction = "outward"
  )
  
  # --- Top gene labels ---
  top_genes <- character(0)
  if (label_top_genes) {
    top_genes <- get_top_genes_for_group(
      dt = dt,
      group_name = group_name,
      top_n = top_n_genes,
      rank_col = gene_rank_col,
      keep_mode = keep_mode
    )
  }
  
  label_mask <- rownames(Xmat) %in% top_genes
  
  if (is.null(main_title)) {
    main_title <- paste("Pairs:", group_name, "|", value_col)
  }
  
  # --- Panel function ---
  panel_fun <- function(x, y, ...) {
    
    panel_scatter_flag_labels(
      x, y,
      flag = flag_vec,
      labels = rownames(Xmat),
      label_mask = label_mask,
      pch = pch, cex = cex
    )
    
    # quadrant lines
    abline(h = 0, v = 0, lty = 2, col = "grey80")
    
    # threshold lines
    # thr <- 2.57
    # abline(h = c(-thr, thr), v = c(-thr, thr),
    #        lty = 2, col = "red")
    
    # quadrant counts
    q1 <- sum(x > 0 & y > 0)
    q2 <- sum(x < 0 & y > 0)
    q3 <- sum(x < 0 & y < 0)
    q4 <- sum(x > 0 & y < 0)
    
    usr <- par("usr")

    text(usr[2], usr[4], q1, adj = c(1,1), cex = 1.2)
    text(usr[1], usr[4], q2, adj = c(0,1), cex = 1.2)
    text(usr[1], usr[3], q3, adj = c(0,0), cex = 1.2)
    text(usr[2], usr[3], q4, adj = c(1,0), cex = 1.2)
  }
  
  # --- Plot ---
  
  n_contrasts <- ncol(Xmat)
  
  if (n_contrasts < 2L) {
    return(invisible(NULL))
  }
  
  if (n_contrasts <= max_vars) {
    contrast_windows <- list(seq_len(n_contrasts))
  } else {
    contrast_windows <- lapply(
      seq_len(n_contrasts - max_vars + 1L),
      function(i) i:(i + max_vars - 1L)
    )
  }
  
  for (w in seq_along(contrast_windows)) {
    
    cols <- contrast_windows[[w]]
    Xsub <- Xmat[, cols, drop = FALSE]
    
    plot_title <- main_title
    
    if (is.null(plot_title)) {
      plot_title <- paste(
        "Pairs:",
        group_name,
        "|",
        value_col,
        "| contrasts",
        paste(range(cols), collapse = "-")
      )
    } else if (length(contrast_windows) > 1L) {
      plot_title <- paste0(
        plot_title,
        " | contrasts ",
        paste(range(cols), collapse = "-")
      )
    }
    
    pairs(
      Xsub,
      lower.panel = panel_fun,
      upper.panel = NULL,
      diag.panel = NULL,
      gap = 0.5,
      main = plot_title
    )
  }
  
}

# ------------------------------------------------------------
# 11. Save embeddings side by side to one PDF
# Updated: class-colored + cluster-colored + cluster labels
# ------------------------------------------------------------

save_embedding_comparison_pdf <- function(filename, dt, annot_dt,
                                          value_col = "X",
                                          keep_mode = c("interesting_only", "rowwise_only", "all"),
                                          seed = 1,
                                          include_summary = TRUE,
                                          label_top_genes = TRUE,
                                          top_n_genes = 10,
                                          gene_rank_col = "mahalanobis",
                                          include_cluster_plots = TRUE,
                                          cluster_col = "hdbscan_cluster",
                                          width = 14, height = 6) {
  keep_mode <- match.arg(keep_mode)
  group_ids <- unique(dt$group_id)
  
  pdf(filename, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  
  if (include_summary) {
    print(summary_bar)
  }
  
  for (g in group_ids) {
    # color by plot class
    p_tsne_class <- plot_tsne_one_group(
      dt = dt, annot_dt = annot_dt, group_name = g,
      value_col = value_col, keep_mode = keep_mode, seed = seed,
      color_col = "plot_class",
      label_top_genes = label_top_genes,
      top_n_genes = top_n_genes,
      gene_rank_col = gene_rank_col
    )
    
    p_umap_class <- plot_umap_one_group(
      dt = dt, annot_dt = annot_dt, group_name = g,
      value_col = value_col, keep_mode = keep_mode, seed = seed,
      color_col = "plot_class",
      label_top_genes = label_top_genes,
      top_n_genes = top_n_genes,
      gene_rank_col = gene_rank_col
      )
    
    if (!is.null(p_tsne_class) && !is.null(p_umap_class)) {
      gridExtra::grid.arrange(p_tsne_class, p_umap_class, ncol = 2)
    } else if (!is.null(p_tsne_class)) {
      print(p_tsne_class)
    } else if (!is.null(p_umap_class)) {
      print(p_umap_class)
    }
    
    # color by cluster
    if (include_cluster_plots && cluster_col %in% colnames(annot_dt)) {
      p_tsne_cluster <- plot_tsne_one_group(
        dt = dt, annot_dt = annot_dt, group_name = g,
        value_col = value_col, keep_mode = keep_mode, seed = seed,
        color_col = cluster_col,
        label_top_genes = label_top_genes,
        top_n_genes = top_n_genes,
        gene_rank_col = gene_rank_col,
        cluster_col = cluster_col
      )
      
      p_umap_cluster <- plot_umap_one_group(
        dt = dt, annot_dt = annot_dt, group_name = g,
        value_col = value_col, keep_mode = keep_mode, seed = seed,
        color_col = cluster_col,
        label_top_genes = label_top_genes,
        top_n_genes = top_n_genes,
        gene_rank_col = gene_rank_col,
        cluster_col = cluster_col
      )
      
      if (!is.null(p_tsne_cluster) && !is.null(p_umap_cluster)) {
        gridExtra::grid.arrange(p_tsne_cluster, p_umap_cluster, ncol = 2)
      } else if (!is.null(p_tsne_cluster)) {
        print(p_tsne_cluster)
      } else if (!is.null(p_umap_cluster)) {
        print(p_umap_cluster)
      }
    }
  }
}

# ------------------------------------------------------------
# 3. Make Plots
# ------------------------------------------------------------

# change out_dir temporarily while troubleshooting
out_dir = "C:/Users/dobodo/Downloads"

make_plot_pipeline = function(out.list, out_dir){


# per gene tables
tout = out.list$main_summary
gout = out.list$gene_summary
setDT(tout)
setDT(gout)

# diagnostic plots
groups <- unique(gout$group_id)
group_pairs <- split(
  groups,
  ceiling(seq_along(groups) / 2)
)

# make heatmaps
gclusts = list()

pdf(
  file.path(out_dir, "diagnostic_plots.pdf"),
  height = 11,
  width = 15
)

for (grp in group_pairs) {
  
  plot_dt <- gout[group_id %in% grp]
  
  group_class_counts <- unique(
    plot_dt[, .N, by = .(group_id, gene_class)])
  res = plot_summary_bar(group_class_counts)
  print(res$prop_plot)
  
  
  print(
    plot_global_local_summary(
      dt = plot_dt,
      x_var = "p_global",
      y_var = "ssZ",
      pattern_col = "pattern",
      top_n_labels = 10,
      label_rank_var = "ssZ"
    )
  )
  
  print(
    plot_global_local_summary(
      dt = plot_dt,
      x_var = "mahalanobis",
      y_var = "ssZ",
      pattern_col = "pattern",
      top_n_labels = 10,
      label_rank_var = "mahalanobis"
    )
  )
  
  print(
    plot_global_local_summary(
      dt = plot_dt,
      x_var = "eff_n",
      y_var = "ssZ",
      pattern_col = "pattern",
      top_n_labels = 10,
      label_rank_var = "ssZ"
    )
  )
  
  # heatmap of top max 100 genes per group and class
  for (g in grp) {
    res <- make_class_heatmap(
      tout = touttout[group_id == g,],
      group_name = g,
      value_col = "Zres",
      top_n_genes = 100,
      rank_var = "mahalanobis",
      max_rows = 100,
      transform = "signed_log10"
    )
    
    if (!is.null(res)) {
      gclusts[[g]] <- res
    }
  }
  
  for (g in grp) {
    
    dt = tout[group_id == g,]
    plot_pairs(
      dt = dt,
      group_name = g,
      keep_mode = "rowwise_only"
      )
     
    # TODO: Add a generic plotting function to counts for only flagged genes; color pair-specific flagged genes if interesting_only; add option for quadrant specific labeling; maybe enable locator gene labeling 
  }
  
  
}

  
# 3. Close the device
dev.off()

print(res$counts_table)

cluster_tbl <- data.table::rbindlist(
  lapply(gclusts, function(group_res) {
    data.table::rbindlist(
      lapply(group_res, `[[`, "cluster_df"),
      fill = TRUE
    )
  }),
  fill = TRUE
)




# save_embedding_comparison_pdf(
#   filename = file.path(out_dir, "embeddings_tsne_umap_interesting_only.pdf"),
#   dt = tout,
#   annot_dt = gene_annot,
#   value_col = "X",
#   keep_mode = "interesting_only",
#   include_summary = TRUE,
#   label_top_genes = TRUE,
#   top_n_genes = 30,
#   gene_rank_col = "mahalanobis",
#   include_cluster_plots = TRUE,
#   cluster_col = "hdbscan_cluster"
#   )
# 
# save_pairs_pdf(
#   filename = file.path(out_dir, "pairs.pdf"),
#   dt = tout,
#   value_col = "X",
#   label_top_genes = TRUE,
#   top_n_genes = 20,
#   gene_rank_col = "mahalanobis"
# )

}


# ------------------------------------------------------------
# 3. Save Plots
# ------------------------------------------------------------
