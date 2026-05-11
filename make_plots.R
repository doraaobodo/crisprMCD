

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
    x_var = c("p_global", "mahalanobis", 'eff_n'),
    y_var = c("ssZ", "mahalanobis"),
    pattern_col = NULL,
    x_cutoff = NULL,
    facet_scales = "free",
    top_n_labels = 5,
    label_rank_var = NULL,
    cap_x = NULL,
    cap_y = NULL,
    point_alpha = 0.5,
    xlabel = NULL,
    ylabel = NULL
) {

  dt <- as.data.table(copy(dt))
  
  x_var <- match.arg(x_var)
  y_var <- match.arg(y_var)
  
  if (is.null(pattern_col)) {
    pattern_col <- if (
      "pattern" %in% names(dt)
      ) "pattern" else if (
        "pattern_flag" %in% names(dt)
        ) "pattern_flag" else stop(
          "Could not find pattern or pattern_flag")
  }
  
  if (!(pattern_col %in% names(dt))) stop("pattern_col not found in dt")
  if (!(x_var %in% names(dt))) stop("x_var not found in dt")
  if (!(y_var %in% names(dt))) stop("y_var not found in dt")
  if (!("group_id" %in% names(dt))) stop("group_id column not found in dt")
  if (!("local_dominant_sign_meaning" %in% names(dt))) stop("local_dominant_sign_meaning column not found in dt")
  if (!("gene" %in% names(dt))) stop("gene column not found in dt")
  if (is.null(label_rank_var)) label_rank_var <- y_var
  if (!(label_rank_var %in% names(dt))) stop("label_rank_var not found in dt")
  
  plot_dt <- dt[get(pattern_col) != "none"]
  if (nrow(plot_dt) == 0L) stop("No rows left after filtering out 'none'")
  
  # cap values for plotting only
  plot_dt[, x_plot := get(x_var)]
  plot_dt[, y_plot := get(y_var)]
  
  if (!is.null(cap_x)) plot_dt[, x_plot := pmin(x_plot, cap_x)]
  if (!is.null(cap_y)) plot_dt[, y_plot := pmin(y_plot, cap_y)]
  
  # choose top genes per facet
  label_dt <- plot_dt[!is.na(get(label_rank_var))]
  label_dt <- plot_dt[, ]
  
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
      x = x_plot,
      y = y_plot,
      color = .data[[pattern_col]],
      shape = .data[["local_dominant_sign_meaning"]]
    )
  ) +
    geom_point(alpha = point_alpha) +
    facet_wrap(
      stats::as.formula("group_id ~ local_dominant_sign_meaning"),
      scales = facet_scales
    ) +
    theme_bw() +
    labs(
      x = ifelse(!is.null(xlabel),xlabel, x_var),
      y = ifelse(!is.null(ylabel),ylabel, y_var),
      color = pattern_col,
      shape = "local_dominant_sign_meaning"
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
make_class_heatmap <- function(tout, group_name, class_name,
                               value_col = c("X", "Zres"),
                               top_n_genes = NULL,
                               rank_var = NULL,
                               max_rows = 200,
                               rows_per_cluster = 10,
                               row_km_repeats = 50,
                               cap_value = NULL,
                               seed = 123,
                               use_row_annotation = TRUE) {
  
  value_col <- match.arg(value_col)
  
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  
  tout <- as.data.table(copy(tout))
  
  # stable class order
  tout[, gene_class := factor(
    gene_class,
    levels = c(
      "None",
      "Local Negative", "Local Positive",
      "Global Negative", "Global Mixed", "Global Positive"
    )
  )]
  
  # genes in this group/class
  keep_genes <- tout[
    group_id == group_name & gene_class == class_name,
    unique(gene)
  ]
  
  # optional filtering to top genes by rank_var
  if (!is.null(top_n_genes) && !is.null(rank_var)) {
    if (!(rank_var %in% names(tout))) {
      stop("rank_var not found in tout")
    }
    
    top_dt <- unique(
      tout[
        group_id == group_name &
          gene_class == class_name &
          !is.na(get(rank_var)),
        .(gene, rank_value = get(rank_var))
      ]
    )
    
    if (nrow(top_dt) > 0L) {
      setorder(top_dt, -rank_value)
      keep_genes <- head(top_dt$gene, top_n_genes)
    }
  }
  
  if (length(keep_genes) == 0L) return(NULL)
  
  # long -> wide
  sub <- unique(tout[
    group_id == group_name & gene %in% keep_genes,
    .(gene, contrast, value = get(value_col))
  ])
  
  if (is.null(cap_value)) {cap_value = quantile(abs(sub$value), 0.99)} 
  sub[, value := pmin(pmax(value, -cap_value), cap_value)]
  if (nrow(sub) == 0L) return(NULL)
  
  wide <- dcast(sub, gene ~ contrast, value.var = "value")
  if (ncol(wide) < 2L) return(NULL)
  
  mat <- as.matrix(wide[, -1])
  rownames(mat) <- wide$gene
  
  # remove incomplete rows
  keep <- complete.cases(mat)
  mat <- mat[keep, , drop = FALSE]
  
  if (nrow(mat) == 0L) return(NULL)
  
  # cap displayed rows for readability
  # if rank_var available, use it; otherwise use row variance
  if (nrow(mat) > max_rows) {
    if (!is.null(rank_var) && rank_var %in% names(tout)) {
      rank_dt <- unique(
        tout[
          group_id == group_name &
            gene_class == class_name &
            gene %in% rownames(mat) &
            !is.na(get(rank_var)),
          .(gene, rank_value = get(rank_var))
        ]
      )
      
      if (nrow(rank_dt) > 0L) {
        setorder(rank_dt, -rank_value)
        keep2 <- intersect(rank_dt$gene, rownames(mat))
        keep2 <- head(keep2, max_rows)
        mat <- mat[keep2, , drop = FALSE]
      } else {
        rv <- apply(mat, 1, var, na.rm = TRUE)
        keep2 <- names(sort(rv, decreasing = TRUE))[seq_len(max_rows)]
        mat <- mat[keep2, , drop = FALSE]
      }
    } else {
      rv <- apply(mat, 1, var, na.rm = TRUE)
      keep2 <- names(sort(rv, decreasing = TRUE))[seq_len(max_rows)]
      mat <- mat[keep2, , drop = FALSE]
    }
  }
  
  if (nrow(mat) == 0L) return(NULL)
  
  # row annotation data aligned to mat rows
  row_annot_dt <- unique(
    tout[
      group_id == group_name &
        gene_class == class_name &
        gene %in% rownames(mat),
      .(gene, hdbscan_cluster)
    ]
  )
  
  row_annot_dt <- row_annot_dt[match(rownames(mat), gene)]
  cluster_chr <- as.character(row_annot_dt$hdbscan_cluster)
  
  # optional replacement for missing annotation labels
  # cluster_chr[is.na(cluster_chr)] <- "NA"
  
  valid_clusters <- sort(unique(cluster_chr[!is.na(cluster_chr)]))
  
  row_ha <- NULL
  if (use_row_annotation && length(valid_clusters) > 0L) {
    cluster_cols <- structure(
      grDevices::hcl.colors(length(valid_clusters), "Dark 3"),
      names = valid_clusters
    )
    
    row_ha <- ComplexHeatmap::rowAnnotation(
      cluster = cluster_chr,
      col = list(cluster = cluster_cols),
      show_annotation_name = TRUE,
      annotation_name_side = "top"
    )
  }
  
  if (!is.null(cap_value)) {
    lim <- cap_value
  } else {
    rng <- range(mat, na.rm = TRUE)
    lim <- max(abs(rng))
  }
  
  col_fun <- circlize::colorRamp2(
    c(-lim, 0, lim),
    c("navy", "white", "firebrick")
  )
  
  # scale row_km by 10 rows per cluster
  n_clusters <- max(1L, ceiling(nrow(mat) / rows_per_cluster))
  n_clusters <- min(n_clusters, ceiling(max_rows / rows_per_cluster))
  
  set.seed(seed)
  
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = value_col,
    col = col_fun,
    left_annotation = row_ha,
    row_km = 3,
    row_km_repeats = row_km_repeats,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    border = TRUE,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 6),
    column_names_rot = 45,
    column_title = paste(group_name, "|", class_name, "|", value_col)
  )
  
  ht_drawn <- ComplexHeatmap::draw(ht)
  
  # extract genes per kmeans cluster/slice
  row_clusters <- ComplexHeatmap::row_order(ht_drawn)
  
  genes_per_cluster <- lapply(row_clusters, function(idx) {
    rownames(mat)[idx]
  })
  names(genes_per_cluster) <- paste0("cluster_", seq_along(genes_per_cluster))
  
  cluster_df <- data.table::rbindlist(
    lapply(seq_along(genes_per_cluster), function(i) {
      data.table(
        km_cluster = paste0("cluster_", i),
        gene = genes_per_cluster[[i]]
      )
    }),
    use.names = TRUE,
    fill = TRUE
  )
  
  list(
    heatmap_drawn = ht_drawn,
    cluster_df = cluster_df
  )
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
    counts_table = counts_table
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
# 10. Pairs plots
# Updated: supports cellwise or rowwise, adds counts + top gene labels
# ------------------------------------------------------------

panel_scatter_flag <- function(x, y, flag, pch = 16, cex = 0.5,
                               col_false = "grey70", col_true = "red") {
  cols <- ifelse(flag, col_true, col_false)
  cols <- grDevices::adjustcolor(cols, alpha.f = 0.6)
  points(x, y, col = cols, pch = pch, cex = cex)
}

panel_scatter_flag_labels <- function(x, y, flag, labels = NULL, label_mask = NULL,
                                      pch = 16, cex = 0.5,
                                      col_false = "grey70", col_true = "red",
                                      text_cex = 0.55) {
  cols <- ifelse(flag, col_true, col_false)
  cols <- grDevices::adjustcolor(cols, alpha.f = 0.6)
  points(x, y, col = cols, pch = pch, cex = cex)
  
  if (!is.null(labels) && !is.null(label_mask) && any(label_mask)) {
    text(
      x = x[label_mask],
      y = y[label_mask],
      labels = labels[label_mask],
      cex = text_cex,
      pos = 3,
      offset = 0.3
    )
  }
}

plot_pairs <- function(dt, group_name, value_col = "X",
                       keep_mode = c("interesting_only", "rowwise_only", "all"),
                       mode = c("cellwise", "rowwise"),
                       main_title = NULL,
                       pch = 16, cex = 0.5,
                       label_top_genes = TRUE,
                       top_n_genes = 10,
                       gene_rank_col = "mahalanobis",
                       lower_only = TRUE) {
  
  keep_mode <- match.arg(keep_mode)
  mode <- match.arg(mode)
  
  sub <- dt[group_id == group_name]
  
  # Wide matrix for all genes in the group
  wide <- prep_embedding_matrix(dt, group_name, value_col = value_col, keep_mode = "all")
  if (is.null(wide)) return(invisible(NULL))
  
  # Wide matrix for selected subset, used only to define flag coloring
  sub_wide <- prep_embedding_matrix(dt, group_name, value_col = value_col, keep_mode = keep_mode)
  if (is.null(sub_wide)) return(invisible(NULL))
  
  keep_genes <- rownames(sub_wide)
  flag_vec <- rownames(wide) %in% keep_genes
  
  Xmat <- wide
  
  # top genes to label
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
    main_title <- paste("Pairs plot -", mode, "outliers:", group_name, "|", value_col, "|", keep_mode)
  }
  
  if (lower_only) {
    pairs(
      Xmat,
      lower.panel = function(x, y, ...) {
        panel_scatter_flag_labels(
          x, y,
          flag = flag_vec,
          labels = rownames(Xmat),
          label_mask = label_mask,
          pch = pch, cex = cex
        )
      },
      upper.panel = function(x, y, ...) {},
      diag.panel = NULL,
      gap = 0.5,
      main = main_title
    )
  } else {
    pairs(
      Xmat,
      lower.panel = function(x, y, ...) {
        panel_scatter_flag_labels(
          x, y,
          flag = flag_vec,
          labels = rownames(Xmat),
          label_mask = label_mask,
          pch = pch, cex = cex
        )
      },
      upper.panel = function(x, y, ...) {
        panel_scatter_flag(
          x, y,
          flag = flag_vec,
          pch = pch, cex = cex
        )
      },
      diag.panel = NULL,
      gap = 0.5,
      main = main_title
    )
  }
  
  txt <- get_flag_summary_text(dt, group_name, mode = mode)
  mtext(txt, side = 3, line = 0.2, adj = 0, cex = 0.9)
  
  legend(
    "topright",
    legend = c("Not highlighted", "Highlighted subset"),
    col = c(
      grDevices::adjustcolor("grey70", alpha.f = 0.6),
      grDevices::adjustcolor("red", alpha.f = 0.6)
    ),
    pch = 16,
    bty = "n"
  )
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


save_pairs_pdf <- function(filename, dt,
                           value_col = "X",
                           width = 11, height = 11,
                           label_top_genes = TRUE,
                           top_n_genes = 10,
                           gene_rank_col = "mahalanobis") {
  group_ids <- unique(dt$group_id)
  
  pdf(filename, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  
  for (g in group_ids) {
    plot_pairs(
      dt = dt,
      group_name = g,
      value_col = value_col,
      keep_mode = "interesting_only",
      mode = "cellwise",
      label_top_genes = label_top_genes,
      top_n_genes = top_n_genes,
      gene_rank_col = gene_rank_col
    )
    
    plot_pairs(
      dt = dt,
      group_name = g,
      value_col = value_col,
      keep_mode = "rowwise_only",
      mode = "rowwise",
      label_top_genes = label_top_genes,
      top_n_genes = top_n_genes,
      gene_rank_col = gene_rank_col
    )
  }
}

# ------------------------------------------------------------
# 3. Save Plots
# ------------------------------------------------------------

# change out_dir temporarily while troubleshooting
# out_dir = "C:/Users/dobodo/Downloads"

# per gene tables

library(data.table)
library(ggplot2)
library(Rtsne)
library(uwot)
library(GGally)
library(gridExtra)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)

tout = out.list$main_summary
gout = out.list$gene_summary
setDT(tout)
setDT(gout)


group_class_counts <- unique(gout[, .N, by = .(group_id, gene_class)])










p1 <- plot_global_local_summary(
  dt = gout,
  x_var = "p_global",
  y_var = "ssZ",
  pattern_col = "pattern",
  x_cutoff = 0.5,
  cap_y = 250,
  top_n_labels = 10,
  label_rank_var = "ssZ",
  xlabel = 'P(Global)',
  ylabel = 'ssZ'
)

p2 <- plot_global_local_summary(
  dt = gout,
  x_var = "mahalanobis",
  y_var = "ssZ",
  pattern_col = "pattern",
  x_cutoff = 0.6,
  cap_x = 250,
  cap_y = 250,
  top_n_labels = 10,
  label_rank_var = "mahalanobis"
)


p3 <- plot_global_local_summary(
  dt = gout,
  x_var = "eff_n",
  y_var = "ssZ",
  pattern_col = "pattern",
  x_cutoff = 1.5,
  cap_x = NULL,
  cap_y = 250,
  top_n_labels = 10,
  label_rank_var = "ssZ",
  xlabel = 'Effective n',
  ylabel = 'ssZ'
)

# make heatmaps
groups <- unique(tout$group_id)
gclusts = list()


# 1. Open the PDF device
pdf(file.path(out_dir,"temo_diagnostic_plots.pdf"),
    height = 11, width = 15)

# Or use a loop for a list of plots
p3;p2;p1
# lapply(gclusts, `[[`, 1)

for (g in groups) {
  
  gclusts[[g]] <- make_class_heatmap(
    tout = tout,
    group_name = g,
    class_name = "Local Positive",
    value_col = "Zres",
    top_n_genes = 100,
    rank_var = "mahalanobis",
    use_row_annotation=FALSE
  )

  gclusts[[g]] <- make_class_heatmap(
    tout = tout,
    group_name = g,
    class_name = "Local Negative",
    value_col = "Zres",
    top_n_genes = 100,
    rank_var = "mahalanobis",
    use_row_annotation=FALSE
  )
  
  gclusts[[g]] <- make_class_heatmap(
    tout = tout,
    group_name = g,
    class_name = "Global Mixed",
    value_col = "Zres",
    top_n_genes = 100,
    rank_var = "mahalanobis",
    use_row_annotation=FALSE
  )
  
  # gclusts[[g]] <- make_class_heatmap(
  #   tout = tout,
  #   group_name = g,
  #   class_name = "Global Positive",
  #   value_col = "Zres",
  #   top_n_genes = 100,
  #   rank_var = "mahalanobis",
  #   use_row_annotation=FALSE
  # )
  # 
  
  gclusts[[g]] <- make_class_heatmap(
    tout = tout,
    group_name = g,
    class_name = "Global Negative",
    value_col = "Zres",
    top_n_genes = 100,
    rank_var = "mahalanobis",
    use_row_annotation=FALSE
  )

}

res = plot_summary_bar(group_class_counts)
print(res$counts_table)

# 3. Close the device
dev.off()




save_embedding_comparison_pdf(
  filename = file.path(out_dir, "embeddings_tsne_umap_interesting_only.pdf"),
  dt = tout,
  annot_dt = gene_annot,
  value_col = "X",
  keep_mode = "interesting_only",
  include_summary = TRUE,
  label_top_genes = TRUE,
  top_n_genes = 30,
  gene_rank_col = "mahalanobis",
  include_cluster_plots = TRUE,
  cluster_col = "hdbscan_cluster"
  )

save_pairs_pdf(
  filename = file.path(out_dir, "pairs.pdf"),
  dt = tout,
  value_col = "X",
  label_top_genes = TRUE,
  top_n_genes = 20,
  gene_rank_col = "mahalanobis"
)
