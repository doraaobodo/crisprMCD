# ------------------------------------------------------------
# Shared filtering
# ------------------------------------------------------------
filter_interesting_genes <- function(gout, tout, ssZ_quantile = NULL) {
  gout <- as.data.table(copy(gout))
  tout <- as.data.table(copy(tout))
  
  gkeep <- gout[pattern != "none"]
  
  if (!is.null(ssZ_quantile)) {
    qdt <- gkeep[, .(ssZ_cut = quantile(ssZ, ssZ_quantile, na.rm = TRUE)), by = group_id]
    gkeep <- qdt[gkeep, on = "group_id"][ssZ >= ssZ_cut]
  }
  
  keep <- unique(gkeep[, .(group_id, gene)])
  
  tout2 <- keep[tout, on = .(group_id, gene), nomatch = 0]
  gout2 <- keep[gout, on = .(group_id, gene), nomatch = 0]
  
  list(gout = gout2, tout = tout2)
}

# ------------------------------------------------------------
# Plot summary bar plots
# ------------------------------------------------------------


plot_grouping_summary_bar <- function(
    gout,
    grouping_var = c("gene_class", "embedding_cluster"),
    group_name = NULL,
    drop_values = c(NA, "None"),
    title = NULL
) {
  grouping_var <- match.arg(grouping_var)
  
  dt <- data.table::as.data.table(data.table::copy(gout))
  
  if (!is.null(group_name)) {
    dt <- dt[group_id == group_name]
  }
  
  dt <- unique(dt[, .(
    group_id,
    gene,
    grouping = as.character(get(grouping_var))
  )])
  
  dt <- dt[!is.na(grouping)]
  dt <- dt[!(grouping %in% drop_values)]
  
  if (nrow(dt) == 0L) return(NULL)
  
  count_dt <- dt[, .N, by = .(group_id, grouping)]
  
  count_dt[
    ,
    prop := N / sum(N),
    by = group_id
  ]
  
  count_dt[
    ,
    label := fifelse(
      prop >= 0.03,
      sprintf("%d\n%.1f%%", N, 100 * prop),
      ""
    )
  ]
  
  if (grouping_var == "gene_class") {
    class_levels <- c(
      "Local Negative",
      "Local Positive",
      "Global Negative",
      "Global Mixed",
      "Global Positive"
    )
    
    count_dt[, grouping := factor(grouping, levels = class_levels)]
    fill_lab <- "Gene class"
    plot_title <- "Gene class proportions"
  } else {
    count_dt[, grouping := factor(grouping, levels = sort(unique(grouping)))]
    fill_lab <- "Cluster"
    plot_title <- "Cluster proportions"
  }
  
  if (!is.null(title)) {
    plot_title <- title
  }
  
  if (!is.null(group_name)) {
    plot_title <- paste(plot_title, "|", group_name)
  }
  
  ggplot(count_dt, aes(
    x = group_id,
    y = prop,
    fill = grouping
  )) +
    geom_col(color = "black", width = 0.85) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 3
    ) +
    coord_flip() +
    scale_y_continuous(
      limits = c(0, 1),
      expand = c(0, 0),
      labels = function(x) paste0(round(100 * x), "%")
    ) +
    labs(
      title = plot_title,
      x = NULL,
      y = "Proportion of significant genes",
      fill = fill_lab
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# ------------------------------------------------------------
# Add embedding clusters per group
# ------------------------------------------------------------
make_cellmcd_pc <- function(tout, 
                            cellmcd_out, 
                            group_name,
                            n_pcs = 50) {
  dt <- data.table::as.data.table(tout)
  
  xwide <- data.table::dcast(
    dt[group_id == group_name],
    gene ~ contrast,
    value.var = "X"
  )
  
  xmat <- as.matrix(xwide[, -1])
  rownames(xmat) <- xwide$gene
  xmat <- xmat[complete.cases(xmat), , drop = FALSE]
  
  S <- cellmcd_out$S
  mu <- as.numeric(cellmcd_out$mu)
  
  common_cols <- intersect(colnames(xmat), colnames(S))
  xmat <- xmat[, common_cols, drop = FALSE]
  S <- S[common_cols, common_cols, drop = FALSE]
  
  
  eig <- eigen(S)
  
  k <- min(n_pcs, ncol(xmat), ncol(eig$vectors))
  
  scores <- sweep(xmat, 2, mu, "-") %*% eig$vectors[, seq_len(k), drop = FALSE]
  colnames(scores) <- paste0("cellMCD_PC", seq_len(k))
  rownames(scores) <- rownames(xmat)

  # list(
  #   scores = scores,
  #   xmat = xmat,
  #   genes = rownames(scores)
  # )
  scores
}


add_hclust_clusters <- function(
    gout,
    scores,
    k = NULL,
    k_range = 4:10,
    hclust_method = "ward.D2",
    cluster_col = NULL
) {
  
  
  D = dist(scores)
  hc <- hclust(D, method = hclust_method)
  
  # automatically choose k using average silhouette width
  if (is.null(k)) {
    if (!requireNamespace("cluster", quietly = TRUE)) {
      stop("Package 'cluster' is needed for automatic k selection.")
    }
    
    n <- attr(D, "Size")
    k_range <- k_range[k_range >= 2 & k_range < n]
    
    if (length(k_range) == 0L) {
      k <- 4
    } else {
      sil_scores <- sapply(k_range, function(kk) {
        cl <- cutree(hc, k = kk)
        mean(cluster::silhouette(cl, D)[, "sil_width"])
      })
      
      k <- k_range[which.max(sil_scores)]
    }
  }
  
  if (is.null(cluster_col)) {
    cluster_col <- paste0("hclust_cluster")
  }
  
  cl <- cutree(hc, k = k)
  
  cl_dt <- data.table(
    gene = names(cl),
    hclust_cluster = factor(cl)
  )
  
  gout[,
       (cluster_col) := cl_dt[.SD, on = "gene", x.hclust_cluster]
  ]
  
  list(
    gout = gout,
    cluster_col = cluster_col
  )
}

# ------------------------------------------------------------
# Add embedding clusters per group
# ------------------------------------------------------------
add_embedding_clusters_all <- function(gout, tout, mcd_res,
                                       hclust_k = NULL,
                                       hclust_k_range = 4:10,
                                       hclust_method = "ward.D2",
                                       cluster_col = NULL) {
  gout <- as.data.table(copy(gout))
  groups <- unique(gout$group_id)
  
  cluster_list <- list()
  
  for (g in groups) {
    
    scores <- make_cellmcd_pc(
      tout = tout[group_id == g],
      cellmcd_out = mcd_res[[g]][[3]],
      group_name = g,
      n_pcs = 50
    )
    
  # compute hclust labels only if requested
    hc_res <- add_hclust_clusters(
      gout[group_id == g],
      scores,
      k = hclust_k,
      k_range = hclust_k_range,
      hclust_method = hclust_method,
      cluster_col = cluster_col
    )

    cluster_col <- hc_res$cluster_col
    
    cluster_list[[g]] <- unique(
      hc_res$gout[, .(
        group_id,
        gene,
        embedding_cluster = as.character(get(cluster_col))
      )]
    )
  }
  
  clust_dt <- rbindlist(cluster_list, fill = TRUE)
  
  gout <- merge(gout, clust_dt, by = c("group_id", "gene"), all.x = TRUE)
  tout <- merge(tout, clust_dt, by = c("group_id", "gene"), all.x = TRUE)
  
  list(gout = gout, tout = tout, clusters = clust_dt)
}

# ------------------------------------------------------------
# t-SNE wrapper for class or cluster
# ------------------------------------------------------------
make_embedding_plots <- function(
    tsub,
    gsub,
    cellmcd_out,
    group_name,
    color_var = c("gene_class", "embedding_cluster"),
    seed = 1,
    perplexity = NULL,
    hclust_k = NULL,
    hclust_k_range = 4:10,
    hclust_method = "ward.D2",
    label_top_genes = TRUE,
    top_n_genes = 10,
    label_rank_var = "ssZ"
    ) {
  
  color_var <- match.arg(color_var)
  
  gsub <- data.table::as.data.table(data.table::copy(gsub))

  scores <- make_cellmcd_pc(
    tout = tsub,
    cellmcd_out = cellmcd_out,
    group_name = group_name,
    n_pcs = 50
  )
  
  
  # compute hclust labels only if requested
  if (color_var == "embedding_cluster") {
    hc_res <- add_hclust_clusters(
      gsub,
      scores,
      k = hclust_k,
      k_range = hclust_k_range,
      hclust_method = hclust_method
    )
    
    gsub <- hc_res$gout
    color_col <- hc_res$cluster_col
  } else {
    color_col <- color_var
  }
  
  
  genes <- rownames(scores)
  
  annot <- gsub[group_id == group_name & gene %in% genes]
  annot <- annot[match(genes, gene)]
  
  n <- length(genes)
  
  if (is.null(perplexity)) {
    perplexity <- max(2, min(30, floor((n - 1) / 3)))
  }
  
  # ----- t-SNE -----------------------------------------------------
  set.seed(seed)
  
  tsne_out <- Rtsne::Rtsne(
    scores,
    # is_distance = TRUE,
    pca = FALSE,
    check_duplicates = FALSE,
    perplexity = perplexity
  )
  
  tsne_dt <- data.table(
    gene = genes,
    dim1 = tsne_out$Y[, 1],
    dim2 = tsne_out$Y[, 2],
    method = "t-SNE"
  )
  
  # ----- PC scores -----------------------------------------------------
  
  signed_log10 <- function(x, pseudocount = 1) {
    sign(x) * log10(abs(x) + pseudocount)
  }
  
  pc_dt <- data.table(
    gene = genes,
    dim1 = signed_log10(scores[, 1]),
    dim2 = signed_log10(scores[, 2]),
    method = "PCA"
  )
  
  
  
  # plot_dt <- rbindlist(list(tsne_dt, umap_dt))
  plot_dt <- rbindlist(list(tsne_dt, pc_dt))
  
  plot_dt <- merge(plot_dt, annot, by = "gene", all.x = TRUE)
  plot_dt = plot_dt[method == 't-SNE']
  if(color_col == 'gene_class'){plot_dt = plot_dt[gene_class != "None"]}
  
  p <- ggplot(
    plot_dt,
    aes(
      x = dim1,
      y = dim2,
      color = .data[[color_col]]
    )
  ) +
    geom_point(alpha = 0.75, size = 2) +
    # facet_wrap(~method, scales = "free") +
    theme_bw() +
    labs(
      x = "t-SNE 1",
      y = "t-SNE 2",
      color = color_col
    )
  
  # --------------- top gene labels------------------------
  
  if (label_top_genes && label_rank_var %in% names(plot_dt)) {
    label_genes <- unique(
      plot_dt[
        !is.na(get(label_rank_var)),
        .(gene, color_col = get(color_col),rank_val = get(label_rank_var))
      ]
    )
    
    label_genes <- label_genes[order(-rank_val)]
    label_genes <- label_genes[, head(.SD, top_n_genes), .SDcol = "gene", by = color_col]
    label_dt <- plot_dt[gene %in% label_genes$gene]
    
    p <- p +
      ggrepel::geom_text_repel(
        data = label_dt, fontface = "bold",
        ggplot2::aes(label = gene),
        size = 5, force_pull = 10,
        max.overlaps = Inf,
        show.legend = FALSE
      ) + labs(title = group_name,
                  subtitle = "Singificant genes classified based on Global/Local structure and dominant sign amongst samples.",
                  caption = paste("Top", top_n_genes, "genes per class labeled based on", label_rank_var))
  }  
  
  p
}

# ------------------------------------------------------------
# Mahalanobis vs eff_n or p_global
# ------------------------------------------------------------

plot_mahalanobis_summary <- function(
    dt,
    x_var = c("p_global", "eff_n"),
    y_var = c("ssZ", "mahalanobis"),
    grouping_var,
    pattern_col = NULL,
    x_cutoff = NULL,
    facet_scales = "fixed",
    top_n_labels = 10,
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
  
  required_cols <- c(pattern_col, x_var, y_var, grouping_var, 
                     "group_id", "local_dominant_sign_meaning","gene")
  
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
    facet_grid(
      stats::as.formula(paste0("local_dominant_sign_meaning ~ ", grouping_var)),
      scales = facet_scales
    ) +
    theme_bw() +
    labs(
      title = unique(dt$group_id)[1],
      x = ifelse(!is.null(xlabel), xlabel, paste0("log10(", x_var, ")")),
      y = ifelse(!is.null(ylabel), ylabel, y_var),
      color = "Dominant Pattern",
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
# Centroid heatmap
# ------------------------------------------------------------
make_centroid_heatmap <- function(tout, group_name, grouping_var,
                                  value_col = "Zres",
                                  transform = "signed_log10") {
  dt <- as.data.table(copy(tout))
  
  signed_log10 <- function(x) sign(x) * log10(abs(x) + 1)
  
  sub <- dt[
    group_id == group_name &
      pattern != "none" &
      !is.na(get(grouping_var)) &
      is.finite(get(value_col))
  ]
  
  if (nrow(sub) == 0L) return(NULL)
  
  cent <- sub[
    ,
    .(value = mean(get(value_col), na.rm = TRUE)),
    by = .(grouping = get(grouping_var), contrast)
  ]
  
  if (transform == "signed_log10") {
    cent[, value := signed_log10(value)]
  }
  
  wide <- dcast(cent, grouping ~ contrast, value.var = "value")
  if (nrow(wide) < 1L || ncol(wide) < 3L) return(NULL)
  
  mat <- as.matrix(wide[, -1, with = FALSE])
  rownames(mat) <- wide$grouping
  
  lim <- max(abs(mat), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) return(NULL)
  
  col_fun <- circlize::colorRamp2(
    c(-lim, 0, lim),
    c("navy", "white", "firebrick")
  )
  
  ComplexHeatmap::Heatmap(
    mat,
    name = paste0("Median ", value_col),
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    border = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 7),
    column_names_rot = 45,
    column_title = paste("Centroids:", group_name, "|", grouping_var)
  )
}

# ------------------------------------------------------------
# Gene heatmap without built-in clustering
# ------------------------------------------------------------
make_gene_heatmap_by_grouping <- function(tout, group_name, grouping_var,
                                          grouping_value,
                                          value_col = "Zres",
                                          top_n_genes = 100,
                                          rank_var = "mahalanobis",
                                          transform = "signed_log10") {
  dt <- as.data.table(copy(tout))
  
  signed_log10 <- function(x) sign(x) * log10(abs(x) + 1)
  
  sub <- dt[
    group_id == group_name &
      get(grouping_var) == grouping_value &
      pattern != "none" &
      is.finite(get(value_col))
  ]
  
  if (nrow(sub) == 0L) return(NULL)
  
  if (!is.null(top_n_genes) && rank_var %in% names(sub)) {
    keep <- sub[
      is.finite(get(rank_var)),
      .(rank_value = max(get(rank_var), na.rm = TRUE)),
      by = gene
    ][order(-rank_value)][seq_len(min(.N, top_n_genes)), gene]
    
    sub <- sub[gene %in% keep]
  }
  
  mat_dt <- sub[, .(value = mean(get(value_col), na.rm = TRUE)), by = .(gene, contrast)]
  
  if (transform == "signed_log10") {
    mat_dt[, value := signed_log10(value)]
  }
  
  wide <- dcast(mat_dt, gene ~ contrast, value.var = "value")
  if (nrow(wide) < 1L || ncol(wide) < 3L) return(NULL)
  
  mat <- as.matrix(wide[, -1, with = FALSE])
  rownames(mat) <- wide$gene
  mat <- mat[complete.cases(mat), , drop = FALSE]
  
  if (nrow(mat) < 1L || ncol(mat) < 2L) return(NULL)
  
  lim <- max(abs(mat), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) return(NULL)
  
  col_fun <- circlize::colorRamp2(
    c(-lim, 0, lim),
    c("navy", "white", "firebrick")
  )
  
  ComplexHeatmap::Heatmap(
    mat,
    name = value_col,
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    border = TRUE,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 7),
    column_names_rot = 45,
    column_title = paste(group_name, "|", 
                         grouping_var, "=", 
                         grouping_value, "| Top", 
                         top_n_genes, "Genes")
  )
}

# ------------------------------------------------------------
# Zres vs gene index, one page per cluster/class
# ------------------------------------------------------------
plot_zres_by_gene_index <- function(tout, group_name, grouping_var,
                                    grouping_value,
                                    top_n_each_sign = 10) {
  dt <- data.table::as.data.table(data.table::copy(tout))
  
  sub <- dt[
    group_id == group_name &
      get(grouping_var) == grouping_value &
      pattern != "none" &
      is.finite(Zres)
  ]
  
  if (nrow(sub) == 0L) return(NULL)
  
  gene_rank <- sub[
    ,
    .(
      z_value = mean(Zres, na.rm = TRUE)
    ),
    by = .(gene, contrast)
  ]
  
  gene_rank <- gene_rank[
    order(contrast, z_value)
  ]
  
  gene_rank[
    ,
    gene_index := seq_len(.N),
    by = contrast
  ]
  
  pos_genes <- gene_rank[
    z_value > 0
  ][
    order(contrast, -z_value),
    head(.SD, top_n_each_sign),
    by = contrast
  ][
    , .(gene, contrast, label_sign = "positive")
  ]
  
  neg_genes <- gene_rank[
    z_value < 0
  ][
    order(contrast, z_value),
    head(.SD, top_n_each_sign),
    by = contrast
  ][
    , .(gene, contrast, label_sign = "negative")
  ]
  
  top_genes <- data.table::rbindlist(
    list(pos_genes, neg_genes),
    fill = TRUE
  )
  
  plot_dt <- merge(
    gene_rank,
    top_genes,
    by = c("gene", "contrast"),
    all.x = TRUE
  )
  
  plot_dt[, label_gene := data.table::fifelse(
    !is.na(label_sign),
    gene,
    NA_character_
  )]
  
  plot_dt[, is_top_gene := !is.na(label_sign)]
  
  
  signed_log10 <- function(x) sign(x) * log10(abs(x) + 1)
  
  ggplot(plot_dt, aes(
    x = gene_index,
    y = signed_log10(z_value)
  )) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(
      aes(color = is_top_gene),
      alpha = 0.75,
      size = 1.4
    ) +
    geom_line(color = "grey70", linewidth = 0.4) +
    ggrepel::geom_text_repel(
      data = plot_dt[!is.na(label_gene)],
      aes(label = label_gene),
      size = 2.6,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    facet_wrap(~contrast, scales = "free") +
    theme_bw() +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
    labs(
      title = paste(group_name, "|", grouping_var, "=", grouping_value),
      x = "Gene index",
      y = "Zres",
      color = "Top gene"
    ) +
    theme(legend.position = "none")
}
# ------------------------------------------------------------
# One full PDF: either clusters or classes
# ------------------------------------------------------------
make_grouping_pdf <- function(gout, tout, mcd_res, out_file,
                              grouping_var = "gene_class",
                              x_summary = c("eff_n", "p_global"),
                              top_n_genes = 10) {
  x_summary <- match.arg(x_summary)
  
  groups <- unique(gout$group_id)
  
  pdf(out_file, width = 15, height = 11)
  
  # first: one summary plot across all groups
  p_all <- plot_grouping_summary_bar(
    gout = gout,
    grouping_var = grouping_var
  )
  
  if (!is.null(p_all)) print(p_all)
  
  for (g in groups) {
    p_group <- plot_grouping_summary_bar(
      gout = gout,
      grouping_var = grouping_var,
      group_name = g
    )
    
    if (!is.null(p_group)) print(p_group)
  }
  
  for (g in groups) {
    
    color_var <- if (grouping_var == "gene_class") "gene_class" else "embedding_cluster"
    
    p <- make_embedding_plots(
      tsub = tout[group_id == g],
      gsub = gout[group_id == g],
      cellmcd_out = mcd_res[[g]][[3]],
      group_name = g,
      color_var = color_var,
      label_top_genes = TRUE,
      top_n_genes = top_n_genes,
      label_rank_var = "mahalanobis"
    )
    
    print(p)
  }
  
  for (g in groups) {

    p <- plot_mahalanobis_summary(
      gout[group_id == g], 
      x_var = x_summary,
      y_var = "mahalanobis",
      grouping_var = grouping_var
    )
    
    print(p)

  }
  for (g in groups) {
    ht_centroid <- make_centroid_heatmap(
      tout = tout,
      group_name = g,
      grouping_var = grouping_var
    )
    
    if (!is.null(ht_centroid)) {
      ComplexHeatmap::draw(ht_centroid)
    }
  }
  
  for (g in groups) {
    vals <- unique(
      na.omit(
        tout[group_id == g, get(grouping_var)]
      )
    )
    
    for (v in vals) {
      ht_gene <- make_gene_heatmap_by_grouping(
        tout = tout,
        group_name = g,
        grouping_var = grouping_var,
        grouping_value = v
      )
      
      if (!is.null(ht_gene)) {
        ComplexHeatmap::draw(ht_gene)
      }
    }
  }
  
  for (g in groups) {
    vals <- unique(na.omit(tout[group_id == g, get(grouping_var)]))
    
    for (v in vals) {
      p <- plot_zres_by_gene_index(
        tout = tout,
        group_name = g,
        grouping_var = grouping_var,
        grouping_value = v
      )
      
      if (!is.null(p)) print(p)
    }
  }
  
  dev.off()
}

# ------------------------------------------------------------
# Pair plot PDF
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


make_pair_pdf <- function(tout, out_file) {
  groups <- unique(tout$group_id)
  
  pdf(out_file, width = 11, height = 11)
  
  for (g in groups) {
    plot_pairs(
      dt = tout,
      group_name = g,
      keep_mode = "rowwise_only",
      value_col = "Zres",
      top_n_genes = 10,
      gene_rank_col = "mahalanobis"
    )
  }
  
  dev.off()
}

# ------------------------------------------------------------
# XLSX export: one worksheet per group_id
# ------------------------------------------------------------

sign_col <- function(x, cutoff = 0.25) {
  fifelse(
    is.na(x), NA_character_,
    fifelse(x >= cutoff, "positive",
            fifelse(x <= -cutoff, "negative", "neutral"))
  )
}

export_significant_gene_table <- function(gout, out_file,
                                          contrast_cols = NULL,
                                          sign_cutoff = 0.25) {
  wb <- openxlsx::createWorkbook()
  
  gout <- data.table::as.data.table(data.table::copy(gout))
  
  if (is.null(contrast_cols)) {
    contrast_cols <- grep("^contrast_|^coef_|^B_|^W_", names(gout), value = TRUE)
  }
  
  export_cols <- intersect(
    c(
      "group_id", "embedding_cluster", "gene_class", "gene",
      "pattern", "mahalanobis", "p_value", "ssZ",
      contrast_cols
    ),
    names(gout)
  )
  
  pos_style <- openxlsx::createStyle(fontColour = "#9C0006", fgFill = "#FFC7CE")
  neg_style <- openxlsx::createStyle(fontColour = "#006100", fgFill = "#C6EFCE")
  
  for (g in unique(gout$group_id)) {
    sheet <- substr(gsub("[^A-Za-z0-9_]", "_", g), 1, 31)
    openxlsx::addWorksheet(wb, sheet)
    
    dt <- unique(gout[group_id == g, ..export_cols])
    data.table::setorder(dt, embedding_cluster, gene_class, -mahalanobis)
    
    openxlsx::writeDataTable(wb, sheet, dt)
    openxlsx::freezePane(wb, sheet, firstRow = TRUE)
    openxlsx::setColWidths(wb, sheet, cols = seq_along(dt), widths = "auto")
    
    for (cc in contrast_cols) {
      if (!cc %in% names(dt)) next
      
      col_idx <- match(cc, names(dt))
      vals <- dt[[cc]]
      
      pos_rows <- which(vals >= sign_cutoff)
      neg_rows <- which(vals <= -sign_cutoff)
      
      if (length(pos_rows) > 0) {
        openxlsx::addStyle(
          wb, sheet, pos_style,
          rows = pos_rows + 1,
          cols = col_idx,
          gridExpand = TRUE,
          stack = TRUE
        )
      }
      
      if (length(neg_rows) > 0) {
        openxlsx::addStyle(
          wb, sheet, neg_style,
          rows = neg_rows + 1,
          cols = col_idx,
          gridExpand = TRUE,
          stack = TRUE
        )
      }
    }
  }
  
  openxlsx::saveWorkbook(wb, out_file, overwrite = TRUE)
}

export_tout_zres_table <- function(tout, gout, out_file) {
  wb <- openxlsx::createWorkbook()
  
  tout <- data.table::as.data.table(data.table::copy(tout))
  gout <- data.table::as.data.table(data.table::copy(gout))
  
  gene_annot_cols <- intersect(
    c(
      "group_id", "gene",
      "embedding_cluster", "gene_class",
      "mahalanobis", "ssZ", "p_global", "eff_n"
    ),
    names(gout)
  )
  
  gene_annot <- unique(gout[, ..gene_annot_cols])
  
  tout2 <- merge(
    tout,
    gene_annot,
    by = intersect(c("group_id", "gene"), names(tout)),
    all.x = TRUE,
    suffixes = c("", "_gene")
  )
  
  export_cols <- intersect(
    c(
      "group_id", "gene", "contrast",
      "treatment", "cell_line", "timepoint",
      "Zres", "W",
      "embedding_cluster", "gene_class",
      "mahalanobis", "ssZ", "p_global", "eff_n"
    ),
    names(tout2)
  )
  
  for (g in unique(tout2$group_id)) {
    sheet <- substr(gsub("[^A-Za-z0-9_]", "_", g), 1, 31)
    openxlsx::addWorksheet(wb, sheet)
    
    dt <- tout2[group_id == g, ..export_cols]
    data.table::setorder(dt, embedding_cluster, gene_class, gene, contrast)
    
    openxlsx::writeDataTable(wb, sheet, dt)
    openxlsx::freezePane(wb, sheet, firstRow = TRUE)
    openxlsx::setColWidths(wb, sheet, cols = seq_along(dt), widths = "auto")
  }
  
  openxlsx::saveWorkbook(wb, out_file, overwrite = TRUE)
}
# ------------------------------------------------------------
# Main pipeline
# ------------------------------------------------------------
make_plot_pipeline <- function(out.list, mcd_res, out_dir,
                               ssZ_quantile = NULL,
                               x_summary = c("eff_n", "p_global"),
                               hclust_k = NULL,
                               top_n_genes=10) {
  x_summary <- match.arg(x_summary)
  
  tout <- as.data.table(out.list$main_summary)
  gout <- as.data.table(out.list$gene_summary)
  
  filt <- filter_interesting_genes(
    gout = gout,
    tout = tout,
    ssZ_quantile = ssZ_quantile
  )
  
  gout <- filt$gout
  tout <- filt$tout
  
  clust_res <- add_embedding_clusters_all(
    gout = gout,
    tout = tout,
    mcd_res = mcd_res,
    hclust_k = hclust_k
  )
  
  gout <- clust_res$gout
  tout <- clust_res$tout
  
  make_grouping_pdf(
    gout = gout,
    tout = tout,
    mcd_res = mcd_res,
    out_file = file.path(out_dir, "clustering_results.pdf"),
    grouping_var = "embedding_cluster",
    x_summary = x_summary,
    top_n_genes = top_n_genes
  )
  
  make_grouping_pdf(
    gout = gout,
    tout = tout,
    mcd_res = mcd_res,
    out_file = file.path(out_dir, "class_results.pdf"),
    grouping_var = "gene_class",
    x_summary = x_summary
  )
  
  make_pair_pdf(
    tout = tout,
    out_file = file.path(out_dir, "pair_plots.pdf")
  )
  
  export_significant_gene_table(
    gout = gout,
    out_file = file.path(out_dir, "significant_gene_table.xlsx"),
    contrast_cols = grep("_vs_", names(gout), value = TRUE),
    sign_cutoff = 0.25
  )
  
  export_tout_zres_table(
    tout = tout,
    gout = gout,
    out_file = file.path(out_dir, "contrasts_table.xlsx")
  )
  

  invisible(list(
    gout = gout,
    tout = tout,
    embedding_clusters = clust_res$clusters
  ))
}