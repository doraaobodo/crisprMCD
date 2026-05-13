###############################################################################
# main.R
# CLI for CRISPR Drug Screen MCD analysis
###############################################################################

###############################################################################
# Dependencies and setup
###############################################################################

source("run_mcd.R")
source("make_plots.R")


check_required_packages = function(pkgs) {
  missing = pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) == 0) {
    return(invisible(TRUE))
  }

  show_header("Missing Required Packages")
  cat("The following packages are required but not installed:\n")
  for (pkg in missing) cat(" - ", pkg, "\n", sep = "")

  cat("\nInstalling them in R with:\n")
  cat(sprintf(
    'install.packages(c(%s))\n',
    paste(sprintf('"%s"', missing), collapse = ", ")
  ))
  cat("\n")
  
  install.packages(missing)

  # invisible(read_cli_input("Press [Enter] to exit..."))
  # quit(status = 1)
}

###############################################################################
# Console and CLI utilities
###############################################################################

read_cli_input = function(prompt_text, lowercase = FALSE) {
  ans = if (interactive()) {
    readline(prompt = prompt_text)
  } else {
    cat(prompt_text)
    readLines("stdin", n = 1)
  }

  ans = trimws(ans)
  if (lowercase) {
    ans = tolower(ans)
  }
  ans
}

show_header = function(title) {
  cat("\n")
  cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat(title, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")
}

run_or_exit = function(expr, title = "Error") {
  tryCatch(
    expr,
    error = function(e) {
      show_header(title)
      cat(conditionMessage(e), "\n")
      invisible(read_cli_input("Press Enter to exit..."))
      quit(status = 1)
    }
  )
}

ask_yes_no = function(prompt, default = TRUE) {
  suffix = if (default) " [Y/n]: " else " [y/N]: "

  repeat {
    ans = read_cli_input(paste0(prompt, suffix), lowercase = TRUE)
    if (ans == "") return(default)
    if (ans %in% c("y", "yes")) return(TRUE)
    if (ans %in% c("n", "no")) return(FALSE)
    cat("Please enter y or n.\n")
  }
}

###############################################################################
# File and path utilities
###############################################################################

choose_input_source_mode = function() {
  show_header("Select Input Source")
  cat("How would you like to provide the input?\n\n")
  cat("  1) Single file\n")
  cat("  2) Multiple files\n")
  cat("  3) Folder containing files\n\n")

  repeat {
    ans = read_cli_input("Enter choice number: ")
    if (ans == "1") return("single")
    if (ans == "2") return("multiple")
    if (ans == "3") return("directory")
    cat("Invalid choice. Please enter 1, 2, or 3.\n")
  }
}

choose_input_files = function(
    prompt_text = "Enter full paths separated by commas: ",
    caption = "Select input files"
) {
  repeat {
    selected = character(0)

    if (.Platform$OS.type == "windows") {
      selected = tryCatch(
        utils::choose.files(caption = caption, multi = TRUE),
        error = function(e) character(0)
      )
    }

    if (length(selected) > 0 && all(file.exists(selected))) {
      return(normalizePath(selected, winslash = "/", mustWork = TRUE))
    }

    cat("File chooser unavailable or no files selected.\n")
    ans = read_cli_input(prompt_text)
    parts = trimws(strsplit(ans, ",", fixed = TRUE)[[1]])
    parts = path.expand(parts[nzchar(parts)])

    if (length(parts) > 0 && all(file.exists(parts))) {
      return(normalizePath(parts, winslash = "/", mustWork = TRUE))
    }

    cat("One or more files were not found. Please try again.\n")
  }
}

choose_output_dir <- function(dirname = "results", base_dir = getwd(), create = TRUE) {
  out_dir <- file.path(base_dir, dirname)
  
  if (dir.exists(out_dir)) {
    return(normalizePath(out_dir, winslash = "/", mustWork = TRUE))
  }
  
  if (create) {
    ok <- dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (ok && dir.exists(out_dir)) {
      return(normalizePath(out_dir, winslash = "/", mustWork = TRUE))
    }
  }
  
  stop(
    "Output directory does not exist and could not be created: ",
    out_dir,
    call. = FALSE
  )
}

collect_input_paths = function() {
  mode = choose_input_source_mode()

  if (mode == "single") {
    return(list(
      mode = mode,
      paths = choose_input_files(
        prompt_text = "Enter full path to input file: ",
        caption = "Select input file"
      )
    ))
  }

  if (mode == "multiple") {
    return(list(
      mode = mode,
      paths = choose_input_files(
        prompt_text = "Enter full paths to input files, separated by commas: ",
        caption = "Select input files"
      )
    ))
  }

  dir_path = choose_output_dir(
    prompt_text = "Enter folder containing input files: ",
    caption = "Select input folder",
    create = FALSE
  )

  files = list.files(dir_path, full.names = TRUE)
  files = files[file.info(files)$isdir %in% FALSE]

  if (length(files) == 0) {
    stop("No files were found in the selected directory.")
  }

  list(mode = mode, paths = normalizePath(files, winslash = "/", mustWork = TRUE))
}

load_table_data = function(file_path) {
  ext = tolower(tools::file_ext(file_path))

  switch(
    ext,
    "csv" = read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE),
    "tsv" = read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE),
    "txt" = read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE),
    stop("Unsupported file type: ", ext, ". Please use csv, tsv, or txt.")
  )
}

detect_input_type_from_paths = function(paths) {
  types = vapply(paths, function(path) {
    dat = load_table_data(path)
    nms = colnames(dat)
    
    if ("neg|lfc" %in% nms) {
      return("mageck_rra")
    }
    
    if (any(grepl("\\|z$", nms)) || any(grepl("\\|beta$", nms))) {
      return("mageck_mle")
    }
    
    "raw_counts"
    
  }, character(1))

  unique_types = unique(types)
  if (length(unique_types) != 1) {
    stop(
      "Selected files appear to contain mixed input types: ",
      paste(sprintf("%s (%s)", basename(paths), types), collapse = "; ")
    )
  }

  unique_types
}

###############################################################################
# Import sheet helpers
###############################################################################

normalize_yes_no = function(x, default = "yes") {
  x = trimws(as.character(x))
  x[x == ""] = default
  x_low = tolower(x)
  x_low[x_low %in% c("y", "yes", "true", "1")] = "yes"
  x_low[x_low %in% c("n", "no", "false", "0")] = "no"
  x_low
}


guess_column_role = function(column_name, dat) {
  low = tolower(column_name)
  is_num = is.numeric(dat[[column_name]])
  
  if (grepl("guide|grna|sgrna|barcode|construct", low)) return("guide")
  if (grepl("gene|symbol|hgnc|target|^id$", low)) return("gene")
  if (is_num) return("sample")
  "ignore"
}


parse_sample_name = function(sample_name) {
  tokens = unlist(strsplit(sample_name, "[-_.]+", perl = TRUE))
  tokens = tokens[nzchar(tokens)]
  tok_upper = toupper(tokens)
  
  out = list(
    sample = sample_name,
    cell_line = if (length(tokens) > 0) tokens[1] else NA_character_,
    treatment = NA_character_,
    treatment_type = NA_character_,
    rep = NA_integer_,
    timepoint = NA_character_
  )
  
  if (length(tokens) == 0) {
    return(out)
  }
  
  time_idx = grep("^(T|TP|D|DAY|TIME|TIMEPOINT).*$", tok_upper)
  if (length(time_idx) > 0) {
    out$timepoint = tokens[time_idx[1]]
    time_idx = time_idx[1]
  }
  
  rep_idx = grep("^(REP|R|REPLICATE).*$", tok_upper)
  if (length(rep_idx) > 0) {
    out$rep = suppressWarnings(as.integer(sub("^(REP|R|REPLICATE)", "", tok_upper[rep_idx[1]])))
    rep_idx = rep_idx[1]
  }
  
  remaining_idx = seq_along(tokens)
  remaining_idx = setdiff(remaining_idx, c(1, time_idx, rep_idx))
  
  if (length(remaining_idx) > 0) {
    trt_tok = tokens[remaining_idx[length(remaining_idx)]]
    trt_core = sub("\\d+$", "", trt_tok)
    out$treatment = toupper(trt_core)
    
    if (is.na(out$rep) && grepl("\\d+$", trt_tok)) {
      out$rep = suppressWarnings(as.integer(sub("^.*?(\\d+)$", "\\1", trt_tok)))
    }
  }
  
  trt = toupper(out$treatment)
  out$treatment_type = if (
    !is.na(trt) && nzchar(trt) && trt %in% c("DMSO", "CTRL", "CONTROL", "UNTREATED", "VEHICLE", "PBS")
  ) {
    "control"
  } else if (!is.na(trt) && nzchar(trt)) {
    "treatment"
  } else {
    NA_character_
  }
  
  out
}


clean_sample_name = function(x) {
  x = trimws(as.character(x))
  x[x == ""] = NA_character_
  x
}


build_counts_import_sheet = function(paths) {
  rows = lapply(paths, function(path) {
    dat = load_table_data(path)
    nms = colnames(dat)
    
    
    out = data.frame(
      file = normalizePath(path, winslash = "/", mustWork = TRUE),
      file_label = tools::file_path_sans_ext(basename(path)),
      column_name = nms,
      role = vapply(nms, function(x) guess_column_role(x, dat), character(1)),
      sample = NA_character_,
      cell_line = NA_character_,
      treatment = NA_character_,
      treatment_type = NA_character_,
      rep = NA_integer_,
      timepoint = NA_character_,
      include_sample = NA_character_,
      Group_By = NA_character_,
      stringsAsFactors = FALSE
    )
    
    sample_idx = which(out$role == "sample")
    if (length(sample_idx) > 0) {
      parsed = lapply(out$column_name[sample_idx], parse_sample_name)
      out$sample[sample_idx] = vapply(parsed, function(x) x$sample, character(1))
      out$cell_line[sample_idx] = vapply(parsed, function(x) x$cell_line, character(1))
      out$treatment[sample_idx] = vapply(parsed, function(x) x$treatment, character(1))
      out$treatment_type[sample_idx] = vapply(parsed, function(x) x$treatment_type, character(1))
      out$rep[sample_idx] = vapply(parsed, function(x) x$rep, integer(1))
      out$timepoint[sample_idx] = vapply(parsed, function(x) x$timepoint, character(1))
      out$include_sample[sample_idx] = "yes"
    }
    out
    
    })

  do.call(rbind, rows)
}

build_standardized_import_sheet = function(paths, input_type) {
  rows = lapply(paths, function(path) {
    dat = load_table_data(path)

    if (input_type == "mageck_rra") {
      sample_name = tools::file_path_sans_ext(basename(path))
      sample_name = gsub(".gene_summary|_gene_summary", "", sample_name)
      sample_name = gsub(".rra|_rra", "", sample_name)
      parsed = parse_sample_name(sample_name)

      data.frame(
        file = path,
        column_name = sample_name,
        role = "sample",
        sample = sample_name,
        cell_line = parsed$cell_line,
        treatment = parsed$treatment,
        treatment_type = parsed$treatment_type,
        rep = if (is.na(parsed$rep)) NA_integer_ else parsed$rep,
        timepoint = parsed$timepoint,
        include_sample = "yes",
        Group_By = parsed$Group_By,
        stringsAsFactors = FALSE
      )
    } else {
      metric_cols = grep("\\|(z|beta)$", colnames(dat), value = TRUE)
      sample_names = sub("\\|(z|beta)$", "", metric_cols)
      guessed = lapply(sample_names, parse_sample_name)

      data.frame(
        file = path,
        column_name = metric_cols,
        role = "sample",
        sample = sample_names,
        cell_line = vapply(guessed, function(x) x$cell_line, character(1)),
        treatment = vapply(guessed, function(x) x$treatment, character(1)),
        treatment_type = vapply(guessed, function(x) x$treatment_type, character(1)),
        rep = vapply(guessed, function(x) if (is.na(x$rep)) NA_integer_ else x$rep, integer(1)),
        timepoint = vapply(guessed, function(x) x$timepoint, character(1)),
        include_sample = "yes",
        Group_By = vapply(guessed, function(x) x$Group_By, character(1)),
        stringsAsFactors = FALSE
      )
    }
  })

  do.call(rbind, rows)
}

build_import_sheet = function(paths, input_type) {
  if (input_type == "raw_counts") {
    return(build_counts_import_sheet(paths))
  }

  build_standardized_import_sheet(paths, input_type)
}

edit_import_sheet = function(import_sheet, input_type) {
  show_header("Import Sheet Preview")
  cat("Review and edit the import sheet.\n")
  cat("\nImportant columns:\n")
  cat(" - include_sample: yes / no\n")
  cat(" - Group_By: grouping label to use later for MCD analysis\n")
  if (input_type == "raw_counts") {
    cat(" - role: gene / guide / sample / ignore\n")
  }
  cat("\n")
  print(import_sheet)
  cat("\n")

  cat("Opening the import sheet editor... \n")
  edited = utils::edit(import_sheet)

  if (is.null(edited)) {
    stop("Import sheet editing was cancelled.")
  }

  as.data.frame(edited, stringsAsFactors = FALSE)
}

validate_import_sheet_common = function(import_sheet) {
  required_cols = c(
    "file", "column_name", "role", "sample", "cell_line", "treatment",
    "treatment_type", "rep", "timepoint", "include_sample", "Group_By"
  )
  missing_cols = setdiff(required_cols, colnames(import_sheet))
  if (length(missing_cols) > 0) {
    stop("Import sheet is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  import_sheet$include_sample = normalize_yes_no(import_sheet$include_sample, default = "yes")
  bad_include = !(is.na(import_sheet$include_sample) | import_sheet$include_sample %in% c("yes", "no"))
  if (any(bad_include)) {
    stop("include_sample must contain only yes or no.")
  }

  import_sheet$sample = clean_sample_name(import_sheet$sample)
  import_sheet$cell_line = clean_sample_name(import_sheet$cell_line)
  import_sheet$treatment = clean_sample_name(import_sheet$treatment)
  import_sheet$treatment_type = clean_sample_name(import_sheet$treatment_type)
  import_sheet$timepoint = clean_sample_name(import_sheet$timepoint)
  import_sheet$Group_By = clean_sample_name(import_sheet$Group_By)

  import_sheet$rep = suppressWarnings(as.integer(import_sheet$rep))
  import_sheet
}

validate_counts_import_sheet = function(import_sheet) {
  import_sheet = validate_import_sheet_common(import_sheet)

  for (path in unique(import_sheet$file)) {
    one = import_sheet[import_sheet$file == path, , drop = FALSE]
    gene_n = sum(one$role == "gene", na.rm = TRUE)
    if (gene_n != 1) {
      stop("Each counts file must have exactly one column marked as 'gene': ", basename(path))
    }

    guide_n = sum(one$role == "guide", na.rm = TRUE)
    if (guide_n > 1) {
      stop("Each counts file can have at most one column marked as 'guide': ", basename(path))
    }
    
    sample_bin = one$role == "sample" & one$include_sample == "yes"
    sample_rows = one[sample_bin, , drop = FALSE]
    if (nrow(sample_rows) == 0) {
      stop("Each counts file must have at least one column marked as 'sample': ", basename(path))
    }

    required_sample_cols = c("sample", "cell_line", "treatment", "treatment_type", "include_sample")
    for (col in required_sample_cols) {
      bad = is.na(sample_rows[[col]]) | trimws(as.character(sample_rows[[col]])) == ""
      if (any(bad)) {
        stop("Counts sample rows must have non-empty values for '", col, "' in file: ", basename(path))
      }
    }
  }
  
  # for group by
  if(!exists('Group_By', import_sheet)){
    stop('Please do not remove Group_By. Leave Blank if no grouping variable required')
  } else if (all(is.na(import_sheet$Group_By))){
    cat("No Grouping Variable Provided. Setting Group By to None")
    import_sheet$Group_By = ifelse(import_sheet$include_sample == 'yes', "none", NA)
  } else {
    ng = paste(import_sheet$Group_By[which(!is.na(import_sheet$Group_By))], collapse = ", ")
    import_sheet$Group_By = ifelse(import_sheet$include_sample == 'yes', ng, NA)
  }

  import_sheet
}

validate_standardized_import_sheet = function(import_sheet) {
  import_sheet = validate_import_sheet_common(import_sheet)

  sample_rows = import_sheet[import_sheet$role == "sample", , drop = FALSE]
  if (nrow(sample_rows) == 0) {
    stop("No sample rows were found in the import sheet.")
  }

  required_sample_cols = c("sample", "cell_line", "treatment", "treatment_type", "include_sample")
  for (col in required_sample_cols) {
    bad = is.na(sample_rows[[col]]) | trimws(as.character(sample_rows[[col]])) == ""
    if (any(bad)) {
      stop("Sample rows must have non-empty values for '", col, "'.")
    }
  }
  
  # for group by
  if(!exists('Group_By', import_sheet)){
    stop('Please do not remove Group_By. Leave Blank if no grouping variable required')
  } else if (all(is.na(import_sheet$Group_By))){
    cat("No Grouping Variable Provided. Setting Group By to None")
    import_sheet$Group_By = ifelse(import_sheet$include_sample == 'yes', "none", NA)
  } else {
    ng = paste(import_sheet$Group_By[which(!is.na(import_sheet$Group_By))], collapse = ", ")
    import_sheet$Group_By = ifelse(import_sheet$include_sample == 'yes', ng, NA)
  }
  

  import_sheet
}

###############################################################################
# Parsing helpers
###############################################################################

merge_gene_tables = function(tbls, by_cols = "gene") {
  Reduce(function(x, y) merge(x, y, by = by_cols, all = TRUE), tbls)
}

make_unique_sample_names = function(x) {
  make.unique(as.character(x), sep = "__")
}

parse_counts_from_import_sheet = function(paths, import_sheet) {
  parsed_files = lapply(paths, function(path) {
    dat = load_table_data(path)
    spec = import_sheet[import_sheet$file == path, , drop = FALSE]

    gene_col = spec$column_name[spec$role == "gene"][1]
    guide_col = if (any(spec$role == "guide")) spec$column_name[spec$role == "guide"][1] else NA_character_
    sample_spec = spec[spec$role == "sample" & spec$include_sample == "yes" , , drop = FALSE]

    if (!all(sample_spec$column_name %in% colnames(dat))) {
      stop("Some selected sample columns were not found in file: ", basename(path))
    }

    annotation = data.frame(
      gene = as.character(dat[[gene_col]]),
      guide = if (!is.na(guide_col)) as.character(dat[[guide_col]]) else as.character(dat[[gene_col]]),
      stringsAsFactors = FALSE
    )

    counts = dat[, sample_spec$column_name, drop = FALSE]
    counts = as.matrix(counts)


    list(
      annotation = annotation,
      counts = counts,
      sample_map = sample_spec,
      source_info = data.frame(
        file = path,
        gene_column = gene_col,
        guide_column = ifelse(is.na(guide_col), "", guide_col),
        stringsAsFactors = FALSE
      )
    )
  })

  combined = parsed_files[[1]]

  if (length(parsed_files) > 1) {
    for (i in 2:length(parsed_files)) {
      nxt = parsed_files[[i]]
      x = cbind(combined$annotation, combined$counts, stringsAsFactors = FALSE)
      y = cbind(nxt$annotation, nxt$counts, stringsAsFactors = FALSE)
      merged = merge(x, y, by = c("gene", "guide"), all = TRUE)

      annotation = merged[, c("gene", "guide"), drop = FALSE]
      counts = merged[, setdiff(colnames(merged), c("gene", "guide")), drop = FALSE]

      combined$annotation = annotation
      combined$counts = counts
      combined$sample_map = rbind(combined$sample_map, nxt$sample_map)
      combined$source_info = rbind(combined$source_info, nxt$source_info)
    }
  }

  list(
    input_type = "raw_counts",
    annotation = combined$annotation,
    counts_matrix = combined$counts,
    sample_names = colnames(combined$counts),
    source_info = combined$source_info,
    sample_map = combined$sample_map,
    import_sheet = import_sheet
  )
}

parse_one_mageck_rra_file = function(path, metric_col = "neg|lfc", gene_col = "id") {
  dat = load_table_data(path)

  required_cols = c(gene_col, metric_col)
  missing_cols = setdiff(required_cols, colnames(dat))
  if (length(missing_cols) > 0) {
    stop(
      basename(path),
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  out = dat[, c(gene_col, metric_col), drop = FALSE]
  colnames(out) = c("gene", "score")
  out$gene = as.character(out$gene)
  out
}

parse_mageck_rra_from_import_sheet = function(paths, import_sheet, metric_col = "neg|lfc") {
  score_list = vector("list", length(paths))
  sample_rows = vector("list", length(paths))

  for (i in seq_along(paths)) {
    path = paths[i]
    one = parse_one_mageck_rra_file(path, metric_col = metric_col)
    spec = import_sheet[import_sheet$file == path & import_sheet$role == "sample", , drop = FALSE]

    if (nrow(spec) != 1) {
      stop("Each RRA file must contribute exactly one sample row in the import sheet: ", basename(path))
    }

    sample_name = spec$sample[1]
    colnames(one)[2] = sample_name
    score_list[[i]] = one
    sample_rows[[i]] = spec
  }

  merged = merge_gene_tables(score_list, by_cols = "gene")
  rownames(merged) = merged$gene
  sample_names = unlist(lapply(sample_rows, function(x) x$sample), use.names = FALSE)
  score_matrix = merged[, sample_names, drop = FALSE]

  list(
    input_type = "mageck_rra",
    gene_data = data.frame(gene = merged$gene, stringsAsFactors = FALSE),
    score_matrix = score_matrix,
    sample_names = colnames(score_matrix),
    source_info = data.frame(file = paths, stringsAsFactors = FALSE),
    selected_metric = metric_col,
    sample_map = do.call(rbind, sample_rows),
    import_sheet = import_sheet
  )
}

parse_one_mageck_mle_file = function(path) {
  dat = load_table_data(path)

  required_cols = "Gene"
  missing_cols = setdiff(required_cols, colnames(dat))
  if (length(missing_cols) > 0) {
    stop(
      basename(path),
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  metric_cols = grep("\\|(z|beta)$", colnames(dat), value = TRUE)
  if (length(metric_cols) == 0) {
    stop(basename(path), " has no columns ending in '|z' or '|beta'.")
  }

  z_cols = grep("\\|z$", metric_cols, value = TRUE)
  keep_cols = if (length(z_cols) > 0) z_cols else metric_cols

  out = dat[, c("Gene", keep_cols), drop = FALSE]
  colnames(out)[1] = "gene"
  out$gene = as.character(out$gene)
  out
}

parse_mageck_mle_from_import_sheet = function(paths, import_sheet) {
  parsed_list = vector("list", length(paths))
  sample_rows = list()

  for (i in seq_along(paths)) {
    path = paths[i]
    dat = parse_one_mageck_mle_file(path)
    original_metric_cols = colnames(dat)[-1]
    sample_keys = sub("\\|(z|beta)$", "", original_metric_cols)

    spec = import_sheet[import_sheet$file == path & import_sheet$role == "sample", , drop = FALSE]
    idx = match(sample_keys, spec$sample)
    if (any(is.na(idx))) {
      stop("MLE import sheet is missing sample metadata for one or more metric columns in ", basename(path))
    }

    spec = spec[idx, , drop = FALSE]
    new_names = make_unique_sample_names(spec$sample)
    colnames(dat)[-1] = new_names
    spec$sample = new_names

    parsed_list[[i]] = dat
    sample_rows[[i]] = spec
  }

  merged = merge_gene_tables(parsed_list, by_cols = "gene")
  rownames(merged) = merged$gene
  score_matrix = merged[, setdiff(colnames(merged), "gene"), drop = FALSE]

  list(
    input_type = "mageck_mle",
    gene_data = data.frame(gene = merged$gene, stringsAsFactors = FALSE),
    score_matrix = score_matrix,
    sample_names = colnames(score_matrix),
    source_info = data.frame(file = paths, stringsAsFactors = FALSE),
    selected_metric = "z_or_beta",
    sample_map = do.call(rbind, sample_rows),
    import_sheet = import_sheet
  )
}

parse_input_from_import_sheet = function(input_type, paths, import_sheet) {
  if (input_type == "raw_counts") {
    return(parse_counts_from_import_sheet(paths, import_sheet))
  }
  if (input_type == "mageck_rra") {
    return(parse_mageck_rra_from_import_sheet(paths, import_sheet))
  }
  if (input_type == "mageck_mle") {
    return(parse_mageck_mle_from_import_sheet(paths, import_sheet))
  }
  stop("Unsupported input type: ", input_type)
}

import_input_pipeline = function(paths) {
  input_type = detect_input_type_from_paths(paths)
  import_sheet = build_import_sheet(paths, input_type)
  import_sheet = edit_import_sheet(import_sheet, input_type)

  if (input_type == "raw_counts") {
    import_sheet = validate_counts_import_sheet(import_sheet)
  } else {
    import_sheet = validate_standardized_import_sheet(import_sheet)
  }

  parsed_input = parse_input_from_import_sheet(input_type, paths, import_sheet)

  list(
    input_type = input_type,
    # import_sheet = import_sheet,
    parsed_input = parsed_input
  )
}

###############################################################################
# Metadata and sample helpers
###############################################################################

make_analysis_sample_name = function(cell_line, treatment, timepoint) {
  parts = c(cell_line, treatment, timepoint)
  parts = trimws(as.character(parts))
  parts = parts[nzchar(parts) & !is.na(parts)]
  paste(parts, collapse = "_")
}

###############################################################################
# Counts preprocessing
###############################################################################

infer_gene_from_guide <- function(guide_values, gene_values) {
  genes <- unique(trimws(as.character(gene_values)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  
  genes_lc <- tolower(genes)
  guides_lc <- tolower(as.character(guide_values))
  
  # Match longer gene symbols first to reduce partial-match errors
  o <- order(nchar(genes_lc), decreasing = TRUE)
  genes <- genes[o]
  genes_lc <- genes_lc[o]
  
  inferred <- vapply(guides_lc, function(g) {
    hit <- which(vapply(genes_lc, function(gene) {
      grepl(gene, g, fixed = TRUE)
    }, logical(1)))
    
    if (length(hit)) genes[hit[1]] else NA_character_
  }, character(1))
  
  inferred
}

aggregate_counts_to_gene_mean = function(annotation, counts_matrix, use_guide_matching = TRUE) {
  annotation <- as.data.frame(annotation, stringsAsFactors = FALSE)
  mat <- as.matrix(counts_matrix)
  
  if (nrow(mat) != nrow(annotation)) {
    stop("nrow(counts_matrix) must match nrow(annotation)")
  }
  
  gene_labels = as.character(annotation$gene)

  if (use_guide_matching && "guide" %in% colnames(annotation)) {
    
    missing_gene <- is.na(gene_labels) | !nzchar(trimws(gene_labels))
    
      if (any(missing_gene)) {
        inferred <- infer_gene_from_guide(
          guide_values = annotation$guide[missing_gene],
          gene_values  = annotation$gene
        )
        gene_labels[missing_gene] <- inferred
      }
    
  }
  
  keep <- !is.na(gene_labels) & nzchar(gene_labels)
  mat <- mat[keep, , drop = FALSE]
  gene_labels <- gene_labels[keep]
  
  gene_sums <- rowsum(mat, group = gene_labels, reorder = FALSE)
  gene_n <- tabulate(match(gene_labels, rownames(gene_sums)), nbins = nrow(gene_sums))
  
  gene_mat = gene_sums / gene_n


  list(
    gene_data = rownames(gene_mat),
    score_matrix = gene_mat
  )
}

collapse_replicate_columns = function(score_matrix, meta) {
  groups = split(seq_len(ncol(score_matrix)), meta$analysis_sample)
  collapsed = lapply(groups, function(idx) {
    block = score_matrix[, idx, drop = FALSE]
    if (ncol(block) == 1) {
      as.numeric(block[, 1])
    } else {
      rowMeans(block, na.rm = TRUE)
    }
  })

  out = do.call(cbind, collapsed)
  out
}

standardize_columns_to_z = function(score_matrix) {
  mat = as.matrix(score_matrix)
  storage.mode(mat) = "numeric"
  
  zmat = apply(mat, 2, function(x) {
    s = stats::sd(x, na.rm = TRUE)
    m = mean(x, na.rm = TRUE)
    
    if (is.na(s) || s == 0) {
      rep(0, length(x))
    } else {
      (x - m) / s
    }
  })
  
  zmat
}


build_gene_score_matrix_from_raw_counts = function(parsed_raw, meta) {

  keep_cols = meta$sample
  counts_subset = parsed_raw$counts_matrix[, keep_cols, drop = FALSE]

  has_user_guide = any(parsed_raw$source_info$guide_column %in% c("", NA) == FALSE)

  aggregated = aggregate_counts_to_gene_mean(
    annotation = parsed_raw$annotation,
    counts_matrix = counts_subset,
    use_guide_matching = has_user_guide
  )

  rep_collapsed = collapse_replicate_columns(aggregated$score_matrix, meta)
  z_scores = standardize_columns_to_z(rep_collapsed)
  
  
  split_meta = split(meta, meta$analysis_sample)
  
  rows = lapply(split_meta, function(df) {
    group_values = unique(df$Group_By)
    group_values = group_values[!is.na(group_values) & nzchar(group_values)]

    data.frame(
      sample = df$analysis_sample[1],
      cell_line = df$cell_line[1],
      treatment = df$treatment[1],
      treatment_type = df$treatment_type[1],
      rep = nrow(df),
      timepoint = df$timepoint[1],
      include_sample = "yes",
      Group_By = if (length(group_values) == 0) NA_character_ else group_values,
      original_samples = paste(df$sample, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  
  collapsed_meta = do.call(rbind, rows)

  list(
    input_type = "raw_counts",
    gene_data = aggregated$gene_data,
    score_matrix = z_scores,
    sample_names = colnames(z_scores),
    source_info = parsed_raw$source_info,
    selected_metric = "gene_level_z_from_raw_counts",
    metadata = collapsed_meta,
    preprocessing = list(
      gene_aggregation = "mean",
      replicate_aggregation = "mean",
      final_transform = "z_score_by_column"
    ),
    import_sheet = parsed_raw$import_sheet
  )
}

###############################################################################
# Final input assembly
###############################################################################

build_final_input = function(parsed_input) {
  
  meta = as.data.frame(parsed_input$sample_map, stringsAsFactors = FALSE)
  meta = meta[meta$role == "sample", , drop = FALSE]
  meta$include_sample = normalize_yes_no(meta$include_sample, default = "yes")
  meta = meta[meta$include_sample == "yes", , drop = FALSE]
  
  # standardize name to "cell line_treatment_timepoint"
  meta$analysis_sample = vapply(
    seq_len(nrow(meta)),
    function(i) make_analysis_sample_name(meta$cell_line[i], 
                                          meta$treatment[i], 
                                          meta$timepoint[i]),
    character(1)
  )
  

  if (nrow(meta) == 0) {
    stop("No samples were marked for inclusion.")
  }

  
  required_col = "sample"
  
  if (!required_col %in% colnames(meta)) {
    stop("Metadata must contain a 'sample' column.")
  }
  
  sample_names = colnames(parsed_input$counts_matrix)
  meta_samples = as.character(meta$sample)
  missing_in_counts = setdiff(meta_samples, sample_names)

  if (length(missing_in_counts) > 0) {
    stop("These metadata samples are not present in the matrix: ", paste(missing_in_counts, collapse = ", "))
  }
  
  if (parsed_input$input_type == "raw_counts") {
    return(build_gene_score_matrix_from_raw_counts(parsed_input, meta))
  }

  meta_aligned = align_metadata_to_samples(parsed_input$sample_names, meta)
  keep_meta = meta_aligned[meta_aligned$include_sample == "yes", , drop = FALSE]
  keep_cols = keep_meta$sample

  final_input = parsed_input
  final_input$score_matrix = final_input$score_matrix[, keep_cols, drop = FALSE]
  final_input$sample_names = colnames(final_input$score_matrix)
  final_input$metadata = keep_meta
  final_input
}


###############################################################################
# Summaries and output helpers
###############################################################################

# source(summary tables R file)

add_value <- function(x, default = NA_character_) {
  if (is.null(x) || length(x) == 0) NA_character_ else as.character(x)
}

add_nrow <- function(x) add_value(if (!is.null(x)) nrow(x) else NULL)
add_ncol <- function(x) add_value(if (!is.null(x)) ncol(x) else NULL)


add_sheet_with_data <- function(wb, sheet_name, x) {
  header_style <- openxlsx::createStyle(textDecoration = "bold")
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet = sheet_name, x = x, withFilter = TRUE)
  if (nrow(x) > 0 && ncol(x) > 0) {
    openxlsx::addStyle(
      wb,
      sheet = sheet_name,
      style = header_style,
      rows = 1,
      cols = seq_len(ncol(x)),
      gridExpand = TRUE
    )
  }
  openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet = sheet_name, cols = 1:ncol(x), widths = "auto")
}


add_score_sheet_with_data <- function(wb, sheet_name, x, highlight_cols = NULL) {
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet = sheet_name, x = x, withFilter = TRUE)
  
  header_style <- openxlsx::createStyle(textDecoration = "bold")
  pos_style <- openxlsx::createStyle(bgFill = "#F4CCCC")  # light red
  neg_style <- openxlsx::createStyle(bgFill = "#CFE2F3")  # light blue
  
  if (nrow(x) > 0 && ncol(x) > 0) {
    # Header styling
    openxlsx::addStyle(
      wb,
      sheet = sheet_name,
      style = header_style,
      rows = 1,
      cols = seq_len(ncol(x)),
      gridExpand = TRUE
    )
    
    # Match only desired columns
    if (!is.null(highlight_cols)) {
      cols_to_format <- which(colnames(x) %in% highlight_cols)
      
      if (length(cols_to_format) > 0) {
        # > 0.5 → light red
        openxlsx::conditionalFormatting(
          wb,
          sheet = sheet_name,
          cols = cols_to_format,
          rows = 2:(nrow(x) + 1),
          rule = ">0.5",
          style = pos_style
        )
        
        # < -0.5 → light blue
        openxlsx::conditionalFormatting(
          wb,
          sheet = sheet_name,
          cols = cols_to_format,
          rows = 2:(nrow(x) + 1),
          rule = "<-0.5",
          style = neg_style
        )
      }
    }
  }
  
  openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet = sheet_name, cols = 1:ncol(x), widths = "auto")
}

save_analysis_outputs <- function(final_input,
                                  res,
                                  out_dir = ".",
                                  file_name = "mcd_analysis_results.xlsx",
                                  alpha = 0.05,
                                  top_n = 25,
                                  overwrite = TRUE) {
 
  out_file <- file.path(out_dir, file_name)
  
  #-----------------------------
  # Extract objects safely
  #-----------------------------
  score_matrix <- final_input$score_matrix
  import_sheet <- final_input$import_sheet
  metadata <- final_input$metadata
  preprocessing <- final_input$preprocessing
  source_info <- final_input$source_info
  
  mcd_tbl <- res$mcd_combined_results
  cw_mcd_tbl <- res$cw_mcd_combined_results
  cw_gene <- lapply(res$cw_mcd_results_by_group, `[[`, 2)
  cw_gene = do.call(rbind, cw_gene)
  rownames(cw_gene) = NULL
  #-----------------------------
  # Build summary sheet
  #-----------------------------
  summary_df <- data.frame(
    field = c(
      "input_type",
      "selected_metric",
      "n_genes",
      "n_samples_score_matrix",
      "score_matrix_nrow",
      "score_matrix_ncol",
      "group_by",
      "contrast_mode",
      "control_treatment",
      "gene_aggregation",
      "replicate_aggregation",
      "final_transform",
      "mcd_combined_results_nrow",
      "cw_mcd_combined_results_nrow"
    ),
    value = c(
      add_value(final_input$input_type),
      add_value(final_input$selected_metric),
      add_value(length(final_input$gene_data), NA_integer_),
      add_value(length(final_input$sample_names), NA_integer_),
      add_nrow(score_matrix),
      add_ncol(score_matrix),
      add_value(paste(res$group_by, collapse = ";")),
      add_value(res$contrast_mode),
      add_value(res$control_treatment),
      add_value(preprocessing$gene_aggregation),
      add_value(preprocessing$replicate_aggregation),
      add_value(preprocessing$final_transform),
      add_nrow(mcd_tbl),
      add_nrow(cw_mcd_tbl)
    ),
    stringsAsFactors = FALSE
  )

  #-----------------------------
  # Convert matrices to data frames for writing
  #-----------------------------

  score_df <- NULL
  if (!is.null(score_matrix)) {
    score_df <- data.frame(
      gene = rownames(score_matrix),
      as.data.frame(score_matrix, check.names = FALSE, stringsAsFactors = FALSE),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  
  score_cols = as.vector(sapply(res$contrast_matrices, colnames))
  
  #-----------------------------
  # Add plotting tables (functions defined in run_mcd.R)
  #-----------------------------
  
  outlier_summary_tbl <- build_outlier_summary_table(res, alpha = 0.05, top_n = 100)
  outlier_summary_tbl = merge(outlier_summary_tbl, cw_gene, by=c('group_id', 'gene'))
  #-----------------------------
  # Add DBSCAN Clusters
  #-----------------------------
  
  db_res <- run_outlier_hdbscan(
    outlier_summary_tbl = outlier_summary_tbl,
    contrast_cols = score_cols,
    alpha = 0.05,
    minPts = 10
  )
  
  
  cw_mcd_tbl = merge(cw_mcd_tbl, 
                     db_res[, setdiff(colnames(db_res), score_cols)], 
                     by = c('gene', 'group_id'))
  

  
  #-----------------------------
  # Create results dictionary
  #-----------------------------
  
  base_dict <- data.frame(
    sheet_name = "mcd_summary",
    column_name = setdiff(colnames(outlier_summary_tbl), score_cols),
    description = c(
      "Group Id as defined in input sheet. Analysis was run separately per group id",
      "Gene symbol.",
      "Mahalanobis distance summarizing multivariate deviation from center across contrasts within the group. \n Larger distance indicate more extreme outliers.",
      "Rank of the Mahalanobis distance within the group. Smaller ranks indicating more extreme outliers.",
      "P-value associated with the Mahalanobis distance.",
      "Logical indicator of whether the gene was flagged as an outlier in the group.",
      # "Standard deviation of absolute contrast effects. Smaller values indicate small variation among contrasts",
      # "Interpretation of effect_indication: 'Specific' if effect_indication > 0.5, otherwise 'Global.'",
      "Number of samples flagged as extreme in the cell-wise analysis.",
      "Proportion of samples flagged as extreme in the cell-wise analysis",
      "Global score = proportion flagged * mean of the absolute Z residual of flagged samples. Global score is is high when many samples are flagged and flagged samples are extreme",
      "Local score = 1-proportion flagged * maximum absolute Z residual. Local Score is high when only a few samples are flagged but at least one is very extreme",
      "P_global = global_score/ (global_score + local_score). P_global is evidence leaning toward widespread outlyingness.",
      "Pattern classifies genes based on p_global. Pattern is global if the the proportion of global to local scores >0.6, meaning that samples indicate a more global pattern of effects within flagged cells.",
      "Dominant effect sign (positive or negative) amongst contrasts. Computed as (# positive contrasts / N) - (# negative contrasts / N).",
      "Interpretation of dominant_sign: Positive if sign_col > 0.5, Negative if sign_col < -0.5, otherwise Mixed.",
      "Dominant effect sign (positive or negative) amongst contrasts. Computed as (# positive contrasts / N) - (# negative contrasts / N).",
      "Interpretation of dominant_sign: Positive if sign_col > 0.5, Negative if sign_col < -0.5, otherwise Mixed.",
      NA,NA,NA,NA,NA, NA,NA,NA, NA
      
    ),
    stringsAsFactors = FALSE
  )
  
  contrast_dict <- data.frame(
    sheet_name = "mcd_summary",
    column_name = score_cols,
    description = paste("Contrast value for", score_cols, " \nRed (+) indicates resistant effect; Blue (-) indicates sensitivity effect."),
    stringsAsFactors = FALSE
  )
  
  results_dictionary <- rbind(contrast_dict,base_dict)
  
  
  #-----------------------------
  # Create workbook
  #-----------------------------
  
  wb <- openxlsx::createWorkbook()
  
  #-----------------------------
  # Save plotting outputs
  #-----------------------------
  
  
  if (!is.null(import_sheet)) {
    add_sheet_with_data(wb, "import_sheet", as.data.frame(import_sheet, stringsAsFactors = FALSE))
  }
  
  # 
  # if (!is.null(outlier_summary_tbl)) {
  #   add_score_sheet_with_data(wb, "mcd_summary", outlier_summary_tbl, score_cols)
  # }
  

  #-----------------------------
  # Raw Outputs
  #-----------------------------
  # if (!is.null(score_df)) {
  #   add_sheet_with_data(wb, "score_matrix", score_df)
  # }
  
  # if (!is.null(mcd_tbl)) {
  #   add_sheet_with_data(wb, "mcd_results", as.data.frame(mcd_tbl, stringsAsFactors = FALSE))
  # }

  if (!is.null(cw_mcd_tbl)) {
    add_score_sheet_with_data(wb, "main_summary", as.data.frame(cw_mcd_tbl, stringsAsFactors = FALSE),
                              c('X', 'Zres'))
  }

  if (!is.null(db_res)) {
    add_score_sheet_with_data(wb, "mcd_summary", db_res, score_cols)  
  }
  
  #-----------------------------
  # Save workbook
  #-----------------------------
  
  
  
  # Add dictionary FIRST
  add_sheet_with_data(wb, "results_dictionary", results_dictionary)
  
  add_sheet_with_data(wb, "run_summary", summary_df)
  
  openxlsx::saveWorkbook(wb, out_file, overwrite = overwrite)
  
  invisible(out_file)
  
  list('out_file'    = out_file, 
        'gene_summary' = db_res,
       'main_summary' = cw_mcd_tbl)
}

###############################################################################
# Main workflow
###############################################################################
# add video demo to github
# add example metadata table.txt to data file 
# update metadata to try and label baseline cols as well.

main = function() {
  
  utils.pkgs = c("tools", "utils", "robustbase", 
                 "cellWise", "fdrtool", "openxlsx")
  
  viz.pkgs = c("data.table","dplyr","ggplot2", 'ggrepel',
               "GGally","gridExtra", "circlize",
               "pheatmap", # "ComlpexHeatmap",
               "Rtsne","dbscan")
  
  
  check_required_packages(c(utils.pkgs, viz.pkgs))
  
  library(tools)
  library(utils)
  library(robustbase)
  library(cellWise)
  library(fdrtool)
  library(openxlsx)
  
  library(ggplot2)
  library(pheatmap)
  library(dbscan)
  library(dplyr)
  
  library(data.table)
  library(Rtsne)
  library(GGally)
  library(gridExtra)
  library(ggrepel)
  # library(ComplexHeatmap)
  library(circlize)
  

  show_header("CRISPR Drug Screen MCD UI")
  cat("Welcome. ")
  cat("This tool will guide you through loading screen data,\n")
  cat("reviewing a single import sheet, and running an MCD analysis.\n\n")
  cat("R version: ", R.version.string, "\n", sep = "")
  cat("OS: ", Sys.info()[["sysname"]], "\n", sep = "")

  input_src = run_or_exit(collect_input_paths(), "Input Selection Error")

  imported = run_or_exit(
    import_input_pipeline(input_src$paths),
    "Import Pipeline Error"
  )

  final_input = run_or_exit(
    build_final_input(imported$parsed_input),
    "Input Processing Error"
  )

  res = run_or_exit(
    run_mcd_pipeline(
      final_input = final_input,
      B=5
      ),
    "MCD Analysis Error"
  )

  out_dir = choose_output_dir()

  # save workbook
  # make sure run_mcd.R is sourced first!
  out.list <- save_analysis_outputs(
    final_input = final_input,
    res = res,
    out_dir = out_dir,
    file_name = 'mcd_analysis_results.xlsx',
    alpha = 0.05,
    top_n = 25
  )
  
  cat("Saved output to:", out.list$out_file, "\n")
  
  # # save RDS inputs for report
  # save_report_inputs(final_input, res, out_dir = out_dir)
  # 
  # # write Rmd
  # write_mcd_report_rmd(
  #   file = file.path(out_dir, "mcd_report.Rmd"),
  #   title = "MCD Analysis Report"
  # )
  
}

main()
