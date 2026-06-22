#!/usr/bin/env Rscript

# Reproducible target-gene expression analysis across clinical and integrated
# mutation subtypes from the latest LSC17 subtype pipeline.

options(stringsAsFactors = FALSE)

message <- function(...) {
  cat(paste0(..., collapse = ""), "\n", sep = "")
}

args <- commandArgs(trailingOnly = TRUE)
target_genes <- c("CD96", "TNFRSF10A")
if (length(args) >= 1 && nzchar(args[1])) {
  target_genes <- unique(trimws(unlist(strsplit(args[1], split = ","))))
  target_genes <- target_genes[nzchar(target_genes)]
}
min_subtype_n <- 10
if (length(args) >= 2 && nzchar(args[2])) {
  min_subtype_n <- as.integer(args[2])
}
if (is.na(min_subtype_n) || min_subtype_n < 2) {
  stop("min_subtype_n must be an integer >= 2")
}

get_script_file <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0) {
    return(NA_character_)
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)
}

script_file <- get_script_file()
if (!is.na(script_file)) {
  project_dir <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = FALSE)
} else {
  project_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
if (!file.exists(file.path(project_dir, "01_data", "Normalized Expression.csv"))) {
  project_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

input_expression <- file.path(project_dir, "01_data", "Normalized Expression.csv")
input_sample_info <- file.path(
  project_dir,
  "03_result", "LSC17", "subtype_surface_latest",
  "lsc17_sample_scores_with_subtypes.csv"
)
input_mutation_long <- file.path(
  project_dir,
  "03_result", "LSC17", "subtype_surface_latest",
  "integrated_mutation_assignments_long.csv"
)
output_dir <- file.path(project_dir, "03_result", "target_gene_expression")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

clinical_subtype_variables <- c(
  "WHO_specificDx",
  "ELN2017",
  "Fusion_subtype",
  "FAB_subtype",
  "Disease_origin",
  "Disease_stage",
  "Specimen_type"
)

message("Target-gene expression subtype analysis")
message("Project directory: ", project_dir)
message("Output directory: ", output_dir)
message("Target genes: ", paste(target_genes, collapse = ", "))
message("Minimum subtype sample size: ", min_subtype_n)
message("Expression input: ", input_expression)
message("Sample subtype input: ", input_sample_info)
message("Integrated mutation input: ", input_mutation_long)

required_files <- c(input_expression, input_sample_info, input_mutation_long)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files: ", paste(missing_files, collapse = "; "))
}

clean_label <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x[tolower(x) %in% c("na", "n/a", "nan", "null")] <- NA_character_
  x
}

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
  message("Wrote: ", path, " (", nrow(x), " rows)")
}

safe_wilcox <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) {
    return(NA_real_)
  }
  if (length(unique(c(x, y))) < 2) {
    return(NA_real_)
  }
  out <- tryCatch(
    stats::wilcox.test(x, y, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  out
}

message("Reading expression matrix")
expr_raw <- utils::read.csv(input_expression, check.names = FALSE)
required_expr_cols <- c("stable_id", "display_label")
if (!all(required_expr_cols %in% names(expr_raw))) {
  stop("Expression matrix must contain columns: ", paste(required_expr_cols, collapse = ", "))
}
sample_cols <- setdiff(names(expr_raw), c("stable_id", "display_label", "description", "biotype"))
message("Expression genes: ", nrow(expr_raw))
message("Expression sample columns: ", length(sample_cols))

message("Reading latest subtype tables")
sample_info <- utils::read.csv(input_sample_info, check.names = FALSE)
mutation_long <- utils::read.csv(input_mutation_long, check.names = FALSE)
if (!"sample_id" %in% names(sample_info)) {
  stop("Sample subtype file must contain sample_id")
}
if (!all(c("sample_id", "mutation") %in% names(mutation_long))) {
  stop("Mutation long file must contain sample_id and mutation")
}

sample_info$sample_id <- clean_label(sample_info$sample_id)
mutation_long$sample_id <- clean_label(mutation_long$sample_id)
mutation_long$mutation <- clean_label(mutation_long$mutation)

common_samples <- intersect(sample_cols, sample_info$sample_id)
message("Samples with both expression and subtype information: ", length(common_samples))
if (length(common_samples) < min_subtype_n) {
  stop("Too few samples with both expression and subtype information")
}

clinical_vars_found <- intersect(clinical_subtype_variables, names(sample_info))
message("Clinical subtype variables used: ", paste(clinical_vars_found, collapse = ", "))
if (length(clinical_vars_found) == 0) {
  warning("No configured clinical subtype variables were found in the sample subtype file")
}

sample_info_aligned <- sample_info[sample_info$sample_id %in% common_samples, , drop = FALSE]
sample_info_aligned <- sample_info_aligned[match(common_samples, sample_info_aligned$sample_id), , drop = FALSE]

message("Extracting target-gene expression")
gene_presence <- data.frame(
  requested_gene = target_genes,
  matched = FALSE,
  matched_display_label = NA_character_,
  matched_stable_id = NA_character_,
  n_matching_rows = 0L,
  stringsAsFactors = FALSE
)

expr_long_list <- list()
display_upper <- toupper(clean_label(expr_raw$display_label))
for (gene in target_genes) {
  idx <- which(display_upper == toupper(gene))
  presence_idx <- match(gene, gene_presence$requested_gene)
  gene_presence$n_matching_rows[presence_idx] <- length(idx)
  if (length(idx) == 0) {
    message("Gene not found in expression matrix: ", gene)
    next
  }
  gene_presence$matched[presence_idx] <- TRUE
  gene_presence$matched_display_label[presence_idx] <- paste(unique(expr_raw$display_label[idx]), collapse = "|")
  gene_presence$matched_stable_id[presence_idx] <- paste(unique(expr_raw$stable_id[idx]), collapse = "|")

  gene_expr_df <- expr_raw[idx, common_samples, drop = FALSE]
  gene_expr_num <- as.data.frame(lapply(gene_expr_df, function(z) suppressWarnings(as.numeric(z))))
  values <- colMeans(gene_expr_num, na.rm = TRUE)
  values[is.nan(values)] <- NA_real_
  expr_long_list[[gene]] <- data.frame(
    sample_id = common_samples,
    gene = gene,
    expression = as.numeric(values),
    stringsAsFactors = FALSE
  )
}
expr_long <- do.call(rbind, expr_long_list)
if (is.null(expr_long) || nrow(expr_long) == 0) {
  write_csv(gene_presence, file.path(output_dir, "target_gene_presence.csv"))
  stop("None of the requested target genes were found in the expression matrix")
}
message("Target-gene expression rows: ", nrow(expr_long))

message("Building clinical subtype assignments")
clinical_assignments <- list()
for (v in clinical_vars_found) {
  subtype_level <- clean_label(sample_info_aligned[[v]])
  clinical_assignments[[v]] <- data.frame(
    sample_id = sample_info_aligned$sample_id,
    subtype_type = "clinical",
    subtype_variable = v,
    subtype_level = subtype_level,
    stringsAsFactors = FALSE
  )
}
clinical_assignments <- if (length(clinical_assignments) > 0) {
  do.call(rbind, clinical_assignments)
} else {
  data.frame(
    sample_id = character(),
    subtype_type = character(),
    subtype_variable = character(),
    subtype_level = character(),
    stringsAsFactors = FALSE
  )
}
clinical_assignments <- clinical_assignments[!is.na(clinical_assignments$subtype_level), , drop = FALSE]

message("Building integrated mutation subtype assignments")
mutation_assignments <- mutation_long[
  mutation_long$sample_id %in% common_samples & !is.na(mutation_long$mutation),
  c("sample_id", "mutation"),
  drop = FALSE
]
mutation_assignments <- unique(mutation_assignments)
mutation_assignments <- data.frame(
  sample_id = mutation_assignments$sample_id,
  subtype_type = "mutation",
  subtype_variable = "integrated_mutation",
  subtype_level = mutation_assignments$mutation,
  stringsAsFactors = FALSE
)

all_assignments <- unique(rbind(clinical_assignments, mutation_assignments))
all_assignments <- all_assignments[all_assignments$sample_id %in% common_samples, , drop = FALSE]
message("Subtype assignment rows: ", nrow(all_assignments))

count_key <- paste(all_assignments$subtype_type, all_assignments$subtype_variable, all_assignments$subtype_level, sep = "\r")
split_samples <- split(all_assignments$sample_id, count_key)
subtype_counts <- do.call(rbind, lapply(names(split_samples), function(k) {
  parts <- strsplit(k, "\r", fixed = TRUE)[[1]]
  data.frame(
    subtype_type = parts[[1]],
    subtype_variable = parts[[2]],
    subtype_level = parts[[3]],
    n_samples = length(unique(split_samples[[k]])),
    stringsAsFactors = FALSE
  )
}))
subtype_counts <- subtype_counts[order(subtype_counts$subtype_type, subtype_counts$subtype_variable, -subtype_counts$n_samples), ]
eligible_groups <- subtype_counts[subtype_counts$n_samples >= min_subtype_n, , drop = FALSE]
message("Eligible subtype groups: ", nrow(eligible_groups))

mutation_list_by_sample <- aggregate(
  subtype_level ~ sample_id,
  data = mutation_assignments,
  FUN = function(x) paste(sort(unique(x)), collapse = "|")
)
names(mutation_list_by_sample)[names(mutation_list_by_sample) == "subtype_level"] <- "integrated_mutations"

sample_keep_cols <- unique(c(
  "sample_id",
  "subject_id",
  "cohort",
  "sex",
  "race",
  "ethnicity",
  "age_at_diagnosis",
  "age_at_specimen",
  clinical_vars_found,
  "lsc17_score"
))
sample_keep_cols <- intersect(sample_keep_cols, names(sample_info_aligned))
sample_annotation <- sample_info_aligned[, sample_keep_cols, drop = FALSE]
sample_annotation <- merge(sample_annotation, mutation_list_by_sample, by = "sample_id", all.x = TRUE, sort = FALSE)
sample_expression_output <- merge(expr_long, sample_annotation, by = "sample_id", all.x = TRUE, sort = FALSE)
sample_expression_output <- sample_expression_output[order(sample_expression_output$gene, sample_expression_output$sample_id), ]

message("Calculating subtype expression summaries and tests")
summary_rows <- list()
test_rows <- list()
all_expression_samples <- unique(expr_long$sample_id[is.finite(expr_long$expression)])

for (gene in unique(expr_long$gene)) {
  gene_expr <- expr_long[expr_long$gene == gene, c("sample_id", "expression"), drop = FALSE]
  gene_expr <- gene_expr[is.finite(gene_expr$expression), , drop = FALSE]
  for (i in seq_len(nrow(eligible_groups))) {
    group <- eligible_groups[i, , drop = FALSE]
    in_ids <- unique(all_assignments$sample_id[
      all_assignments$subtype_type == group$subtype_type &
        all_assignments$subtype_variable == group$subtype_variable &
        all_assignments$subtype_level == group$subtype_level
    ])
    in_ids <- intersect(in_ids, gene_expr$sample_id)
    if (group$subtype_type == "clinical") {
      background_ids <- unique(all_assignments$sample_id[
        all_assignments$subtype_type == "clinical" &
          all_assignments$subtype_variable == group$subtype_variable
      ])
    } else {
      background_ids <- all_expression_samples
    }
    background_ids <- intersect(background_ids, gene_expr$sample_id)
    out_ids <- setdiff(background_ids, in_ids)

    in_expr <- gene_expr$expression[match(in_ids, gene_expr$sample_id)]
    out_expr <- gene_expr$expression[match(out_ids, gene_expr$sample_id)]
    in_expr <- in_expr[is.finite(in_expr)]
    out_expr <- out_expr[is.finite(out_expr)]

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      gene = gene,
      subtype_type = group$subtype_type,
      subtype_variable = group$subtype_variable,
      subtype_level = group$subtype_level,
      n_samples = length(in_expr),
      mean_expression = if (length(in_expr) > 0) mean(in_expr) else NA_real_,
      median_expression = if (length(in_expr) > 0) stats::median(in_expr) else NA_real_,
      sd_expression = if (length(in_expr) > 1) stats::sd(in_expr) else NA_real_,
      min_expression = if (length(in_expr) > 0) min(in_expr) else NA_real_,
      max_expression = if (length(in_expr) > 0) max(in_expr) else NA_real_,
      percent_expression_gt_0 = if (length(in_expr) > 0) mean(in_expr > 0) * 100 else NA_real_,
      stringsAsFactors = FALSE
    )

    p_value <- safe_wilcox(in_expr, out_expr)
    test_rows[[length(test_rows) + 1L]] <- data.frame(
      gene = gene,
      subtype_type = group$subtype_type,
      subtype_variable = group$subtype_variable,
      subtype_level = group$subtype_level,
      n_in = length(in_expr),
      n_out = length(out_expr),
      mean_in = if (length(in_expr) > 0) mean(in_expr) else NA_real_,
      mean_out = if (length(out_expr) > 0) mean(out_expr) else NA_real_,
      median_in = if (length(in_expr) > 0) stats::median(in_expr) else NA_real_,
      median_out = if (length(out_expr) > 0) stats::median(out_expr) else NA_real_,
      mean_difference = if (length(in_expr) > 0 && length(out_expr) > 0) mean(in_expr) - mean(out_expr) else NA_real_,
      median_difference = if (length(in_expr) > 0 && length(out_expr) > 0) stats::median(in_expr) - stats::median(out_expr) else NA_real_,
      wilcox_p = p_value,
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, summary_rows)
test_df <- do.call(rbind, test_rows)
test_df$fdr_bh_global <- NA_real_
test_df$fdr_bh_by_gene <- NA_real_
test_df$fdr_bh_by_gene_variable <- NA_real_
valid_p <- is.finite(test_df$wilcox_p)
test_df$fdr_bh_global[valid_p] <- stats::p.adjust(test_df$wilcox_p[valid_p], method = "BH")
for (g in unique(test_df$gene)) {
  idx <- which(test_df$gene == g & is.finite(test_df$wilcox_p))
  test_df$fdr_bh_by_gene[idx] <- stats::p.adjust(test_df$wilcox_p[idx], method = "BH")
}
gene_var_key <- paste(test_df$gene, test_df$subtype_type, test_df$subtype_variable, sep = "\r")
for (k in unique(gene_var_key)) {
  idx <- which(gene_var_key == k & is.finite(test_df$wilcox_p))
  test_df$fdr_bh_by_gene_variable[idx] <- stats::p.adjust(test_df$wilcox_p[idx], method = "BH")
}
test_df$is_higher_median <- test_df$median_difference > 0
test_df$is_higher_fdr05_by_gene_variable <- test_df$is_higher_median & test_df$fdr_bh_by_gene_variable < 0.05
test_df <- test_df[order(test_df$gene, test_df$fdr_bh_by_gene_variable, -test_df$median_difference), ]

high_df <- test_df[
  isTRUE(nrow(test_df) > 0) &
    !is.na(test_df$is_higher_fdr05_by_gene_variable) &
    test_df$is_higher_fdr05_by_gene_variable,
  ,
  drop = FALSE
]
high_df <- high_df[order(high_df$gene, high_df$fdr_bh_by_gene_variable, -high_df$median_difference), ]

message("Creating plot data")
eligible_key <- paste(eligible_groups$subtype_type, eligible_groups$subtype_variable, eligible_groups$subtype_level, sep = "\r")
assign_key <- paste(all_assignments$subtype_type, all_assignments$subtype_variable, all_assignments$subtype_level, sep = "\r")
plot_assignments <- all_assignments[assign_key %in% eligible_key, , drop = FALSE]
plot_df <- merge(expr_long, plot_assignments, by = "sample_id", all.x = FALSE, all.y = FALSE, sort = FALSE)
plot_df <- plot_df[is.finite(plot_df$expression), , drop = FALSE]

has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

make_clinical_pdf <- function(plot_df, out_pdf) {
  clinical_plot <- plot_df[plot_df$subtype_type == "clinical", , drop = FALSE]
  if (nrow(clinical_plot) == 0) {
    return(FALSE)
  }
  grDevices::pdf(out_pdf, width = 9.5, height = 6.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (has_ggplot2) {
    for (v in unique(clinical_plot$subtype_variable)) {
      d <- clinical_plot[clinical_plot$subtype_variable == v, , drop = FALSE]
      med_order <- aggregate(expression ~ subtype_level, data = d, FUN = median)
      med_order <- med_order[order(med_order$expression), ]
      d$subtype_level <- factor(d$subtype_level, levels = med_order$subtype_level)
      p <- ggplot2::ggplot(d, ggplot2::aes(x = subtype_level, y = expression)) +
        ggplot2::geom_boxplot(outlier.shape = NA, fill = "#D9EAF7", color = "#2D3E50", linewidth = 0.25) +
        ggplot2::geom_jitter(width = 0.16, height = 0, alpha = 0.38, size = 0.7, color = "#C44E52") +
        ggplot2::coord_flip() +
        ggplot2::facet_wrap(stats::as.formula("~ gene"), scales = "free_x") +
        ggplot2::labs(
          title = paste("Target-gene expression by", v),
          x = NULL,
          y = "Normalized expression"
        ) +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0),
          panel.grid.major.y = ggplot2::element_blank(),
          strip.background = ggplot2::element_rect(fill = "#F0F0F0", color = "#BDBDBD")
        )
      print(p)
    }
  } else {
    for (v in unique(clinical_plot$subtype_variable)) {
      for (g in unique(clinical_plot$gene)) {
        d <- clinical_plot[clinical_plot$subtype_variable == v & clinical_plot$gene == g, , drop = FALSE]
        graphics::boxplot(
          expression ~ subtype_level,
          data = d,
          horizontal = TRUE,
          las = 1,
          main = paste(g, "by", v),
          xlab = "Normalized expression"
        )
      }
    }
  }
  TRUE
}

make_mutation_pdf <- function(plot_df, out_pdf) {
  mutation_plot <- plot_df[plot_df$subtype_type == "mutation", , drop = FALSE]
  if (nrow(mutation_plot) == 0) {
    return(FALSE)
  }
  grDevices::pdf(out_pdf, width = 9.5, height = 7.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (has_ggplot2) {
    for (g in unique(mutation_plot$gene)) {
      d <- mutation_plot[mutation_plot$gene == g, , drop = FALSE]
      med_order <- aggregate(expression ~ subtype_level, data = d, FUN = median)
      med_order <- med_order[order(med_order$expression), ]
      d$subtype_level <- factor(d$subtype_level, levels = med_order$subtype_level)
      p <- ggplot2::ggplot(d, ggplot2::aes(x = subtype_level, y = expression)) +
        ggplot2::geom_boxplot(outlier.shape = NA, fill = "#E6D7B9", color = "#2D3E50", linewidth = 0.25) +
        ggplot2::geom_jitter(width = 0.16, height = 0, alpha = 0.38, size = 0.7, color = "#4C72B0") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = paste(g, "expression by integrated mutation subtype"),
          x = NULL,
          y = "Normalized expression"
        ) +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0),
          panel.grid.major.y = ggplot2::element_blank()
        )
      print(p)
    }
  } else {
    for (g in unique(mutation_plot$gene)) {
      d <- mutation_plot[mutation_plot$gene == g, , drop = FALSE]
      graphics::boxplot(
        expression ~ subtype_level,
        data = d,
        horizontal = TRUE,
        las = 1,
        main = paste(g, "by integrated mutation subtype"),
        xlab = "Normalized expression"
      )
    }
  }
  TRUE
}

make_dotplot_pdf <- function(test_df, out_pdf) {
  plot_tests <- test_df[
    is.finite(test_df$median_difference) &
      is.finite(test_df$fdr_bh_by_gene_variable) &
      test_df$n_in >= min_subtype_n &
      test_df$n_out >= min_subtype_n,
    ,
    drop = FALSE
  ]
  if (nrow(plot_tests) == 0) {
    return(FALSE)
  }
  plot_tests <- do.call(rbind, lapply(split(plot_tests, plot_tests$gene), function(d) {
    d <- d[order(-d$median_difference, d$fdr_bh_by_gene_variable), , drop = FALSE]
    head(d, 30)
  }))
  plot_tests$subtype_label <- paste(plot_tests$subtype_variable, plot_tests$subtype_level, sep = ": ")
  grDevices::pdf(out_pdf, width = 9.5, height = 7.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (has_ggplot2) {
    for (g in unique(plot_tests$gene)) {
      d <- plot_tests[plot_tests$gene == g, , drop = FALSE]
      d <- d[order(d$median_difference), , drop = FALSE]
      d$subtype_label <- factor(d$subtype_label, levels = unique(d$subtype_label))
      p <- ggplot2::ggplot(
        d,
        ggplot2::aes(
          x = median_difference,
          y = subtype_label,
          size = n_in,
          color = -log10(pmax(fdr_bh_by_gene_variable, .Machine$double.xmin))
        )
      ) +
        ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", color = "#777777") +
        ggplot2::geom_point(alpha = 0.86) +
        ggplot2::scale_color_gradient(low = "#4C72B0", high = "#C44E52", name = "-log10 FDR") +
        ggplot2::scale_size_continuous(name = "n") +
        ggplot2::labs(
          title = paste(g, "subtypes ranked by higher median expression"),
          x = "Median expression difference vs comparator",
          y = NULL
        ) +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0),
          panel.grid.major.y = ggplot2::element_blank()
        )
      print(p)
    }
  } else {
    for (g in unique(plot_tests$gene)) {
      d <- plot_tests[plot_tests$gene == g, , drop = FALSE]
      d <- d[order(d$median_difference), , drop = FALSE]
      graphics::dotchart(
        d$median_difference,
        labels = d$subtype_label,
        main = paste(g, "higher-expression subtype ranking"),
        xlab = "Median expression difference vs comparator"
      )
      graphics::abline(v = 0, lty = 2, col = "gray50")
    }
  }
  TRUE
}

clinical_pdf <- file.path(output_dir, "target_gene_expression_clinical_subtype_boxplots.pdf")
mutation_pdf <- file.path(output_dir, "target_gene_expression_mutation_subtype_boxplots.pdf")
dotplot_pdf <- file.path(output_dir, "target_gene_expression_high_subtype_dotplot.pdf")
clinical_pdf_made <- make_clinical_pdf(plot_df, clinical_pdf)
mutation_pdf_made <- make_mutation_pdf(plot_df, mutation_pdf)
dotplot_pdf_made <- make_dotplot_pdf(test_df, dotplot_pdf)
message("Clinical subtype PDF created: ", clinical_pdf_made)
message("Mutation subtype PDF created: ", mutation_pdf_made)
message("High-subtype dotplot PDF created: ", dotplot_pdf_made)

pipeline_info <- data.frame(
  item = c(
    "project_dir",
    "expression_input",
    "sample_subtype_input",
    "mutation_input",
    "output_dir",
    "target_genes",
    "min_subtype_n",
    "expression_sample_columns",
    "samples_with_expression_and_subtype",
    "clinical_variables_used",
    "eligible_subtype_groups",
    "ggplot2_available"
  ),
  value = c(
    project_dir,
    input_expression,
    input_sample_info,
    input_mutation_long,
    output_dir,
    paste(target_genes, collapse = ","),
    as.character(min_subtype_n),
    as.character(length(sample_cols)),
    as.character(length(common_samples)),
    paste(clinical_vars_found, collapse = ","),
    as.character(nrow(eligible_groups)),
    as.character(has_ggplot2)
  ),
  stringsAsFactors = FALSE
)

write_csv(gene_presence, file.path(output_dir, "target_gene_presence.csv"))
write_csv(sample_expression_output, file.path(output_dir, "target_gene_expression_sample_level.csv"))
write_csv(subtype_counts, file.path(output_dir, "target_gene_expression_subtype_sample_counts.csv"))
write_csv(eligible_groups, file.path(output_dir, "target_gene_expression_eligible_subtypes.csv"))
write_csv(summary_df, file.path(output_dir, "target_gene_expression_by_subtype_summary.csv"))
write_csv(test_df, file.path(output_dir, "target_gene_expression_subtype_tests.csv"))
write_csv(high_df, file.path(output_dir, "target_gene_expression_high_subtypes_FDR05.csv"))
write_csv(pipeline_info, file.path(output_dir, "pipeline_run_info.csv"))

message("Top higher-expression subtype results by gene")
if (nrow(high_df) == 0) {
  message("No subtype passed median_difference > 0 and FDR < 0.05 within gene/subtype-variable families.")
} else {
  for (g in unique(high_df$gene)) {
    message("Gene: ", g)
    d <- high_df[high_df$gene == g, , drop = FALSE]
    d <- head(d[order(d$fdr_bh_by_gene_variable, -d$median_difference), ], 10)
    printable <- d[, c(
      "subtype_type",
      "subtype_variable",
      "subtype_level",
      "n_in",
      "median_in",
      "median_out",
      "median_difference",
      "fdr_bh_by_gene_variable"
    ), drop = FALSE]
    print(printable, row.names = FALSE)
  }
}

message("Session information")
print(utils::sessionInfo())
message("Analysis complete")
