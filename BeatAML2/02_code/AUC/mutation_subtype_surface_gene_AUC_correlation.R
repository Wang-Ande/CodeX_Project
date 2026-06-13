# Mutation-subtype surface gene expression correlation with inhibitor AUC.
# Mutation subtypes and LSC17 correlations come from subtype_surface_latest.
# Outputs are written to 03_result/LSC17_AUC/mutation_subtype_surface_gene_auc.

suppressWarnings(suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(matrixStats)
  library(ggplot2)
}))

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

min_n_mutation_drug_auc <- 10
significance_fdr <- 0.05
top_n_auc_genes_per_pair <- 20

input_auc <- file.path(project_root, "01_data", "Inhibitor AUC values.csv")
input_expr <- file.path(project_root, "01_data", "Normalized Expression.csv")
input_lsc17_score <- file.path(project_root, "03_result", "LSC17", "LSC17_sample_score_df.csv")
input_mutation_assignments <- file.path(project_root, "03_result", "LSC17", "subtype_surface_latest", "integrated_mutation_assignments_long.csv")
input_surface_genes <- file.path(project_root, "03_result", "LSC17", "subtype_surface_latest", "surface_gene_list_used.csv")
input_lsc17_subtype_cor <- file.path(project_root, "03_result", "LSC17", "subtype_surface_latest", "subtype_surface_correlation_all.csv")
output_dir <- file.path(project_root, "03_result", "LSC17_AUC", "mutation_subtype_surface_gene_auc")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(paste0(...), "\n", sep = "")
}

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path, call. = FALSE)
  }
}

for (path in c(
  input_auc,
  input_expr,
  input_lsc17_score,
  input_mutation_assignments,
  input_surface_genes,
  input_lsc17_subtype_cor
)) {
  stop_if_missing(path)
}

spearman_matrix_vs_vector <- function(mat, y) {
  keep <- !is.na(y)
  mat <- mat[, keep, drop = FALSE]
  y <- y[keep]
  n <- length(y)

  expr_rank <- matrixStats::rowRanks(mat, ties.method = "average", preserveShape = TRUE)
  y_rank <- rank(y, ties.method = "average")
  rho <- as.vector(cor(t(expr_rank), y_rank, method = "pearson"))

  denom <- pmax(1 - rho^2, .Machine$double.eps)
  t_stat <- rho * sqrt((n - 2) / denom)
  pvalue <- 2 * pt(-abs(t_stat), df = n - 2)
  pvalue[is.na(rho)] <- NA_real_

  data.frame(rho = rho, pvalue = pvalue, stringsAsFactors = FALSE)
}

log_msg("Reading subtype, score, and gene-set inputs...")
mutation_assignments <- data.table::fread(input_mutation_assignments)
setnames(mutation_assignments, old = c("sample_id", "mutation"), new = c("Sample_id", "mutation"), skip_absent = TRUE)
mutation_assignments <- unique(mutation_assignments[!is.na(Sample_id) & Sample_id != "" & !is.na(mutation) & mutation != "", .(Sample_id, mutation)])

lsc17_score <- data.table::fread(input_lsc17_score)
lsc17_score <- unique(lsc17_score[!is.na(Sample_id) & Sample_id != "" & !is.na(LSC17_score), .(Sample_id, LSC17_score)], by = "Sample_id")

surface_genes <- data.table::fread(input_surface_genes)
if (!"gene" %in% names(surface_genes)) {
  stop("surface_gene_list_used.csv must contain a gene column.", call. = FALSE)
}
surface_gene_list <- unique(na.omit(as.character(surface_genes$gene)))

lsc17_subtype_cor <- data.table::fread(input_lsc17_subtype_cor)
required_lsc17_cols <- c("subtype_type", "subtype_variable", "subtype_level", "gene", "rho", "pvalue", "FDR_subtype")
missing_lsc17_cols <- setdiff(required_lsc17_cols, names(lsc17_subtype_cor))
if (length(missing_lsc17_cols) > 0) {
  stop("LSC17 subtype correlation file is missing columns: ", paste(missing_lsc17_cols, collapse = ", "), call. = FALSE)
}
lsc17_mutation_cor <- lsc17_subtype_cor[
  subtype_type == "mutation",
  .(
    mutation = subtype_level,
    gene,
    n_samples_lsc17 = n_samples,
    rho_lsc17 = rho,
    pvalue_lsc17 = pvalue,
    FDR_lsc17 = FDR_subtype
  )
]
lsc17_mutation_cor <- lsc17_mutation_cor[gene %in% surface_gene_list]

log_msg("Reading AUC table...")
auc_dt <- data.table::fread(input_auc, nThread = 2)
auc_clean <- auc_dt[
  !is.na(dbgap_rnaseq_sample) &
    dbgap_rnaseq_sample != "" &
    !is.na(inhibitor) &
    inhibitor != "" &
    !is.na(auc),
  .(
    Sample_id = dbgap_rnaseq_sample,
    inhibitor,
    auc,
    dbgap_subject_id,
    type,
    status,
    paper_inclusion,
    curve_type
  )
]
if (auc_clean[, any(duplicated(paste(Sample_id, inhibitor, sep = "||")))]) {
  log_msg("Duplicate sample-drug AUC rows detected; using mean AUC per Sample_id/inhibitor.")
  auc_clean <- auc_clean[
    ,
    .(
      auc = mean(auc, na.rm = TRUE),
      dbgap_subject_id = dplyr::first(dbgap_subject_id),
      type = dplyr::first(type),
      status = dplyr::first(status),
      paper_inclusion = dplyr::first(paper_inclusion),
      curve_type = dplyr::first(curve_type)
    ),
    by = .(Sample_id, inhibitor)
  ]
}

log_msg("Reading normalized expression matrix and selecting surface genes...")
expr_dt <- data.table::fread(input_expr, nThread = 2)
if (!"display_label" %in% names(expr_dt)) {
  stop("Expression file is missing display_label.", call. = FALSE)
}
gene_names <- expr_dt[["display_label"]]
if (any(is.na(gene_names) | gene_names == "")) {
  stop("Expression file contains missing display_label values.", call. = FALSE)
}
if (anyDuplicated(gene_names) > 0) {
  stop("Expression file contains duplicated display_label values.", call. = FALSE)
}

surface_gene_list <- intersect(surface_gene_list, gene_names)
if (length(surface_gene_list) == 0) {
  stop("No surface genes overlap with expression matrix.", call. = FALSE)
}
expr_dt <- expr_dt[display_label %in% surface_gene_list]
expr_sample_cols <- setdiff(names(expr_dt), c("stable_id", "display_label", "description", "biotype"))
expr_mat <- as.matrix(expr_dt[, ..expr_sample_cols])
storage.mode(expr_mat) <- "numeric"
rownames(expr_mat) <- expr_dt[["display_label"]]
rm(expr_dt)

common_lsc17_expr_samples <- intersect(lsc17_score$Sample_id, colnames(expr_mat))
mutation_assignments <- mutation_assignments[Sample_id %in% common_lsc17_expr_samples]
auc_clean <- auc_clean[Sample_id %in% common_lsc17_expr_samples]

alignment_summary <- data.frame(
  item = c(
    "mutation_assignment_rows_after_lsc17_expression_filter",
    "unique_mutation_samples_after_filter",
    "unique_mutation_subtypes_after_filter",
    "auc_rows_after_lsc17_expression_filter",
    "unique_auc_samples_after_filter",
    "surface_genes_used"
  ),
  value = c(
    nrow(mutation_assignments),
    uniqueN(mutation_assignments$Sample_id),
    uniqueN(mutation_assignments$mutation),
    nrow(auc_clean),
    uniqueN(auc_clean$Sample_id),
    length(surface_gene_list)
  )
)
data.table::fwrite(alignment_summary, file.path(output_dir, "sample_alignment_summary.csv"))

log_msg("Counting eligible mutation-drug pairs...")
mutation_auc_sample_level <- merge(
  mutation_assignments,
  auc_clean,
  by = "Sample_id",
  allow.cartesian = TRUE
)
mutation_auc_sample_level <- merge(
  mutation_auc_sample_level,
  lsc17_score,
  by = "Sample_id",
  all.x = TRUE
)
setorder(mutation_auc_sample_level, mutation, inhibitor, Sample_id)
data.table::fwrite(mutation_auc_sample_level, file.path(output_dir, "mutation_drug_auc_sample_level.csv"))

pair_counts <- mutation_auc_sample_level[
  ,
  .(
    n_samples = uniqueN(Sample_id),
    auc_median = median(auc, na.rm = TRUE),
    auc_mean = mean(auc, na.rm = TRUE),
    lsc17_median = median(LSC17_score, na.rm = TRUE),
    lsc17_mean = mean(LSC17_score, na.rm = TRUE)
  ),
  by = .(mutation, inhibitor)
]
pair_counts[, eligible_auc_correlation := n_samples >= min_n_mutation_drug_auc]
setorder(pair_counts, -n_samples, mutation, inhibitor)
data.table::fwrite(pair_counts, file.path(output_dir, "mutation_drug_auc_sample_counts.csv"))

eligible_pairs <- pair_counts[eligible_auc_correlation == TRUE, .(mutation, inhibitor, n_samples)]
log_msg("Eligible mutation-drug pairs with n >= ", min_n_mutation_drug_auc, ": ", nrow(eligible_pairs))

log_msg("Computing mutation-drug surface gene vs AUC Spearman correlations...")
auc_cor_list <- vector("list", nrow(eligible_pairs))
for (i in seq_len(nrow(eligible_pairs))) {
  current_mutation <- eligible_pairs$mutation[i]
  current_drug <- eligible_pairs$inhibitor[i]
  if (i %% 100 == 0 || i == 1 || i == nrow(eligible_pairs)) {
    log_msg("AUC correlation pair ", i, "/", nrow(eligible_pairs), ": ", current_mutation, " | ", current_drug)
  }

  current_samples <- mutation_auc_sample_level[
    mutation == current_mutation & inhibitor == current_drug,
    .(Sample_id, auc)
  ]
  current_samples <- current_samples[Sample_id %in% colnames(expr_mat)]
  setorder(current_samples, Sample_id)

  mat <- expr_mat[, current_samples$Sample_id, drop = FALSE]
  res <- spearman_matrix_vs_vector(mat, current_samples$auc)
  res$gene <- rownames(mat)
  res$mutation <- current_mutation
  res$inhibitor <- current_drug
  res$n_samples_auc <- nrow(current_samples)
  res$FDR_auc <- p.adjust(res$pvalue, method = "BH")
  res <- res[, c("mutation", "inhibitor", "n_samples_auc", "gene", "rho", "pvalue", "FDR_auc")]
  setnames(res, old = c("rho", "pvalue"), new = c("rho_auc", "pvalue_auc"))
  auc_cor_list[[i]] <- res
}

auc_cor <- data.table::rbindlist(auc_cor_list, use.names = TRUE)
auc_cor[, abs_rho_auc := abs(rho_auc)]
setorder(auc_cor, mutation, inhibitor, pvalue_auc, -abs_rho_auc)
auc_cor[, abs_rho_auc := NULL]
data.table::fwrite(auc_cor, file.path(output_dir, "mutation_drug_surface_gene_auc_spearman_all.csv"))

auc_cor_sig <- auc_cor[FDR_auc < significance_fdr]
data.table::fwrite(auc_cor_sig, file.path(output_dir, "mutation_drug_surface_gene_auc_spearman_significant_FDR05.csv"))

top_auc_by_pair <- auc_cor[
  ,
  .SD[order(-abs(rho_auc), pvalue_auc)][seq_len(min(.N, top_n_auc_genes_per_pair))],
  by = .(mutation, inhibitor)
]
data.table::fwrite(top_auc_by_pair, file.path(output_dir, "mutation_drug_surface_gene_auc_top20_by_abs_rho.csv"))

log_msg("Joining AUC correlations with mutation-subtype LSC17 correlations...")
combined <- merge(
  auc_cor,
  lsc17_mutation_cor,
  by = c("mutation", "gene"),
  all.x = TRUE
)
combined[, lsc17_significant := FDR_lsc17 < significance_fdr]
combined[, auc_significant := FDR_auc < significance_fdr]
combined[, dual_significant := lsc17_significant & auc_significant]
combined[, dual_positive := dual_significant & rho_lsc17 > 0 & rho_auc > 0]
combined[, dual_negative := dual_significant & rho_lsc17 < 0 & rho_auc < 0]
combined[, discordant_significant := dual_significant & sign(rho_lsc17) != sign(rho_auc)]
combined[, combined_direction := fifelse(
  dual_positive,
  "dual_positive",
  fifelse(
    dual_negative,
    "dual_negative",
    fifelse(discordant_significant, "discordant", "not_dual_significant")
  )
)]
combined[, rank_score := abs(rho_lsc17) * abs(rho_auc)]
setorder(combined, mutation, inhibitor, -dual_positive, -dual_significant, -rank_score)
data.table::fwrite(combined, file.path(output_dir, "mutation_drug_surface_gene_auc_lsc17_combined_all.csv"))

dual_sig <- combined[dual_significant == TRUE]
data.table::fwrite(dual_sig, file.path(output_dir, "mutation_drug_surface_gene_dual_significant_FDR05.csv"))

dual_positive <- combined[dual_positive == TRUE]
data.table::fwrite(dual_positive, file.path(output_dir, "mutation_drug_surface_gene_dual_positive_FDR05.csv"))

dual_summary <- combined[
  ,
  .(
    n_samples_auc = dplyr::first(n_samples_auc),
    n_surface_genes_tested = .N,
    n_auc_FDR05 = sum(auc_significant, na.rm = TRUE),
    n_lsc17_FDR05 = sum(lsc17_significant, na.rm = TRUE),
    n_dual_significant = sum(dual_significant, na.rm = TRUE),
    n_dual_positive = sum(dual_positive, na.rm = TRUE),
    n_dual_negative = sum(dual_negative, na.rm = TRUE),
    n_discordant_significant = sum(discordant_significant, na.rm = TRUE),
    top_dual_positive_gene = {
      x <- .SD[dual_positive == TRUE]
      if (nrow(x) == 0) NA_character_ else x$gene[which.max(x$rank_score)]
    },
    top_dual_positive_score = {
      x <- .SD[dual_positive == TRUE]
      if (nrow(x) == 0) NA_real_ else max(x$rank_score, na.rm = TRUE)
    }
  ),
  by = .(mutation, inhibitor)
]
setorder(dual_summary, -n_dual_positive, -n_dual_significant, mutation, inhibitor)
data.table::fwrite(dual_summary, file.path(output_dir, "mutation_drug_dual_correlation_summary.csv"))

mutation_summary <- dual_summary[
  ,
  .(
    n_drugs_tested = .N,
    n_drugs_with_dual_positive = sum(n_dual_positive > 0),
    total_dual_positive_gene_drug_hits = sum(n_dual_positive),
    total_dual_significant_gene_drug_hits = sum(n_dual_significant)
  ),
  by = mutation
]
setorder(mutation_summary, -total_dual_positive_gene_drug_hits, mutation)
data.table::fwrite(mutation_summary, file.path(output_dir, "mutation_dual_correlation_summary.csv"))

drug_summary <- dual_summary[
  ,
  .(
    n_mutations_tested = .N,
    n_mutations_with_dual_positive = sum(n_dual_positive > 0),
    total_dual_positive_gene_mutation_hits = sum(n_dual_positive),
    total_dual_significant_gene_mutation_hits = sum(n_dual_significant)
  ),
  by = inhibitor
]
setorder(drug_summary, -total_dual_positive_gene_mutation_hits, inhibitor)
data.table::fwrite(drug_summary, file.path(output_dir, "drug_dual_correlation_summary.csv"))

log_msg("Creating summary PDF figures...")
top_pairs_plot <- dual_summary[n_dual_positive > 0][1:min(.N, 50)]
if (nrow(top_pairs_plot) > 0) {
  top_pairs_plot[, pair_label := paste(mutation, inhibitor, sep = " | ")]
  top_pairs_plot[, pair_label := factor(pair_label, levels = rev(pair_label))]
  p_pair <- ggplot(top_pairs_plot, aes(x = pair_label, y = n_dual_positive)) +
    geom_col(fill = "#4E79A7", width = 0.72) +
    coord_flip() +
    labs(
      x = NULL,
      y = "Dual-positive surface genes",
      title = "Top mutation-drug pairs with dual-positive LSC17 and AUC correlations"
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.text.y = element_text(size = 5.5)
    )
  ggsave(
    filename = file.path(output_dir, "top50_mutation_drug_dual_positive_counts.pdf"),
    plot = p_pair,
    width = 9,
    height = 8,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
  )
}

top_mutation_plot <- mutation_summary[1:min(.N, 30)]
if (nrow(top_mutation_plot) > 0) {
  top_mutation_plot[, mutation := factor(mutation, levels = rev(mutation))]
  p_mut <- ggplot(top_mutation_plot, aes(x = mutation, y = total_dual_positive_gene_drug_hits)) +
    geom_col(fill = "#B07AA1", width = 0.72) +
    coord_flip() +
    labs(
      x = NULL,
      y = "Total dual-positive gene-drug hits",
      title = "Dual-positive candidates summarized by mutation subtype"
    ) +
    theme_classic(base_size = 8) +
    theme(plot.title = element_text(face = "bold", size = 10))
  ggsave(
    filename = file.path(output_dir, "mutation_dual_positive_summary.pdf"),
    plot = p_mut,
    width = 7,
    height = 5,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
  )
}

run_info <- data.frame(
  item = c(
    "project_root",
    "auc_file",
    "expression_file",
    "lsc17_score_file",
    "mutation_assignment_file",
    "surface_gene_file",
    "lsc17_subtype_correlation_file",
    "output_dir",
    "min_n_mutation_drug_auc",
    "surface_genes_used",
    "mutation_subtypes_after_filter",
    "eligible_mutation_drug_pairs",
    "auc_correlation_rows",
    "combined_rows",
    "dual_significant_rows",
    "dual_positive_rows"
  ),
  value = c(
    project_root,
    input_auc,
    input_expr,
    input_lsc17_score,
    input_mutation_assignments,
    input_surface_genes,
    input_lsc17_subtype_cor,
    output_dir,
    min_n_mutation_drug_auc,
    length(surface_gene_list),
    uniqueN(mutation_assignments$mutation),
    nrow(eligible_pairs),
    nrow(auc_cor),
    nrow(combined),
    nrow(dual_sig),
    nrow(dual_positive)
  )
)
data.table::fwrite(run_info, file.path(output_dir, "pipeline_run_info.csv"))

log_msg("Done.")
log_msg("Outputs written to: ", output_dir)
log_msg("Session info:")
print(sessionInfo())
