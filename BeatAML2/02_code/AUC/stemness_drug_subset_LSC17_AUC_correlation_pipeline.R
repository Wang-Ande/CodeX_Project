# LSC17 and gene-expression correlation with inhibitor AUC.
# Restricted to drugs reported as resistance-associated with stemness skewing.
# Outputs are written to 03_result/LSC17_AUC/stemness_drug_subset_global.

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

min_n_lsc17_auc <- 10
min_n_gene_auc <- 100
top_n_genes_per_drug <- 50
significance_fdr <- 0.05

target_drugs <- c(
  "Panobinostat",
  "SNS-032 (BMS-387032)",
  "AT7519",
  "Selumetinib (AZD6244)",
  "Trametinib (GSK1120212)",
  "17-AAG (Tanespimycin)",
  "Idelalisib",
  "CI-1040 (PD184352)",
  "GDC-0941",
  "Bortezomib (Velcade)",
  "PHT-427",
  "Linifanib (ABT-869)",
  "INK-128",
  "Axitinib (AG-013736)",
  "Doramapimod (BIRB 796)",
  "Flavopiridol",
  "MK-2206",
  "OTX-015",
  "Rapamycin",
  "Cediranib (AZD2171)",
  "Nilotinib",
  "Dasatinib",
  "Motesanib (AMG-706)"
)
target_drugs_dt <- data.table::data.table(
  inhibitor = target_drugs,
  target_order = seq_along(target_drugs)
)

input_auc <- file.path(project_root, "01_data", "Inhibitor AUC values.csv")
input_expr <- file.path(project_root, "01_data", "Normalized Expression.csv")
input_lsc17_score <- file.path(project_root, "03_result", "LSC17", "LSC17_sample_score_df.csv")
output_dir <- file.path(project_root, "03_result", "LSC17_AUC", "stemness_drug_subset_global")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(paste0(...), "\n", sep = "")
}

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path, call. = FALSE)
  }
}

stop_if_missing(input_auc)
stop_if_missing(input_expr)
stop_if_missing(input_lsc17_score)

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

log_msg("Reading LSC17 sample score table...")
lsc17_score <- data.table::fread(input_lsc17_score)
required_score_cols <- c("Sample_id", "LSC17_score")
missing_score_cols <- setdiff(required_score_cols, names(lsc17_score))
if (length(missing_score_cols) > 0) {
  stop("LSC17 score table is missing columns: ", paste(missing_score_cols, collapse = ", "), call. = FALSE)
}
lsc17_score <- lsc17_score[, .(Sample_id, LSC17_score)]
lsc17_score <- lsc17_score[!is.na(Sample_id) & Sample_id != "" & !is.na(LSC17_score)]
lsc17_score <- unique(lsc17_score, by = "Sample_id")

log_msg("Reading inhibitor AUC table...")
auc_dt <- data.table::fread(input_auc, nThread = 2)
required_auc_cols <- c("dbgap_rnaseq_sample", "inhibitor", "auc")
missing_auc_cols <- setdiff(required_auc_cols, names(auc_dt))
if (length(missing_auc_cols) > 0) {
  stop("AUC table is missing columns: ", paste(missing_auc_cols, collapse = ", "), call. = FALSE)
}

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

target_presence <- auc_clean[
  inhibitor %in% target_drugs,
  .(
    n_auc_samples = uniqueN(Sample_id),
    n_auc_rows = .N
  ),
  by = inhibitor
][target_drugs_dt, on = "inhibitor"]
target_presence[, in_auc := !is.na(n_auc_rows)]
target_presence[is.na(n_auc_samples), n_auc_samples := 0L]
target_presence[is.na(n_auc_rows), n_auc_rows := 0L]
setcolorder(target_presence, c("target_order", "inhibitor", "in_auc", "n_auc_samples", "n_auc_rows"))
data.table::fwrite(target_presence, file.path(output_dir, "target_drug_auc_presence.csv"))

missing_target_drugs <- target_presence[in_auc == FALSE, inhibitor]
if (length(missing_target_drugs) > 0) {
  log_msg("Target drugs missing from AUC table: ", paste(missing_target_drugs, collapse = "; "))
}

auc_clean <- auc_clean[inhibitor %in% target_drugs]
if (nrow(auc_clean) == 0) {
  stop("No target drugs are present in the cleaned AUC table.", call. = FALSE)
}

auc_lsc17 <- merge(auc_clean, lsc17_score, by = "Sample_id", all = FALSE)

drug_counts <- auc_clean[
  ,
  .(
    n_auc_samples = uniqueN(Sample_id),
    n_auc_rows = .N
  ),
  by = inhibitor
][
  auc_lsc17[, .(n_auc_lsc17_overlap = uniqueN(Sample_id)), by = inhibitor],
  on = "inhibitor"
]
drug_counts[is.na(n_auc_lsc17_overlap), n_auc_lsc17_overlap := 0L]
drug_counts[, eligible_lsc17_auc := n_auc_lsc17_overlap >= min_n_lsc17_auc]
drug_counts[, eligible_gene_auc := n_auc_lsc17_overlap >= min_n_gene_auc]
setorder(drug_counts, -n_auc_lsc17_overlap, inhibitor)
data.table::fwrite(drug_counts, file.path(output_dir, "drug_auc_lsc17_overlap_sample_counts.csv"))

eligible_lsc17_drugs <- drug_counts[eligible_lsc17_auc == TRUE, inhibitor]
eligible_gene_drugs <- drug_counts[eligible_gene_auc == TRUE, inhibitor]

log_msg("Eligible drugs for LSC17-AUC correlation (n >= ", min_n_lsc17_auc, "): ", length(eligible_lsc17_drugs))
log_msg("Eligible drugs for gene-AUC correlation (n >= ", min_n_gene_auc, "): ", length(eligible_gene_drugs))

sample_level_lsc17 <- auc_lsc17[inhibitor %in% eligible_lsc17_drugs]
setorder(sample_level_lsc17, inhibitor, Sample_id)
data.table::fwrite(sample_level_lsc17, file.path(output_dir, "drug_auc_lsc17_sample_level_min10.csv"))

log_msg("Computing per-drug LSC17 score vs AUC Spearman correlations...")
lsc17_auc_cor <- sample_level_lsc17[
  ,
  {
    test <- suppressWarnings(cor.test(LSC17_score, auc, method = "spearman", exact = FALSE))
    .(
      n_samples = uniqueN(Sample_id),
      rho = unname(test$estimate),
      pvalue = test$p.value,
      auc_median = median(auc, na.rm = TRUE),
      auc_mean = mean(auc, na.rm = TRUE),
      lsc17_median = median(LSC17_score, na.rm = TRUE),
      lsc17_mean = mean(LSC17_score, na.rm = TRUE)
    )
  },
  by = inhibitor
]
lsc17_auc_cor[, FDR := p.adjust(pvalue, method = "BH")]
lsc17_auc_cor[, abs_rho := abs(rho)]
setorder(lsc17_auc_cor, pvalue, -abs_rho)
lsc17_auc_cor[, abs_rho := NULL]
data.table::fwrite(lsc17_auc_cor, file.path(output_dir, "drug_lsc17_auc_spearman_min10.csv"))

log_msg("Reading normalized expression matrix...")
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

needed_samples <- intersect(unique(auc_lsc17[inhibitor %in% eligible_gene_drugs, Sample_id]), names(expr_dt))
if (length(needed_samples) == 0) {
  stop("No expression samples overlap with eligible gene-AUC drug subsets.", call. = FALSE)
}
expr_mat <- as.matrix(expr_dt[, ..needed_samples])
storage.mode(expr_mat) <- "numeric"
rownames(expr_mat) <- gene_names
rm(expr_dt)

gene_drug_sample_counts <- auc_lsc17[
  inhibitor %in% eligible_gene_drugs,
  .(n_auc_lsc17_samples = uniqueN(Sample_id)),
  by = inhibitor
]
gene_drug_sample_counts[
  ,
  n_auc_lsc17_expression_samples := vapply(
    inhibitor,
    function(drug) {
      length(intersect(auc_lsc17[inhibitor == drug, Sample_id], colnames(expr_mat)))
    },
    integer(1)
  )
]
gene_drugs_final <- gene_drug_sample_counts[n_auc_lsc17_expression_samples >= min_n_gene_auc, inhibitor]
data.table::fwrite(gene_drug_sample_counts, file.path(output_dir, "drug_gene_auc_expression_overlap_sample_counts.csv"))

log_msg("Final drugs for gene-AUC correlation after expression overlap: ", length(gene_drugs_final))

gene_cor_list <- vector("list", length(gene_drugs_final))
names(gene_cor_list) <- gene_drugs_final

for (i in seq_along(gene_drugs_final)) {
  drug <- gene_drugs_final[i]
  log_msg("Gene-AUC correlation: ", i, "/", length(gene_drugs_final), " ", drug)
  drug_auc <- auc_lsc17[inhibitor == drug, .(Sample_id, auc)]
  drug_auc <- drug_auc[Sample_id %in% colnames(expr_mat)]
  setorder(drug_auc, Sample_id)
  samples <- drug_auc$Sample_id
  mat <- expr_mat[, samples, drop = FALSE]
  res <- spearman_matrix_vs_vector(mat, drug_auc$auc)
  res$gene <- rownames(mat)
  res$inhibitor <- drug
  res$n_samples <- length(samples)
  res$FDR <- p.adjust(res$pvalue, method = "BH")
  res <- res[, c("inhibitor", "n_samples", "gene", "rho", "pvalue", "FDR")]
  gene_cor_list[[i]] <- res
}

gene_auc_cor <- data.table::rbindlist(gene_cor_list, use.names = TRUE)
gene_auc_cor[, abs_rho := abs(rho)]
setorder(gene_auc_cor, inhibitor, pvalue, -abs_rho)
gene_auc_cor[, abs_rho := NULL]
data.table::fwrite(gene_auc_cor, file.path(output_dir, "drug_gene_auc_spearman_min100_all.csv"))

gene_auc_sig <- gene_auc_cor[FDR < significance_fdr]
data.table::fwrite(gene_auc_sig, file.path(output_dir, "drug_gene_auc_spearman_min100_significant_FDR05.csv"))

top_gene_auc <- gene_auc_cor[
  ,
  .SD[order(-abs(rho), pvalue)][seq_len(min(.N, top_n_genes_per_drug))],
  by = inhibitor
]
data.table::fwrite(top_gene_auc, file.path(output_dir, "drug_gene_auc_top50_by_abs_rho_min100.csv"))

gene_auc_summary <- gene_auc_cor[
  ,
  .(
    n_samples = dplyr::first(n_samples),
    n_genes_tested = .N,
    n_FDR05 = sum(FDR < significance_fdr, na.rm = TRUE),
    n_FDR05_positive = sum(FDR < significance_fdr & rho > 0, na.rm = TRUE),
    n_FDR05_negative = sum(FDR < significance_fdr & rho < 0, na.rm = TRUE),
    top_abs_gene = gene[which.max(abs(rho))],
    top_abs_rho = rho[which.max(abs(rho))],
    top_abs_FDR = FDR[which.max(abs(rho))]
  ),
  by = inhibitor
]
gene_auc_summary[, abs_top_abs_rho := abs(top_abs_rho)]
setorder(gene_auc_summary, -n_FDR05, -abs_top_abs_rho, inhibitor)
gene_auc_summary[, abs_top_abs_rho := NULL]
data.table::fwrite(gene_auc_summary, file.path(output_dir, "drug_gene_auc_summary_min100.csv"))

log_msg("Creating summary PDF figures...")

if (nrow(lsc17_auc_cor) > 0) {
  plot_df <- copy(lsc17_auc_cor)
  plot_df[, inhibitor := factor(inhibitor, levels = inhibitor[order(rho)])]
  p_lsc17 <- ggplot(plot_df, aes(x = inhibitor, y = rho, color = FDR < significance_fdr)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_point(size = 1.5) +
    coord_flip() +
    scale_color_manual(values = c("FALSE" = "grey45", "TRUE" = "#C94C4C"), name = "FDR < 0.05") +
    labs(
      x = NULL,
      y = "Spearman rho: LSC17 score vs AUC",
      title = "Per-drug LSC17 score association with inhibitor AUC"
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.text.y = element_text(size = 5)
    )
  ggsave(
    filename = file.path(output_dir, "drug_lsc17_auc_spearman_min10.pdf"),
    plot = p_lsc17,
    width = 8,
    height = max(6, 0.12 * nrow(plot_df)),
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
  )
}

if (nrow(gene_auc_summary) > 0) {
  summary_plot <- copy(gene_auc_summary)
  summary_plot[, inhibitor := factor(inhibitor, levels = inhibitor[order(n_FDR05)])]
  p_gene_summary <- ggplot(summary_plot, aes(x = inhibitor, y = n_FDR05)) +
    geom_col(fill = "#4E79A7", width = 0.72) +
    coord_flip() +
    labs(
      x = NULL,
      y = "Genes with FDR < 0.05",
      title = "Gene expression associations with inhibitor AUC by drug"
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.text.y = element_text(size = 5)
    )
  ggsave(
    filename = file.path(output_dir, "drug_gene_auc_significant_gene_counts_min100.pdf"),
    plot = p_gene_summary,
    width = 8,
    height = max(5, 0.18 * nrow(summary_plot)),
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
    "output_dir",
    "target_drugs_requested",
    "target_drugs_found_in_auc",
    "target_drugs_missing_in_auc",
    "target_drug_list",
    "min_n_lsc17_auc",
    "min_n_gene_auc",
    "n_lsc17_score_samples",
    "n_auc_rows_nonmissing",
    "n_unique_drugs_with_auc",
    "n_eligible_lsc17_auc_drugs",
    "n_eligible_gene_auc_drugs",
    "n_final_gene_auc_drugs",
    "n_genes_tested"
  ),
  value = c(
    project_root,
    input_auc,
    input_expr,
    input_lsc17_score,
    output_dir,
    length(target_drugs),
    sum(target_presence$in_auc),
    length(missing_target_drugs),
    paste(target_drugs, collapse = "; "),
    min_n_lsc17_auc,
    min_n_gene_auc,
    nrow(lsc17_score),
    nrow(auc_clean),
    uniqueN(auc_clean$inhibitor),
    length(eligible_lsc17_drugs),
    length(eligible_gene_drugs),
    length(gene_drugs_final),
    length(gene_names)
  )
)
data.table::fwrite(run_info, file.path(output_dir, "pipeline_run_info.csv"))

log_msg("Done.")
log_msg("Outputs written to: ", output_dir)
log_msg("Session info:")
print(sessionInfo())
