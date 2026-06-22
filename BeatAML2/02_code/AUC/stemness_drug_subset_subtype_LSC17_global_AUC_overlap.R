# Combine subtype-specific LSC17 surface-gene correlations with global drug AUC correlations.
# AUC correlations are not stratified by subtype and are restricted to target stemness-skewing resistance drugs.
# Outputs are written to 03_result/LSC17_AUC/stemness_drug_subset_subtype_lsc17_global_auc_overlap.

suppressWarnings(suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
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

significance_fdr <- 0.05
top_n_per_subtype_drug <- 20

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

input_lsc17_subtype_cor <- file.path(project_root, "03_result", "LSC17", "subtype_surface_latest", "subtype_surface_correlation_all.csv")
input_surface_genes <- file.path(project_root, "03_result", "LSC17", "subtype_surface_latest", "surface_gene_list_used.csv")
input_global_auc_cor <- file.path(project_root, "03_result", "LSC17_AUC", "stemness_drug_subset_global", "drug_gene_auc_spearman_min100_all.csv")
input_global_auc_counts <- file.path(project_root, "03_result", "LSC17_AUC", "stemness_drug_subset_global", "drug_gene_auc_expression_overlap_sample_counts.csv")
output_dir <- file.path(project_root, "03_result", "LSC17_AUC", "stemness_drug_subset_subtype_lsc17_global_auc_overlap")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(paste0(...), "\n", sep = "")
}

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path, call. = FALSE)
  }
}

for (path in c(input_lsc17_subtype_cor, input_surface_genes, input_global_auc_cor, input_global_auc_counts)) {
  stop_if_missing(path)
}

log_msg("Reading surface gene list...")
surface_genes <- data.table::fread(input_surface_genes)
if (!"gene" %in% names(surface_genes)) {
  stop("surface_gene_list_used.csv must contain a gene column.", call. = FALSE)
}
surface_gene_list <- unique(na.omit(as.character(surface_genes$gene)))

log_msg("Reading subtype-specific LSC17 correlations...")
lsc17_cor <- data.table::fread(input_lsc17_subtype_cor, nThread = 2)
required_lsc17_cols <- c("subtype_type", "subtype_variable", "subtype_level", "n_samples", "gene", "rho", "pvalue", "FDR_subtype")
missing_lsc17_cols <- setdiff(required_lsc17_cols, names(lsc17_cor))
if (length(missing_lsc17_cols) > 0) {
  stop("LSC17 subtype correlation file is missing columns: ", paste(missing_lsc17_cols, collapse = ", "), call. = FALSE)
}
lsc17_cor <- lsc17_cor[
  gene %in% surface_gene_list,
  .(
    gene,
    subtype_type,
    subtype_variable,
    subtype_level,
    n_samples_lsc17 = n_samples,
    rho_lsc17 = rho,
    pvalue_lsc17 = pvalue,
    FDR_lsc17 = FDR_subtype
  )
]

log_msg("Reading global gene-AUC correlations and subsetting to surface genes...")
auc_cor <- data.table::fread(
  input_global_auc_cor,
  nThread = 2,
  select = c("inhibitor", "n_samples", "gene", "rho", "pvalue", "FDR")
)
required_auc_cols <- c("inhibitor", "n_samples", "gene", "rho", "pvalue", "FDR")
missing_auc_cols <- setdiff(required_auc_cols, names(auc_cor))
if (length(missing_auc_cols) > 0) {
  stop("Global AUC correlation file is missing columns: ", paste(missing_auc_cols, collapse = ", "), call. = FALSE)
}
auc_cor <- auc_cor[
  inhibitor %in% target_drugs & gene %in% surface_gene_list,
  .(
    gene,
    inhibitor,
    n_samples_auc = n_samples,
    rho_auc = rho,
    pvalue_auc = pvalue,
    FDR_auc_all_genes = FDR
  )
]
auc_cor[, FDR_auc_surface := p.adjust(pvalue_auc, method = "BH"), by = inhibitor]
setorder(auc_cor, inhibitor, pvalue_auc)
data.table::fwrite(auc_cor, file.path(output_dir, "surface_gene_global_auc_correlations_min100.csv"))

target_drug_status <- auc_cor[
  ,
  .(n_surface_auc_rows = .N, n_samples_auc = dplyr::first(n_samples_auc)),
  by = inhibitor
][target_drugs_dt, on = "inhibitor"]
target_drug_status[, included_in_global_auc_overlap := !is.na(n_surface_auc_rows)]
target_drug_status[is.na(n_surface_auc_rows), n_surface_auc_rows := 0L]
target_drug_status[is.na(n_samples_auc), n_samples_auc := 0L]
data.table::fwrite(target_drug_status, file.path(output_dir, "target_drug_global_auc_status.csv"))

auc_counts <- data.table::fread(input_global_auc_counts)

log_msg("Combining subtype LSC17 correlations with global AUC correlations by gene...")
combined <- merge(lsc17_cor, auc_cor, by = "gene", allow.cartesian = TRUE)
combined[, lsc17_significant := FDR_lsc17 < significance_fdr]
combined[, auc_significant_all_genes := FDR_auc_all_genes < significance_fdr]
combined[, auc_significant_surface := FDR_auc_surface < significance_fdr]
combined[, dual_significant_all_gene_FDR := lsc17_significant & auc_significant_all_genes]
combined[, dual_positive_all_gene_FDR := dual_significant_all_gene_FDR & rho_lsc17 > 0 & rho_auc > 0]
combined[, dual_negative_all_gene_FDR := dual_significant_all_gene_FDR & rho_lsc17 < 0 & rho_auc < 0]
combined[, discordant_all_gene_FDR := dual_significant_all_gene_FDR & sign(rho_lsc17) != sign(rho_auc)]
combined[, dual_significant_surface_FDR := lsc17_significant & auc_significant_surface]
combined[, dual_positive_surface_FDR := dual_significant_surface_FDR & rho_lsc17 > 0 & rho_auc > 0]
combined[, rank_score := abs(rho_lsc17) * abs(rho_auc)]
combined[, combined_direction_all_gene_FDR := fifelse(
  dual_positive_all_gene_FDR,
  "dual_positive",
  fifelse(
    dual_negative_all_gene_FDR,
    "dual_negative",
    fifelse(discordant_all_gene_FDR, "discordant", "not_dual_significant")
  )
)]
setcolorder(
  combined,
  c(
    "subtype_type",
    "subtype_variable",
    "subtype_level",
    "inhibitor",
    "gene",
    "n_samples_lsc17",
    "n_samples_auc",
    "rho_lsc17",
    "pvalue_lsc17",
    "FDR_lsc17",
    "rho_auc",
    "pvalue_auc",
    "FDR_auc_all_genes",
    "FDR_auc_surface",
    "lsc17_significant",
    "auc_significant_all_genes",
    "auc_significant_surface",
    "dual_significant_all_gene_FDR",
    "dual_positive_all_gene_FDR",
    "dual_negative_all_gene_FDR",
    "discordant_all_gene_FDR",
    "dual_significant_surface_FDR",
    "dual_positive_surface_FDR",
    "combined_direction_all_gene_FDR",
    "rank_score"
  )
)
setorder(combined, subtype_type, subtype_variable, subtype_level, inhibitor, -dual_positive_all_gene_FDR, -dual_significant_all_gene_FDR, -rank_score)
data.table::fwrite(combined, file.path(output_dir, "subtype_lsc17_global_auc_combined_all.csv"))

dual_sig_all <- combined[dual_significant_all_gene_FDR == TRUE]
data.table::fwrite(dual_sig_all, file.path(output_dir, "subtype_lsc17_global_auc_dual_significant_all_gene_FDR05.csv"))

dual_positive_all <- combined[dual_positive_all_gene_FDR == TRUE]
data.table::fwrite(dual_positive_all, file.path(output_dir, "subtype_lsc17_global_auc_dual_positive_all_gene_FDR05.csv"))

dual_positive_surface <- combined[dual_positive_surface_FDR == TRUE]
data.table::fwrite(dual_positive_surface, file.path(output_dir, "subtype_lsc17_global_auc_dual_positive_surface_FDR05.csv"))

top_dual_positive <- dual_positive_all[
  ,
  .SD[order(-rank_score, FDR_lsc17, FDR_auc_all_genes)][seq_len(min(.N, top_n_per_subtype_drug))],
  by = .(subtype_type, subtype_variable, subtype_level, inhibitor)
]
data.table::fwrite(top_dual_positive, file.path(output_dir, "subtype_lsc17_global_auc_top20_dual_positive_by_pair.csv"))

pair_summary <- combined[
  ,
  .(
    n_samples_lsc17 = dplyr::first(n_samples_lsc17),
    n_samples_auc = dplyr::first(n_samples_auc),
    n_surface_genes_tested = .N,
    n_lsc17_FDR05 = sum(lsc17_significant, na.rm = TRUE),
    n_auc_FDR05_all_gene = sum(auc_significant_all_genes, na.rm = TRUE),
    n_auc_FDR05_surface = sum(auc_significant_surface, na.rm = TRUE),
    n_dual_significant_all_gene_FDR = sum(dual_significant_all_gene_FDR, na.rm = TRUE),
    n_dual_positive_all_gene_FDR = sum(dual_positive_all_gene_FDR, na.rm = TRUE),
    n_dual_negative_all_gene_FDR = sum(dual_negative_all_gene_FDR, na.rm = TRUE),
    n_discordant_all_gene_FDR = sum(discordant_all_gene_FDR, na.rm = TRUE),
    n_dual_positive_surface_FDR = sum(dual_positive_surface_FDR, na.rm = TRUE),
    top_dual_positive_gene = {
      x <- .SD[dual_positive_all_gene_FDR == TRUE]
      if (nrow(x) == 0) NA_character_ else x$gene[which.max(x$rank_score)]
    },
    top_dual_positive_score = {
      x <- .SD[dual_positive_all_gene_FDR == TRUE]
      if (nrow(x) == 0) NA_real_ else max(x$rank_score, na.rm = TRUE)
    }
  ),
  by = .(subtype_type, subtype_variable, subtype_level, inhibitor)
]
setorder(pair_summary, -n_dual_positive_all_gene_FDR, -n_dual_significant_all_gene_FDR, subtype_type, subtype_variable, subtype_level, inhibitor)
data.table::fwrite(pair_summary, file.path(output_dir, "subtype_drug_dual_correlation_summary.csv"))

subtype_summary <- pair_summary[
  ,
  .(
    n_drugs_tested = .N,
    n_drugs_with_dual_positive = sum(n_dual_positive_all_gene_FDR > 0),
    total_dual_positive_gene_drug_hits = sum(n_dual_positive_all_gene_FDR),
    total_dual_significant_gene_drug_hits = sum(n_dual_significant_all_gene_FDR)
  ),
  by = .(subtype_type, subtype_variable, subtype_level)
]
setorder(subtype_summary, -total_dual_positive_gene_drug_hits, subtype_type, subtype_variable, subtype_level)
data.table::fwrite(subtype_summary, file.path(output_dir, "subtype_dual_correlation_summary.csv"))

drug_summary <- pair_summary[
  ,
  .(
    n_subtypes_tested = .N,
    n_subtypes_with_dual_positive = sum(n_dual_positive_all_gene_FDR > 0),
    total_dual_positive_gene_subtype_hits = sum(n_dual_positive_all_gene_FDR),
    total_dual_significant_gene_subtype_hits = sum(n_dual_significant_all_gene_FDR)
  ),
  by = inhibitor
]
setorder(drug_summary, -total_dual_positive_gene_subtype_hits, inhibitor)
data.table::fwrite(drug_summary, file.path(output_dir, "drug_dual_correlation_summary.csv"))

gene_summary <- dual_positive_all[
  ,
  .(
    n_subtype_drug_hits = .N,
    n_subtypes = uniqueN(paste(subtype_type, subtype_variable, subtype_level, sep = "||")),
    n_drugs = uniqueN(inhibitor),
    max_rank_score = max(rank_score, na.rm = TRUE),
    max_rho_lsc17 = rho_lsc17[which.max(rank_score)],
    max_rho_auc = rho_auc[which.max(rank_score)]
  ),
  by = gene
]
setorder(gene_summary, -n_subtype_drug_hits, -max_rank_score, gene)
data.table::fwrite(gene_summary, file.path(output_dir, "gene_dual_positive_summary.csv"))

log_msg("Creating summary PDF figures...")
top_pair_plot <- pair_summary[n_dual_positive_all_gene_FDR > 0][1:min(.N, 50)]
if (nrow(top_pair_plot) > 0) {
  top_pair_plot[, pair_label := paste(subtype_variable, subtype_level, inhibitor, sep = " | ")]
  top_pair_plot[, pair_label := factor(pair_label, levels = rev(pair_label))]
  p_pair <- ggplot(top_pair_plot, aes(x = pair_label, y = n_dual_positive_all_gene_FDR, fill = subtype_type)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = c(clinical = "#4E79A7", mutation = "#B07AA1")) +
    coord_flip() +
    labs(
      x = NULL,
      y = "Dual-positive surface genes",
      title = "Top subtype-drug pairs using subtype LSC17 correlation and global AUC correlation"
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.text.y = element_text(size = 5),
      legend.position = "top"
    )
  ggsave(
    filename = file.path(output_dir, "top50_subtype_drug_dual_positive_counts.pdf"),
    plot = p_pair,
    width = 9,
    height = 8,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
  )
}

top_subtype_plot <- subtype_summary[1:min(.N, 40)]
if (nrow(top_subtype_plot) > 0) {
  top_subtype_plot[, subtype_label := paste(subtype_variable, subtype_level, sep = ": ")]
  top_subtype_plot[, subtype_label := factor(subtype_label, levels = rev(subtype_label))]
  p_subtype <- ggplot(top_subtype_plot, aes(x = subtype_label, y = total_dual_positive_gene_drug_hits, fill = subtype_type)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = c(clinical = "#4E79A7", mutation = "#B07AA1")) +
    coord_flip() +
    labs(
      x = NULL,
      y = "Total dual-positive gene-drug hits",
      title = "Dual-positive candidates summarized by subtype"
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.text.y = element_text(size = 5.5),
      legend.position = "top"
    )
  ggsave(
    filename = file.path(output_dir, "subtype_dual_positive_summary.pdf"),
    plot = p_subtype,
    width = 8,
    height = 7,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
  )
}

top_drug_plot <- drug_summary[1:min(.N, 40)]
if (nrow(top_drug_plot) > 0) {
  top_drug_plot[, inhibitor := factor(inhibitor, levels = rev(inhibitor))]
  p_drug <- ggplot(top_drug_plot, aes(x = inhibitor, y = total_dual_positive_gene_subtype_hits)) +
    geom_col(fill = "#59A14F", width = 0.72) +
    coord_flip() +
    labs(
      x = NULL,
      y = "Total dual-positive gene-subtype hits",
      title = "Dual-positive candidates summarized by drug"
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.text.y = element_text(size = 5.5)
    )
  ggsave(
    filename = file.path(output_dir, "drug_dual_positive_summary.pdf"),
    plot = p_drug,
    width = 8,
    height = 7,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
  )
}

run_info <- data.frame(
  item = c(
    "project_root",
    "lsc17_subtype_correlation_file",
    "surface_gene_file",
    "global_auc_correlation_file",
    "global_auc_sample_count_file",
    "output_dir",
    "target_drugs_requested",
    "target_drugs_included",
    "target_drugs_missing_from_global_auc_input",
    "target_drug_list",
    "subtype_groups",
    "drugs",
    "surface_genes",
    "lsc17_correlation_rows",
    "global_auc_surface_rows",
    "combined_rows",
    "dual_significant_all_gene_FDR_rows",
    "dual_positive_all_gene_FDR_rows",
    "dual_positive_surface_FDR_rows"
  ),
  value = c(
    project_root,
    input_lsc17_subtype_cor,
    input_surface_genes,
    input_global_auc_cor,
    input_global_auc_counts,
    output_dir,
    length(target_drugs),
    sum(target_drug_status$included_in_global_auc_overlap),
    sum(!target_drug_status$included_in_global_auc_overlap),
    paste(target_drugs, collapse = "; "),
    uniqueN(paste(lsc17_cor$subtype_type, lsc17_cor$subtype_variable, lsc17_cor$subtype_level, sep = "||")),
    uniqueN(auc_cor$inhibitor),
    length(surface_gene_list),
    nrow(lsc17_cor),
    nrow(auc_cor),
    nrow(combined),
    nrow(dual_sig_all),
    nrow(dual_positive_all),
    nrow(dual_positive_surface)
  )
)
data.table::fwrite(run_info, file.path(output_dir, "pipeline_run_info.csv"))

input_summary <- data.frame(
  item = c(
    "interpretation",
    "auc_scope",
    "lsc17_scope",
    "dual_positive_definition",
    "primary_auc_FDR",
    "target_drug_scope"
  ),
  value = c(
    "Subtype-specific LSC17 association combined with target-drug global, non-subtype-stratified AUC association.",
    "Global AUC correlations from stemness_drug_subset_global/drug_gene_auc_spearman_min100_all.csv; not stratified by subtype.",
    "Subtype-specific LSC17 correlations from subtype_surface_latest; clinical and mutation subtypes included.",
    "FDR_lsc17 < 0.05, FDR_auc_all_genes < 0.05, rho_lsc17 > 0, rho_auc > 0.",
    "FDR_auc_all_genes is used as the primary significance threshold; FDR_auc_surface is also reported.",
    paste(target_drugs, collapse = "; ")
  )
)
data.table::fwrite(input_summary, file.path(output_dir, "analysis_definition.csv"))

log_msg("Done.")
log_msg("Outputs written to: ", output_dir)
log_msg("Session info:")
print(sessionInfo())
