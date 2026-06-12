# LSC17 subtype-stratified surface gene correlation pipeline
# Outputs are written to 03_result/LSC17/subtype_surface.

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(matrixStats)
})

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

min_subtype_n <- 20
top_n_per_group <- 20
top_n_heatmap_genes <- 50
top_n_scatter_pairs <- 8
significance_fdr <- 0.05

input_clinical <- file.path(project_root, "01_data", "Clinical Summary.xlsx")
input_expression <- file.path(project_root, "01_data", "Normalized Expression.csv")
input_surface_gene_source <- file.path(project_root, "03_result", "LSC17", "LASA_LSC17_cor_df.csv")
output_dir <- file.path(project_root, "03_result", "LSC17", "subtype_surface")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path, call. = FALSE)
  }
}

stop_if_missing(input_clinical)
stop_if_missing(input_expression)
stop_if_missing(input_surface_gene_source)

clean_label <- function(x, drop_unknown = TRUE) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  if (drop_unknown) {
    x[tolower(x) %in% c("unknown", "n/a", "na", "missing", "<na>")] <- NA_character_
  }
  x
}

mutation_status <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  ifelse(is.na(x) | x == "" | tolower(x) %in% c("na", "n/a", "<na>"),
         "not_detected",
         "mutated")
}

positive_negative_status <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[tolower(x) %in% c("", "na", "n/a", "<na>", "unknown")] <- NA_character_
  x
}

disease_origin <- function(is_denovo, is_transformed) {
  is_denovo <- tolower(as.character(is_denovo))
  is_transformed <- tolower(as.character(is_transformed))
  out <- rep(NA_character_, length(is_denovo))
  out[is_denovo == "true"] <- "de_novo"
  out[is_transformed == "true"] <- "transformed"
  out[is_denovo == "false" & is_transformed == "false"] <- "other_or_treated"
  out
}

message("Reading clinical data...")
clinical <- openxlsx::read.xlsx(input_clinical)

aml_dx <- "ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS"
sample_info <- clinical %>%
  filter(
    .data[["dxAtInclusion"]] == aml_dx,
    !is.na(.data[["dbgap_rnaseq_sample"]]),
    .data[["dbgap_rnaseq_sample"]] != ""
  ) %>%
  distinct(.data[["dbgap_rnaseq_sample"]], .keep_all = TRUE) %>%
  transmute(
    sample_id = .data[["dbgap_rnaseq_sample"]],
    WHO_specificDx = clean_label(.data[["specificDxAtInclusion"]]),
    ELN2017 = clean_label(.data[["ELN2017"]]),
    Fusion_subtype = clean_label(.data[["consensusAMLFusions"]]),
    FAB_subtype = clean_label(.data[["fabBlastMorphology"]]),
    FLT3_ITD = positive_negative_status(.data[["FLT3-ITD"]]),
    NPM1 = positive_negative_status(.data[["NPM1"]]),
    RUNX1 = mutation_status(.data[["RUNX1"]]),
    ASXL1 = mutation_status(.data[["ASXL1"]]),
    TP53 = mutation_status(.data[["TP53"]]),
    Disease_origin = disease_origin(.data[["isDenovo"]], .data[["isTransformed"]]),
    Disease_stage = clean_label(.data[["diseaseStageAtSpecimenCollection"]]),
    Specimen_type = clean_label(.data[["specimenType"]], drop_unknown = FALSE)
  )

message("Reading expression matrix...")
expr_dt <- data.table::fread(input_expression, nThread = 2)
required_expr_cols <- c("display_label")
missing_expr_cols <- setdiff(required_expr_cols, names(expr_dt))
if (length(missing_expr_cols) > 0) {
  stop("Expression file is missing columns: ", paste(missing_expr_cols, collapse = ", "), call. = FALSE)
}

gene_names <- expr_dt[["display_label"]]
if (any(is.na(gene_names) | gene_names == "")) {
  stop("Expression file contains missing display_label values.", call. = FALSE)
}
if (anyDuplicated(gene_names) > 0) {
  stop("Expression file contains duplicated display_label values. Resolve duplicates before running.", call. = FALSE)
}

sample_cols <- intersect(sample_info$sample_id, names(expr_dt))
if (length(sample_cols) < min_subtype_n) {
  stop("Too few AML RNA-seq samples overlap with the expression matrix.", call. = FALSE)
}

expr_mat <- as.matrix(expr_dt[, ..sample_cols])
storage.mode(expr_mat) <- "numeric"
rownames(expr_mat) <- gene_names
rm(expr_dt)

sample_info <- sample_info %>%
  filter(.data[["sample_id"]] %in% colnames(expr_mat))

expr_mat <- expr_mat[, sample_info$sample_id, drop = FALSE]

weights_LSC17 <- c(
  GPR56 = 0.0501,
  AKR1C3 = -0.0402,
  CD34 = 0.0338,
  NGFRAP1 = 0.0465,
  EMP1 = 0.0146,
  C19orf77 = -0.0226,
  SOCS2 = 0.0271,
  CPXM1 = -0.0258,
  CDK6 = -0.0704,
  KIAA0125 = 0.0196,
  DPYSL3 = 0.0284,
  MMRN1 = 0.0258,
  LAPTM4B = 0.00582,
  ARHGAP22 = -0.0138,
  NYNRIN = 0.00865,
  ZBTB46 = -0.0347,
  DNMT3B = 0.0874
)

missing_lsc17 <- setdiff(names(weights_LSC17), rownames(expr_mat))
if (length(missing_lsc17) > 0) {
  stop("Missing LSC17 genes: ", paste(missing_lsc17, collapse = ", "), call. = FALSE)
}

sample_info$lsc17_score <- colSums(
  sweep(expr_mat[names(weights_LSC17), , drop = FALSE], 1, weights_LSC17, `*`)
)

message("Reading surface gene list from previous project-local LASA result...")
surface_source <- data.table::fread(input_surface_gene_source)
if (!"gene" %in% names(surface_source)) {
  stop("Surface gene source must contain a 'gene' column: ", input_surface_gene_source, call. = FALSE)
}

surface_gene_list <- unique(na.omit(as.character(surface_source$gene)))
surface_gene_list <- surface_gene_list[surface_gene_list %in% rownames(expr_mat)]
if (length(surface_gene_list) == 0) {
  stop("No surface genes overlap with the expression matrix.", call. = FALSE)
}

data.table::fwrite(
  data.frame(gene = surface_gene_list),
  file.path(output_dir, "surface_gene_list_used.csv")
)

sample_scores_out <- sample_info %>%
  arrange(desc(.data[["lsc17_score"]]))
data.table::fwrite(sample_scores_out, file.path(output_dir, "lsc17_sample_scores_with_subtypes.csv"))

subtype_variables <- c(
  "WHO_specificDx",
  "ELN2017",
  "Fusion_subtype",
  "FAB_subtype",
  "FLT3_ITD",
  "NPM1",
  "RUNX1",
  "ASXL1",
  "TP53",
  "Disease_origin",
  "Disease_stage",
  "Specimen_type"
)

subtype_assignments <- bind_rows(lapply(subtype_variables, function(v) {
  data.frame(
    sample_id = sample_info$sample_id,
    subtype_variable = v,
    subtype_level = sample_info[[v]],
    stringsAsFactors = FALSE
  )
})) %>%
  filter(!is.na(.data[["subtype_level"]]))

subtype_counts <- subtype_assignments %>%
  dplyr::count(subtype_variable, subtype_level, name = "n_samples") %>%
  arrange(subtype_variable, desc(n_samples), subtype_level)

eligible_subtypes <- subtype_counts %>%
  filter(n_samples >= min_subtype_n)

skipped_subtypes <- subtype_counts %>%
  filter(n_samples < min_subtype_n)

data.table::fwrite(subtype_counts, file.path(output_dir, "subtype_sample_counts.csv"))
data.table::fwrite(eligible_subtypes, file.path(output_dir, "subtype_eligible_groups.csv"))
data.table::fwrite(skipped_subtypes, file.path(output_dir, "subtype_skipped_groups.csv"))

surface_expr <- expr_mat[surface_gene_list, , drop = FALSE]
rm(expr_mat)

correlate_surface_genes <- function(sample_ids, subtype_variable, subtype_level) {
  sample_ids <- intersect(sample_ids, colnames(surface_expr))
  score <- sample_info$lsc17_score[match(sample_ids, sample_info$sample_id)]
  keep <- !is.na(score)
  sample_ids <- sample_ids[keep]
  score <- score[keep]
  n <- length(sample_ids)

  mat <- surface_expr[, sample_ids, drop = FALSE]
  expr_rank <- matrixStats::rowRanks(mat, ties.method = "average", preserveShape = TRUE)
  score_rank <- rank(score, ties.method = "average")
  rho <- as.vector(cor(t(expr_rank), score_rank, method = "pearson"))

  denom <- pmax(1 - rho^2, .Machine$double.eps)
  t_stat <- rho * sqrt((n - 2) / denom)
  pvalue <- 2 * pt(-abs(t_stat), df = n - 2)
  pvalue[is.na(rho)] <- NA_real_

  data.frame(
    subtype_variable = subtype_variable,
    subtype_level = subtype_level,
    n_samples = n,
    gene = rownames(mat),
    rho = rho,
    pvalue = pvalue,
    FDR_subtype = p.adjust(pvalue, method = "BH"),
    stringsAsFactors = FALSE
  )
}

message("Running subtype-stratified Spearman correlations...")
cor_list <- vector("list", nrow(eligible_subtypes))
for (i in seq_len(nrow(eligible_subtypes))) {
  current_var <- eligible_subtypes$subtype_variable[i]
  current_level <- eligible_subtypes$subtype_level[i]
  current_samples <- subtype_assignments %>%
    filter(
      subtype_variable == current_var,
      subtype_level == current_level
    ) %>%
    pull(sample_id)

  cor_list[[i]] <- correlate_surface_genes(
    sample_ids = current_samples,
    subtype_variable = current_var,
    subtype_level = current_level
  )
}

cor_df <- bind_rows(cor_list) %>%
  mutate(
    FDR_global = p.adjust(pvalue, method = "BH"),
    direction = case_when(
      rho > 0 ~ "positive",
      rho < 0 ~ "negative",
      TRUE ~ "zero"
    )
  ) %>%
  arrange(subtype_variable, subtype_level, desc(abs(rho)))

data.table::fwrite(cor_df, file.path(output_dir, "subtype_surface_correlation_all.csv"))

sig_df <- cor_df %>%
  filter(FDR_subtype < significance_fdr) %>%
  arrange(subtype_variable, subtype_level, FDR_subtype, desc(abs(rho)))

data.table::fwrite(sig_df, file.path(output_dir, "subtype_surface_correlation_significant_FDR05.csv"))

top_by_abs <- cor_df %>%
  group_by(subtype_variable, subtype_level) %>%
  slice_max(order_by = abs(rho), n = top_n_per_group, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(subtype_variable, subtype_level, desc(abs(rho)))

data.table::fwrite(top_by_abs, file.path(output_dir, "subtype_top_surface_genes_by_abs_rho.csv"))

sig_counts <- cor_df %>%
  mutate(significant = FDR_subtype < significance_fdr) %>%
  filter(significant) %>%
  dplyr::count(
    subtype_variable,
    subtype_level,
    direction,
    name = "n_significant_genes"
  )

sig_counts_complete <- eligible_subtypes %>%
  crossing(direction = c("positive", "negative")) %>%
  left_join(sig_counts, by = c("subtype_variable", "subtype_level", "direction")) %>%
  mutate(n_significant_genes = replace_na(n_significant_genes, 0L))

data.table::fwrite(sig_counts_complete, file.path(output_dir, "subtype_significant_gene_counts.csv"))

summary_by_group <- cor_df %>%
  group_by(subtype_variable, subtype_level, n_samples) %>%
  summarise(
    n_surface_genes_tested = n(),
    n_FDR05 = sum(FDR_subtype < significance_fdr, na.rm = TRUE),
    n_FDR05_positive = sum(FDR_subtype < significance_fdr & rho > 0, na.rm = TRUE),
    n_FDR05_negative = sum(FDR_subtype < significance_fdr & rho < 0, na.rm = TRUE),
    max_abs_rho = max(abs(rho), na.rm = TRUE),
    top_abs_gene = gene[which.max(abs(rho))],
    top_abs_rho = rho[which.max(abs(rho))],
    top_abs_FDR = FDR_subtype[which.max(abs(rho))],
    .groups = "drop"
  ) %>%
  arrange(desc(n_FDR05), desc(max_abs_rho))

data.table::fwrite(summary_by_group, file.path(output_dir, "subtype_correlation_summary_by_group.csv"))

run_info <- data.frame(
  item = c(
    "project_root",
    "clinical_file",
    "expression_file",
    "surface_gene_source",
    "output_dir",
    "min_subtype_n",
    "n_aml_rnaseq_samples_in_clinical",
    "n_aml_expression_overlap_samples",
    "n_surface_genes_tested",
    "n_eligible_subtype_groups",
    "n_skipped_subtype_groups"
  ),
  value = c(
    project_root,
    input_clinical,
    input_expression,
    input_surface_gene_source,
    output_dir,
    min_subtype_n,
    nrow(clinical %>% filter(.data[["dxAtInclusion"]] == aml_dx, !is.na(.data[["dbgap_rnaseq_sample"]]))),
    nrow(sample_info),
    length(surface_gene_list),
    nrow(eligible_subtypes),
    nrow(skipped_subtypes)
  )
)
data.table::fwrite(run_info, file.path(output_dir, "pipeline_run_info.csv"))

message("Creating PDF figures...")

plot_group_label <- function(variable, level) {
  paste(variable, level, sep = ": ")
}

summary_plot_df <- summary_by_group %>%
  mutate(group_label = plot_group_label(subtype_variable, subtype_level)) %>%
  arrange(n_FDR05) %>%
  mutate(group_label = factor(group_label, levels = unique(group_label)))

p_counts <- ggplot(summary_plot_df, aes(x = group_label, y = n_FDR05)) +
  geom_col(fill = "#4E79A7", width = 0.72) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Surface genes with subtype FDR < 0.05",
    title = "Significant LSC17-associated surface genes by subtype"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 7)
  )

ggsave(
  filename = file.path(output_dir, "subtype_significant_gene_counts.pdf"),
  plot = p_counts,
  width = 8,
  height = max(5, 0.18 * nrow(summary_plot_df)),
  units = "in"
)

heatmap_genes <- sig_df %>%
  group_by(gene) %>%
  summarise(max_abs_rho = max(abs(rho), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_abs_rho)) %>%
  slice_head(n = top_n_heatmap_genes) %>%
  pull(gene)

if (length(heatmap_genes) == 0) {
  heatmap_genes <- cor_df %>%
    group_by(gene) %>%
    summarise(max_abs_rho = max(abs(rho), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_abs_rho)) %>%
    slice_head(n = top_n_heatmap_genes) %>%
    pull(gene)
}

heatmap_df <- cor_df %>%
  filter(gene %in% heatmap_genes) %>%
  mutate(
    group_label = plot_group_label(subtype_variable, subtype_level),
    gene = factor(gene, levels = rev(heatmap_genes))
  )

group_order <- summary_by_group %>%
  mutate(group_label = plot_group_label(subtype_variable, subtype_level)) %>%
  arrange(subtype_variable, desc(n_samples), subtype_level) %>%
  pull(group_label)

heatmap_df$group_label <- factor(heatmap_df$group_label, levels = group_order)

p_heatmap <- ggplot(heatmap_df, aes(x = group_label, y = gene, fill = rho)) +
  geom_tile(color = "grey88", linewidth = 0.15) +
  scale_fill_gradient2(
    low = "#3B7CB7",
    mid = "white",
    high = "#C94C4C",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman rho"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Top subtype-specific surface gene correlations with LSC17 score"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank()
  )

ggsave(
  filename = file.path(output_dir, "subtype_top_correlation_heatmap.pdf"),
  plot = p_heatmap,
  width = max(10, 0.24 * length(group_order)),
  height = max(7, 0.16 * length(heatmap_genes)),
  units = "in",
  limitsize = FALSE
)

score_long <- subtype_assignments %>%
  left_join(sample_info %>% select(sample_id, lsc17_score), by = "sample_id") %>%
  inner_join(eligible_subtypes, by = c("subtype_variable", "subtype_level")) %>%
  group_by(subtype_variable) %>%
  mutate(subtype_level = reorder(subtype_level, lsc17_score, FUN = median, na.rm = TRUE)) %>%
  ungroup()

p_score <- ggplot(score_long, aes(x = subtype_level, y = lsc17_score)) +
  geom_boxplot(outlier.shape = NA, fill = "#A0CBE8", color = "grey30", width = 0.6) +
  geom_jitter(width = 0.12, size = 0.45, alpha = 0.45, color = "grey25") +
  facet_wrap(~ subtype_variable, scales = "free_x", ncol = 2) +
  labs(
    x = NULL,
    y = "LSC17 score",
    title = "LSC17 score distribution across eligible subtype groups"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    strip.text = element_text(face = "bold", size = 8)
  )

ggsave(
  filename = file.path(output_dir, "subtype_lsc17_score_distribution.pdf"),
  plot = p_score,
  width = 10,
  height = 8,
  units = "in"
)

scatter_pairs <- cor_df %>%
  filter(FDR_subtype < significance_fdr) %>%
  arrange(desc(abs(rho))) %>%
  distinct(subtype_variable, subtype_level, gene, .keep_all = TRUE) %>%
  slice_head(n = top_n_scatter_pairs)

if (nrow(scatter_pairs) > 0) {
  scatter_df <- bind_rows(lapply(seq_len(nrow(scatter_pairs)), function(i) {
    row <- scatter_pairs[i, ]
    samples <- subtype_assignments %>%
      filter(
        subtype_variable == row$subtype_variable,
        subtype_level == row$subtype_level
      ) %>%
      pull(sample_id)
    samples <- intersect(samples, colnames(surface_expr))
    data.frame(
      sample_id = samples,
      gene = row$gene,
      subtype_variable = row$subtype_variable,
      subtype_level = row$subtype_level,
      rho = row$rho,
      FDR_subtype = row$FDR_subtype,
      expression = as.numeric(surface_expr[row$gene, samples]),
      lsc17_score = sample_info$lsc17_score[match(samples, sample_info$sample_id)],
      stringsAsFactors = FALSE
    )
  })) %>%
    mutate(
      facet_label = paste0(
        subtype_variable, ": ", subtype_level, "\n",
        gene, " rho=", sprintf("%.2f", rho), ", FDR=", format(FDR_subtype, digits = 2, scientific = TRUE)
      )
    )

  p_scatter <- ggplot(scatter_df, aes(x = expression, y = lsc17_score)) +
    geom_point(size = 1.4, alpha = 0.75, color = "#4E79A7") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.45, color = "#C94C4C") +
    facet_wrap(~ facet_label, scales = "free_x", ncol = 2) +
    labs(
      x = "Surface gene expression",
      y = "LSC17 score",
      title = "Top subtype-specific gene-LSC17 associations"
    ) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      strip.text = element_text(size = 7)
    )

  ggsave(
    filename = file.path(output_dir, "subtype_top_gene_scatterplots.pdf"),
    plot = p_scatter,
    width = 8,
    height = max(5, 2.2 * ceiling(nrow(scatter_pairs) / 2)),
    units = "in"
  )
}

message("Done. Outputs written to: ", output_dir)
