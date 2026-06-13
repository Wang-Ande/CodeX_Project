# LSC17 subtype-stratified surface gene correlation pipeline.
# Mutation subtypes integrate variantSummary tokens with structured mutation columns.
# Outputs are written to 03_result/LSC17/subtype_surface_latest.

suppressWarnings(suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(matrixStats)
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

min_subtype_n <- 10
top_n_per_group <- 20
top_n_global_heatmap_genes <- 50
top_n_scatter_pairs <- 8
top_n_mutation_heatmap <- 10
top_n_clinical_heatmap <- 20
significance_fdr <- 0.05

input_clinical <- file.path(project_root, "01_data", "Clinical Summary.xlsx")
input_expression <- file.path(project_root, "01_data", "Normalized Expression.csv")
input_surface_gene_source <- file.path(project_root, "03_result", "LSC17", "LASA_LSC17_cor_df.csv")
output_dir <- file.path(project_root, "03_result", "LSC17", "subtype_surface_latest")
output_filtered_clinical <- file.path(project_root, "01_data", "LSC17_filtered_clinical_integrated_mutation.csv")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(paste0(...), "\n", sep = "")
}

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path, call. = FALSE)
  }
}

clean_label <- function(x, drop_unknown = TRUE) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  if (drop_unknown) {
    x[tolower(x) %in% c("unknown", "n/a", "na", "missing", "<na>")] <- NA_character_
  }
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

parse_variant_tokens_one <- function(x) {
  if (length(x) == 0 || is.na(x) || trimws(x) == "") {
    return(character(0))
  }
  tokens <- unlist(strsplit(as.character(x), "\\|", fixed = FALSE), use.names = FALSE)
  tokens <- trimws(tokens)
  tokens <- tokens[tokens != ""]
  tokens <- sub("\\s*\\(.*$", "", tokens)
  tokens <- trimws(tokens)
  tokens <- sub("^([A-Za-z0-9]+(?:-[A-Za-z0-9]+)*).*$", "\\1", tokens, perl = TRUE)
  tokens <- toupper(trimws(tokens))
  tokens <- tokens[!tolower(tokens) %in% c("na", "n/a", "unknown", "none", "<na>")]
  unique(tokens[tokens != ""])
}

is_mutation_positive <- function(x) {
  x <- tolower(trimws(as.character(x)))
  !is.na(x) & !x %in% c("", "na", "n/a", "<na>", "unknown", "negative", "false", "0")
}

stop_if_missing(input_clinical)
stop_if_missing(input_expression)
stop_if_missing(input_surface_gene_source)

log_msg("Reading clinical data...")
clinical <- openxlsx::read.xlsx(input_clinical)
aml_dx <- "ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS"

sample_info <- clinical %>%
  dplyr::filter(
    .data[["dxAtInclusion"]] == aml_dx,
    !is.na(.data[["dbgap_rnaseq_sample"]]),
    .data[["dbgap_rnaseq_sample"]] != ""
  ) %>%
  dplyr::distinct(.data[["dbgap_rnaseq_sample"]], .keep_all = TRUE) %>%
  dplyr::transmute(
    sample_id = .data[["dbgap_rnaseq_sample"]],
    subject_id = .data[["dbgap_subject_id"]],
    dnaseq_sample = .data[["dbgap_dnaseq_sample"]],
    cohort = .data[["cohort"]],
    sex = .data[["consensus_sex"]],
    race = .data[["reportedRace"]],
    ethnicity = .data[["reportedEthnicity"]],
    age_at_diagnosis = .data[["ageAtDiagnosis"]],
    age_at_specimen = .data[["ageAtSpecimenAcquisition"]],
    dx_at_inclusion = .data[["dxAtInclusion"]],
    WHO_specificDx = clean_label(.data[["specificDxAtInclusion"]]),
    ELN2017 = clean_label(.data[["ELN2017"]]),
    Fusion_subtype = clean_label(.data[["consensusAMLFusions"]]),
    FAB_subtype = clean_label(.data[["fabBlastMorphology"]]),
    Disease_origin = disease_origin(.data[["isDenovo"]], .data[["isTransformed"]]),
    is_denovo = .data[["isDenovo"]],
    is_relapse = .data[["isRelapse"]],
    is_transformed = .data[["isTransformed"]],
    Disease_stage = clean_label(.data[["diseaseStageAtSpecimenCollection"]]),
    Specimen_type = clean_label(.data[["specimenType"]], drop_unknown = FALSE),
    response_to_induction = .data[["responseToInductionTx"]],
    vital_status = .data[["vitalStatus"]],
    overall_survival = .data[["overallSurvival"]],
    blasts_bm_percent = .data[["%.Blasts.in.BM"]],
    blasts_pb_percent = .data[["%.Blasts.in.PB"]],
    wbc_count = .data[["wbcCount"]],
    platelet_count = .data[["plateletCount"]],
    hemoglobin = .data[["hemoglobin"]],
    variantSummary_raw = .data[["variantSummary"]],
    FLT3_ITD_structured = .data[["FLT3-ITD"]],
    NPM1_structured = .data[["NPM1"]],
    RUNX1_structured = .data[["RUNX1"]],
    ASXL1_structured = .data[["ASXL1"]],
    TP53_structured = .data[["TP53"]],
    CEBPA_biallelic_structured = .data[["CEBPA_Biallelic"]]
  )

variant_list <- lapply(sample_info$variantSummary_raw, parse_variant_tokens_one)
variant_assignments <- data.frame(
  sample_id = rep(sample_info$sample_id, lengths(variant_list)),
  mutation = unlist(variant_list, use.names = FALSE),
  source = "variantSummary",
  stringsAsFactors = FALSE
) %>%
  dplyr::distinct()

make_structured_assignment <- function(sample_ids, values, mutation_name, source_name, positive_fun = is_mutation_positive) {
  positive <- positive_fun(values)
  data.frame(
    sample_id = sample_ids[positive],
    mutation = mutation_name,
    source = source_name,
    stringsAsFactors = FALSE
  )
}

structured_assignments <- dplyr::bind_rows(
  make_structured_assignment(sample_info$sample_id, sample_info$FLT3_ITD_structured, "FLT3-ITD", "structured_FLT3_ITD"),
  make_structured_assignment(sample_info$sample_id, sample_info$NPM1_structured, "NPM1", "structured_NPM1"),
  make_structured_assignment(sample_info$sample_id, sample_info$RUNX1_structured, "RUNX1", "structured_RUNX1"),
  make_structured_assignment(sample_info$sample_id, sample_info$ASXL1_structured, "ASXL1", "structured_ASXL1"),
  make_structured_assignment(sample_info$sample_id, sample_info$TP53_structured, "TP53", "structured_TP53"),
  make_structured_assignment(
    sample_info$sample_id,
    sample_info$CEBPA_biallelic_structured,
    "CEBPA",
    "structured_CEBPA_Biallelic",
    positive_fun = function(x) tolower(trimws(as.character(x))) == "true"
  )
)

mutation_assignments_with_source <- dplyr::bind_rows(variant_assignments, structured_assignments) %>%
  dplyr::mutate(mutation = toupper(.data[["mutation"]])) %>%
  dplyr::filter(!is.na(.data[["mutation"]]), .data[["mutation"]] != "") %>%
  dplyr::distinct()

mutation_assignments <- mutation_assignments_with_source %>%
  dplyr::distinct(.data[["sample_id"]], .data[["mutation"]])

variant_parsed <- vapply(variant_list, function(x) paste(x, collapse = "|"), character(1))
variant_parsed[variant_parsed == ""] <- NA_character_
sample_info$variantSummary_parsed <- variant_parsed

structured_by_sample <- structured_assignments %>%
  dplyr::distinct(.data[["sample_id"]], .data[["mutation"]]) %>%
  dplyr::arrange(.data[["sample_id"]], .data[["mutation"]]) %>%
  dplyr::group_by(.data[["sample_id"]]) %>%
  dplyr::summarise(structured_mutations = paste(.data[["mutation"]], collapse = "|"), .groups = "drop")

integrated_by_sample <- mutation_assignments %>%
  dplyr::arrange(.data[["sample_id"]], .data[["mutation"]]) %>%
  dplyr::group_by(.data[["sample_id"]]) %>%
  dplyr::summarise(mutation_subtypes_integrated = paste(.data[["mutation"]], collapse = "|"), .groups = "drop")

sample_info <- sample_info %>%
  dplyr::left_join(structured_by_sample, by = "sample_id") %>%
  dplyr::left_join(integrated_by_sample, by = "sample_id")

log_msg("Reading expression matrix...")
expr_dt <- data.table::fread(input_expression, nThread = 2)
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

sample_cols <- intersect(sample_info$sample_id, names(expr_dt))
if (length(sample_cols) < min_subtype_n) {
  stop("Too few AML RNA-seq samples overlap with the expression matrix.", call. = FALSE)
}

expr_mat <- as.matrix(expr_dt[, ..sample_cols])
storage.mode(expr_mat) <- "numeric"
rownames(expr_mat) <- gene_names
rm(expr_dt)

sample_info <- sample_info %>%
  dplyr::filter(.data[["sample_id"]] %in% colnames(expr_mat))
expr_mat <- expr_mat[, sample_info$sample_id, drop = FALSE]
mutation_assignments <- mutation_assignments %>%
  dplyr::filter(.data[["sample_id"]] %in% sample_info$sample_id)
mutation_assignments_with_source <- mutation_assignments_with_source %>%
  dplyr::filter(.data[["sample_id"]] %in% sample_info$sample_id)

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

filtered_clinical <- sample_info %>%
  dplyr::select(
    sample_id,
    subject_id,
    dnaseq_sample,
    cohort,
    sex,
    race,
    ethnicity,
    age_at_diagnosis,
    age_at_specimen,
    dx_at_inclusion,
    WHO_specificDx,
    ELN2017,
    Fusion_subtype,
    FAB_subtype,
    Disease_origin,
    is_denovo,
    is_relapse,
    is_transformed,
    Disease_stage,
    Specimen_type,
    response_to_induction,
    vital_status,
    overall_survival,
    blasts_bm_percent,
    blasts_pb_percent,
    wbc_count,
    platelet_count,
    hemoglobin,
    lsc17_score,
    variantSummary_raw,
    variantSummary_parsed,
    structured_mutations,
    mutation_subtypes_integrated
  )
data.table::fwrite(filtered_clinical, output_filtered_clinical)

log_msg("Reading surface gene list from previous project-local LASA result...")
surface_source <- data.table::fread(input_surface_gene_source)
if (!"gene" %in% names(surface_source)) {
  stop("Surface gene source must contain a gene column: ", input_surface_gene_source, call. = FALSE)
}

surface_gene_list <- unique(na.omit(as.character(surface_source$gene)))
surface_gene_list <- surface_gene_list[surface_gene_list %in% rownames(expr_mat)]
if (length(surface_gene_list) == 0) {
  stop("No surface genes overlap with the expression matrix.", call. = FALSE)
}

data.table::fwrite(data.frame(gene = surface_gene_list), file.path(output_dir, "surface_gene_list_used.csv"))
data.table::fwrite(sample_info, file.path(output_dir, "lsc17_sample_scores_with_subtypes.csv"))
data.table::fwrite(variant_assignments, file.path(output_dir, "variantSummary_mutation_assignments_long.csv"))
data.table::fwrite(structured_assignments, file.path(output_dir, "structured_mutation_assignments_long.csv"))
data.table::fwrite(mutation_assignments_with_source, file.path(output_dir, "integrated_mutation_assignments_with_source_long.csv"))
data.table::fwrite(mutation_assignments, file.path(output_dir, "integrated_mutation_assignments_long.csv"))

clinical_subtype_variables <- c(
  "WHO_specificDx",
  "ELN2017",
  "Fusion_subtype",
  "FAB_subtype",
  "Disease_origin",
  "Disease_stage",
  "Specimen_type"
)

clinical_assignments <- dplyr::bind_rows(lapply(clinical_subtype_variables, function(v) {
  data.frame(
    sample_id = sample_info$sample_id,
    subtype_type = "clinical",
    subtype_variable = v,
    subtype_level = sample_info[[v]],
    stringsAsFactors = FALSE
  )
})) %>%
  dplyr::filter(!is.na(.data[["subtype_level"]]))

mutation_subtype_assignments <- mutation_assignments %>%
  dplyr::transmute(
    sample_id = .data[["sample_id"]],
    subtype_type = "mutation",
    subtype_variable = "integrated_mutation",
    subtype_level = .data[["mutation"]]
  )

subtype_assignments <- dplyr::bind_rows(clinical_assignments, mutation_subtype_assignments)

subtype_counts <- subtype_assignments %>%
  dplyr::count(subtype_type, subtype_variable, subtype_level, name = "n_samples") %>%
  dplyr::arrange(subtype_type, subtype_variable, dplyr::desc(n_samples), subtype_level)

eligible_subtypes <- subtype_counts %>%
  dplyr::filter(n_samples >= min_subtype_n)

skipped_subtypes <- subtype_counts %>%
  dplyr::filter(n_samples < min_subtype_n)

data.table::fwrite(subtype_counts, file.path(output_dir, "subtype_sample_counts.csv"))
data.table::fwrite(eligible_subtypes, file.path(output_dir, "subtype_eligible_groups.csv"))
data.table::fwrite(skipped_subtypes, file.path(output_dir, "subtype_skipped_groups.csv"))

surface_expr <- expr_mat[surface_gene_list, , drop = FALSE]
rm(expr_mat)

correlate_surface_genes <- function(sample_ids, subtype_type, subtype_variable, subtype_level) {
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
    subtype_type = subtype_type,
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

log_msg("Running subtype-stratified Spearman correlations...")
cor_list <- vector("list", nrow(eligible_subtypes))
for (i in seq_len(nrow(eligible_subtypes))) {
  current_type <- eligible_subtypes$subtype_type[i]
  current_var <- eligible_subtypes$subtype_variable[i]
  current_level <- eligible_subtypes$subtype_level[i]
  current_samples <- subtype_assignments %>%
    dplyr::filter(
      .data[["subtype_type"]] == current_type,
      .data[["subtype_variable"]] == current_var,
      .data[["subtype_level"]] == current_level
    ) %>%
    dplyr::pull(.data[["sample_id"]])

  cor_list[[i]] <- correlate_surface_genes(
    sample_ids = current_samples,
    subtype_type = current_type,
    subtype_variable = current_var,
    subtype_level = current_level
  )
}

cor_df <- dplyr::bind_rows(cor_list) %>%
  dplyr::mutate(
    FDR_global = p.adjust(.data[["pvalue"]], method = "BH"),
    direction = dplyr::case_when(
      .data[["rho"]] > 0 ~ "positive",
      .data[["rho"]] < 0 ~ "negative",
      TRUE ~ "zero"
    )
  ) %>%
  dplyr::arrange(.data[["subtype_type"]], .data[["subtype_variable"]], .data[["subtype_level"]], dplyr::desc(abs(.data[["rho"]])))

data.table::fwrite(cor_df, file.path(output_dir, "subtype_surface_correlation_all.csv"))

sig_df <- cor_df %>%
  dplyr::filter(.data[["FDR_subtype"]] < significance_fdr) %>%
  dplyr::arrange(.data[["subtype_type"]], .data[["subtype_variable"]], .data[["subtype_level"]], .data[["FDR_subtype"]], dplyr::desc(abs(.data[["rho"]])))
data.table::fwrite(sig_df, file.path(output_dir, "subtype_surface_correlation_significant_FDR05.csv"))

top_by_abs <- cor_df %>%
  dplyr::filter(!is.na(.data[["rho"]])) %>%
  dplyr::group_by(.data[["subtype_type"]], .data[["subtype_variable"]], .data[["subtype_level"]]) %>%
  dplyr::slice_max(order_by = abs(.data[["rho"]]), n = top_n_per_group, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(.data[["subtype_type"]], .data[["subtype_variable"]], .data[["subtype_level"]], dplyr::desc(abs(.data[["rho"]])))
data.table::fwrite(top_by_abs, file.path(output_dir, "subtype_top_surface_genes_by_abs_rho.csv"))

summary_by_group <- cor_df %>%
  dplyr::group_by(.data[["subtype_type"]], .data[["subtype_variable"]], .data[["subtype_level"]], .data[["n_samples"]]) %>%
  dplyr::summarise(
    n_surface_genes_tested = dplyr::n(),
    n_tested_non_na = sum(!is.na(.data[["rho"]])),
    n_FDR05 = sum(.data[["FDR_subtype"]] < significance_fdr, na.rm = TRUE),
    n_FDR05_positive = sum(.data[["FDR_subtype"]] < significance_fdr & .data[["rho"]] > 0, na.rm = TRUE),
    n_FDR05_negative = sum(.data[["FDR_subtype"]] < significance_fdr & .data[["rho"]] < 0, na.rm = TRUE),
    max_abs_rho = max(abs(.data[["rho"]]), na.rm = TRUE),
    top_abs_gene = .data[["gene"]][which.max(abs(.data[["rho"]]))],
    top_abs_rho = .data[["rho"]][which.max(abs(.data[["rho"]]))],
    top_abs_FDR = .data[["FDR_subtype"]][which.max(abs(.data[["rho"]]))],
    .groups = "drop"
  ) %>%
  dplyr::arrange(.data[["subtype_type"]], dplyr::desc(.data[["n_FDR05"]]), dplyr::desc(.data[["max_abs_rho"]]))
data.table::fwrite(summary_by_group, file.path(output_dir, "subtype_correlation_summary_by_group.csv"))

sig_counts <- cor_df %>%
  dplyr::mutate(significant = .data[["FDR_subtype"]] < significance_fdr) %>%
  dplyr::filter(.data[["significant"]]) %>%
  dplyr::count(.data[["subtype_type"]], .data[["subtype_variable"]], .data[["subtype_level"]], .data[["direction"]], name = "n_significant_genes")

sig_counts_complete <- eligible_subtypes %>%
  tidyr::crossing(direction = c("positive", "negative")) %>%
  dplyr::left_join(sig_counts, by = c("subtype_type", "subtype_variable", "subtype_level", "direction")) %>%
  dplyr::mutate(n_significant_genes = tidyr::replace_na(.data[["n_significant_genes"]], 0L))
data.table::fwrite(sig_counts_complete, file.path(output_dir, "subtype_significant_gene_counts.csv"))

plot_group_label <- function(variable, level) {
  paste(variable, level, sep = ": ")
}

log_msg("Creating PDF figures...")

summary_plot_df <- summary_by_group %>%
  dplyr::mutate(group_label = dplyr::if_else(
    .data[["subtype_variable"]] == "integrated_mutation",
    .data[["subtype_level"]],
    plot_group_label(.data[["subtype_variable"]], .data[["subtype_level"]])
  )) %>%
  dplyr::arrange(.data[["n_FDR05"]]) %>%
  dplyr::mutate(group_label = factor(.data[["group_label"]], levels = unique(.data[["group_label"]])))

p_counts <- ggplot(summary_plot_df, aes(x = group_label, y = n_FDR05, fill = subtype_type)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c(clinical = "#4E79A7", mutation = "#B07AA1")) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Surface genes with subtype FDR < 0.05",
    title = "Significant LSC17-associated surface genes by subtype"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 5.5),
    legend.position = "top"
  )
ggsave(
  filename = file.path(output_dir, "subtype_significant_gene_counts.pdf"),
  plot = p_counts,
  width = 8,
  height = max(6, 0.12 * nrow(summary_plot_df)),
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

heatmap_genes <- sig_df %>%
  dplyr::group_by(.data[["gene"]]) %>%
  dplyr::summarise(max_abs_rho = max(abs(.data[["rho"]]), na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(.data[["max_abs_rho"]]), .data[["gene"]]) %>%
  dplyr::slice_head(n = top_n_global_heatmap_genes) %>%
  dplyr::pull(.data[["gene"]])

heatmap_df <- cor_df %>%
  dplyr::filter(.data[["gene"]] %in% heatmap_genes) %>%
  dplyr::mutate(
    group_label = dplyr::if_else(
      .data[["subtype_variable"]] == "integrated_mutation",
      .data[["subtype_level"]],
      plot_group_label(.data[["subtype_variable"]], .data[["subtype_level"]])
    ),
    gene = factor(.data[["gene"]], levels = rev(heatmap_genes))
  )
group_order <- summary_by_group %>%
  dplyr::mutate(group_label = dplyr::if_else(
    .data[["subtype_variable"]] == "integrated_mutation",
    .data[["subtype_level"]],
    plot_group_label(.data[["subtype_variable"]], .data[["subtype_level"]])
  )) %>%
  dplyr::arrange(.data[["subtype_type"]], .data[["subtype_variable"]], dplyr::desc(.data[["n_samples"]]), .data[["subtype_level"]]) %>%
  dplyr::pull(.data[["group_label"]])
heatmap_df$group_label <- factor(heatmap_df$group_label, levels = group_order)

p_heatmap <- ggplot(heatmap_df, aes(x = group_label, y = gene, fill = rho)) +
  geom_tile(color = "grey88", linewidth = 0.12) +
  facet_grid(. ~ subtype_type, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "#2F6FA3",
    mid = "white",
    high = "#C94C4C",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman rho"
  ) +
  labs(x = NULL, y = NULL, title = "Top global subtype-specific surface gene correlations with LSC17 score") +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 4.5),
    axis.text.y = element_text(size = 5.5),
    strip.text = element_text(face = "bold", size = 7),
    panel.grid = element_blank()
  )
ggsave(
  filename = file.path(output_dir, "subtype_top_correlation_heatmap.pdf"),
  plot = p_heatmap,
  width = 16,
  height = 8,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

score_long <- subtype_assignments %>%
  dplyr::inner_join(eligible_subtypes, by = c("subtype_type", "subtype_variable", "subtype_level")) %>%
  dplyr::left_join(sample_info %>% dplyr::select(sample_id, lsc17_score), by = "sample_id") %>%
  dplyr::mutate(group_label = dplyr::if_else(
    .data[["subtype_variable"]] == "integrated_mutation",
    .data[["subtype_level"]],
    .data[["subtype_level"]]
  ))

p_score <- ggplot(score_long, aes(x = group_label, y = lsc17_score)) +
  geom_boxplot(outlier.shape = NA, fill = "#A0CBE8", color = "grey30", width = 0.55) +
  geom_jitter(width = 0.12, size = 0.35, alpha = 0.35, color = "grey25") +
  facet_wrap(~ subtype_variable, scales = "free_x", ncol = 2) +
  labs(x = NULL, y = "LSC17 score", title = "LSC17 score distribution across eligible subtype groups") +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    strip.text = element_text(face = "bold", size = 7)
  )
ggsave(
  filename = file.path(output_dir, "subtype_lsc17_score_distribution.pdf"),
  plot = p_score,
  width = 11,
  height = 9,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

scatter_pairs <- cor_df %>%
  dplyr::filter(.data[["FDR_subtype"]] < significance_fdr) %>%
  dplyr::arrange(dplyr::desc(abs(.data[["rho"]]))) %>%
  dplyr::distinct(.data[["subtype_type"]], .data[["subtype_variable"]], .data[["subtype_level"]], .data[["gene"]], .keep_all = TRUE) %>%
  dplyr::slice_head(n = top_n_scatter_pairs)

if (nrow(scatter_pairs) > 0) {
  scatter_df <- dplyr::bind_rows(lapply(seq_len(nrow(scatter_pairs)), function(i) {
    row <- scatter_pairs[i, ]
    samples <- subtype_assignments %>%
      dplyr::filter(
        .data[["subtype_type"]] == row$subtype_type,
        .data[["subtype_variable"]] == row$subtype_variable,
        .data[["subtype_level"]] == row$subtype_level
      ) %>%
      dplyr::pull(.data[["sample_id"]])
    samples <- intersect(samples, colnames(surface_expr))
    data.frame(
      sample_id = samples,
      gene = row$gene,
      subtype_type = row$subtype_type,
      subtype_variable = row$subtype_variable,
      subtype_level = row$subtype_level,
      rho = row$rho,
      FDR_subtype = row$FDR_subtype,
      expression = as.numeric(surface_expr[row$gene, samples]),
      lsc17_score = sample_info$lsc17_score[match(samples, sample_info$sample_id)],
      stringsAsFactors = FALSE
    )
  })) %>%
    dplyr::mutate(
      facet_label = paste0(
        subtype_type, " | ", subtype_level, "\n",
        gene, " rho=", sprintf("%.2f", rho), ", FDR=", format(FDR_subtype, digits = 2, scientific = TRUE)
      )
    )

  p_scatter <- ggplot(scatter_df, aes(x = expression, y = lsc17_score)) +
    geom_point(size = 1.2, alpha = 0.75, color = "#4E79A7") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.35, color = "#C94C4C") +
    facet_wrap(~ facet_label, scales = "free_x", ncol = 2) +
    labs(x = "Surface gene expression", y = "LSC17 score", title = "Top subtype-specific gene-LSC17 associations") +
    theme_classic(base_size = 8) +
    theme(plot.title = element_text(face = "bold", size = 10), strip.text = element_text(size = 6))

  ggsave(
    filename = file.path(output_dir, "subtype_top_gene_scatterplots.pdf"),
    plot = p_scatter,
    width = 8,
    height = max(5, 2.1 * ceiling(nrow(scatter_pairs) / 2)),
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
  )
}

make_bubble_heatmap <- function(cor_input, subtype_type_filter, top_n, output_file, title_text) {
  top_pairs <- cor_input %>%
    dplyr::filter(.data[["subtype_type"]] == subtype_type_filter, !is.na(.data[["rho"]])) %>%
    dplyr::group_by(.data[["subtype_variable"]], .data[["subtype_level"]]) %>%
    dplyr::slice_max(order_by = abs(.data[["rho"]]), n = top_n, with_ties = FALSE) %>%
    dplyr::ungroup()

  data.table::fwrite(top_pairs, file.path(output_dir, paste0(subtype_type_filter, "_bubble_heatmap_top_pairs.csv")))

  selected_genes <- top_pairs %>%
    dplyr::group_by(.data[["gene"]]) %>%
    dplyr::summarise(max_abs_rho = max(abs(.data[["rho"]]), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data[["max_abs_rho"]]), .data[["gene"]]) %>%
    dplyr::pull(.data[["gene"]])

  plot_df <- cor_input %>%
    dplyr::filter(.data[["subtype_type"]] == subtype_type_filter, .data[["gene"]] %in% selected_genes) %>%
    dplyr::mutate(
      subtype_label = dplyr::if_else(
        .data[["subtype_variable"]] == "integrated_mutation",
        .data[["subtype_level"]],
        plot_group_label(.data[["subtype_variable"]], .data[["subtype_level"]])
      ),
      neg_log10_fdr = -log10(pmax(.data[["FDR_subtype"]], .Machine$double.xmin)),
      neg_log10_fdr_capped = pmin(.data[["neg_log10_fdr"]], 20),
      significant = .data[["FDR_subtype"]] < significance_fdr
    )

  column_order <- plot_df %>%
    dplyr::distinct(.data[["subtype_variable"]], .data[["subtype_level"]], .data[["subtype_label"]], .data[["n_samples"]]) %>%
    dplyr::arrange(.data[["subtype_variable"]], dplyr::desc(.data[["n_samples"]]), .data[["subtype_level"]]) %>%
    dplyr::pull(.data[["subtype_label"]])

  plot_df$gene <- factor(plot_df$gene, levels = rev(selected_genes))
  plot_df$subtype_label <- factor(plot_df$subtype_label, levels = column_order)

  n_cols <- length(column_order)
  n_rows <- length(selected_genes)
  plot_width <- max(8, min(18, 0.32 * n_cols + 3.5))
  plot_height <- max(6, min(22, 0.13 * n_rows + 2.8))

  p <- ggplot(plot_df, aes(x = subtype_label, y = gene)) +
    geom_point(
      aes(size = neg_log10_fdr_capped, fill = rho, alpha = significant),
      shape = 21,
      color = "grey20",
      stroke = 0.16
    ) +
    scale_fill_gradient2(
      low = "#2F6FA3",
      mid = "white",
      high = "#C94C4C",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Spearman rho"
    ) +
    scale_size_continuous(
      range = c(0.5, 3.5),
      name = "-log10(FDR)",
      breaks = c(2, 5, 10, 20),
      labels = c("2", "5", "10", ">=20")
    ) +
    scale_alpha_manual(values = c("FALSE" = 0.25, "TRUE" = 0.95), guide = "none") +
    labs(x = NULL, y = NULL, title = title_text) +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
      axis.text.y = element_text(size = 5.5),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.18),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6)
    )

  if (subtype_type_filter == "clinical") {
    p <- p + facet_grid(. ~ subtype_variable, scales = "free_x", space = "free_x")
  }

  ggsave(
    filename = output_file,
    plot = p,
    width = plot_width,
    height = plot_height,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
  )
}

make_bubble_heatmap(
  cor_input = cor_df,
  subtype_type_filter = "mutation",
  top_n = top_n_mutation_heatmap,
  output_file = file.path(output_dir, "mutation_subtype_top10_bubble_heatmap.pdf"),
  title_text = "Integrated mutation subtype-specific surface gene correlations with LSC17 score"
)

make_bubble_heatmap(
  cor_input = cor_df,
  subtype_type_filter = "clinical",
  top_n = top_n_clinical_heatmap,
  output_file = file.path(output_dir, "clinical_subtype_top20_bubble_heatmap.pdf"),
  title_text = "Clinical subtype-specific surface gene correlations with LSC17 score"
)

mutation_source_summary <- mutation_assignments_with_source %>%
  dplyr::group_by(.data[["mutation"]]) %>%
  dplyr::summarise(
    n_integrated = dplyr::n_distinct(.data[["sample_id"]]),
    n_variantSummary = dplyr::n_distinct(.data[["sample_id"]][.data[["source"]] == "variantSummary"]),
    n_structured = dplyr::n_distinct(.data[["sample_id"]][.data[["source"]] != "variantSummary"]),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(.data[["n_integrated"]]), .data[["mutation"]])
data.table::fwrite(mutation_source_summary, file.path(output_dir, "mutation_source_summary.csv"))

run_info <- data.frame(
  item = c(
    "project_root",
    "clinical_file",
    "expression_file",
    "surface_gene_source",
    "filtered_clinical_output",
    "output_dir",
    "min_subtype_n",
    "top_n_per_group",
    "top_n_mutation_heatmap",
    "top_n_clinical_heatmap",
    "n_aml_rnaseq_samples_in_clinical",
    "n_aml_expression_overlap_samples",
    "n_surface_genes_tested",
    "n_variantSummary_assignments",
    "n_structured_assignments",
    "n_integrated_mutation_assignments",
    "n_unique_integrated_mutation_tokens",
    "n_eligible_subtype_groups",
    "n_eligible_mutation_groups",
    "n_eligible_clinical_groups",
    "n_skipped_subtype_groups"
  ),
  value = c(
    project_root,
    input_clinical,
    input_expression,
    input_surface_gene_source,
    output_filtered_clinical,
    output_dir,
    min_subtype_n,
    top_n_per_group,
    top_n_mutation_heatmap,
    top_n_clinical_heatmap,
    sum(clinical[["dxAtInclusion"]] == aml_dx & !is.na(clinical[["dbgap_rnaseq_sample"]])),
    nrow(sample_info),
    length(surface_gene_list),
    nrow(variant_assignments),
    nrow(structured_assignments),
    nrow(mutation_assignments),
    length(unique(mutation_assignments$mutation)),
    nrow(eligible_subtypes),
    sum(eligible_subtypes$subtype_type == "mutation"),
    sum(eligible_subtypes$subtype_type == "clinical"),
    nrow(skipped_subtypes)
  )
)
data.table::fwrite(run_info, file.path(output_dir, "pipeline_run_info.csv"))

log_msg("Done.")
log_msg("Filtered clinical file: ", output_filtered_clinical)
log_msg("Outputs written to: ", output_dir)
log_msg("Session info:")
print(sessionInfo())
