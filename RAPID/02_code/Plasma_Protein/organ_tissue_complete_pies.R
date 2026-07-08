library(dplyr)
library(ggplot2)
library(readxl)

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1])
if (is.na(script_file) || !nzchar(script_file)) {
  script_file <- file.path(getwd(), "02_code", "Plasma_Protein", "organ_tissue_complete_pies.R")
}
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
data_dir <- file.path(repo_dir, "01_data", "Plasma_Protein")
out_dir <- file.path(repo_dir, "03_result", "Plasma_Protein", "Organ_tissue_complete_pies")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Organ tissue complete pies pipeline\n")
cat("Run time:", format(Sys.time()), "\n")
cat("Repository directory:", repo_dir, "\n")
cat("Data directory:", data_dir, "\n")
cat("Output directory:", out_dir, "\n")

organs <- c(
  "Brain",
  "Cerebellum",
  "Neuronal",
  "Oligodendrocytes",
  "Retina",
  "Choroid plexus",
  "Adrenal gland",
  "Pituitary gland",
  "Parathyroid gland",
  "Bone marrow",
  "Spleen",
  "Lymphoid tissue",
  "Thymus",
  "Lung",
  "Liver",
  "Stomach",
  "Small intestine",
  "Intestine",
  "Pancreas",
  "Salivary gland",
  "Kidney",
  "Breast",
  "Fallopian tube",
  "Epididymis",
  "Spermatids",
  "Placenta",
  "Skeletal muscle",
  "Heart muscle",
  "Smooth muscle",
  "Epithelium",
  "Squamous epithelium",
  "Skin",
  "Ciliated cells",
  "Connective tissue"
)

tissue_colors <- c(
  "Brain" = "#D8C9EA",
  "Cerebellum" = "#8B315D",
  "Neuronal" = "#6A4C93",
  "Oligodendrocytes" = "#B8897C",
  "Retina" = "#294E7A",
  "Choroid plexus" = "#D36B6B",
  "Adrenal gland" = "#8E6F9E",
  "Pituitary gland" = "#B7A7C8",
  "Parathyroid gland" = "#D6C7A8",
  "Bone marrow" = "#A9C88E",
  "Spleen" = "#3B6B61",
  "Lymphoid tissue" = "#6FA59C",
  "Thymus" = "#1F3F3A",
  "Lung" = "#9AB8C9",
  "Liver" = "#B93A32",
  "Stomach" = "#E2C8B8",
  "Small intestine" = "#C0A0C8",
  "Intestine" = "#5AB7AD",
  "Pancreas" = "#A6B5D8",
  "Salivary gland" = "#E0A15E",
  "Kidney" = "#346B8F",
  "Breast" = "#DFA4A5",
  "Fallopian tube" = "#8A4E78",
  "Epididymis" = "#6E5A8A",
  "Spermatids" = "#C7A7DD",
  "Placenta" = "#D68AA8",
  "Skeletal muscle" = "#B47C72",
  "Heart muscle" = "#5E7892",
  "Smooth muscle" = "#C9B27E",
  "Epithelium" = "#7BAF9E",
  "Squamous epithelium" = "#D9CDEB",
  "Skin" = "#C85C72",
  "Ciliated cells" = "#82A6B1",
  "Connective tissue" = "#A97C50",
  "Others" = "#D7D7D7"
)

normalize_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

normalize_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "", x)
}

format_percent <- function(x, digits = 1) {
  paste0(formatC(x, format = "f", digits = digits), "%")
}

required_files <- c(
  low_abundance = file.path(data_dir, "Low_abundance_15ul_tidy.csv"),
  total_protein = file.path(data_dir, "Plasma_Total_protein_tidy.csv"),
  specific_tissue = file.path(data_dir, "Specific tissue cluster2.xlsx")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files: ", paste(missing_files, collapse = "; "))
}

cat("\nInput files:\n")
print(required_files)

low_abundance <- read.csv(required_files[["low_abundance"]], row.names = 1, check.names = FALSE)
low_abundance <- low_abundance[!is.na(low_abundance$Low_abundance_15ul), , drop = FALSE]
total_protein <- read.csv(required_files[["total_protein"]], row.names = 1, check.names = FALSE)

if (!"Gene" %in% names(low_abundance)) {
  stop("Low_abundance_15ul_tidy.csv is missing column: Gene")
}
if (!"Gene" %in% names(total_protein)) {
  stop("Plasma_Total_protein_tidy.csv is missing column: Gene")
}

neat_genes <- total_protein %>%
  transmute(Gene = as.character(Gene), Gene_key = normalize_gene(Gene)) %>%
  filter(!is.na(Gene_key), nzchar(Gene_key)) %>%
  distinct(Gene_key, .keep_all = TRUE)

enriched_genes <- low_abundance %>%
  transmute(Gene = as.character(Gene), Gene_key = normalize_gene(Gene)) %>%
  filter(!is.na(Gene_key), nzchar(Gene_key)) %>%
  distinct(Gene_key, .keep_all = TRUE)

neat_keys <- neat_genes$Gene_key
enriched_keys <- enriched_genes$Gene_key
enriched_only_keys <- setdiff(enriched_keys, neat_keys)

cat("\nInput gene-set summary:\n")
cat("Neat unique proteins:", length(neat_keys), "\n")
cat("Enriched unique proteins:", length(enriched_keys), "\n")
cat("Enriched-only proteins:", length(enriched_only_keys), "\n")

excel_sheets <- excel_sheets(required_files[["specific_tissue"]])
tissue_all <- bind_rows(lapply(excel_sheets, function(sheet_name) {
  sheet_data <- read_excel(required_files[["specific_tissue"]], sheet = sheet_name)
  sheet_data$Source_sheet <- sheet_name
  sheet_data
}))

required_tissue_cols <- c("Gene", "Tissue", "Function")
missing_tissue_cols <- setdiff(required_tissue_cols, names(tissue_all))
if (length(missing_tissue_cols) > 0) {
  stop("Specific tissue workbook is missing required columns: ", paste(missing_tissue_cols, collapse = ", "))
}

requested_organs <- data.frame(
  Tissue = organs,
  Tissue_norm = normalize_label(organs),
  Tissue_order = seq_along(organs),
  stringsAsFactors = FALSE
)

tissue_lookup <- tissue_all %>%
  transmute(Tissue_original = as.character(Tissue), Tissue_norm = normalize_label(Tissue)) %>%
  distinct()

tissue_name_matching <- requested_organs %>%
  left_join(tissue_lookup, by = "Tissue_norm") %>%
  arrange(Tissue_order, Tissue_original)

missing_tissues <- tissue_name_matching %>%
  filter(is.na(Tissue_original)) %>%
  distinct(Tissue) %>%
  pull(Tissue)

tissue_selected <- tissue_all %>%
  mutate(
    Tissue_original = as.character(Tissue),
    Tissue_norm = normalize_label(Tissue_original)
  ) %>%
  inner_join(requested_organs, by = "Tissue_norm", suffix = c("_workbook", "")) %>%
  mutate(
    Gene = as.character(Gene),
    Gene_key = normalize_gene(Gene),
    Function = as.character(Function),
    Source_sheet = as.character(Source_sheet)
  ) %>%
  filter(!is.na(Gene_key), nzchar(Gene_key)) %>%
  group_by(Gene_key, Tissue, Tissue_order) %>%
  summarise(
    Example_gene_label = dplyr::first(stats::na.omit(Gene)),
    Tissue_original = paste(sort(unique(Tissue_original)), collapse = "; "),
    Function = paste(sort(unique(Function)), collapse = "; "),
    Source_sheets = paste(sort(unique(Source_sheet)), collapse = "; "),
    .groups = "drop"
  )

protein_sets <- bind_rows(
  neat_genes %>% mutate(Protein_set = "Neat"),
  enriched_genes %>% mutate(Protein_set = "Enriched"),
  enriched_genes %>% filter(Gene_key %in% enriched_only_keys) %>% mutate(Protein_set = "Enriched only new")
) %>%
  mutate(
    Protein_set = factor(Protein_set, levels = c("Neat", "Enriched", "Enriched only new"))
  )

matched_assignments <- protein_sets %>%
  inner_join(tissue_selected, by = "Gene_key") %>%
  group_by(Protein_set, Gene_key) %>%
  mutate(
    matched_tissue_count = n_distinct(Tissue),
    fractional_weight = 1 / matched_tissue_count,
    Assignment_type = "Selected tissue"
  ) %>%
  ungroup()

others_assignments <- protein_sets %>%
  anti_join(tissue_selected %>% distinct(Gene_key), by = "Gene_key") %>%
  transmute(
    Protein_set,
    Gene_key,
    Gene,
    Tissue = "Others",
    Tissue_order = length(organs) + 1L,
    Example_gene_label = Gene,
    Tissue_original = "Others",
    Function = "",
    Source_sheets = "",
    matched_tissue_count = 0L,
    fractional_weight = 1,
    Assignment_type = "Others"
  )

complete_assignments <- bind_rows(matched_assignments, others_assignments) %>%
  mutate(
    Tissue = factor(Tissue, levels = c(organs, "Others")),
    Protein_set = factor(as.character(Protein_set), levels = c("Neat", "Enriched", "Enriched only new"))
  ) %>%
  arrange(Protein_set, Tissue, Gene_key)

complete_tissue_summary <- complete_assignments %>%
  group_by(Protein_set, Tissue) %>%
  summarise(
    fractional_count = sum(fractional_weight),
    assigned_row_count = n(),
    unique_proteins = n_distinct(Gene_key),
    .groups = "drop"
  ) %>%
  group_by(Protein_set) %>%
  mutate(
    set_total_fractional_count = sum(fractional_count),
    percent = 100 * fractional_count / set_total_fractional_count
  ) %>%
  ungroup() %>%
  arrange(Protein_set, desc(percent), Tissue)

summary_metrics <- complete_tissue_summary %>%
  group_by(Protein_set) %>%
  summarise(
    total_proteins = sum(fractional_count),
    selected_tissue_fractional_count = sum(fractional_count[as.character(Tissue) != "Others"]),
    others_count = sum(fractional_count[as.character(Tissue) == "Others"]),
    selected_tissue_percent = 100 * selected_tissue_fractional_count / total_proteins,
    others_percent = 100 * others_count / total_proteins,
    .groups = "drop"
  )

color_key <- data.frame(
  Tissue = names(tissue_colors),
  Color_hex = unname(tissue_colors),
  stringsAsFactors = FALSE
)

write.csv(tissue_name_matching, file.path(out_dir, "tissue_name_matching.csv"), row.names = FALSE)
write.csv(tissue_selected, file.path(out_dir, "selected_tissue_gene_annotations.csv"), row.names = FALSE)
write.csv(complete_assignments, file.path(out_dir, "complete_gene_tissue_fractional_assignments.csv"), row.names = FALSE)
write.csv(complete_tissue_summary, file.path(out_dir, "complete_tissue_pie_summary.csv"), row.names = FALSE)
write.csv(summary_metrics, file.path(out_dir, "summary_metrics.csv"), row.names = FALSE)
write.csv(color_key, file.path(out_dir, "tissue_color_key.csv"), row.names = FALSE)

cat("\nMatched requested tissues:", length(setdiff(organs, missing_tissues)), "\n")
cat("Missing requested tissues:", ifelse(length(missing_tissues) > 0, paste(missing_tissues, collapse = " | "), "None"), "\n")
cat("\nComplete-pie summary metrics:\n")
print(summary_metrics)

make_complete_pie <- function(set_name, title_text, filename) {
  plot_data <- complete_tissue_summary %>%
    filter(Protein_set == set_name) %>%
    mutate(
      Tissue = factor(as.character(Tissue), levels = c(organs, "Others")),
      label = ifelse(percent >= 4, format_percent(percent), "")
    )

  total_n <- unique(plot_data$set_total_fractional_count)

  p <- ggplot(plot_data, aes(x = "", y = fractional_count, fill = Tissue)) +
    geom_col(width = 1, color = "white", linewidth = 0.45) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 3.0,
      color = "black"
    ) +
    scale_fill_manual(values = tissue_colors, drop = FALSE) +
    labs(
      title = title_text,
      subtitle = paste0("Complete denominator: ", formatC(total_n, format = "f", digits = 0), " proteins; unmatched proteins are Others"),
      fill = NULL
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 9.2, hjust = 0.5, color = "grey30"),
      legend.position = "right",
      legend.text = element_text(size = 7.8),
      plot.margin = margin(10, 8, 10, 8)
    ) +
    guides(fill = guide_legend(ncol = 2, override.aes = list(linewidth = 0)))

  pdf(file.path(out_dir, filename), width = 9.2, height = 6.2)
  print(p)
  dev.off()
}

make_complete_pie("Neat", "Neat proteins: tissue origin composition", "neat_tissue_complete_pie.pdf")
make_complete_pie("Enriched", "Enriched proteins: tissue origin composition", "enriched_tissue_complete_pie.pdf")
make_complete_pie("Enriched only new", "Enriched-only new proteins: tissue origin composition", "enriched_only_new_tissue_complete_pie.pdf")

get_metric_line <- function(set_name) {
  row <- summary_metrics %>% filter(Protein_set == set_name)
  paste0(
    "- ",
    set_name,
    ": ",
    formatC(row$total_proteins, format = "f", digits = 0),
    " total proteins; selected tissue fraction ",
    format_percent(row$selected_tissue_percent),
    "; Others ",
    format_percent(row$others_percent),
    "."
  )
}

top_lines <- function(set_name, n = 6) {
  complete_tissue_summary %>%
    filter(Protein_set == set_name) %>%
    arrange(desc(percent), Tissue) %>%
    slice_head(n = n) %>%
    transmute(line = paste0("- ", Tissue, ": ", format_percent(percent), " (", formatC(fractional_count, format = "f", digits = 1), ")")) %>%
    pull(line)
}

readme_lines <- c(
  "# Complete tissue-origin pie charts",
  "",
  "This folder contains exactly three complete pie charts requested for Neat, Enriched, and Enriched-only new proteins. Each pie uses all proteins in that set as the denominator.",
  "",
  "## Inputs",
  "",
  "- `Plasma_Total_protein_tidy.csv`: unique `Gene` labels define the Neat set.",
  "- `Low_abundance_15ul_tidy.csv`: unique non-missing `Gene` labels define the Enriched set.",
  "- `Specific tissue cluster2.xlsx`: all sheets are read; only the selected tissue list used in the Sankey analyses is counted as tissue source.",
  "",
  "## Counting rules",
  "",
  "- Gene labels are trimmed and matched case-insensitively.",
  "- `Enriched only new` means proteins in Enriched but not in Neat.",
  "- Each protein contributes total weight 1 to its pie.",
  "- If a protein matches multiple selected tissues, its weight is split equally across those tissues.",
  "- If a protein does not match any selected tissue, its full weight is assigned to `Others`.",
  "- Tissue colors are copied from `03_result/Plasma_Protein/Sankey_reference_palette_organ_flow_colors/01_cell_fraction_reference_organ_flow_colors/palette_key.csv`; `Others` is neutral grey.",
  "",
  "## Main metrics",
  "",
  get_metric_line("Neat"),
  get_metric_line("Enriched"),
  get_metric_line("Enriched only new"),
  "",
  "## PDF outputs",
  "",
  "- `neat_tissue_complete_pie.pdf`: complete tissue-origin pie for all Neat proteins; non-tissue/unmatched proteins are `Others`.",
  "- `enriched_tissue_complete_pie.pdf`: complete tissue-origin pie for all Enriched proteins; non-tissue/unmatched proteins are `Others`.",
  "- `enriched_only_new_tissue_complete_pie.pdf`: complete tissue-origin pie for Enriched-only new proteins; non-tissue/unmatched proteins are `Others`.",
  "",
  "## Table outputs",
  "",
  "- `complete_tissue_pie_summary.csv`: fractional counts and percentages used by the three pie charts.",
  "- `complete_gene_tissue_fractional_assignments.csv`: per-protein tissue assignment with fractional weights; unmatched proteins are assigned to `Others`.",
  "- `summary_metrics.csv`: set-level selected-tissue and Others percentages.",
  "- `selected_tissue_gene_annotations.csv`: selected tissue annotations from `Specific tissue cluster2.xlsx` after tissue-name matching.",
  "- `tissue_name_matching.csv`: requested tissue names and matched workbook tissue labels.",
  "- `tissue_color_key.csv`: exact tissue color hex codes used in the pie charts.",
  "- `report.txt`: console log from the pipeline run.",
  "",
  "## Largest slices",
  "",
  "Neat:",
  top_lines("Neat"),
  "",
  "Enriched:",
  top_lines("Enriched"),
  "",
  "Enriched only new:",
  top_lines("Enriched only new"),
  "",
  "Generated by `02_code/Plasma_Protein/organ_tissue_complete_pies.R`."
)
writeLines(readme_lines, file.path(out_dir, "README.md"))

cat("\nOutput files written under:\n")
print(out_dir)
cat("\nPipeline completed successfully.\n")
