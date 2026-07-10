library(dplyr)
library(ggplot2)
library(readxl)

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1])
if (is.na(script_file) || !nzchar(script_file)) {
  script_file <- file.path(getwd(), "02_code", "Plasma_Protein", "enriched_tissue_detection_pie.R")
}
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
data_dir <- file.path(repo_dir, "01_data", "Plasma_Protein")
out_dir <- file.path(repo_dir, "03_result", "Plasma_Protein", "Organ_tissue_complete_pies")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Enriched tissue detection pie pipeline\n")
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

cat("\nInput gene-set summary:\n")
cat("Neat unique proteins:", length(neat_keys), "\n")
cat("Enriched unique proteins:", length(enriched_keys), "\n")
cat("Enriched-only proteins:", length(setdiff(enriched_keys, neat_keys)), "\n")

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

enriched_tissue_gene_status <- enriched_genes %>%
  inner_join(
    tissue_selected %>%
      group_by(Gene_key) %>%
      summarise(
        matched_tissue_count = n_distinct(Tissue),
        matched_tissues = paste(sort(unique(Tissue)), collapse = "; "),
        .groups = "drop"
      ),
    by = "Gene_key"
  ) %>%
  mutate(
    Detection_status = ifelse(Gene_key %in% neat_keys, "Detected in Neat", "New in Enriched"),
    Detection_status = factor(Detection_status, levels = c("Detected in Neat", "New in Enriched"))
  ) %>%
  arrange(Detection_status, Gene_key)

summary_table <- enriched_tissue_gene_status %>%
  count(Detection_status, name = "protein_count") %>%
  mutate(
    total_enriched_tissue_source_proteins = sum(protein_count),
    percent = 100 * protein_count / total_enriched_tissue_source_proteins
  )

write.csv(
  enriched_tissue_gene_status,
  file.path(out_dir, "enriched_tissue_detected_gene_status.csv"),
  row.names = FALSE
)
write.csv(
  summary_table,
  file.path(out_dir, "enriched_tissue_detected_in_neat_vs_new_summary.csv"),
  row.names = FALSE
)

cat("\nEnriched tissue-source detection summary:\n")
print(summary_table)

status_colors <- c(
  "Detected in Neat" = "#5E7892",
  "New in Enriched" = "#D68AA8"
)

plot_data <- summary_table %>%
  mutate(
    label = paste0(protein_count, "\n", format_percent(percent))
  )

total_n <- unique(plot_data$total_enriched_tissue_source_proteins)

pie_plot <- ggplot(plot_data, aes(x = "", y = protein_count, fill = Detection_status)) +
  geom_col(width = 1, color = "white", linewidth = 0.9) +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4.2,
    lineheight = 0.9,
    color = "black"
  ) +
  scale_fill_manual(values = status_colors, drop = FALSE) +
  labs(
    title = "Enriched tissue-source proteins",
    subtitle = paste0("Denominator: ", total_n, " Enriched proteins matched to selected tissues"),
    fill = NULL
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    plot.margin = margin(12, 12, 12, 12)
  )

pdf(file.path(out_dir, "enriched_tissue_detected_in_neat_vs_new_pie.pdf"), width = 6.2, height = 5.4)
print(pie_plot)
dev.off()

detected_row <- summary_table %>% filter(Detection_status == "Detected in Neat")
new_row <- summary_table %>% filter(Detection_status == "New in Enriched")

readme_lines <- c(
  "# Enriched tissue-source detected in Neat vs new",
  "",
  "This supplementary result adds one pie chart under `Organ_tissue_complete_pies`.",
  "",
  "## Question",
  "",
  "Among Enriched proteins that match at least one selected tissue, how many were already detected in Neat and how many are newly detected only in Enriched?",
  "",
  "## Counting rules",
  "",
  "- The denominator is unique Enriched proteins matched to at least one selected tissue in `Specific tissue cluster2.xlsx`.",
  "- `Detected in Neat` means the same normalized `Gene` label is present in `Plasma_Total_protein_tidy.csv`.",
  "- `New in Enriched` means the normalized `Gene` label is present in `Low_abundance_15ul_tidy.csv` but absent from `Plasma_Total_protein_tidy.csv`.",
  "- Protein labels are trimmed and matched case-insensitively.",
  "",
  "## Result",
  "",
  paste0("- Total Enriched tissue-source proteins: ", total_n),
  paste0("- Detected in Neat: ", detected_row$protein_count, " (", format_percent(detected_row$percent), ")"),
  paste0("- New in Enriched: ", new_row$protein_count, " (", format_percent(new_row$percent), ")"),
  "",
  "## Outputs",
  "",
  "- `enriched_tissue_detected_in_neat_vs_new_pie.pdf`: pie chart requested here.",
  "- `enriched_tissue_detected_in_neat_vs_new_summary.csv`: counts and percentages used by the pie chart.",
  "- `enriched_tissue_detected_gene_status.csv`: per-protein status table with matched tissue list.",
  "- `report.txt`: console log from this supplementary pipeline run.",
  "",
  "Generated by `02_code/Plasma_Protein/enriched_tissue_detection_pie.R`."
)
writeLines(readme_lines, file.path(out_dir, "README_enriched_tissue_detection_pie.md"))

cat("\nOutput files written under:\n")
print(out_dir)
cat("\nPipeline completed successfully.\n")
