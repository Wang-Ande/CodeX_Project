library(dplyr)
library(ggplot2)
library(ggalluvial)
library(readxl)

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1])
if (is.na(script_file) || !nzchar(script_file)) {
  script_file <- file.path(getwd(), "02_code", "Plasma_Protein", "organ_to_neat_sankey.R")
}
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
data_dir <- file.path(repo_dir, "01_data", "Plasma_Protein")
out_dir <- file.path(repo_dir, "03_result", "Plasma_Protein", "Organ_neat_sankey")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Organ-to-neat Sankey pipeline\n")
cat("Run time:", format(Sys.time()), "\n")
cat("Working directory:", getwd(), "\n")
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

normalize_label <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "", trimws(as.character(x))))
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
numeric_cols <- names(total_protein)[vapply(total_protein, is.numeric, logical(1))]
cat("\nNumeric columns used for neat/full mean abundance:\n")
print(numeric_cols)

total_protein <- total_protein %>%
  mutate(
    Full_mean = rowMeans(across(all_of(numeric_cols)), na.rm = TRUE),
    Full_mean = ifelse(is.nan(Full_mean), NA, Full_mean)
  )

cat("\nInput dimensions:\n")
cat("Low abundance rows:", nrow(low_abundance), "columns:", ncol(low_abundance), "\n")
cat("Total protein rows:", nrow(total_protein), "columns:", ncol(total_protein), "\n")
cat("Duplicated low abundance genes:", sum(duplicated(low_abundance$Gene)), "\n")
cat("Duplicated total protein genes:", sum(duplicated(total_protein$Gene)), "\n")

full_genes <- unique(total_protein$Gene)
low_genes <- unique(low_abundance$Gene)
overlap_genes <- intersect(full_genes, low_genes)

overlap_rank <- data.frame(
  Gene = overlap_genes,
  Full_log2 = log2(total_protein$Full_mean[match(overlap_genes, total_protein$Gene)] + 1),
  Low_log2 = log2(low_abundance$Low_abundance_15ul[match(overlap_genes, low_abundance$Gene)] + 1),
  stringsAsFactors = FALSE
)

overlap_rank <- overlap_rank %>%
  mutate(
    Full_quantile = ntile(Full_log2, 4),
    Low_quantile = ntile(Low_log2, 4)
  )

full_only <- setdiff(full_genes, low_genes)
low_only <- setdiff(low_genes, full_genes)

overlap_data <- overlap_rank %>%
  mutate(
    Full_quantile = factor(Full_quantile, levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4")),
    Low_quantile = factor(Low_quantile, levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4")),
    Type = "Overlap"
  ) %>%
  select(Gene, Full_quantile, Low_quantile, Type)

full_only_data <- total_protein %>%
  filter(Gene %in% full_only) %>%
  mutate(
    Full_quantile = ntile(log2(Full_mean + 1), 4),
    Full_quantile = factor(Full_quantile, levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4")),
    Low_quantile = factor("Not detected", levels = c("Q4", "Q3", "Q2", "Q1", "Not detected")),
    Type = "Full_only"
  ) %>%
  select(Gene, Full_quantile, Low_quantile, Type)

low_only_data <- low_abundance %>%
  filter(Gene %in% low_only) %>%
  mutate(
    Low_quantile = ntile(log2(Low_abundance_15ul + 1), 4),
    Low_quantile = factor(Low_quantile, levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4")),
    Full_quantile = factor("Not detected", levels = c("Q4", "Q3", "Q2", "Q1", "Not detected")),
    Type = "Low_only"
  ) %>%
  select(Gene, Full_quantile, Low_quantile, Type)

protein_quantile_map <- bind_rows(overlap_data, full_only_data, low_only_data) %>%
  mutate(
    Neat_quantile = factor(
      as.character(Full_quantile),
      levels = c("Q4", "Q3", "Q2", "Q1", "Not detected")
    ),
    Low_quantile = factor(
      as.character(Low_quantile),
      levels = c("Q4", "Q3", "Q2", "Q1", "Not detected")
    )
  ) %>%
  select(Gene, Type, Neat_quantile, Low_quantile)

cat("\nProtein universe summary:\n")
cat("Full/neat genes:", length(full_genes), "\n")
cat("Low abundance genes:", length(low_genes), "\n")
cat("Overlap genes:", length(overlap_genes), "\n")
cat("Full/neat only genes:", length(full_only), "\n")
cat("Low only genes, not detected in neat/full:", length(low_only), "\n")

excel_sheets <- excel_sheets(required_files[["specific_tissue"]])
cat("\nSpecific tissue workbook sheets:\n")
print(excel_sheets)

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
  Organ = organs,
  Organ_norm = normalize_label(organs),
  Organ_order = seq_along(organs),
  stringsAsFactors = FALSE
)

tissue_lookup <- tissue_all %>%
  transmute(Tissue_original = as.character(Tissue), Organ_norm = normalize_label(Tissue)) %>%
  distinct()

organ_name_matching <- requested_organs %>%
  left_join(tissue_lookup, by = "Organ_norm") %>%
  arrange(Organ_order, Tissue_original)

write.csv(
  organ_name_matching,
  file.path(out_dir, "organ_name_matching.csv"),
  row.names = FALSE
)

missing_organs <- organ_name_matching %>%
  filter(is.na(Tissue_original)) %>%
  distinct(Organ) %>%
  pull(Organ)

cat("\nRequested organ matching:\n")
cat("Requested organs:", length(organs), "\n")
cat("Matched requested organs:", length(setdiff(organs, missing_organs)), "\n")
if (length(missing_organs) > 0) {
  cat("Missing requested organs after normalized matching:", paste(missing_organs, collapse = " | "), "\n")
}

tissue_selected <- tissue_all %>%
  mutate(
    Tissue_original = as.character(Tissue),
    Organ_norm = normalize_label(Tissue_original)
  ) %>%
  inner_join(requested_organs, by = "Organ_norm") %>%
  mutate(
    Gene = as.character(Gene),
    Function = as.character(Function),
    Source_sheet = as.character(Source_sheet)
  ) %>%
  group_by(Gene, Organ, Organ_order) %>%
  summarise(
    Tissue_original = paste(sort(unique(Tissue_original)), collapse = "; "),
    Function = paste(sort(unique(Function)), collapse = "; "),
    Source_sheets = paste(sort(unique(Source_sheet)), collapse = "; "),
    .groups = "drop"
  )

gene_assignments <- tissue_selected %>%
  inner_join(protein_quantile_map, by = "Gene") %>%
  mutate(
    Organ = factor(Organ, levels = organs),
    Neat_quantile = factor(as.character(Neat_quantile), levels = c("Q4", "Q3", "Q2", "Q1", "Not detected")),
    Low_quantile = factor(as.character(Low_quantile), levels = c("Q4", "Q3", "Q2", "Q1", "Not detected"))
  ) %>%
  arrange(Organ_order, Neat_quantile, Low_quantile, Gene)

unmatched_tissue_genes <- setdiff(unique(tissue_selected$Gene), unique(protein_quantile_map$Gene))

cat("\nSpecific tissue selection summary:\n")
cat("Selected requested organ gene-organ assignments:", nrow(tissue_selected), "\n")
cat("Selected requested organ unique genes:", n_distinct(tissue_selected$Gene), "\n")
cat("Assignments matched to neat/low protein universe:", nrow(gene_assignments), "\n")
cat("Unique selected tissue genes matched to protein universe:", n_distinct(gene_assignments$Gene), "\n")
cat("Selected tissue genes absent from protein universe:", length(unmatched_tissue_genes), "\n")

write.csv(
  gene_assignments,
  file.path(out_dir, "selected_organ_neat_low_gene_assignments.csv"),
  row.names = FALSE
)

organ_neat_counts <- gene_assignments %>%
  count(Organ, Organ_order, Neat_quantile, name = "gene_organ_assignments") %>%
  group_by(Organ, Organ_order, Neat_quantile) %>%
  mutate(unique_genes = gene_organ_assignments) %>%
  ungroup() %>%
  arrange(Organ_order, Neat_quantile)

organ_neat_low_counts <- gene_assignments %>%
  count(Organ, Organ_order, Neat_quantile, Low_quantile, name = "gene_organ_assignments") %>%
  group_by(Organ, Organ_order, Neat_quantile, Low_quantile) %>%
  mutate(unique_genes = gene_organ_assignments) %>%
  ungroup() %>%
  arrange(Organ_order, Neat_quantile, Low_quantile)

summary_metrics <- data.frame(
  metric = c(
    "full_neat_genes",
    "low_abundance_genes",
    "overlap_genes",
    "full_neat_only_genes",
    "low_only_not_detected_in_neat_genes",
    "requested_organs",
    "matched_requested_organs",
    "selected_requested_organ_gene_assignments",
    "selected_requested_organ_unique_genes",
    "matched_gene_organ_assignments",
    "matched_unique_genes",
    "selected_tissue_genes_absent_from_protein_universe"
  ),
  value = c(
    length(full_genes),
    length(low_genes),
    length(overlap_genes),
    length(full_only),
    length(low_only),
    length(organs),
    length(setdiff(organs, missing_organs)),
    nrow(tissue_selected),
    n_distinct(tissue_selected$Gene),
    nrow(gene_assignments),
    n_distinct(gene_assignments$Gene),
    length(unmatched_tissue_genes)
  )
)

write.csv(organ_neat_counts, file.path(out_dir, "organ_to_neat_counts.csv"), row.names = FALSE)
write.csv(organ_neat_low_counts, file.path(out_dir, "organ_neat_low_counts.csv"), row.names = FALSE)
write.csv(summary_metrics, file.path(out_dir, "summary_metrics.csv"), row.names = FALSE)

cat("\nTop organ-to-neat counts:\n")
print(head(organ_neat_counts[order(-organ_neat_counts$gene_organ_assignments), ], 20))

flow_colors <- c(
  "Q4" = "#1F4E6D",
  "Q3" = "#4F86A5",
  "Q2" = "#90B4C8",
  "Q1" = "#B0C8D8",
  "Not detected" = "#BDBDBD"
)

plot_theme <- theme_bw(base_size = 11) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  )

organ_neat_plot <- ggplot(
  organ_neat_counts,
  aes(axis1 = Organ, axis2 = Neat_quantile, y = gene_organ_assignments)
) +
  geom_alluvium(aes(fill = Neat_quantile), alpha = 0.8, width = 1 / 12) +
  geom_stratum(width = 1 / 12, fill = "grey92", color = "black", linewidth = 0.25) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.1) +
  scale_x_discrete(limits = c("Organ category", "Neat quantile"), expand = c(0.08, 0.08)) +
  scale_fill_manual(values = flow_colors, drop = FALSE) +
  labs(
    title = "Organ category to neat protein abundance category",
    subtitle = "Counts are unique gene-organ assignments from the requested organ list",
    x = NULL
  ) +
  plot_theme

pdf(file.path(out_dir, "organ_to_neat_sankey.pdf"), width = 8.5, height = 16)
print(organ_neat_plot)
dev.off()

combined_plot <- ggplot(
  organ_neat_low_counts,
  aes(
    axis1 = Organ,
    axis2 = Neat_quantile,
    axis3 = Low_quantile,
    y = gene_organ_assignments
  )
) +
  geom_alluvium(aes(fill = Neat_quantile), alpha = 0.78, width = 1 / 12) +
  geom_stratum(width = 1 / 12, fill = "grey92", color = "black", linewidth = 0.25) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 1.95) +
  scale_x_discrete(
    limits = c("Organ category", "Neat quantile", "Low quantile"),
    expand = c(0.06, 0.06)
  ) +
  scale_fill_manual(values = flow_colors, drop = FALSE) +
  labs(
    title = "Organ category to neat and low-abundance protein categories",
    subtitle = "Left panel uses only the requested organ list from Specific tissue cluster2.xlsx",
    x = NULL
  ) +
  plot_theme

pdf(file.path(out_dir, "combined_organ_neat_low_sankey.pdf"), width = 11, height = 16)
print(combined_plot)
dev.off()

readme_lines <- c(
  "# Organ to neat protein Sankey analysis",
  "",
  "This folder contains the added Sankey analysis for the plasma protein data.",
  "",
  "## Inputs",
  "",
  "- `Low_abundance_15ul_tidy.csv`: low-abundance protein table used to define the low-abundance quantile and the proteins not detected in the neat/full protein table.",
  "- `Plasma_Total_protein_tidy.csv`: neat/full plasma protein table. Numeric sample columns are averaged per gene, then used to define the current neat-side `Q1`-`Q4` categories.",
  "- `Specific tissue cluster2.xlsx`: tissue-specific gene table. The workbook is treated as a long table with `Gene`, `Tissue`, and `Function` columns. Only the requested organ list is used. Organ matching ignores spaces and punctuation so workbook labels such as `Bonemarrow` can match requested labels such as `Bone marrow`.",
  "",
  "## Outputs",
  "",
  "- `organ_name_matching.csv`: how each requested organ name matched the tissue names in `Specific tissue cluster2.xlsx`.",
  "- `selected_organ_neat_low_gene_assignments.csv`: one row per matched gene-organ assignment, with source tissue annotation, neat quantile, low-abundance quantile, and overlap/full-only/low-only type.",
  "- `organ_to_neat_counts.csv`: count table for the standalone `Organ category -> Neat quantile` Sankey plot.",
  "- `organ_neat_low_counts.csv`: count table for the combined `Organ category -> Neat quantile -> Low quantile` Sankey plot.",
  "- `summary_metrics.csv`: compact metrics describing input sizes, matched organs, and matched genes.",
  "- `organ_to_neat_sankey.pdf`: standalone Sankey plot showing the requested organ categories flowing into the existing neat-side categories `Q4`, `Q3`, `Q2`, `Q1`, and `Not detected`. This PDF is generated from `Specific tissue cluster2.xlsx` restricted to the requested organs, joined to the neat-side categories derived from `Plasma_Total_protein_tidy.csv` and `Low_abundance_15ul_tidy.csv`.",
  "- `combined_organ_neat_low_sankey.pdf`: three-axis Sankey plot adding the organ category panel to the existing neat-to-low relationship. It uses the same requested organ-specific gene assignments, then maps each assignment to the current neat-side category and low-abundance category.",
  "- `report.txt`: appended console log from pipeline runs.",
  "",
  "## Run summary",
  "",
  paste0("- Requested organs: ", length(organs), ". Matched organs in `Specific tissue cluster2.xlsx`: ", length(setdiff(organs, missing_organs)), "."),
  paste0("- Missing requested organs: ", ifelse(length(missing_organs) > 0, paste(missing_organs, collapse = ", "), "None"), "."),
  paste0("- Matched gene-organ assignments used for plotting: ", nrow(gene_assignments), "."),
  paste0("- Unique matched genes used for plotting: ", n_distinct(gene_assignments$Gene), "."),
  "",
  "## Interpretation notes",
  "",
  "- Counts represent gene-organ assignments, not mutually exclusive proteins. If a gene is annotated to multiple requested organs, it contributes once to each organ.",
  "- `Not detected` on the neat side means the gene appears in the low-abundance table but not in the neat/full total protein table, matching the left-side category logic in the original Sankey template.",
  "- The quantile construction follows the original `Sankey plot.R` template: overlap, neat-only, and low-only groups are categorized using the same group-specific `ntile` approach as the template.",
  "",
  "Generated by `02_code/Plasma_Protein/organ_to_neat_sankey.R`."
)

writeLines(readme_lines, file.path(out_dir, "README.md"))

cat("\nOutput files written:\n")
print(list.files(out_dir, full.names = FALSE))
cat("\nPipeline completed successfully.\n")
