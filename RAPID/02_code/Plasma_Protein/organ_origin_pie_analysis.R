library(dplyr)
library(ggplot2)
library(readxl)

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1])
if (is.na(script_file) || !nzchar(script_file)) {
  script_file <- file.path(getwd(), "02_code", "Plasma_Protein", "organ_origin_pie_analysis.R")
}
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
data_dir <- file.path(repo_dir, "01_data", "Plasma_Protein")
out_dir <- file.path(repo_dir, "03_result", "Plasma_Protein", "Organ_origin_pie_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Organ origin pie analysis pipeline\n")
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
overlap_keys <- intersect(neat_keys, enriched_keys)
enriched_only_keys <- setdiff(enriched_keys, neat_keys)

cat("\nInput gene-set summary:\n")
cat("Neat unique proteins:", length(neat_keys), "\n")
cat("Enriched unique proteins:", length(enriched_keys), "\n")
cat("Overlap proteins:", length(overlap_keys), "\n")
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

missing_organs <- organ_name_matching %>%
  filter(is.na(Tissue_original)) %>%
  distinct(Organ) %>%
  pull(Organ)

tissue_selected <- tissue_all %>%
  mutate(
    Tissue_original = as.character(Tissue),
    Organ_norm = normalize_label(Tissue_original)
  ) %>%
  inner_join(requested_organs, by = "Organ_norm") %>%
  mutate(
    Gene = as.character(Gene),
    Gene_key = normalize_gene(Gene),
    Function = as.character(Function),
    Source_sheet = as.character(Source_sheet)
  ) %>%
  filter(!is.na(Gene_key), nzchar(Gene_key)) %>%
  group_by(Gene_key, Organ, Organ_order) %>%
  summarise(
    Example_gene_label = dplyr::first(stats::na.omit(Gene)),
    Tissue_original = paste(sort(unique(Tissue_original)), collapse = "; "),
    Function = paste(sort(unique(Function)), collapse = "; "),
    Source_sheets = paste(sort(unique(Source_sheet)), collapse = "; "),
    .groups = "drop"
  )

all_set_genes <- bind_rows(
  neat_genes %>% mutate(Protein_set = "Neat"),
  enriched_genes %>% mutate(Protein_set = "Enriched"),
  enriched_genes %>% filter(Gene_key %in% enriched_only_keys) %>% mutate(Protein_set = "Enriched only new")
) %>%
  mutate(
    Protein_set = factor(
      Protein_set,
      levels = c("Enriched only new", "Neat", "Enriched")
    )
  )

gene_organ_summary <- tissue_selected %>%
  group_by(Gene_key) %>%
  summarise(
    matched_organ_count = n_distinct(Organ),
    matched_organs = paste(sort(unique(Organ)), collapse = "; "),
    .groups = "drop"
  )

gene_origin_status <- all_set_genes %>%
  left_join(gene_organ_summary, by = "Gene_key") %>%
  mutate(
    matched_organ_count = ifelse(is.na(matched_organ_count), 0L, matched_organ_count),
    matched_organs = ifelse(is.na(matched_organs), "", matched_organs),
    Organ_origin_status = ifelse(matched_organ_count > 0, "Matched selected organ", "No selected organ match")
  ) %>%
  arrange(Protein_set, Gene_key)

origin_summary <- gene_origin_status %>%
  count(Protein_set, Organ_origin_status, name = "protein_count") %>%
  group_by(Protein_set) %>%
  mutate(
    set_total = sum(protein_count),
    percent = 100 * protein_count / set_total
  ) %>%
  ungroup() %>%
  arrange(Protein_set, desc(Organ_origin_status))

summary_metrics <- origin_summary %>%
  filter(Organ_origin_status == "Matched selected organ") %>%
  transmute(
    Protein_set = as.character(Protein_set),
    total_proteins = set_total,
    organ_origin_proteins = protein_count,
    organ_origin_percent = percent
  ) %>%
  arrange(factor(Protein_set, levels = c("Enriched only new", "Neat", "Enriched")))

organ_assignment_by_set <- all_set_genes %>%
  inner_join(tissue_selected, by = "Gene_key") %>%
  group_by(Protein_set, Gene_key) %>%
  mutate(fractional_weight = 1 / n_distinct(Organ)) %>%
  ungroup() %>%
  mutate(
    Organ = factor(Organ, levels = organs),
    Protein_set = factor(as.character(Protein_set), levels = c("Enriched only new", "Neat", "Enriched"))
  ) %>%
  arrange(Protein_set, Organ, Gene_key)

organ_distribution <- organ_assignment_by_set %>%
  group_by(Protein_set, Organ) %>%
  summarise(
    fractional_count = sum(fractional_weight),
    unique_proteins_with_organ = n_distinct(Gene_key),
    .groups = "drop"
  ) %>%
  group_by(Protein_set) %>%
  mutate(
    organ_derived_fractional_total = sum(fractional_count),
    percent_of_organ_derived = 100 * fractional_count / organ_derived_fractional_total
  ) %>%
  ungroup() %>%
  arrange(Protein_set, desc(percent_of_organ_derived), Organ)

top_organs <- organ_distribution %>%
  filter(Protein_set %in% c("Neat", "Enriched")) %>%
  group_by(Organ) %>%
  summarise(combined_fractional_count = sum(fractional_count), .groups = "drop") %>%
  arrange(desc(combined_fractional_count), Organ) %>%
  slice_head(n = 12) %>%
  pull(Organ) %>%
  as.character()

organ_distribution_top <- organ_distribution %>%
  filter(Protein_set %in% c("Neat", "Enriched")) %>%
  mutate(Organ_plot = ifelse(as.character(Organ) %in% top_organs, as.character(Organ), "Other selected organs")) %>%
  group_by(Protein_set, Organ_plot) %>%
  summarise(
    fractional_count = sum(fractional_count),
    unique_proteins_with_organ = sum(unique_proteins_with_organ),
    .groups = "drop"
  ) %>%
  group_by(Protein_set) %>%
  mutate(
    organ_derived_fractional_total = sum(fractional_count),
    percent_of_organ_derived = 100 * fractional_count / organ_derived_fractional_total
  ) %>%
  ungroup() %>%
  arrange(Protein_set, desc(percent_of_organ_derived), Organ_plot)

write.csv(organ_name_matching, file.path(out_dir, "organ_name_matching.csv"), row.names = FALSE)
write.csv(gene_origin_status, file.path(out_dir, "gene_organ_origin_status.csv"), row.names = FALSE)
write.csv(tissue_selected, file.path(out_dir, "selected_organ_gene_annotations.csv"), row.names = FALSE)
write.csv(organ_assignment_by_set, file.path(out_dir, "gene_organ_fractional_assignments.csv"), row.names = FALSE)
write.csv(origin_summary, file.path(out_dir, "organ_origin_binary_summary.csv"), row.names = FALSE)
write.csv(summary_metrics, file.path(out_dir, "summary_metrics.csv"), row.names = FALSE)
write.csv(organ_distribution, file.path(out_dir, "organ_distribution_fractional_all_sets.csv"), row.names = FALSE)
write.csv(organ_distribution_top, file.path(out_dir, "neat_enriched_organ_distribution_top12.csv"), row.names = FALSE)

cat("\nMatched requested organs:", length(setdiff(organs, missing_organs)), "\n")
cat("Missing requested organs:", ifelse(length(missing_organs) > 0, paste(missing_organs, collapse = " | "), "None"), "\n")
cat("\nOrgan-origin summary:\n")
print(summary_metrics)

origin_colors <- c(
  "Matched selected organ" = "#5DA5A4",
  "No selected organ match" = "#D8D8D8"
)

origin_plot_data <- origin_summary %>%
  mutate(
    Organ_origin_status = factor(
      Organ_origin_status,
      levels = c("Matched selected organ", "No selected organ match")
    ),
    label = paste0(protein_count, "\n", format_percent(percent))
  )

origin_pie <- ggplot(origin_plot_data, aes(x = "", y = protein_count, fill = Organ_origin_status)) +
  geom_col(width = 1, color = "white", linewidth = 0.8) +
  coord_polar(theta = "y") +
  facet_wrap(~Protein_set, nrow = 1) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.4,
    color = "black",
    lineheight = 0.92
  ) +
  scale_fill_manual(values = origin_colors, drop = FALSE) +
  labs(
    title = "Organ-origin coverage",
    subtitle = "Protein counts and percentages are based on unique protein labels",
    fill = NULL
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    plot.margin = margin(10, 10, 12, 10)
  )

pdf(file.path(out_dir, "organ_origin_coverage_pies.pdf"), width = 9.2, height = 4.6)
print(origin_pie)
dev.off()

organ_palette <- c(
  "Brain" = "#EDA66E",
  "Cerebellum" = "#93CEDA",
  "Neuronal" = "#5678B8",
  "Oligodendrocytes" = "#E8B186",
  "Retina" = "#A7D8E2",
  "Choroid plexus" = "#6E8DC2",
  "Adrenal gland" = "#E59662",
  "Pituitary gland" = "#7FC2D0",
  "Parathyroid gland" = "#466BA8",
  "Bone marrow" = "#F0C19A",
  "Spleen" = "#B8E0E7",
  "Lymphoid tissue" = "#7E9BC9",
  "Thymus" = "#D98555",
  "Lung" = "#69B8C9",
  "Liver" = "#3E5F9A",
  "Stomach" = "#F4CEAD",
  "Small intestine" = "#CBE9EE",
  "Intestine" = "#90ABD3",
  "Pancreas" = "#DD9A70",
  "Salivary gland" = "#84CFDA",
  "Kidney" = "#5D83BC",
  "Breast" = "#EFAE7B",
  "Fallopian tube" = "#9DD4DF",
  "Epididymis" = "#6A89BD",
  "Spermatids" = "#D97650",
  "Placenta" = "#5FAEC2",
  "Skeletal muscle" = "#486CA1",
  "Heart muscle" = "#F2BA8F",
  "Smooth muscle" = "#B1DDE6",
  "Epithelium" = "#7894C7",
  "Squamous epithelium" = "#ECA16C",
  "Skin" = "#8CC8D5",
  "Ciliated cells" = "#5A75AC",
  "Connective tissue" = "#F6D6BD",
  "Other selected organs" = "#C9C9C9"
)

top_levels <- c(top_organs, "Other selected organs")
organ_top_plot_data <- organ_distribution_top %>%
  mutate(
    Organ_plot = factor(Organ_plot, levels = top_levels),
    label = ifelse(
      percent_of_organ_derived >= 5,
      paste0(format_percent(percent_of_organ_derived)),
      ""
    )
  )

organ_top_pie <- ggplot(organ_top_plot_data, aes(x = "", y = fractional_count, fill = Organ_plot)) +
  geom_col(width = 1, color = "white", linewidth = 0.55) +
  coord_polar(theta = "y") +
  facet_wrap(~Protein_set, nrow = 1) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.1,
    color = "black"
  ) +
  scale_fill_manual(values = organ_palette, drop = FALSE) +
  labs(
    title = "Organ distribution among organ-derived proteins",
    subtitle = "Neat and Enriched; top 12 organs by combined fractional count, remaining organs collapsed",
    fill = NULL
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 9.5, hjust = 0.5, color = "grey30"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.text = element_text(size = 8.5),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(linewidth = 0)))

pdf(file.path(out_dir, "neat_enriched_organ_distribution_top12_pies.pdf"), width = 10.8, height = 5.5)
print(organ_top_pie)
dev.off()

organ_all_plot_data <- organ_distribution %>%
  filter(Protein_set %in% c("Neat", "Enriched")) %>%
  mutate(
    Organ = factor(as.character(Organ), levels = organs),
    label = ifelse(
      percent_of_organ_derived >= 6,
      paste0(format_percent(percent_of_organ_derived)),
      ""
    )
  )

organ_all_pie <- ggplot(organ_all_plot_data, aes(x = "", y = fractional_count, fill = Organ)) +
  geom_col(width = 1, color = "white", linewidth = 0.45) +
  coord_polar(theta = "y") +
  facet_wrap(~Protein_set, nrow = 1) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 2.7,
    color = "black"
  ) +
  scale_fill_manual(values = organ_palette[organs], drop = FALSE) +
  labs(
    title = "All selected-organ distribution among organ-derived proteins",
    subtitle = "Fractional counting: one protein split equally across matched organs",
    fill = NULL
  ) +
  theme_void(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 9.2, hjust = 0.5, color = "grey30"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.text = element_text(size = 7.8),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(ncol = 2, override.aes = list(linewidth = 0)))

pdf(file.path(out_dir, "neat_enriched_organ_distribution_all_organs_pies.pdf"), width = 11.5, height = 6.4)
print(organ_all_pie)
dev.off()

enriched_only_row <- summary_metrics %>%
  filter(Protein_set == "Enriched only new")

neat_row <- summary_metrics %>%
  filter(Protein_set == "Neat")

enriched_row <- summary_metrics %>%
  filter(Protein_set == "Enriched")

top_neat <- organ_distribution %>%
  filter(Protein_set == "Neat") %>%
  arrange(desc(percent_of_organ_derived)) %>%
  slice_head(n = 5)

top_enriched <- organ_distribution %>%
  filter(Protein_set == "Enriched") %>%
  arrange(desc(percent_of_organ_derived)) %>%
  slice_head(n = 5)

readme_lines <- c(
  "# Organ origin pie analysis",
  "",
  "This folder summarizes whether proteins are matched to the selected organs in `Specific tissue cluster2.xlsx` and visualizes organ-origin proportions as pie charts.",
  "",
  "## Inputs",
  "",
  "- `Low_abundance_15ul_tidy.csv`: enriched low-abundance protein table. Unique non-missing `Gene` labels define the Enriched set.",
  "- `Plasma_Total_protein_tidy.csv`: neat plasma protein table. Unique `Gene` labels define the Neat set.",
  "- `Specific tissue cluster2.xlsx`: tissue-specific annotation table. All sheets are read; only the selected organs used in the Sankey analyses are counted.",
  "",
  "## Counting rules",
  "",
  "- Protein labels are trimmed and matched case-insensitively using normalized `Gene` labels.",
  "- `Enriched only new` means proteins in Enriched but not in Neat.",
  "- A protein is counted as organ-derived if it matches at least one selected organ.",
  "- For organ distribution pies, each protein contributes total weight 1. If it matches multiple organs, the weight is split equally across those organs.",
  "",
  "## Main result",
  "",
  paste0(
    "- Enriched-only new proteins: ",
    enriched_only_row$total_proteins,
    " total; ",
    enriched_only_row$organ_origin_proteins,
    " matched selected organs (",
    format_percent(enriched_only_row$organ_origin_percent),
    ")."
  ),
  paste0(
    "- Neat proteins: ",
    neat_row$total_proteins,
    " total; ",
    neat_row$organ_origin_proteins,
    " matched selected organs (",
    format_percent(neat_row$organ_origin_percent),
    ")."
  ),
  paste0(
    "- Enriched proteins: ",
    enriched_row$total_proteins,
    " total; ",
    enriched_row$organ_origin_proteins,
    " matched selected organs (",
    format_percent(enriched_row$organ_origin_percent),
    ")."
  ),
  "",
  "## Outputs",
  "",
  "- `organ_origin_coverage_pies.pdf`: binary pie charts showing matched selected organ vs no selected organ match for Enriched-only new, Neat, and Enriched proteins.",
  "- `neat_enriched_organ_distribution_top12_pies.pdf`: Neat and Enriched organ-distribution pies using fractional counting; top 12 organs by combined fractional count are shown, with the rest collapsed as `Other selected organs`.",
  "- `neat_enriched_organ_distribution_all_organs_pies.pdf`: Neat and Enriched organ-distribution pies using all selected organs. Small slices are shown mainly through the legend.",
  "- `summary_metrics.csv`: total proteins, organ-derived proteins, and organ-derived percentages for Enriched-only new, Neat, and Enriched.",
  "- `organ_origin_binary_summary.csv`: full binary counts and percentages used by `organ_origin_coverage_pies.pdf`.",
  "- `organ_distribution_fractional_all_sets.csv`: fractional organ distribution for Enriched-only new, Neat, and Enriched.",
  "- `neat_enriched_organ_distribution_top12.csv`: collapsed table used by `neat_enriched_organ_distribution_top12_pies.pdf`.",
  "- `gene_organ_origin_status.csv`: one row per protein per set with matched organ status and matched organ list.",
  "- `gene_organ_fractional_assignments.csv`: one row per protein-organ assignment with fractional weight.",
  "- `selected_organ_gene_annotations.csv`: selected organ annotations from `Specific tissue cluster2.xlsx` after organ-name matching.",
  "- `organ_name_matching.csv`: requested organ names and matched tissue labels from the workbook.",
  "- `report.txt`: console log from the pipeline run.",
  "",
  "## Top organ distribution",
  "",
  "Top Neat organs among organ-derived proteins:",
  paste0("- ", as.character(top_neat$Organ), ": ", format_percent(top_neat$percent_of_organ_derived)),
  "",
  "Top Enriched organs among organ-derived proteins:",
  paste0("- ", as.character(top_enriched$Organ), ": ", format_percent(top_enriched$percent_of_organ_derived)),
  "",
  "Generated by `02_code/Plasma_Protein/organ_origin_pie_analysis.R`."
)
writeLines(readme_lines, file.path(out_dir, "README.md"))

cat("\nOutput files written under:\n")
print(out_dir)
cat("\nPipeline completed successfully.\n")
