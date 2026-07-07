library(dplyr)
library(ggplot2)
library(ggalluvial)

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1])
if (is.na(script_file) || !nzchar(script_file)) {
  script_file <- file.path(getwd(), "02_code", "Plasma_Protein", "abundance_sankey_color_optimized.R")
}
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
data_dir <- file.path(repo_dir, "01_data", "Plasma_Protein")
out_dir <- file.path(repo_dir, "03_result", "Plasma_Protein", "Sankey_color_optimized")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Abundance Sankey color-optimized pipeline\n")
cat("Run time:", format(Sys.time()), "\n")
cat("Repository directory:", repo_dir, "\n")
cat("Data directory:", data_dir, "\n")
cat("Output directory:", out_dir, "\n")

required_files <- c(
  low_abundance = file.path(data_dir, "Low_abundance_15ul_tidy.csv"),
  total_protein = file.path(data_dir, "Plasma_Total_protein_tidy.csv")
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

full_genes <- unique(total_protein$Gene)
low_genes <- unique(low_abundance$Gene)
overlap_genes <- intersect(full_genes, low_genes)
full_only <- setdiff(full_genes, low_genes)
low_only <- setdiff(low_genes, full_genes)

overlap_rank <- data.frame(
  Gene = overlap_genes,
  Full_log2 = log2(total_protein$Full_mean[match(overlap_genes, total_protein$Gene)] + 1),
  Low_log2 = log2(low_abundance$Low_abundance_15ul[match(overlap_genes, low_abundance$Gene)] + 1),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Full_quantile = ntile(Full_log2, 4),
    Low_quantile = ntile(Low_log2, 4)
  )

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
    Full_quantile = factor(as.character(Full_quantile), levels = c("Q4", "Q3", "Q2", "Q1", "Not detected")),
    Low_quantile = factor(as.character(Low_quantile), levels = c("Q4", "Q3", "Q2", "Q1", "Not detected")),
    Flow_category = ifelse(Type == "Overlap", as.character(Full_quantile), "Not detected"),
    Flow_category = factor(Flow_category, levels = c("Q4", "Q3", "Q2", "Q1", "Not detected"))
  )

abundance_counts <- protein_quantile_map %>%
  count(Full_quantile, Low_quantile, Type, Flow_category, name = "protein_count") %>%
  arrange(Full_quantile, Low_quantile, Type)

axis_totals <- bind_rows(
  protein_quantile_map %>%
    count(Full_quantile, name = "axis_total") %>%
    transmute(axis = "Neat quantile", category = as.character(Full_quantile), axis_total),
  protein_quantile_map %>%
    count(Low_quantile, name = "axis_total") %>%
    transmute(axis = "Low quantile", category = as.character(Low_quantile), axis_total)
) %>%
  mutate(category = factor(category, levels = c("Q4", "Q3", "Q2", "Q1", "Not detected"))) %>%
  arrange(axis, category)

alignment_check <- axis_totals %>%
  filter(axis == "Neat quantile") %>%
  transmute(
    Neat_quantile = category,
    abundance_neat_axis_total = axis_total
  )

summary_metrics <- data.frame(
  metric = c(
    "full_neat_genes",
    "low_abundance_genes",
    "overlap_genes",
    "full_neat_only_genes",
    "low_only_not_detected_in_neat_genes",
    "abundance_sankey_total_height"
  ),
  value = c(
    length(full_genes),
    length(low_genes),
    length(overlap_genes),
    length(full_only),
    length(low_only),
    sum(abundance_counts$protein_count)
  )
)

write.csv(protein_quantile_map, file.path(out_dir, "abundance_gene_quantile_map.csv"), row.names = FALSE)
write.csv(abundance_counts, file.path(out_dir, "abundance_sankey_counts.csv"), row.names = FALSE)
write.csv(axis_totals, file.path(out_dir, "abundance_axis_totals.csv"), row.names = FALSE)
write.csv(alignment_check, file.path(out_dir, "abundance_neat_axis_alignment.csv"), row.names = FALSE)
write.csv(summary_metrics, file.path(out_dir, "abundance_summary_metrics.csv"), row.names = FALSE)

cat("\nProtein universe summary:\n")
cat("Full/neat genes:", length(full_genes), "\n")
cat("Low abundance genes:", length(low_genes), "\n")
cat("Overlap genes:", length(overlap_genes), "\n")
cat("Full/neat only genes:", length(full_only), "\n")
cat("Low only genes, not detected in neat/full:", length(low_only), "\n")
cat("\nAxis totals:\n")
print(axis_totals)

plot_data <- abundance_counts %>%
  mutate(alluvium = row_number()) %>%
  select(alluvium, Full_quantile, Low_quantile, Type, Flow_category, protein_count)

lodes <- to_lodes_form(
  plot_data,
  axes = c("Full_quantile", "Low_quantile"),
  key = "axis",
  id = "alluvium"
) %>%
  mutate(
    axis_x = ifelse(axis == "Full_quantile", 1, 2),
    stratum = factor(as.character(stratum), levels = c("Q4", "Q3", "Q2", "Q1", "Not detected"))
  )

matched_colors <- c(
  "Q4" = "#3B2F70",
  "Q3" = "#2F6F9E",
  "Q2" = "#5AA6A6",
  "Q1" = "#B7D7C2",
  "Not detected" = "#D1D1D1"
)
flow_alpha <- c(
  "Q4" = 0.74,
  "Q3" = 0.72,
  "Q2" = 0.68,
  "Q1" = 0.64,
  "Not detected" = 0.46
)
label_map <- c(
  "Q4" = "Q4",
  "Q3" = "Q3",
  "Q2" = "Q2",
  "Q1" = "Q1",
  "Not detected" = "NA"
)

base_plot <- ggplot(
  lodes,
  aes(
    x = axis_x,
    stratum = stratum,
    alluvium = alluvium,
    y = protein_count
  )
) +
  geom_alluvium(
    aes(fill = Flow_category, alpha = Flow_category),
    width = 0.13,
    knot.pos = 0.35,
    curve_type = "xspline",
    show.legend = FALSE
  ) +
  geom_stratum(
    aes(fill = after_stat(stratum)),
    width = 0.15,
    color = "black",
    linewidth = 0.32,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = matched_colors, drop = FALSE) +
  scale_alpha_manual(values = flow_alpha, drop = FALSE, guide = "none") +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Global proteins", "Enriched proteins"),
    limits = c(0.84, 2.24),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 12, color = "black", margin = margin(t = 4)),
    plot.margin = margin(8, 70, 18, 14)
  )

strata_coords <- ggplot_build(base_plot)$data[[2]]
if (!all(c("x", "y", "stratum") %in% names(strata_coords))) {
  stop("Could not recover stratum coordinates from ggplot_build().")
}

left_labels <- strata_coords %>%
  filter(abs(x - 1) < 1e-6) %>%
  mutate(label = unname(label_map[as.character(stratum)]), x_label = 0.92)

right_labels <- strata_coords %>%
  filter(abs(x - 2) < 1e-6) %>%
  mutate(label = unname(label_map[as.character(stratum)]), x_label = 2.08)

abundance_plot <- base_plot +
  geom_text(
    data = left_labels,
    inherit.aes = FALSE,
    aes(x = x_label, y = y, label = label),
    hjust = 1,
    size = 3.8,
    color = "black"
  ) +
  geom_text(
    data = right_labels,
    inherit.aes = FALSE,
    aes(x = x_label, y = y, label = label),
    hjust = 0,
    size = 3.8,
    color = "black"
  )

pdf(file.path(out_dir, "abundance_sankey_color_optimized.pdf"), width = 5.2, height = 8.8)
print(abundance_plot)
dev.off()

readme_lines <- c(
  "# Sankey color optimized",
  "",
  "This folder contains color-optimized Sankey panels for the tissue-derived protein, global protein, and enriched protein composition figure.",
  "",
  "## Inputs",
  "",
  "- `Low_abundance_15ul_tidy.csv`: low-abundance protein table.",
  "- `Plasma_Total_protein_tidy.csv`: neat/full plasma protein table.",
  "- `Specific tissue cluster2.xlsx`: tissue-specific gene table used by the left organ panel.",
  "",
  "## Outputs",
  "",
  "- `organ_to_neat_aligned_panel.pdf`: left panel, showing tissue-derived organ categories flowing into aligned Global protein Q4/Q3/Q2/Q1/NA categories. It is based on `Specific tissue cluster2.xlsx`, `Plasma_Total_protein_tidy.csv`, and the abundance target proportions derived from the full abundance Sankey universe.",
  "- `organ_to_neat_aligned_counts.csv`: organ-to-global category raw counts plus `aligned_display_weight`, which is the plotted value used to make the Global protein axis align with the abundance panel.",
  "- `organ_system_color_key.csv`: organ-to-system grouping and left-module color key.",
  "- `selected_organ_neat_gene_assignments.csv`: gene-level organ assignments used for the left panel.",
  "- `organ_name_matching.csv`: requested organ matching against `Specific tissue cluster2.xlsx`.",
  "- `alignment_check.csv`: confirms the organ panel's Global protein axis totals match the abundance-panel Global protein axis totals.",
  "- `organ_summary_metrics.csv`: compact metrics for the organ-side panel.",
  "- `abundance_sankey_color_optimized.pdf`: abundance Sankey with the same PDF height and the same Q4/Q3/Q2/Q1/NA palette as the color-optimized organ panel.",
  "- `abundance_sankey_counts.csv`: aggregated flow counts used for plotting.",
  "- `abundance_axis_totals.csv`: left and right axis totals.",
  "- `abundance_neat_axis_alignment.csv`: left neat-axis totals to compare with the organ aligned panel.",
  "- `abundance_gene_quantile_map.csv`: gene-level category map.",
  "- `abundance_summary_metrics.csv`: compact metrics for the abundance panel.",
  "- `README_organ_panel.md`: organ-side panel notes written by the organ script.",
  "- `report.txt`: console log from the color-optimized organ and abundance pipeline runs.",
  "",
  "## Style Notes",
  "",
  "- PDF height is 8.8 inches, matching the aligned organ panel.",
  "- Colors use the optimized abundance palette: Q4 `#3B2F70`, Q3 `#2F6F9E`, Q2 `#5AA6A6`, Q1 `#B7D7C2`, Not detected/NA `#D1D1D1`.",
  "- The plotted label `NA` is used for the original `Not detected` category to avoid label clipping during figure assembly.",
  "- Q4/Q3/Q2/Q1 flows use alpha 0.74/0.72/0.68/0.64. NA flows use alpha 0.46 so the large undetected component remains visible on a white page without overpowering the quantile flows.",
  "- Organ modules are colored by organ system with low-saturation colors, while flows and Global/Enriched modules use the abundance palette.",
  "- Overlap flows are colored by the source Global protein quantile. Full-only and low-only flows are rendered in the same neutral grey used for `Not detected`.",
  "",
  "Generated by `02_code/Plasma_Protein/abundance_sankey_color_optimized.R`."
)

writeLines(readme_lines, file.path(out_dir, "README.md"))

cat("\nOutput files written:\n")
print(list.files(out_dir, full.names = FALSE))
cat("\nPipeline completed successfully.\n")
