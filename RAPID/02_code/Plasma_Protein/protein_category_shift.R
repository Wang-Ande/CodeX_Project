suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1])
if (is.na(script_file) || !nzchar(script_file)) {
  script_file <- file.path(getwd(), "02_code", "Plasma_Protein", "protein_category_shift.R")
}

script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
data_dir <- file.path(repo_dir, "01_data", "Plasma_Protein")
out_dir <- file.path(repo_dir, "03_result", "Plasma_Protein", "Protein_category_shift")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

old_outputs <- list.files(
  out_dir,
  pattern = "^(protein_category_|set_category_summary|summary_metrics|README).*[.](csv|pdf|md)$",
  full.names = TRUE
)
if (length(old_outputs) > 0) {
  unlink(old_outputs)
}

cat("Protein category shift pipeline\n")
cat("Run time:", format(Sys.time()), "\n")
cat("Repository directory:", repo_dir, "\n")
cat("Data directory:", data_dir, "\n")
cat("Output directory:", out_dir, "\n")

required_files <- c(
  neat = file.path(data_dir, "Plasma_Total_protein_tidy.csv"),
  enriched = file.path(data_dir, "Low_abundance_15ul_tidy.csv"),
  category = file.path(data_dir, "Protein category.xlsx")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files: ", paste(missing_files, collapse = "; "))
}

cat("\nInput files:\n")
print(required_files)

normalize_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  toupper(x)
}

percent_label <- function(x) {
  paste0(round(x * 100), "%")
}

pp_label <- function(x) {
  paste0(ifelse(x > 0, "+", ""), sprintf("%.2f", x))
}

neat_raw <- read.csv(required_files[["neat"]], row.names = 1, check.names = FALSE)
enriched_raw <- read.csv(required_files[["enriched"]], row.names = 1, check.names = FALSE)
enriched_raw <- enriched_raw[!is.na(enriched_raw$Low_abundance_15ul), , drop = FALSE]

if (!"Gene" %in% names(neat_raw)) {
  stop("Column `Gene` was not found in Plasma_Total_protein_tidy.csv.")
}
if (!"Gene" %in% names(enriched_raw)) {
  stop("Column `Gene` was not found in Low_abundance_15ul_tidy.csv.")
}

neat_genes <- data.frame(
  Set = "Neat",
  Gene = as.character(neat_raw$Gene),
  Gene_clean = normalize_gene(neat_raw$Gene),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Gene_clean)) %>%
  distinct(Set, Gene_clean, .keep_all = TRUE)

enriched_genes <- data.frame(
  Set = "Enriched",
  Gene = as.character(enriched_raw$Gene),
  Gene_clean = normalize_gene(enriched_raw$Gene),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Gene_clean)) %>%
  distinct(Set, Gene_clean, .keep_all = TRUE)

set_levels <- c("Neat", "Enriched")
gene_sets <- bind_rows(neat_genes, enriched_genes) %>%
  mutate(Set = factor(Set, levels = set_levels))

category_raw <- read_excel(required_files[["category"]], sheet = "Sheet1")
required_category_cols <- c("Gene", "Protein_Class")
missing_category_cols <- setdiff(required_category_cols, names(category_raw))
if (length(missing_category_cols) > 0) {
  stop("Missing required category columns: ", paste(missing_category_cols, collapse = ", "))
}

category_map <- category_raw %>%
  transmute(
    Gene = as.character(Gene),
    Gene_clean = normalize_gene(Gene),
    Protein_Class = trimws(as.character(Protein_Class)),
    Database = if ("Database" %in% names(category_raw)) trimws(as.character(Database)) else NA_character_
  ) %>%
  filter(!is.na(Gene_clean), !is.na(Protein_Class), Protein_Class != "") %>%
  distinct(Gene_clean, Protein_Class, .keep_all = TRUE)

annotation_raw <- gene_sets %>%
  left_join(
    category_map %>% select(Gene_clean, Protein_Class),
    by = "Gene_clean",
    relationship = "many-to-many"
  ) %>%
  mutate(
    Protein_Class = ifelse(is.na(Protein_Class) | Protein_Class == "", "Unclassified", Protein_Class)
  ) %>%
  distinct(Set, Gene, Gene_clean, Protein_Class)

gene_category_counts <- annotation_raw %>%
  group_by(Set, Gene, Gene_clean) %>%
  summarise(n_categories = n_distinct(Protein_Class), .groups = "drop")

gene_category_assignments <- annotation_raw %>%
  left_join(gene_category_counts, by = c("Set", "Gene", "Gene_clean")) %>%
  mutate(fractional_weight = 1 / n_categories) %>%
  arrange(Set, Gene_clean, Protein_Class)

total_gene_counts <- gene_sets %>%
  group_by(Set) %>%
  summarise(total_genes = n_distinct(Gene_clean), .groups = "drop")

all_categories <- sort(unique(gene_category_assignments$Protein_Class))
category_grid <- expand.grid(
  Set = set_levels,
  Protein_Class = all_categories,
  stringsAsFactors = FALSE
)

category_summary_raw <- gene_category_assignments %>%
  group_by(Set, Protein_Class) %>%
  summarise(
    fractional_gene_count = sum(fractional_weight),
    unique_gene_count = n_distinct(Gene_clean),
    .groups = "drop"
  )

category_summary <- category_grid %>%
  left_join(category_summary_raw, by = c("Set", "Protein_Class")) %>%
  left_join(total_gene_counts, by = "Set") %>%
  mutate(
    fractional_gene_count = ifelse(is.na(fractional_gene_count), 0, fractional_gene_count),
    unique_gene_count = ifelse(is.na(unique_gene_count), 0, unique_gene_count),
    percentage = 100 * fractional_gene_count / total_genes
  ) %>%
  arrange(Set, desc(percentage), Protein_Class)

neat_summary <- category_summary %>%
  filter(Set == "Neat") %>%
  transmute(
    Protein_Class,
    neat_fractional_gene_count = fractional_gene_count,
    neat_unique_gene_count = unique_gene_count,
    neat_total_genes = total_genes,
    neat_percentage = percentage
  )

enriched_summary <- category_summary %>%
  filter(Set == "Enriched") %>%
  transmute(
    Protein_Class,
    enriched_fractional_gene_count = fractional_gene_count,
    enriched_unique_gene_count = unique_gene_count,
    enriched_total_genes = total_genes,
    enriched_percentage = percentage
  )

delta_summary <- category_grid %>%
  select(Protein_Class) %>%
  distinct() %>%
  left_join(neat_summary, by = "Protein_Class") %>%
  left_join(enriched_summary, by = "Protein_Class") %>%
  mutate(
    across(where(is.numeric), ~ ifelse(is.na(.x), 0, .x)),
    delta_percentage_points = enriched_percentage - neat_percentage,
    ratio_enriched_to_neat = ifelse(neat_percentage > 0, enriched_percentage / neat_percentage, NA_real_),
    direction = ifelse(delta_percentage_points >= 0, "Higher in Enriched", "Higher in Neat")
  ) %>%
  arrange(desc(abs(delta_percentage_points)), Protein_Class)

category_gene_diagnostics <- gene_category_assignments %>%
  group_by(Set, Gene, Gene_clean) %>%
  summarise(
    categories = paste(sort(unique(Protein_Class)), collapse = "; "),
    n_categories = first(n_categories),
    .groups = "drop"
  ) %>%
  arrange(Set, desc(n_categories), Gene_clean)

set_diagnostics <- gene_sets %>%
  mutate(has_protein_category = Gene_clean %in% category_map$Gene_clean) %>%
  group_by(Set) %>%
  summarise(
    total_genes = n_distinct(Gene_clean),
    category_matched_genes = n_distinct(Gene_clean[has_protein_category]),
    unclassified_genes = n_distinct(Gene_clean[!has_protein_category]),
    .groups = "drop"
  )

multi_category_diagnostics <- category_gene_diagnostics %>%
  group_by(Set) %>%
  summarise(
    multi_category_genes = sum(n_categories > 1),
    max_categories_per_gene = max(n_categories),
    .groups = "drop"
  )

summary_metrics <- set_diagnostics %>%
  left_join(multi_category_diagnostics, by = "Set")

relationship_metrics <- data.frame(
  metric = c(
    "neat_genes",
    "enriched_genes",
    "overlap_genes",
    "neat_only_genes",
    "enriched_only_genes",
    "protein_category_rows",
    "protein_category_unique_genes",
    "protein_category_classes"
  ),
  value = c(
    n_distinct(neat_genes$Gene_clean),
    n_distinct(enriched_genes$Gene_clean),
    length(intersect(neat_genes$Gene_clean, enriched_genes$Gene_clean)),
    length(setdiff(neat_genes$Gene_clean, enriched_genes$Gene_clean)),
    length(setdiff(enriched_genes$Gene_clean, neat_genes$Gene_clean)),
    nrow(category_raw),
    n_distinct(category_map$Gene_clean),
    n_distinct(category_map$Protein_Class)
  )
)

write.csv(gene_category_assignments, file.path(out_dir, "protein_category_gene_assignments.csv"), row.names = FALSE)
write.csv(category_summary, file.path(out_dir, "protein_category_fractional_composition.csv"), row.names = FALSE)
write.csv(delta_summary, file.path(out_dir, "protein_category_delta.csv"), row.names = FALSE)
write.csv(category_gene_diagnostics, file.path(out_dir, "protein_category_gene_diagnostics.csv"), row.names = FALSE)
write.csv(summary_metrics, file.path(out_dir, "set_category_summary.csv"), row.names = FALSE)
write.csv(relationship_metrics, file.path(out_dir, "summary_metrics.csv"), row.names = FALSE)

cat("\nSet diagnostics:\n")
print(summary_metrics)
cat("\nRelationship metrics:\n")
print(relationship_metrics)
cat("\nTop category shifts by absolute delta percentage points:\n")
print(head(delta_summary, 15))

top_n <- min(14, sum(delta_summary$Protein_Class != "Unclassified"))
top_categories <- delta_summary %>%
  filter(Protein_Class != "Unclassified") %>%
  slice_head(n = top_n) %>%
  pull(Protein_Class)

plot_composition <- category_summary %>%
  mutate(
    Plot_category = case_when(
      Protein_Class == "Unclassified" ~ "Unclassified",
      Protein_Class %in% top_categories ~ Protein_Class,
      TRUE ~ "Other annotated categories"
    )
  ) %>%
  group_by(Set, Plot_category) %>%
  summarise(
    fractional_gene_count = sum(fractional_gene_count),
    unique_gene_count = sum(unique_gene_count),
    total_genes = first(total_genes),
    percentage = sum(percentage),
    .groups = "drop"
  )

plot_delta <- plot_composition %>%
  select(Set, Plot_category, percentage) %>%
  reshape(
    idvar = "Plot_category",
    timevar = "Set",
    direction = "wide"
  )
names(plot_delta) <- sub("^percentage\\.", "", names(plot_delta))
if (!"Neat" %in% names(plot_delta)) {
  plot_delta$Neat <- 0
}
if (!"Enriched" %in% names(plot_delta)) {
  plot_delta$Enriched <- 0
}
plot_delta <- plot_delta %>%
  mutate(delta_percentage_points = Enriched - Neat) %>%
  arrange(desc(ifelse(Plot_category == "Unclassified", Inf, abs(delta_percentage_points))))

stack_levels <- plot_delta %>%
  arrange(desc(Enriched), desc(Neat), Plot_category) %>%
  pull(Plot_category)
stack_levels <- c(setdiff(stack_levels, c("Other annotated categories", "Unclassified")), "Other annotated categories", "Unclassified")
stack_levels <- stack_levels[stack_levels %in% unique(plot_composition$Plot_category)]

plot_composition <- plot_composition %>%
  mutate(
    Set = factor(Set, levels = set_levels),
    Plot_category = factor(Plot_category, levels = rev(stack_levels))
  )

write.csv(plot_composition, file.path(out_dir, "protein_category_plot_composition.csv"), row.names = FALSE)

base_palette <- c(
  "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#86BCB6",
  "#2F4B7C", "#A05195", "#D45087", "#F95D6A", "#FFA600",
  "#665191", "#003F5C", "#7A5195", "#EF5675", "#6A994E",
  "#577590", "#BC5090", "#8CD17D", "#B6992D", "#499894"
)
named_categories <- setdiff(stack_levels, c("Other annotated categories", "Unclassified"))
category_colors <- setNames(rep(base_palette, length.out = length(named_categories)), named_categories)
category_colors <- c(
  category_colors,
  "Other annotated categories" = "#E0E0E0",
  "Unclassified" = "#BDBDBD"
)
category_colors <- category_colors[names(category_colors) %in% stack_levels]

stacked_plot <- ggplot(
  plot_composition,
  aes(x = Set, y = percentage / 100, fill = Plot_category)
) +
  geom_col(width = 0.58, color = "white", linewidth = 0.18) +
  scale_y_continuous(labels = percent_label, expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = category_colors, drop = FALSE) +
  labs(
    x = NULL,
    y = "Percentage of proteins",
    fill = "Protein category",
    title = "Protein category composition",
    subtitle = "Each protein contributes total weight 1; multi-category proteins are fractionally weighted."
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5, color = "grey35"),
    axis.text.x = element_text(size = 10.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.35, "cm"),
    legend.key.width = unit(0.35, "cm"),
    legend.position = "right",
    plot.margin = margin(8, 10, 8, 8)
  ) +
  guides(fill = guide_legend(reverse = TRUE, ncol = 1))

pdf(file.path(out_dir, "protein_category_stacked_bar.pdf"), width = 7.4, height = 5.2)
print(stacked_plot)
dev.off()

delta_plot_n <- min(22, sum(delta_summary$Protein_Class != "Unclassified"))
delta_plot_data <- delta_summary %>%
  filter(Protein_Class != "Unclassified") %>%
  slice_head(n = delta_plot_n) %>%
  arrange(delta_percentage_points) %>%
  mutate(
    Protein_Class = factor(Protein_Class, levels = Protein_Class),
    direction = factor(direction, levels = c("Higher in Enriched", "Higher in Neat"))
  )

delta_plot <- ggplot(
  delta_plot_data,
  aes(x = delta_percentage_points, y = Protein_Class, fill = direction)
) +
  geom_vline(xintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_col(width = 0.68, color = "white", linewidth = 0.15) +
  geom_text(
    aes(label = pp_label(delta_percentage_points)),
    hjust = ifelse(delta_plot_data$delta_percentage_points >= 0, -0.12, 1.12),
    size = 2.7,
    color = "black"
  ) +
  scale_fill_manual(
    values = c("Higher in Enriched" = "#2C7FB8", "Higher in Neat" = "#A67C52"),
    drop = FALSE
  ) +
  labs(
    x = "Delta percentage points (Enriched - Neat)",
    y = NULL,
    fill = NULL,
    title = "Protein category shift",
    subtitle = "Positive values indicate higher fractional category percentage in enriched proteins."
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5, color = "grey35"),
    axis.text.y = element_text(size = 8.2, color = "black"),
    axis.text.x = element_text(color = "black"),
    legend.position = "top",
    plot.margin = margin(8, 28, 8, 8)
  )

pdf(file.path(out_dir, "protein_category_delta_plot.pdf"), width = 7.8, height = max(4.8, 0.22 * nrow(delta_plot_data) + 2.1))
print(delta_plot)
dev.off()

all_categories_for_plot <- delta_summary %>%
  arrange(desc(enriched_percentage), desc(neat_percentage), Protein_Class) %>%
  pull(Protein_Class)

all_plot <- category_summary %>%
  mutate(
    Set = factor(Set, levels = set_levels),
    Protein_Class = factor(Protein_Class, levels = rev(all_categories_for_plot))
  )

all_colors <- setNames(rep(base_palette, length.out = length(all_categories_for_plot)), all_categories_for_plot)
all_colors["Unclassified"] <- "#BDBDBD"

all_stacked_plot <- ggplot(
  all_plot,
  aes(x = Set, y = percentage / 100, fill = Protein_Class)
) +
  geom_col(width = 0.58, color = "white", linewidth = 0.12) +
  scale_y_continuous(labels = percent_label, expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = all_colors, drop = FALSE) +
  labs(
    x = NULL,
    y = "Percentage of proteins",
    fill = "Protein category",
    title = "Protein category composition - all categories"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 7.2),
    legend.key.height = unit(0.28, "cm"),
    legend.key.width = unit(0.28, "cm"),
    legend.position = "right",
    plot.margin = margin(8, 8, 8, 8)
  ) +
  guides(fill = guide_legend(reverse = TRUE, ncol = max(1, ceiling(length(all_categories_for_plot) / 28))))

pdf(file.path(out_dir, "protein_category_stacked_bar_all_categories.pdf"), width = 8.2, height = max(5.2, 0.14 * length(all_categories_for_plot) + 4.0))
print(all_stacked_plot)
dev.off()

top_readme <- delta_summary %>%
  slice_head(n = min(8, nrow(delta_summary))) %>%
  transmute(
    line = paste0(
      "- `", Protein_Class, "`: Neat ",
      sprintf("%.2f", neat_percentage), "%, Enriched ",
      sprintf("%.2f", enriched_percentage), "%, delta ",
      pp_label(delta_percentage_points), " pp"
    )
  ) %>%
  pull(line)

readme_lines <- c(
  "# Protein category shift",
  "",
  "This folder compares the protein category composition of neat plasma proteins and enriched low-abundance proteins.",
  "",
  "## Inputs",
  "",
  "- `Plasma_Total_protein_tidy.csv`: neat protein table. Unique `Gene` labels define the Neat set.",
  "- `Low_abundance_15ul_tidy.csv`: enriched low-abundance protein table. Unique non-missing `Gene` labels define the Enriched set.",
  "- `Protein category.xlsx`: category annotation table. `Sheet1` columns `Gene` and `Protein_Class` are used.",
  "",
  "## Counting rule",
  "",
  "Gene labels are trimmed and matched case-insensitively. Each protein contributes total weight 1. If a protein has multiple `Protein_Class` annotations, its weight is split equally across those categories. Proteins without a matching category annotation are assigned to `Unclassified`. This makes each stacked bar sum to 100%.",
  "",
  "## Outputs",
  "",
  "- `protein_category_stacked_bar.pdf`: readable stacked bar plot using the largest changing categories plus `Other annotated categories` and `Unclassified`.",
  "- `protein_category_delta_plot.pdf`: delta plot of the largest category shifts, shown as Enriched minus Neat percentage points.",
  "- `protein_category_stacked_bar_all_categories.pdf`: stacked bar plot using all protein categories.",
  "- `protein_category_fractional_composition.csv`: all category percentages and fractional counts for Neat and Enriched.",
  "- `protein_category_delta.csv`: Neat percentage, Enriched percentage, and delta percentage points for every category.",
  "- `protein_category_plot_composition.csv`: collapsed category table used by `protein_category_stacked_bar.pdf`.",
  "- `protein_category_gene_assignments.csv`: gene-category assignments and fractional weights.",
  "- `protein_category_gene_diagnostics.csv`: per-gene category list and number of assigned categories.",
  "- `set_category_summary.csv`: set-level annotation coverage and multi-category diagnostics.",
  "- `summary_metrics.csv`: compact input-set relationship metrics.",
  "- `report.txt`: console log from the pipeline run.",
  "",
  "## Main observed shifts",
  "",
  top_readme,
  "",
  "Interpretation: positive delta values mean that a category has a higher fractional percentage in the Enriched set than in the Neat set. These are descriptive composition shifts rather than a formal independence test, because the two protein sets can overlap.",
  "",
  "Generated by `02_code/Plasma_Protein/protein_category_shift.R`."
)

writeLines(readme_lines, file.path(out_dir, "README.md"))

cat("\nOutput files written:\n")
print(list.files(out_dir, full.names = FALSE))
cat("\nPipeline completed successfully.\n")
