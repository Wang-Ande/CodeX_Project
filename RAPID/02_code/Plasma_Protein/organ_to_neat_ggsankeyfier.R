library(dplyr)
library(ggplot2)
library(ggsankeyfier)
library(readxl)

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1])
if (is.na(script_file) || !nzchar(script_file)) {
  script_file <- file.path(getwd(), "02_code", "Plasma_Protein", "organ_to_neat_ggsankeyfier.R")
}
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
data_dir <- file.path(repo_dir, "01_data", "Plasma_Protein")
out_dir <- file.path(repo_dir, "03_result", "Plasma_Protein", "Organ_neat_ggsankeyfier")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Organ-to-neat ggsankeyfier pipeline\n")
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

organ_neat_counts <- gene_assignments %>%
  count(Organ, Organ_order, Neat_quantile, name = "gene_organ_assignments") %>%
  mutate(
    Organ = factor(as.character(Organ), levels = organs),
    Neat_quantile = factor(as.character(Neat_quantile), levels = c("Q4", "Q3", "Q2", "Q1", "Not detected")),
    Flow_category = as.character(Neat_quantile)
  ) %>%
  arrange(Organ_order, Neat_quantile)

summary_metrics <- data.frame(
  metric = c(
    "requested_organs",
    "matched_requested_organs",
    "matched_gene_organ_assignments",
    "matched_unique_genes",
    "missing_requested_organs"
  ),
  value = c(
    length(organs),
    length(setdiff(organs, missing_organs)),
    nrow(gene_assignments),
    n_distinct(gene_assignments$Gene),
    paste(missing_organs, collapse = "; ")
  )
)

write.csv(organ_name_matching, file.path(out_dir, "organ_name_matching.csv"), row.names = FALSE)
write.csv(gene_assignments, file.path(out_dir, "selected_organ_neat_gene_assignments.csv"), row.names = FALSE)
write.csv(organ_neat_counts, file.path(out_dir, "organ_to_neat_counts.csv"), row.names = FALSE)
write.csv(summary_metrics, file.path(out_dir, "summary_metrics.csv"), row.names = FALSE)

cat("\nRequested organ matching:\n")
cat("Requested organs:", length(organs), "\n")
cat("Matched requested organs:", length(setdiff(organs, missing_organs)), "\n")
cat("Missing requested organs:", ifelse(length(missing_organs) > 0, paste(missing_organs, collapse = " | "), "None"), "\n")
cat("\nPlotting summary:\n")
cat("Matched gene-organ assignments:", nrow(gene_assignments), "\n")
cat("Unique matched genes:", n_distinct(gene_assignments$Gene), "\n")

plot_counts <- organ_neat_counts %>%
  transmute(
    Organ = Organ,
    `Neat quantile` = Neat_quantile,
    Count = gene_organ_assignments,
    Flow_category = Flow_category
  )

sankey_long <- pivot_stages_longer(
  data = plot_counts,
  stages_from = c("Organ", "Neat quantile"),
  values_from = "Count",
  additional_aes_from = "Flow_category"
) %>%
  mutate(
    stage = factor(as.character(stage), levels = c("Organ", "Neat quantile")),
    node = factor(as.character(node), levels = c(rev(organs), "Not detected", "Q1", "Q2", "Q3", "Q4")),
    Node_fill = as.character(node),
    Label = as.character(node)
  )

neat_colors <- c(
  "Q4" = "#1F4E6D",
  "Q3" = "#4F86A5",
  "Q2" = "#90B4C8",
  "Q1" = "#B0C8D8",
  "Not detected" = "#BDBDBD"
)

present_organs <- levels(droplevels(sankey_long$node[sankey_long$stage == "Organ"]))
present_organs <- present_organs[present_organs %in% organs]
set.seed(20260625)
organ_colors <- grDevices::hcl.colors(length(present_organs), palette = "Dark 3")
organ_colors <- sample(organ_colors, length(organ_colors))
names(organ_colors) <- present_organs
all_fill_colors <- c(organ_colors, neat_colors)

pos <- position_sankey(
  width = 0.10,
  align = "justify",
  order = "as_is",
  h_space = 0.34,
  v_space = 11
)
pos_organ_label <- position_sankey(
  width = 0.10,
  align = "justify",
  order = "as_is",
  h_space = 0.34,
  v_space = 11,
  nudge_x = -0.14
)
pos_neat_label <- position_sankey(
  width = 0.10,
  align = "justify",
  order = "as_is",
  h_space = 0.34,
  v_space = 11,
  nudge_x = 0.14
)

organ_plot <- ggplot(
  sankey_long,
  aes(
    x = stage,
    y = Count,
    group = node,
    connector = connector,
    edge_id = edge_id
  )
) +
  geom_sankeyedge(
    aes(fill = Flow_category),
    position = pos,
    alpha = 0.42,
    slope = 0.72,
    color = NA,
    show.legend = FALSE
  ) +
  geom_sankeynode(
    aes(fill = Node_fill),
    position = pos,
    color = "white",
    linewidth = 0.18,
    show.legend = FALSE
  ) +
  geom_text(
    data = subset(sankey_long, stage == "Organ"),
    aes(label = Label),
    stat = "sankeynode",
    position = pos_organ_label,
    hjust = 1,
    size = 2.55,
    color = "#1F2933",
    show.legend = FALSE
  ) +
  geom_text(
    data = subset(sankey_long, stage == "Neat quantile"),
    aes(label = Label),
    stat = "sankeynode",
    position = pos_neat_label,
    hjust = 0,
    size = 3.1,
    color = "#111827",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = all_fill_colors, drop = FALSE) +
  scale_x_discrete(
    labels = c("Organ category", "Neat quantile"),
    expand = expansion(add = c(0.95, 0.75))
  ) +
  labs(
    title = "Organ category to neat protein abundance category",
    subtitle = "Generated with ggsankeyfier; organ labels are placed left of the modules",
    x = NULL,
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "plain", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 10.5, color = "#4B5563", hjust = 0, margin = margin(b = 8)),
    axis.text.x = element_text(size = 10.5, color = "#374151", margin = margin(t = 6)),
    plot.margin = margin(14, 96, 18, 118)
  )

pdf(file.path(out_dir, "organ_to_neat_ggsankeyfier.pdf"), width = 12.5, height = 16)
print(organ_plot)
dev.off()

readme_lines <- c(
  "# Organ to neat protein Sankey analysis with ggsankeyfier",
  "",
  "This folder contains the standalone organ-side Sankey plot requested for manual assembly with the existing neat-to-low Sankey plot.",
  "",
  "## Inputs",
  "",
  "- `Low_abundance_15ul_tidy.csv`: low-abundance protein table used to identify proteins not detected in the neat/full protein table.",
  "- `Plasma_Total_protein_tidy.csv`: neat/full plasma protein table. Numeric columns are averaged per gene, then converted into neat-side `Q1`-`Q4` categories.",
  "- `Specific tissue cluster2.xlsx`: tissue-specific gene table. It is read as a long table with `Gene`, `Tissue`, and `Function`; only the requested organ list is used.",
  "",
  "## Outputs",
  "",
  "- `organ_to_neat_ggsankeyfier.pdf`: standalone `Organ category -> Neat quantile` Sankey plot. The PDF is generated from `Specific tissue cluster2.xlsx` restricted to the requested organs, joined to neat categories derived from `Plasma_Total_protein_tidy.csv` and `Low_abundance_15ul_tidy.csv`. Organ modules use a coordinated multi-color palette; flows and neat modules keep the original neat-category colors.",
  "- `organ_to_neat_counts.csv`: count table used for the PDF.",
  "- `selected_organ_neat_gene_assignments.csv`: matched gene-organ assignments with neat quantile, low quantile, and type.",
  "- `organ_name_matching.csv`: requested organ names matched to tissue names in `Specific tissue cluster2.xlsx` after normalized matching.",
  "- `summary_metrics.csv`: compact run metrics.",
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
  "- Counts represent gene-organ assignments, not mutually exclusive proteins.",
  "- `Not detected` on the neat side means the gene appears in the low-abundance table but not in the neat/full total protein table.",
  "- `ggsankeyfier` is used for this plot so vertical spacing between modules can be controlled with `v_space`.",
  "",
  "Generated by `02_code/Plasma_Protein/organ_to_neat_ggsankeyfier.R`."
)

writeLines(readme_lines, file.path(out_dir, "README.md"))

cat("\nOutput files written:\n")
print(list.files(out_dir, full.names = FALSE))
cat("\nPipeline completed successfully.\n")
