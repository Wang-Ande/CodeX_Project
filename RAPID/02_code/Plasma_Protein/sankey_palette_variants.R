library(dplyr)
library(ggplot2)
library(ggalluvial)
library(ggrepel)
library(readxl)

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1])
if (is.na(script_file) || !nzchar(script_file)) {
  script_file <- file.path(getwd(), "02_code", "Plasma_Protein", "sankey_palette_variants.R")
}
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
data_dir <- file.path(repo_dir, "01_data", "Plasma_Protein")
out_dir <- file.path(repo_dir, "03_result", "Plasma_Protein", "Sankey_palette_variants")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Sankey palette variants pipeline\n")
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

organ_system_map <- data.frame(
  Organ = organs,
  Organ_system = c(
    rep("Nervous / sensory", 6),
    rep("Endocrine", 3),
    rep("Immune / hematopoietic", 4),
    "Respiratory",
    rep("Digestive / metabolic", 6),
    rep("Urogenital / reproductive", 6),
    rep("Muscle / connective", 3),
    rep("Epithelial / skin", 4),
    "Muscle / connective"
  ),
  stringsAsFactors = FALSE
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
    Full_quantile = factor(as.character(Full_quantile), levels = c("Q4", "Q3", "Q2", "Q1", "Not detected")),
    Low_quantile = factor(as.character(Low_quantile), levels = c("Q4", "Q3", "Q2", "Q1", "Not detected")),
    Neat_quantile = factor(as.character(Full_quantile), levels = c("Q4", "Q3", "Q2", "Q1", "Not detected")),
    Flow_category = ifelse(Type == "Overlap", as.character(Full_quantile), "Not detected"),
    Flow_category = factor(Flow_category, levels = c("Q4", "Q3", "Q2", "Q1", "Not detected"))
  ) %>%
  select(Gene, Type, Full_quantile, Low_quantile, Neat_quantile, Flow_category)

abundance_counts <- protein_quantile_map %>%
  count(Full_quantile, Low_quantile, Type, Flow_category, name = "protein_count") %>%
  arrange(Full_quantile, Low_quantile, Type)

abundance_neat_targets <- protein_quantile_map %>%
  count(Neat_quantile, name = "abundance_neat_axis_total") %>%
  mutate(Neat_quantile = factor(as.character(Neat_quantile), levels = c("Q4", "Q3", "Q2", "Q1", "Not detected")))

axis_totals <- bind_rows(
  protein_quantile_map %>%
    count(Full_quantile, name = "axis_total") %>%
    transmute(axis = "Global proteins", category = as.character(Full_quantile), axis_total),
  protein_quantile_map %>%
    count(Low_quantile, name = "axis_total") %>%
    transmute(axis = "Enriched proteins", category = as.character(Low_quantile), axis_total)
) %>%
  mutate(category = factor(category, levels = c("Q4", "Q3", "Q2", "Q1", "Not detected"))) %>%
  arrange(axis, category)

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
  inner_join(protein_quantile_map %>% select(Gene, Type, Neat_quantile, Low_quantile), by = "Gene") %>%
  mutate(
    Organ = factor(Organ, levels = organs),
    Neat_quantile = factor(as.character(Neat_quantile), levels = c("Q4", "Q3", "Q2", "Q1", "Not detected")),
    Low_quantile = factor(as.character(Low_quantile), levels = c("Q4", "Q3", "Q2", "Q1", "Not detected"))
  ) %>%
  arrange(Organ_order, Neat_quantile, Low_quantile, Gene)

min_flow_basis <- 2
organ_neat_counts <- gene_assignments %>%
  count(Organ, Organ_order, Neat_quantile, name = "gene_organ_assignments") %>%
  mutate(
    Organ = factor(as.character(Organ), levels = organs),
    Neat_quantile = factor(as.character(Neat_quantile), levels = c("Q4", "Q3", "Q2", "Q1", "Not detected")),
    Flow_category = as.character(Neat_quantile)
  ) %>%
  arrange(Organ_order, Neat_quantile) %>%
  left_join(abundance_neat_targets, by = "Neat_quantile") %>%
  group_by(Neat_quantile) %>%
  mutate(
    display_basis = pmax(gene_organ_assignments, min_flow_basis),
    organ_neat_raw_total = sum(gene_organ_assignments),
    display_basis_total = sum(display_basis),
    alignment_scale = abundance_neat_axis_total / display_basis_total,
    aligned_display_weight = display_basis * alignment_scale
  ) %>%
  ungroup()

alignment_check <- organ_neat_counts %>%
  group_by(Neat_quantile, abundance_neat_axis_total) %>%
  summarise(aligned_neat_total = sum(aligned_display_weight), .groups = "drop") %>%
  mutate(delta = aligned_neat_total - abundance_neat_axis_total)

write.csv(organ_name_matching, file.path(out_dir, "organ_name_matching.csv"), row.names = FALSE)
write.csv(gene_assignments, file.path(out_dir, "selected_organ_neat_gene_assignments.csv"), row.names = FALSE)
write.csv(organ_neat_counts, file.path(out_dir, "organ_to_neat_aligned_counts.csv"), row.names = FALSE)
write.csv(abundance_neat_targets, file.path(out_dir, "abundance_neat_axis_targets.csv"), row.names = FALSE)
write.csv(alignment_check, file.path(out_dir, "alignment_check.csv"), row.names = FALSE)
write.csv(protein_quantile_map, file.path(out_dir, "abundance_gene_quantile_map.csv"), row.names = FALSE)
write.csv(abundance_counts, file.path(out_dir, "abundance_sankey_counts.csv"), row.names = FALSE)
write.csv(axis_totals, file.path(out_dir, "abundance_axis_totals.csv"), row.names = FALSE)

cat("\nData summary:\n")
cat("Full/neat genes:", length(full_genes), "\n")
cat("Low abundance genes:", length(low_genes), "\n")
cat("Overlap genes:", length(overlap_genes), "\n")
cat("Matched requested organs:", length(setdiff(organs, missing_organs)), "\n")
cat("Missing requested organs:", ifelse(length(missing_organs) > 0, paste(missing_organs, collapse = " | "), "None"), "\n")
cat("Matched gene-organ assignments:", nrow(gene_assignments), "\n")
cat("Aligned display total:", sum(organ_neat_counts$aligned_display_weight), "\n")
cat("\nAlignment check:\n")
print(alignment_check)

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

palette_variants <- list(
  list(
    id = "01_sapphire_coral",
    name = "Sapphire coral",
    note = "High-contrast cool abundance colors with warm organ accents.",
    abundance = c("Q4" = "#1B1B4D", "Q3" = "#2364AA", "Q2" = "#3DA5D9", "Q1" = "#C2EABD", "Not detected" = "#D4D4D4"),
    organs = c(
      "Nervous / sensory" = "#7B61A7",
      "Endocrine" = "#F4A261",
      "Immune / hematopoietic" = "#2A9D8F",
      "Respiratory" = "#5DADE2",
      "Digestive / metabolic" = "#43AA8B",
      "Urogenital / reproductive" = "#E76F9A",
      "Muscle / connective" = "#B87935",
      "Epithelial / skin" = "#F08A8A"
    )
  ),
  list(
    id = "02_royal_sunset",
    name = "Royal sunset",
    note = "Dramatic purple-to-coral abundance palette, warmer and less blue.",
    abundance = c("Q4" = "#2D1E59", "Q3" = "#7A1E7A", "Q2" = "#D1495B", "Q1" = "#F4A261", "Not detected" = "#D8D2C4"),
    organs = c(
      "Nervous / sensory" = "#5E4FA2",
      "Endocrine" = "#FDB863",
      "Immune / hematopoietic" = "#1B9E77",
      "Respiratory" = "#3288BD",
      "Digestive / metabolic" = "#66C2A5",
      "Urogenital / reproductive" = "#E7298A",
      "Muscle / connective" = "#A6761D",
      "Epithelial / skin" = "#FC8D62"
    )
  ),
  list(
    id = "03_viridis_gold",
    name = "Viridis gold",
    note = "Colorblind-friendly viridis-like abundance ramp with bright Q1.",
    abundance = c("Q4" = "#440154", "Q3" = "#31688E", "Q2" = "#35B779", "Q1" = "#FDE725", "Not detected" = "#D0D0D0"),
    organs = c(
      "Nervous / sensory" = "#6A3D9A",
      "Endocrine" = "#FFB000",
      "Immune / hematopoietic" = "#33A02C",
      "Respiratory" = "#1F78B4",
      "Digestive / metabolic" = "#00A087",
      "Urogenital / reproductive" = "#C51B7D",
      "Muscle / connective" = "#8C510A",
      "Epithelial / skin" = "#FB9A99"
    )
  ),
  list(
    id = "04_teal_ruby",
    name = "Teal ruby",
    note = "Crisp teal abundance scale with ruby and ochre organ accents.",
    abundance = c("Q4" = "#003049", "Q3" = "#0A9396", "Q2" = "#94D2BD", "Q1" = "#E9D8A6", "Not detected" = "#D6D6D6"),
    organs = c(
      "Nervous / sensory" = "#6D597A",
      "Endocrine" = "#F77F00",
      "Immune / hematopoietic" = "#2E7D32",
      "Respiratory" = "#277DA1",
      "Digestive / metabolic" = "#00A896",
      "Urogenital / reproductive" = "#B5179E",
      "Muscle / connective" = "#9C6644",
      "Epithelial / skin" = "#EF476F"
    )
  ),
  list(
    id = "05_ink_lotus",
    name = "Ink lotus",
    note = "Deep ink abundance colors with a soft lotus Q1 and vivid organ groups.",
    abundance = c("Q4" = "#0B132B", "Q3" = "#3A506B", "Q2" = "#5BC0BE", "Q1" = "#F2D7EE", "Not detected" = "#D7D7D7"),
    organs = c(
      "Nervous / sensory" = "#845EC2",
      "Endocrine" = "#FFC75F",
      "Immune / hematopoietic" = "#4D9F0C",
      "Respiratory" = "#0081CF",
      "Digestive / metabolic" = "#00C9A7",
      "Urogenital / reproductive" = "#D65DB1",
      "Muscle / connective" = "#B08968",
      "Epithelial / skin" = "#FF8066"
    )
  )
)

make_organ_plot <- function(variant) {
  organ_plot_data <- organ_neat_counts %>%
    mutate(alluvium = row_number()) %>%
    select(alluvium, Organ, Neat_quantile, Flow_category, aligned_display_weight)

  lodes <- to_lodes_form(
    organ_plot_data,
    axes = c("Organ", "Neat_quantile"),
    key = "axis",
    id = "alluvium"
  ) %>%
    mutate(
      axis_x = ifelse(axis == "Organ", 1, 2),
      stratum = factor(as.character(stratum), levels = c(organs, "Q4", "Q3", "Q2", "Q1", "Not detected"))
    )

  present_organs <- organs[organs %in% unique(as.character(organ_neat_counts$Organ))]
  organ_colors <- organ_system_map %>%
    filter(.data$Organ %in% present_organs) %>%
    mutate(color = variant$organs[.data$Organ_system]) %>%
    select("Organ", "color")
  organ_colors <- setNames(organ_colors$color, organ_colors$Organ)
  fill_colors <- c(organ_colors, variant$abundance)

  base_plot <- ggplot(
    lodes,
    aes(
      x = axis_x,
      stratum = stratum,
      alluvium = alluvium,
      y = aligned_display_weight
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
    scale_fill_manual(values = fill_colors, drop = FALSE) +
    scale_alpha_manual(values = flow_alpha, drop = FALSE, guide = "none") +
    scale_x_continuous(
      breaks = c(1, 2),
      labels = c("Tissue-derived proteins", "Global proteins"),
      limits = c(0.30, 2.24),
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 12, color = "black", margin = margin(t = 4)),
      plot.margin = margin(8, 70, 18, 118)
    )

  strata_coords <- ggplot_build(base_plot)$data[[2]]
  if (!all(c("x", "y", "stratum") %in% names(strata_coords))) {
    stop("Could not recover organ stratum coordinates from ggplot_build().")
  }

  left_labels <- strata_coords %>%
    filter(abs(x - 1) < 1e-6) %>%
    mutate(label = as.character(stratum), x_anchor = 0.93)

  right_labels <- strata_coords %>%
    filter(abs(x - 2) < 1e-6) %>%
    mutate(label = unname(label_map[as.character(stratum)]), x_label = 2.08)

  base_plot +
    geom_text_repel(
      data = left_labels,
      inherit.aes = FALSE,
      aes(x = x_anchor, y = y, label = label),
      nudge_x = -0.34,
      direction = "y",
      hjust = 1,
      size = 3.6,
      color = "black",
      segment.color = "grey45",
      segment.size = 0.28,
      segment.alpha = 0.9,
      min.segment.length = 0,
      box.padding = 0.06,
      point.padding = 0.03,
      force = 1.1,
      force_pull = 0.12,
      max.overlaps = Inf,
      seed = 20260625
    ) +
    geom_text(
      data = right_labels,
      inherit.aes = FALSE,
      aes(x = x_label, y = y, label = label),
      hjust = 0,
      size = 3.8,
      color = "black"
    )
}

make_abundance_plot <- function(variant) {
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
    scale_fill_manual(values = variant$abundance, drop = FALSE) +
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
    stop("Could not recover abundance stratum coordinates from ggplot_build().")
  }

  left_labels <- strata_coords %>%
    filter(abs(x - 1) < 1e-6) %>%
    mutate(label = unname(label_map[as.character(stratum)]), x_label = 0.92)

  right_labels <- strata_coords %>%
    filter(abs(x - 2) < 1e-6) %>%
    mutate(label = unname(label_map[as.character(stratum)]), x_label = 2.08)

  base_plot +
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
}

write_palette_key <- function(variant, variant_dir) {
  abundance_key <- data.frame(
    Palette = variant$id,
    Palette_name = variant$name,
    Color_role = names(variant$abundance),
    Color_hex = unname(variant$abundance),
    stringsAsFactors = FALSE
  )

  organ_key <- data.frame(
    Palette = variant$id,
    Palette_name = variant$name,
    Color_role = names(variant$organs),
    Color_hex = unname(variant$organs),
    stringsAsFactors = FALSE
  )

  write.csv(
    bind_rows(
      abundance_key %>% mutate(Color_group = "Abundance / flow"),
      organ_key %>% mutate(Color_group = "Organ system")
    ),
    file.path(variant_dir, "palette_key.csv"),
    row.names = FALSE
  )
}

palette_summary <- bind_rows(lapply(palette_variants, function(variant) {
  data.frame(
    palette_id = variant$id,
    palette_name = variant$name,
    note = variant$note,
    Q4 = variant$abundance[["Q4"]],
    Q3 = variant$abundance[["Q3"]],
    Q2 = variant$abundance[["Q2"]],
    Q1 = variant$abundance[["Q1"]],
    NA_color = variant$abundance[["Not detected"]],
    stringsAsFactors = FALSE
  )
}))
write.csv(palette_summary, file.path(out_dir, "palette_summary.csv"), row.names = FALSE)

for (variant in palette_variants) {
  variant_dir <- file.path(out_dir, variant$id)
  dir.create(variant_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\nRendering palette:", variant$id, "-", variant$name, "\n")

  organ_plot <- make_organ_plot(variant)
  abundance_plot <- make_abundance_plot(variant)

  pdf(file.path(variant_dir, paste0(variant$id, "_organ_to_global.pdf")), width = 8.2, height = 8.8)
  print(organ_plot)
  dev.off()

  pdf(file.path(variant_dir, paste0(variant$id, "_global_to_enriched.pdf")), width = 5.2, height = 8.8)
  print(abundance_plot)
  dev.off()

  write_palette_key(variant, variant_dir)
}

readme_lines <- c(
  "# Sankey palette variants",
  "",
  "This folder contains alternative color effects for the latest Sankey layout. The layout, data processing, axis alignment, labels, and flow alpha settings are unchanged from the final color-optimized version; only color assignments differ.",
  "",
  "## Inputs",
  "",
  "- `Low_abundance_15ul_tidy.csv`: low-abundance/enriched protein table.",
  "- `Plasma_Total_protein_tidy.csv`: neat/global protein table.",
  "- `Specific tissue cluster2.xlsx`: tissue-specific organ annotation table.",
  "",
  "## Shared data outputs",
  "",
  "- `organ_to_neat_aligned_counts.csv`: organ-to-global category counts and aligned display weights used by all organ panels.",
  "- `abundance_sankey_counts.csv`: global-to-enriched flow counts used by all abundance panels.",
  "- `abundance_axis_totals.csv`: global and enriched axis totals.",
  "- `alignment_check.csv`: confirms the organ panel Global protein axis matches the abundance panel Global protein axis.",
  "- `palette_summary.csv`: compact list of the five palette variants and their Q4/Q3/Q2/Q1/NA colors.",
  "- `report.txt`: console log from this pipeline run.",
  "",
  "## Palette folders",
  "",
  "- `01_sapphire_coral`: high-contrast cool abundance colors with warm organ accents.",
  "- `02_royal_sunset`: dramatic purple-to-coral abundance palette, warmer and less blue.",
  "- `03_viridis_gold`: colorblind-friendly viridis-like abundance ramp with bright Q1.",
  "- `04_teal_ruby`: crisp teal abundance scale with ruby and ochre organ accents.",
  "- `05_ink_lotus`: deep ink abundance colors with a soft lotus Q1 and vivid organ groups.",
  "",
  "Each palette folder contains:",
  "",
  "- `<palette>_organ_to_global.pdf`: tissue-derived proteins to Global proteins.",
  "- `<palette>_global_to_enriched.pdf`: Global proteins to Enriched proteins.",
  "- `palette_key.csv`: exact abundance and organ-system color hex codes.",
  "",
  "## Plotting notes",
  "",
  "- Q4/Q3/Q2/Q1 flows use alpha 0.74/0.72/0.68/0.64.",
  "- NA flows use alpha 0.46 so the large undetected component remains visible on white backgrounds without dominating the colored Q flows.",
  "- Organ modules are colored by organ system; flows and Global/Enriched modules are colored by the abundance category.",
  "- The plotted label `NA` corresponds to the original `Not detected` category.",
  "",
  "Generated by `02_code/Plasma_Protein/sankey_palette_variants.R`."
)
writeLines(readme_lines, file.path(out_dir, "README.md"))

cat("\nOutput files written under:\n")
print(out_dir)
cat("\nPipeline completed successfully.\n")
