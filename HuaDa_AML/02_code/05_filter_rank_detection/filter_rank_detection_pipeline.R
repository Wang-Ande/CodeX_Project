#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  idx <- grep(file_arg, args, fixed = TRUE)
  if (length(idx) == 0) {
    return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  }
  normalizePath(sub(file_arg, "", args[idx[1]], fixed = TRUE),
                winslash = "/", mustWork = TRUE)
}

script_path <- get_script_path()
project_dir <- normalizePath(file.path(dirname(script_path), "..", ".."),
                              winslash = "/", mustWork = TRUE)
input_file <- file.path(project_dir, "01_data", "Total_pg_Report.csv")
out_dir <- file.path(project_dir, "03_result", "05_filter_rank_detection")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

report_file <- file.path(out_dir, "report.txt")
report_con <- file(report_file, open = "wt", encoding = "UTF-8")
sink(report_con, split = TRUE)
sink(report_con, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(report_con)
}, add = TRUE)

cat("Pipeline: filter samples, protein rank, target detection\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Project directory:", project_dir, "\n")
cat("Input file:", input_file, "\n")
cat("Output directory:", out_dir, "\n\n")

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

na_values <- c("", "NA", "NaN", "NAN", "nan", "NULL", "null")
raw_df <- read.csv(input_file,
                   check.names = FALSE,
                   stringsAsFactors = FALSE,
                   na.strings = na_values)

if (ncol(raw_df) < 2) {
  stop("Input file must contain one protein column and at least one sample column.")
}

protein <- as.character(raw_df[[1]])
keep_protein <- !is.na(protein) & nzchar(protein)
raw_df <- raw_df[keep_protein, , drop = FALSE]
protein <- protein[keep_protein]
protein_key <- make.unique(protein, sep = "__dup")

expr_df <- raw_df[, -1, drop = FALSE]
expr_df[] <- lapply(expr_df, function(x) suppressWarnings(as.numeric(x)))
expr_mat <- as.matrix(expr_df)
rownames(expr_mat) <- protein_key

non_finite_count <- sum(!is.finite(expr_mat), na.rm = TRUE)
non_positive_count <- sum(is.finite(expr_mat) & !is.na(expr_mat) & expr_mat <= 0)
expr_mat[!is.finite(expr_mat)] <- NA_real_
expr_mat[expr_mat <= 0] <- NA_real_

cat("Input proteins:", nrow(expr_mat), "\n")
cat("Input samples:", ncol(expr_mat), "\n")
cat("Non-finite values converted to missing:", non_finite_count, "\n")
cat("Non-positive values converted to missing:", non_positive_count, "\n\n")

is_missing <- is.na(expr_mat)
sample_missing <- data.frame(
  Sample = colnames(expr_mat),
  total_proteins = nrow(expr_mat),
  detected_proteins = colSums(!is_missing),
  missing_proteins = colSums(is_missing),
  missing_rate = colMeans(is_missing),
  stringsAsFactors = FALSE
)
sample_missing$status <- ifelse(sample_missing$missing_rate > 0.5,
                                "removed_missing_rate_gt_50pct",
                                "kept")

sample_missing_file <- file.path(out_dir, "sample_missing_rate.csv")
write.csv(sample_missing, sample_missing_file, row.names = FALSE, na = "")

removed_samples <- sample_missing$Sample[sample_missing$status != "kept"]
kept_samples <- sample_missing$Sample[sample_missing$status == "kept"]

cat("Samples removed for missing rate > 50%:", length(removed_samples), "\n")
if (length(removed_samples) > 0) {
  cat("Removed samples:", paste(removed_samples, collapse = ", "), "\n")
}
cat("Samples kept:", length(kept_samples), "\n\n")

filtered_mat <- expr_mat[, kept_samples, drop = FALSE]

filtered_file <- file.path(out_dir, "Total_pg_Report_sample_missing_le50pct.csv")
filtered_df <- data.frame(Protein = protein, filtered_mat, check.names = FALSE)
write.csv(filtered_df, filtered_file, row.names = FALSE, na = "")

sample_pdf <- file.path(out_dir, "sample_missing_rate.pdf")
pdf(sample_pdf, width = 8, height = 10)
sample_plot_df <- sample_missing[order(sample_missing$missing_rate, decreasing = TRUE), ]
bar_cols <- ifelse(sample_plot_df$status == "kept", "#2f7f5f", "#b13b3b")
par(mar = c(5, 8, 4, 2))
barplot(sample_plot_df$missing_rate * 100,
        names.arg = sample_plot_df$Sample,
        horiz = TRUE,
        las = 1,
        col = bar_cols,
        border = NA,
        xlab = "Missing proteins (%)",
        main = "Sample missing rate")
abline(v = 50, col = "#b13b3b", lty = 2, lwd = 1.5)
legend("bottomright",
       legend = c("Kept", "Removed (>50%)"),
       fill = c("#2f7f5f", "#b13b3b"),
       border = NA,
       bty = "n")
dev.off()

protein_detected <- rowSums(!is.na(filtered_mat))
protein_missing <- rowSums(is.na(filtered_mat))
protein_mean <- rowMeans(filtered_mat, na.rm = TRUE)
protein_mean[protein_detected == 0] <- NA_real_
protein_median <- apply(filtered_mat, 1, median, na.rm = TRUE)
protein_median[protein_detected == 0] <- NA_real_
protein_rank <- rank(-protein_mean, ties.method = "min", na.last = "keep")

rank_df <- data.frame(
  Protein = protein,
  mean_abundance = protein_mean,
  median_abundance = protein_median,
  rank_by_mean_desc = protein_rank,
  detected_samples = protein_detected,
  missing_samples = protein_missing,
  detection_rate = protein_detected / ncol(filtered_mat),
  missing_rate = protein_missing / ncol(filtered_mat),
  stringsAsFactors = FALSE
)
rank_df <- rank_df[order(rank_df$rank_by_mean_desc, na.last = TRUE), ]

rank_file <- file.path(out_dir, "protein_mean_rank.csv")
write.csv(rank_df, rank_file, row.names = FALSE, na = "")

target_genes <- c("CD96", "TNFRSF10A")
target_idx <- which(toupper(protein) %in% toupper(target_genes))

rank_pdf <- file.path(out_dir, "protein_mean_rank.pdf")
pdf(rank_pdf, width = 8, height = 5)
plot(rank_df$rank_by_mean_desc,
     log10(rank_df$mean_abundance),
     pch = 16,
     cex = 0.35,
     col = "#315c7c",
     xlab = "Protein rank by mean abundance",
     ylab = "log10(mean abundance)",
     main = "Protein mean abundance rank")
grid(col = "#dddddd")
if (length(target_idx) > 0) {
  target_rows <- rank_df[toupper(rank_df$Protein) %in% toupper(target_genes), , drop = FALSE]
  target_rows <- target_rows[!is.na(target_rows$rank_by_mean_desc), , drop = FALSE]
  if (nrow(target_rows) > 0) {
    points(target_rows$rank_by_mean_desc,
           log10(target_rows$mean_abundance),
           pch = 19,
           col = "#c43b3b",
           cex = 1.1)
    text(target_rows$rank_by_mean_desc,
         log10(target_rows$mean_abundance),
         labels = target_rows$Protein,
         pos = 4,
         cex = 0.75,
         col = "#8f2525")
  }
}
dev.off()

target_summary_list <- list()
target_long_list <- list()

for (target in target_genes) {
  idx <- which(toupper(protein) == toupper(target))
  if (length(idx) == 0) {
    target_summary_list[[target]] <- data.frame(
      Target = target,
      Protein = NA_character_,
      matched = FALSE,
      samples_after_filter = ncol(filtered_mat),
      detected_samples = NA_integer_,
      missing_samples = NA_integer_,
      detection_rate = NA_real_,
      missing_rate = NA_real_,
      mean_detected_abundance = NA_real_,
      median_detected_abundance = NA_real_,
      rank_by_mean_desc = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }

  for (i in idx) {
    values <- as.numeric(filtered_mat[i, ])
    detected <- !is.na(values)
    rank_match <- rank_df$rank_by_mean_desc[match(protein[i], rank_df$Protein)]

    target_summary_list[[paste(target, i, sep = "_")]] <- data.frame(
      Target = target,
      Protein = protein[i],
      matched = TRUE,
      samples_after_filter = ncol(filtered_mat),
      detected_samples = sum(detected),
      missing_samples = sum(!detected),
      detection_rate = mean(detected),
      missing_rate = mean(!detected),
      mean_detected_abundance = ifelse(any(detected), mean(values[detected]), NA_real_),
      median_detected_abundance = ifelse(any(detected), median(values[detected]), NA_real_),
      rank_by_mean_desc = rank_match,
      stringsAsFactors = FALSE
    )

    target_long_list[[paste(target, i, sep = "_")]] <- data.frame(
      Target = target,
      Protein = protein[i],
      Sample = colnames(filtered_mat),
      Abundance = values,
      Detected = detected,
      stringsAsFactors = FALSE
    )
  }
}

target_summary <- do.call(rbind, target_summary_list)
rownames(target_summary) <- NULL
target_summary_file <- file.path(out_dir, "CD96_TNFRSF10A_detection_summary.csv")
write.csv(target_summary, target_summary_file, row.names = FALSE, na = "")

if (length(target_long_list) > 0) {
  target_long <- do.call(rbind, target_long_list)
  rownames(target_long) <- NULL
} else {
  target_long <- data.frame(
    Target = character(),
    Protein = character(),
    Sample = character(),
    Abundance = numeric(),
    Detected = logical()
  )
}
target_by_sample_file <- file.path(out_dir, "CD96_TNFRSF10A_detection_by_sample.csv")
write.csv(target_long, target_by_sample_file, row.names = FALSE, na = "")

target_pdf <- file.path(out_dir, "CD96_TNFRSF10A_detection.pdf")
pdf(target_pdf, width = 7, height = 5)
matched_summary <- target_summary[target_summary$matched, , drop = FALSE]
if (nrow(matched_summary) > 0) {
  rate <- matched_summary$detection_rate * 100
  names(rate) <- matched_summary$Protein
  par(mar = c(5, 5, 4, 2))
  mids <- barplot(rate,
                  ylim = c(0, 100),
                  col = c("#4f7cac", "#d08a36")[seq_along(rate)],
                  border = NA,
                  ylab = "Detected samples (%)",
                  main = "CD96 and TNFRSF10A detection rate")
  text(mids,
       pmin(rate + 5, 96),
       labels = paste0(matched_summary$detected_samples, "/",
                       matched_summary$samples_after_filter),
       cex = 0.85)
  abline(h = seq(0, 100, 25), col = "#eeeeee", lty = 1)

  detected_long <- target_long[target_long$Detected, , drop = FALSE]
  if (nrow(detected_long) > 0) {
    boxplot(log10(Abundance) ~ Protein,
            data = detected_long,
            col = c("#b8d6f0", "#f0c894"),
            border = "#555555",
            ylab = "log10(abundance)",
            main = "Detected abundance distribution")
    stripchart(log10(Abundance) ~ Protein,
               data = detected_long,
               method = "jitter",
               vertical = TRUE,
               pch = 16,
               cex = 0.5,
               col = "#333333",
               add = TRUE)
  }
} else {
  plot.new()
  text(0.5, 0.5, "CD96 and TNFRSF10A were not found in the matrix.")
}
dev.off()

top10 <- head(rank_df[!is.na(rank_df$rank_by_mean_desc), ], 10)
summary_lines <- c(
  "Brief result summary",
  paste0("Input matrix: ", nrow(expr_mat), " proteins x ", ncol(expr_mat), " samples."),
  paste0("Sample filter: removed ", length(removed_samples),
         " samples with missing rate > 50%; kept ", length(kept_samples), " samples."),
  paste0("Filtered matrix: ", nrow(filtered_mat), " proteins x ", ncol(filtered_mat), " samples."),
  paste0("Top ranked proteins by mean abundance: ",
         paste(top10$Protein, collapse = ", ")),
  "Target detection:"
)

for (i in seq_len(nrow(target_summary))) {
  row <- target_summary[i, ]
  if (!isTRUE(row$matched)) {
    summary_lines <- c(summary_lines,
                       paste0("  - ", row$Target, ": not found in the protein matrix."))
  } else {
    summary_lines <- c(
      summary_lines,
      paste0("  - ", row$Protein, ": detected in ", row$detected_samples, "/",
             row$samples_after_filter, " samples (",
             sprintf("%.1f", row$detection_rate * 100), "%), rank ",
             row$rank_by_mean_desc, ".")
    )
  }
}

summary_file <- file.path(out_dir, "analysis_summary.txt")
writeLines(summary_lines, con = summary_file, useBytes = TRUE)

cat("Output files:\n")
cat("- Sample missing table:", sample_missing_file, "\n")
cat("- Sample missing plot:", sample_pdf, "\n")
cat("- Filtered matrix:", filtered_file, "\n")
cat("- Protein rank table:", rank_file, "\n")
cat("- Protein rank plot:", rank_pdf, "\n")
cat("- Target detection summary:", target_summary_file, "\n")
cat("- Target detection by sample:", target_by_sample_file, "\n")
cat("- Target detection plot:", target_pdf, "\n")
cat("- Brief summary:", summary_file, "\n\n")

cat(paste(summary_lines, collapse = "\n"), "\n")
cat("\nFinished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
