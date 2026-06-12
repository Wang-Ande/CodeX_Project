#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  idx <- grep("--file=", args, fixed = TRUE)
  if (length(idx) == 0) {
    return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  }
  normalizePath(sub("--file=", "", args[idx[1]], fixed = TRUE),
                winslash = "/", mustWork = TRUE)
}

message_line <- function(...) {
  cat(..., "\n", sep = "")
}

script_path <- get_script_path()
project_dir <- normalizePath(file.path(dirname(script_path), "..", ".."),
                              winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_dir, "03_result", "05_filter_rank_detection")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

report_file <- file.path(out_dir, "report.txt")
report_con <- file(report_file, open = "at", encoding = "UTF-8")
sink(report_con, split = TRUE)
sink(report_con, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(report_con)
}, add = TRUE)

message_line("")
message_line("============================================================")
message_line("Pipeline: add CD33 and IL3RA to target detection and rank plot")
message_line("Started at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
message_line("Project directory: ", project_dir)
message_line("Output directory: ", out_dir)

input_file <- file.path(project_dir, "01_data", "Total_pg_Report.csv")
filtered_file <- file.path(out_dir, "Total_pg_Report_sample_missing_le50pct.csv")
rank_file <- file.path(out_dir, "protein_mean_rank.csv")
rank_pdf <- file.path(out_dir, "protein_mean_rank.pdf")

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

read_expression_matrix <- function(path) {
  na_values <- c("", "NA", "NaN", "NAN", "nan", "NULL", "null")
  df <- read.csv(path,
                 check.names = FALSE,
                 stringsAsFactors = FALSE,
                 na.strings = na_values)
  if (ncol(df) < 2) {
    stop("Expression file must contain one protein column and sample columns: ", path)
  }

  protein <- as.character(df[[1]])
  keep <- !is.na(protein) & nzchar(protein)
  df <- df[keep, , drop = FALSE]
  protein <- protein[keep]

  expr <- df[, -1, drop = FALSE]
  expr[] <- lapply(expr, function(x) suppressWarnings(as.numeric(x)))
  mat <- as.matrix(expr)
  mat[!is.finite(mat)] <- NA_real_
  mat[mat <= 0] <- NA_real_
  rownames(mat) <- make.unique(protein, sep = "__dup")

  list(protein = protein, matrix = mat)
}

original <- read_expression_matrix(input_file)
expr_mat <- original$matrix
protein <- original$protein

sample_missing <- data.frame(
  Sample = colnames(expr_mat),
  total_proteins = nrow(expr_mat),
  detected_proteins = colSums(!is.na(expr_mat)),
  missing_proteins = colSums(is.na(expr_mat)),
  missing_rate = colMeans(is.na(expr_mat)),
  stringsAsFactors = FALSE
)
sample_missing$status <- ifelse(sample_missing$missing_rate > 0.5,
                                "removed_missing_rate_gt_50pct",
                                "kept")
kept_samples <- sample_missing$Sample[sample_missing$status == "kept"]
removed_samples <- sample_missing$Sample[sample_missing$status != "kept"]

filtered_mat <- expr_mat[, kept_samples, drop = FALSE]
filtered_df <- data.frame(Protein = protein, filtered_mat, check.names = FALSE)
write.csv(filtered_df, filtered_file, row.names = FALSE, na = "")
write.csv(sample_missing,
          file.path(out_dir, "sample_missing_rate.csv"),
          row.names = FALSE,
          na = "")

message_line("Input matrix: ", nrow(expr_mat), " proteins x ", ncol(expr_mat), " samples.")
message_line("Samples removed for missing rate > 50%: ", length(removed_samples))
if (length(removed_samples) > 0) {
  message_line("Removed samples: ", paste(removed_samples, collapse = ", "))
}
message_line("Filtered matrix: ", nrow(filtered_mat), " proteins x ", ncol(filtered_mat), " samples.")

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
write.csv(rank_df, rank_file, row.names = FALSE, na = "")

target_genes <- c("CD96", "TNFRSF10A", "CD33", "IL3RA")
target_colors <- c(
  CD96 = "#c43b3b",
  TNFRSF10A = "#7448a6",
  CD33 = "#1b8a6b",
  IL3RA = "#d08a36"
)

target_rows <- rank_df[toupper(rank_df$Protein) %in% toupper(target_genes), , drop = FALSE]
target_rows <- target_rows[!is.na(target_rows$rank_by_mean_desc), , drop = FALSE]

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
if (nrow(target_rows) > 0) {
  for (i in seq_len(nrow(target_rows))) {
    target_name <- target_rows$Protein[i]
    target_color <- target_colors[[toupper(target_name)]]
    if (is.null(target_color)) {
      target_color <- "#c43b3b"
    }
    points(target_rows$rank_by_mean_desc[i],
           log10(target_rows$mean_abundance[i]),
           pch = 19,
           col = target_color,
           cex = 1.1)
    text(target_rows$rank_by_mean_desc[i],
         log10(target_rows$mean_abundance[i]),
         labels = target_name,
         pos = ifelse(i %% 2 == 0, 2, 4),
         cex = 0.75,
         col = target_color)
  }
  legend("topright",
         legend = target_rows$Protein,
         col = target_colors[toupper(target_rows$Protein)],
         pch = 19,
         bty = "n",
         cex = 0.75)
}
dev.off()

build_target_detection <- function(targets, protein, mat, rank_df) {
  summary_list <- list()
  sample_list <- list()

  for (target in targets) {
    idx <- which(toupper(protein) == toupper(target))
    if (length(idx) == 0) {
      summary_list[[target]] <- data.frame(
        Target = target,
        Protein = NA_character_,
        matched = FALSE,
        samples_after_filter = ncol(mat),
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
      values <- as.numeric(mat[i, ])
      detected <- !is.na(values)
      rank_match <- rank_df$rank_by_mean_desc[match(protein[i], rank_df$Protein)]
      summary_list[[paste(target, i, sep = "_")]] <- data.frame(
        Target = target,
        Protein = protein[i],
        matched = TRUE,
        samples_after_filter = ncol(mat),
        detected_samples = sum(detected),
        missing_samples = sum(!detected),
        detection_rate = mean(detected),
        missing_rate = mean(!detected),
        mean_detected_abundance = ifelse(any(detected), mean(values[detected]), NA_real_),
        median_detected_abundance = ifelse(any(detected), median(values[detected]), NA_real_),
        rank_by_mean_desc = rank_match,
        stringsAsFactors = FALSE
      )
      sample_list[[paste(target, i, sep = "_")]] <- data.frame(
        Target = target,
        Protein = protein[i],
        Sample = colnames(mat),
        Abundance = values,
        Detected = detected,
        stringsAsFactors = FALSE
      )
    }
  }

  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  if (length(sample_list) > 0) {
    sample_df <- do.call(rbind, sample_list)
    rownames(sample_df) <- NULL
  } else {
    sample_df <- data.frame(Target = character(),
                            Protein = character(),
                            Sample = character(),
                            Abundance = numeric(),
                            Detected = logical())
  }

  list(summary = summary_df, by_sample = sample_df)
}

target_detection <- build_target_detection(target_genes, protein, filtered_mat, rank_df)
target_summary <- target_detection$summary
target_by_sample <- target_detection$by_sample

all_target_summary_file <- file.path(out_dir, "target_protein_detection_summary.csv")
all_target_by_sample_file <- file.path(out_dir, "target_protein_detection_by_sample.csv")
new_target_summary_file <- file.path(out_dir, "CD33_IL3RA_detection_summary.csv")
new_target_by_sample_file <- file.path(out_dir, "CD33_IL3RA_detection_by_sample.csv")

write.csv(target_summary, all_target_summary_file, row.names = FALSE, na = "")
write.csv(target_by_sample, all_target_by_sample_file, row.names = FALSE, na = "")
write.csv(target_summary[target_summary$Target %in% c("CD33", "IL3RA"), , drop = FALSE],
          new_target_summary_file,
          row.names = FALSE,
          na = "")
write.csv(target_by_sample[target_by_sample$Target %in% c("CD33", "IL3RA"), , drop = FALSE],
          new_target_by_sample_file,
          row.names = FALSE,
          na = "")

target_pdf <- file.path(out_dir, "target_protein_detection.pdf")
new_target_pdf <- file.path(out_dir, "CD33_IL3RA_detection.pdf")

plot_detection_pdf <- function(summary_df, sample_df, pdf_file, title_prefix) {
  pdf(pdf_file, width = 7, height = 5)
  matched <- summary_df[summary_df$matched, , drop = FALSE]
  if (nrow(matched) == 0) {
    plot.new()
    text(0.5, 0.5, paste0(title_prefix, ": no target proteins were found."))
    dev.off()
    return(invisible(NULL))
  }

  rate <- matched$detection_rate * 100
  names(rate) <- matched$Protein
  colors <- target_colors[toupper(matched$Protein)]
  colors[is.na(colors)] <- "#4f7cac"

  par(mar = c(5, 5, 4, 2))
  mids <- barplot(rate,
                  ylim = c(0, 100),
                  col = colors,
                  border = NA,
                  ylab = "Detected samples (%)",
                  main = paste0(title_prefix, " detection rate"))
  text(mids,
       pmin(rate + 5, 96),
       labels = paste0(matched$detected_samples, "/",
                       matched$samples_after_filter),
       cex = 0.85)
  abline(h = seq(0, 100, 25), col = "#eeeeee", lty = 1)

  detected_sample <- sample_df[sample_df$Detected, , drop = FALSE]
  if (nrow(detected_sample) > 0) {
    boxplot(log10(Abundance) ~ Protein,
            data = detected_sample,
            col = "#d6e3f0",
            border = "#555555",
            ylab = "log10(abundance)",
            main = paste0(title_prefix, " detected abundance"))
    stripchart(log10(Abundance) ~ Protein,
               data = detected_sample,
               method = "jitter",
               vertical = TRUE,
               pch = 16,
               cex = 0.5,
               col = "#333333",
               add = TRUE)
  }
  dev.off()
}

plot_detection_pdf(target_summary, target_by_sample, target_pdf, "Target proteins")
plot_detection_pdf(target_summary[target_summary$Target %in% c("CD33", "IL3RA"), , drop = FALSE],
                   target_by_sample[target_by_sample$Target %in% c("CD33", "IL3RA"), , drop = FALSE],
                   new_target_pdf,
                   "CD33 and IL3RA")

summary_lines <- c(
  "Brief result summary",
  paste0("Input matrix: ", nrow(expr_mat), " proteins x ", ncol(expr_mat), " samples."),
  paste0("Sample filter: removed ", length(removed_samples),
         " samples with missing rate > 50%; kept ", length(kept_samples), " samples."),
  paste0("Rank plot was regenerated with labels for: ",
         paste(target_genes, collapse = ", "), "."),
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

summary_file <- file.path(out_dir, "target_protein_analysis_summary.txt")
writeLines(summary_lines, con = summary_file, useBytes = TRUE)

message_line("Output files:")
message_line("- Updated rank table: ", rank_file)
message_line("- Updated rank plot: ", rank_pdf)
message_line("- All target detection summary: ", all_target_summary_file)
message_line("- All target detection by sample: ", all_target_by_sample_file)
message_line("- CD33/IL3RA detection summary: ", new_target_summary_file)
message_line("- CD33/IL3RA detection by sample: ", new_target_by_sample_file)
message_line("- All target detection plot: ", target_pdf)
message_line("- CD33/IL3RA detection plot: ", new_target_pdf)
message_line("- Brief summary: ", summary_file)
message_line("")
message_line(paste(summary_lines, collapse = "\n"))
message_line("Finished at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
