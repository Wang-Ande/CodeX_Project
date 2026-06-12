library(dplyr)
library(openxlsx)
library(readr)

# gene expr matrix ----
# 行是基因名，列是样本名
expr_mat <- read.csv("./01_data/Normalized Expression.csv", row.names = 1)
rownames(expr_mat) <- expr_mat$display_label
expr_mat <- expr_mat[,-c(1,2,3)]

# clinical info ----
clinical_summary <- read.xlsx("./01_data/Clinical Summary.xlsx")
table(clinical_summary$dxAtInclusion)
table(clinical_summary$specimenType)
# 只保留 BM来源的，且为AML的样本
aml_bm <- clinical_summary %>%
  filter(
    dxAtInclusion == "ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS",
    # !specimenType == "Bone Marrow Aspirate",
    !is.na(dbgap_rnaseq_sample)
  )
# 查看筛选后的数量
table(aml_bm$dxAtInclusion)
table(aml_bm$specimenType)

# select aml_bm ----
expr_mat_filter <- expr_mat[,colnames(expr_mat)%in%aml_bm$dbgap_rnaseq_sample]
#expr_mat_filter <- expr_mat
# weights_LSC17 ----
weights_LSC17 <- c(
  GPR56    =  0.0501,
  AKR1C3   = -0.0402,
  CD34     =  0.0338,
  NGFRAP1  =  0.0465,
  EMP1     =  0.0146,
  C19orf77   = -0.0226,  # SMIM24 别名
  SOCS2    =  0.0271,
  CPXM1    = -0.0258,
  CDK6     = -0.0704,
  KIAA0125 =  0.0196,
  DPYSL3   =  0.0284,
  MMRN1    =  0.0258,
  LAPTM4B  =  0.00582,
  ARHGAP22 = -0.0138,
  NYNRIN   =  0.00865,
  ZBTB46   = -0.0347,
  DNMT3B   =  0.0874
)

missing_genes <- setdiff(names(weights_LSC17), rownames(expr_mat_filter))
if (length(missing_genes) > 0) {
  stop("Missing genes: ", paste(missing_genes, collapse = ", "))
}

LSC17_score <- colSums(
  sweep(expr_mat_filter[names(weights_LSC17), , drop = FALSE], 1, weights_LSC17, `*`)
)

library(matrixStats)
expr_rank <- rowRanks(as.matrix(expr_mat_filter), ties.method = "average")
score_rank <- rank(LSC17_score, ties.method = "average")
# Spearman 本质上就是对秩做 Pearson 相关
gene_rho <- as.vector(cor(t(expr_rank), score_rank, method = "pearson"))

n <- ncol(expr_mat_filter)
t_stat <- gene_rho * sqrt((n - 2) / (1 - gene_rho^2))
gene_p <- 2 * pt(-abs(t_stat), df = n - 2)

gene_cor_df <- data.frame(
  gene = rownames(expr_mat_filter),
  rho = gene_rho,
  pvalue = gene_p,
  FDR = p.adjust(gene_p, method = "BH")
)

gene_cor_df <- gene_cor_df[order(gene_cor_df$rho, decreasing = TRUE), ]
gene_cor_df[gene_cor_df$gene == "GSN",]
gene_rank <- gene_cor_df$rho
names(gene_rank) <- gene_cor_df$gene
write.csv(gene_cor_df, file = "./03_result/LSC17/Gene_cor_df.csv")

# 膜蛋白子集 ----
library(dplyr)
library(tibble)
set_tag <- "LASA"
aml_surface_genes <- read.delim("D:/R/R-Project/Database/Leucegene/LASA/percenage_presence_table.tsv")

# aml_surface_genes 是你下载/整理好的 AML 表面膜蛋白基因列表
aml_surface_cor_df <- gene_cor_df %>%
  filter(gene %in% aml_surface_genes$Gene.names) %>%
  arrange(desc(abs(rho)))

# 可选：只在 AML surface genes 内重新校正 FDR
aml_surface_cor_df <- aml_surface_cor_df %>%
  mutate(FDR_surface = p.adjust(pvalue, method = "BH"))
write.csv(aml_surface_cor_df, file = paste0("./03_result/LSC17/",set_tag,"_LSC17_cor_df.csv"))

# rank plot
library(dplyr)
library(ggplot2)
library(scales)
library(ggrepel)
plot_df <- aml_surface_cor_df %>%
  arrange(desc(rho)) %>%
  mutate(
    rank = row_number(),
    corr_t = rho * sqrt((n - 2) / (1 - rho^2)),
    corr_t_plot = pmax(pmin(corr_t, 20), -20)
  )
genes_to_label <- c("SCARF1", "CALCRL", "CD109", "F2RL1", "JAM3", "F2R", "CD34", "SV2A", "SPINT2", "PEAR1",
                    "SIRPB2", "INSR", "P2RY2", "FUT4", "ITGA7","IL31RA", "LRPAP1", "ATP9B", "SLC2A13", "P2RX4")

label_df <- plot_df %>%
  filter(gene %in% genes_to_label)

p_rank <- ggplot(plot_df, aes(x = rank, y = rho, color = corr_t_plot)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray55") +
  geom_point(size = 1.2, alpha = 0.9) +
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    color = "black",
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.color = "black"
  ) +
  scale_color_gradient2(
    low = "#6baed6",
    mid = "white",
    high = "#fb3b2f",
    midpoint = 0,
    limits = c(-20, 20),
    name = "LSC17 Corr.t",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(2.9, "cm"),
      barheight = unit(0.35, "cm")
  )) +
  scale_x_continuous(labels = comma) +
  coord_cartesian(ylim = c(-0.6, 0.65)) +
  labs(
    x = "Rank",
    y = "LSC17 Correlation"
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 15),
    legend.position = c(0.72, 0.88),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
p_rank
ggsave(
  filename = paste0("./03_result/LSC17/",set_tag,"_LSC17_AML_surface_rank_plot.pdf"),
  plot = p_rank,
  width = 5,
  height = 5,
  units = "in"
)
dev.off()
