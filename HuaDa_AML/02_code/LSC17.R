lsc17_genes <- c(
  "DNMT3B", "ZBTB46", "NYNRIN", "ARHGAP22", "LAPTM4B",
  "MMRN1", "DPYSL3", "FAM30A", "CDK6", "CPXM1",
  "SOCS2", "SMIM24", "EMP1", "BEX3", "CD34",
  "AKR1C3", "ADGRG1"
)
lsc_df <- as.data.frame(lsc17_genes)
# 行是gene symbol，列是AML样本
# values是log2转换后的蛋白丰度，例如log2 LFQ / log2 intensity / normalized TMT abundance
prot_mat <- read.csv("./01_data/Total_pg_norm_filled_normed.csv", row.names = 1)
# 检查实际测到了多少个
lsc17_in_data <- intersect(lsc17_genes, rownames(prot_mat))
lsc17_missing <- setdiff(lsc17_genes, rownames(prot_mat))

lsc17_in_data
lsc17_missing
length(lsc17_in_data)

library(dplyr)
library(ggplot2)
library(ggpubr)
out_dir <- "./03_result/04_LSC17_score/"
dir.create(out_dir)

# 1. 对每个蛋白跨样本做z-score
prot_mat <- log2(prot_mat)
prot_z <- t(scale(t(prot_mat)))

# 2. 计算LSC17-like score
# LSC17-like score = mean(z-scored abundance of detected LSC17 proteins)
# 先把每个LSC17相关蛋白在样本间标准化成 z-score，再对每个样本的这些蛋白 z-score 取平均。
# 这个方法在蛋白组里是可以用的，因为蛋白组经常无法完整检测到所有17个基因对应蛋白，也不一定适合直接套用RNA层面的回归系数。
# 但它不是原始 LSC17 score，而是“基于LSC17基因集的蛋白组干性评分”。
lsc17_score <- colMeans(prot_z[lsc17_in_data, , drop = FALSE], na.rm = TRUE)

# 3. 提取目标蛋白表达
gene_interest <- "GSN"

plot_df <- data.frame(
  Sample = colnames(prot_mat),
  LSC17_score = as.numeric(lsc17_score),
  Gene_expression = as.numeric(prot_mat[gene_interest, ])
)

# 4. 去除缺失值
plot_df <- plot_df %>%
  filter(!is.na(LSC17_score), !is.na(Gene_expression))

cor_res <- cor.test(
  plot_df$LSC17_score,
  plot_df$Gene_expression,
  method = "pearson"
)

r_value <- unname(cor_res$estimate)
p_value <- cor_res$p.value

label_text <- paste0(
  "R = ", round(r_value, 2),
  ", P ", ifelse(p_value < 0.001, "< 0.001", paste0("= ", signif(p_value, 2)))
)

label_text

# corr plot
p <- ggplot(plot_df, aes(x = LSC17_score, y = Gene_expression)) +
  geom_point(
    size = 3,
    alpha = 0.45,
    color = "#C95F5F"
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    linetype = "dashed",
    linewidth = 0.6,
    fill = "grey80",
    alpha = 0.35
  ) +
  annotate(
    "text",
    x = min(plot_df$LSC17_score, na.rm = TRUE),
    y = min(plot_df$Gene_expression, na.rm = TRUE),
    label = label_text,
    hjust = 0,
    vjust = -0.5,
    size = 5,
    fontface = "italic"
  ) +
  labs(
    x = "LSC17 score",
    y = paste0(gene_interest, " expression")
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 18, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.line = element_line(linewidth = 0.6, color = "grey40"),
    plot.margin = margin(10, 10, 10, 10)
  )

p

ggsave(
  filename = paste0(out_dir, gene_interest, "_LSC17_correlation.pdf"),
  plot = p,
  width = 4.5,
  height = 4.2,
  useDingbats = FALSE
)
dev.off()
