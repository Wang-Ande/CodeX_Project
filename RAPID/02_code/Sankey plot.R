library(dplyr)
library(ggplot2)
library(ggalluvial)

# ---------------------------
# 1. 准备数据
# ---------------------------
low_abundance <- read.csv("./01_data/Plasma_Protein/Low_abundance_15ul_tidy.csv", row.names = 1)
low_abundance <- low_abundance[!is.na(low_abundance$Low_abundance_15ul),]
total_protein <- read.csv("./01_data/Plasma_Protein/Plasma_Total_protein_tidy.csv", row.names = 1)
total_protein <- total_protein %>%
  mutate(
    Full_mean = rowMeans(across(where(is.numeric)), na.rm = TRUE),
    Full_mean = ifelse(is.nan(Full_mean), NA, Full_mean)
  )
full_genes <- unique(total_protein$Gene)
low_genes  <- unique(low_abundance$Gene)
overlap_genes <- intersect(full_genes, low_genes)

overlap_rank <- data.frame(
  Gene = overlap_genes,
  Full_log2 = log2(total_protein$Full_mean[match(overlap_genes, total_protein$Gene)] + 1),
  Low_log2  = log2(low_abundance$Low_abundance_15ul[match(overlap_genes, low_abundance$Gene)] + 1)
)

overlap_rank <- overlap_rank %>%
  mutate(
    Full_quantile = ntile(Full_log2, 4),
    Low_quantile  = ntile(Low_log2, 4)
  )


full_genes <- total_protein$Gene
lap_genes  <- low_abundance$Gene
overlap_genes <- intersect(full_genes, lap_genes)
full_only <- setdiff(full_genes, lap_genes)
lap_only  <- setdiff(lap_genes, full_genes)

# 处理 overlap 蛋白
overlap_data <- overlap_rank %>%
  mutate(
    Full_quantile = factor(Full_quantile, levels=1:4, labels=c("Q1","Q2","Q3","Q4")),
    Low_quantile  = factor(Low_quantile,  levels=1:4, labels=c("Q1","Q2","Q3","Q4")),
    Type = "Overlap"
  )

# 处理 Full only 蛋白
full_only_data <- total_protein %>%
  filter(Gene %in% full_only) %>%
  mutate(
    Full_quantile = ntile(log2(Full_mean + 1), 4),
    Full_quantile = factor(Full_quantile, levels=1:4, labels=c("Q1","Q2","Q3","Q4")),
    Low_quantile = factor("Not detected", levels=c("Q4","Q3","Q2","Q1","Not detected")),
    Type = "Full_only"
  )

# 处理 LAP only 蛋白
lap_only_data <- low_abundance %>%
  filter(Gene %in% lap_only) %>%
  mutate(
    Low_quantile = ntile(log2(Low_abundance_15ul + 1), 4),
    Low_quantile = factor(Low_quantile, levels = 1:4, labels = c("Q1","Q2","Q3","Q4")),
    Full_quantile = factor("Not detected", levels = c("Q4","Q3","Q2","Q1","Not detected")),
    Type = "LAP_only"
  )
# 合并数据
alluvial_data_full <- bind_rows(overlap_data, full_only_data, lap_only_data)
# 调整 Low_quantile 顺序
alluvial_data_full <- alluvial_data_full %>%
  mutate(
    Full_quantile = factor(Full_quantile, levels = c("Q4", "Q3", "Q2", "Q1", "Not detected")),
    Low_quantile  = factor(Low_quantile,  levels = c("Q4", "Q3", "Q2", "Q1", "Not detected"))
  )
# ---------------------------
# 2. 绘制 Alluvial 图
# ---------------------------
# 确保 Protein_ID 列存在
alluvial_data_full <- alluvial_data_full %>%
  mutate(Protein_ID = Gene)
alluvial_data_full <- alluvial_data_full[,-c(7:11)]
# 1. 给 overlap 蛋白指定 Flow 颜色
alluvial_data_full <- alluvial_data_full %>%
  mutate(
    Flow_Color = case_when(
      Type == "Overlap" & Full_quantile == "Q4" ~ "#1F4E6D",  # 深蓝
      Type == "Overlap" & Full_quantile == "Q3" ~ "#4F86A5",
      Type == "Overlap" & Full_quantile == "Q2" ~ "#90B4C8",
      Type == "Overlap" & Full_quantile == "Q1" ~ "#B0C8D8",  # 浅蓝加深
      Type == "Full_only" ~ "#c0c0c0",
      Type == "LAP_only"  ~ "#b0b0b0",
      TRUE ~ "grey50"
    )
  )
# 2. 转换为 lodes form
lodes <- to_lodes_form(alluvial_data_full,
                       key = "axis",
                       axes = c("Full_quantile", "Low_quantile"),
                       id = "Protein_ID") %>%
  rename(alluvium = Protein_ID) # ggalluvial 要求列名叫 alluvium


# 绘图
pdf("./03_result/Plasma_Protein/Intersect Sankey plot.pdf", width=6, height=5.5)
ggplot(lodes,
       aes(x = axis, stratum = stratum, alluvium = alluvium, fill = Flow_Color, y = 1)) +
  geom_alluvium(alpha = 0.8, width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_fill_identity() +   # 直接使用 Flow_Color 中的颜色
  labs(x = "Proteome", y = "Number of proteins") +
  theme_bw(base_size = 13) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
dev.off()

