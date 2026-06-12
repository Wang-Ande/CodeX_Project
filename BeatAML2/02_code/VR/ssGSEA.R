# library packages
library(GSVA)
library(msigdbr)
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(openxlsx)

# set output dir
dir_ssGSEA <- "./03_result/VR/ssGSEA/"
if (!dir.exists(dir_ssGSEA)) {
  dir.create(dir_ssGSEA)
  print(paste("Folder", dir_ssGSEA, "created."))
} else {
  print(paste("Folder", dir_ssGSEA, "already exists."))
}


# data input
expr_matrix <- read.csv("./01_data/Venetoclax/Venetoclax 20propsubgroup Norm data.csv", row.names = 1)
head(expr_matrix)

# 转 numeric matrix
expr_matrix <- as.matrix(expr_matrix)

# 过滤低表达蛋白

dim(expr_matrix)
summary(expr_matrix[, 1])

# group input
sample_info <- read.xlsx("./01_data/Venetoclax/Venetoclax 20propsubgroup Response.xlsx")
sample_info <- sample_info[,-2]
sample_info <- sample_info[order(sample_info$group, decreasing = TRUE), ]
head(sample_info)

# gene pathway set
# 从MSigDB提取通路，也可以是GO KEGG 自定义基因集
msig <- msigdbr(species = "Homo sapiens")
a <- msig[grep("SERINE", msig$gs_name),] # 大写

pathway_list <- list(
  Serine1 = msig %>%
    filter(gs_name == "GOBP_L_SERINE_BIOSYNTHETIC_PROCESS") %>%
    pull(gene_symbol),

  Serine2 = msig %>%
    filter(gs_name == "GOBP_L_SERINE_CATABOLIC_PROCESS") %>%
    pull(gene_symbol),

  Serine3 = msig %>%
    filter(gs_name == "GOBP_SERINE_FAMILY_AMINO_ACID_CATABOLIC_PROCESS") %>%
    pull(gene_symbol),

  Serine4 = msig %>%
    filter(gs_name == "GOBP_SERINE_TRANSPORT") %>%
    pull(gene_symbol),

  Serine5 = msig %>%
    filter(gs_name == "REACTOME_SERINE_METABOLISM") %>%
    pull(gene_symbol),
  
  Serine6 = msig %>%
    filter(gs_name == "KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM") %>%
    pull(gene_symbol),
  
  Serine7 = msig %>%
    filter(gs_name == "WP_SERINE_METABOLISM") %>%
    pull(gene_symbol)%>%
    unique()
)
summary(pathway_list)
save(pathway_list, file = paste0(dir_ssGSEA, "Merge_pathway_list.RData"))
pathway_long <- enframe(pathway_list, name = "pathway", value = "gene_list") %>%
  unnest(cols = c(gene_list))
write.csv(pathway_long, paste0(dir_ssGSEA, "Merge_pathway_list.csv"), row.names = FALSE)

# 已有signature
signature <- read.xlsx("D://R/R-Project/tools/Human/Signature/AML Hierarchies signature Cell 2019.xlsx")
# signature <- signature[c(1:30),]
pathway_list <- lapply(signature, function(x) {
  x <- na.omit(x)
  x <- trimws(as.character(x))
  x <- x[x != ""]
  unique(x)
})

# ssGSEA scores
# 先定义ssGSEA分析的参数
ssGSEA_param <- ssgseaParam(
  exprData = expr_matrix,
  geneSets = pathway_list,
  alpha = 0.25,
  normalize = TRUE,
  checkNA = "yes",
  use = "na.rm",
  verbose = TRUE
)

# 将设置好的参数传给gsva函数进行ssGSEA分析
ssGSEA_mat <- gsva(
  ssGSEA_param,            # 设置好的分析参数
  verbose = TRUE           # 将给出分析过程的信息
)


# pathway score heatmap
# annotation_col <- sample_info %>%
#   column_to_rownames("id")
# 
# pdf(paste0(dir_ssGSEA, "Merge activity heatmap.pdf"), width = 8, height = 3)
# pheatmap(
#   ssGSEA_mat,
#   scale = "row",cluster_cols = FALSE,
#   clustering_method = "complete",
#   clustering_distance_rows = "correlation",
#   clustering_distance_cols = "correlation",
#   annotation_col = annotation_col,
#   border_color = NA,
#   show_colnames = TRUE,     # 样本多时建议关掉
#   fontsize_row = 10,
#   fontsize_col = 9,
#   main = "Pathway activity (ssGSEA)"
# )
# dev.off()

library(ComplexHeatmap)
library(circlize)
library(grid)

# 1) 等价于 pheatmap(scale="row")
mat <- ssGSEA_mat
mat_z <- t(scale(t(mat)))   # row-wise z-score
colnames(mat_z)
mat_z <- mat_z[,sample_info$sample_id]
# 调整排列顺序（不适用行聚类）
rownames(mat_z)
mat_z <- mat_z[c("HSC-like", "Progenitor-like", "GMP-like", "Promono-like", "Monocyte-like", "cDC-like"), ]

# 2) 颜色（和你图里类似的蓝-白-红）
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#327eba", "white", "#e06663"))
signature <- "Cell_status_Cell2019"

group <- sample_info$group
names(group) <- colnames(mat_z)
group <- factor(group, levels = c("Sensitivity", "Resistance"))
group_col <- c(
  "Sensitivity" = "#4F7C82",
  "Resistance" = "#9B5F6D"
)

ord <- order(group)
mat_z2 <- mat_z[, ord]
group2 <- group[ord]

ha_col <- HeatmapAnnotation(
  Group = group2,
  col = list(Group = group_col),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10),
  simple_anno_size = unit(4, "mm")
)

pdf(paste0(dir_ssGSEA, signature, "_activity heatmap.pdf"),
    width = 5.5, height = 3.5)

ht <- Heatmap(
  mat_z2,
  name = "Z-score",
  col = col_fun,
  
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  
  column_split = group2,
  gap = unit(2, "mm"),
  column_title = NULL,
  
  show_row_names = TRUE,
  show_column_names = FALSE,
  
  top_annotation = ha_col,
  row_names_gp = gpar(fontsize = 10),
  
  border = FALSE,
  
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter"
  )
)

draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legend = TRUE
)

dev.off()

# boxplot by groups
df_long <- as.data.frame(t(ssGSEA_mat)) %>%
  rownames_to_column("sample_id") %>%
  left_join(sample_info, by = "sample_id") %>%
  pivot_longer(
    cols = -c(sample_id, group),
    names_to = "Pathway",
    values_to = "Score"
  )
table(df_long$group)
# plot
df_long$group <- factor(df_long$group, levels= c("Sensitivity", "Resistance"))
my_cols <- c(
  Sensitivity = "#4DA3FF",        # 蓝
  Resistance = "#F8766D"  # 绿
)
my_comparisons <- list(
  c("Sensitivity", "Resistance")
)


p1 <- ggplot(df_long, aes(x = group, y = Score, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 1) +
  facet_wrap(~ Pathway, scales = "free_y", ncol = 3) +
  stat_compare_means(
    method = "t.test",
    comparisons = my_comparisons,
    label = "p.signif",
    label.y.npc = 0.95
  ) +
  scale_fill_manual(values = my_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(p1)
ggsave(paste0(signature, "_activity boxplot.pdf"), plot = p1, device = "pdf",
       path = dir_ssGSEA, width = 5.5, height = 4.5)
dev.off()

# FDR compute
stat_res <- df_long %>%
  group_by(Pathway) %>%
  wilcox_test(Score ~ Group) %>%
  adjust_pvalue(method = "BH") %>%
  mutate(p.adj = p.adj)

stat_res

# 
# library(ggplot2)
# library(ggpubr)
# library(dplyr)
# 
# # 建议固定组顺序（避免 facet 内顺序乱）
# df_long <- df_long %>%
#   mutate(group = factor(group, levels = c("WT", "GilR300nm", "GilR3000nm")))

# 和你热图一致的配色（蓝/红/绿的顺序你按实际想要的改）
my_cols <- c(
  WT = "#4DA3FF",        # 蓝
  GilR300nm = "#00BA38",  # 红
  GilR3000nm = "#F8766D"  # 绿
)

# p1 <- ggplot(df_long, aes(x = group, y = Score, fill = group)) +
#   geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.5, color = "black") +
#   geom_point(
#     position = position_jitter(width = 0.12, height = 0),
#     size = 1.6, alpha = 0.8
#   ) +
#   facet_wrap(~ Pathway, scales = "free_y", ncol = 3) +
#   stat_compare_means(
#     method = "wilcox.test",
#     comparisons = my_comparisons,
#     label = "p.signif",
#     tip.length = 0.01,     # 线条末端更短更精致
#     label.y.npc = "top",  # 统一显示高度
#     size = 3.6
#   ) +
#   scale_fill_manual(values = my_cols) +
#   labs(x = NULL, y = "ssGSEA score") +
#   theme_bw(base_size = 11) +
#   theme(
#     legend.position = "none",
#     strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
#     strip.text = element_text(size = 11, face = "bold"),
#     panel.grid = element_blank(),
#     panel.border = element_rect(color = "black", linewidth = 0.6),
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#     axis.title.y = element_text(margin = margin(r = 6)),
#     plot.margin = margin(6, 10, 6, 6),
#     panel.spacing = unit(0.9, "lines")
#   )
# 
# p1


## =========================
## 0) packages
## =========================
suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggalluvial)
})

## =========================
## 1) 输入：ssGSEA_mat
##    行 = pathway, 列 = sample
## =========================
# ssGSEA_mat <- ...  # 你已经有了

stopifnot(is.matrix(ssGSEA_mat) || is.data.frame(ssGSEA_mat))
ss <- as.matrix(ssGSEA_mat)

## =========================
## 2) 构建样本分组信息 meta
##    依据列名匹配 WT / GilR300nm / GilR3000nm
## =========================
samples <- colnames(ss)

group_map <- case_when(
  str_detect(samples, "^WT_") ~ "WT",
  str_detect(samples, "^GilR300nm_") ~ "ER",
  str_detect(samples, "^GilR3000nm_") ~ "LR",
  TRUE ~ NA_character_
)

meta <- data.frame(
  sample = samples,
  group  = factor(group_map, levels = c("WT","ER","LR")),
  stringsAsFactors = FALSE
)

if (any(is.na(meta$group))) {
  stop("有些样本列名没有被分到组里：\n",
       paste(meta$sample[is.na(meta$group)], collapse = ", "),
       "\n请检查列名规则或修改 group_map。")
}

## =========================
## 3) limma：在通路分数层面做差异
##    Early vs WT ; Late vs Early
## =========================
design <- model.matrix(~0 + group, data = meta)
colnames(design) <- levels(meta$group)

fit <- lmFit(ss, design)
contr <- makeContrasts(
  ER_vs_WT   = ER - WT,
  LR_vs_ER = LR  - ER,
  levels = design
)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)

res_early <- topTable(fit2, coef="ER_vs_WT", number=Inf, sort.by="none")
res_late  <- topTable(fit2, coef="LR_vs_ER", number=Inf, sort.by="none")

## =========================
## 4) 打方向标签：Up/Down/NS
## =========================
tag_dir <- function(logFC, P.Value, cutoff=0.05) {
  ifelse(P.Value < cutoff & logFC > 0, "↑",
         ifelse(P.Value < cutoff & logFC < 0, "↓", "≈"))
}

df <- data.frame(
  pathway     = rownames(ss),
  early_logFC = res_early$logFC,
  early_FDR   = res_early$P.Value,
  late_logFC  = res_late$logFC,
  late_FDR    = res_late$P.Value,
  stringsAsFactors = FALSE
) %>%
  mutate(
    ER_dir = tag_dir(early_logFC, early_FDR, cutoff = 0.05),
    LR_dir  = tag_dir(late_logFC,  late_FDR,  cutoff = 0.05)
  )

## =========================
## 5) 定义 Category（左列）
##    你目前 6-7 条通路可以先用“手工映射”
##    以后通路多了也可以扩展规则
## =========================
df <- df %>%
  mutate(Category = case_when(
    pathway %in% c("OXPHOS", "Mito_translation", "Cristae_MICOS") ~ "Mito structure / OXPHOS",
    pathway %in% c("TCA_cycle") ~ "TCA / Metabolism",
    pathway %in% c("Glycolysis") ~ "Glycolysis",
    pathway %in% c("Mitophagy") ~ "Stress / Mitophagy",
    pathway %in% c("RAS_MAPK") ~ "RAS / MAPK",
    TRUE ~ "Other"
  ))

## =========================
## 6) 生成 alluvial 数据：Category -> Early_dir -> Late_dir
##    y = 1（每条通路一条“丝带”），
##    通路多了就会更饱满；现在少也没问题
## =========================
allu <- df %>%
  transmute(
    pathway,
    axis1 = Category,
    axis2 = paste0("ER", ER_dir),
    axis3 = paste0("LR",  LR_dir),
    y     = 1
  )

## 为了让方向顺序固定（↑在上/或你想要的顺序）
allu$axis2 <- factor(allu$axis2, levels = c("ER ↑","ER ≈","ER ↓"))
allu$axis3 <- factor(allu$axis3, levels = c("LR ↑","LR ≈","LR ↓"))

## =========================
## 7) 画 “标书风” Alluvial（ggalluvial）
## =========================
p_alluvial <- ggplot(allu,
                     aes(axis1 = axis1, axis2 = axis2, axis3 = axis3, y = y)
) +
  geom_alluvium(aes(fill = axis1), width = 1/10, alpha = 0.85) +
  geom_stratum(width = 1/7, color = "grey30") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
  scale_x_discrete(limits = c("Category", "ER vs WT", "LR vs Early"),
                   expand = c(0.08, 0.08)) +
  labs(title = "Mechanistic rewiring (ssGSEA)",
       y = "Number of pathways") +
  theme_classic(base_size = 13) +
  theme(
    axis.title.x = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )

print(p_alluvial)

## =========================
## 8) 可选：净变化条形图（放在 alluvial 旁边作为统计支撑）
##    这里用各 Category 的平均 logFC
## =========================
net <- df %>%
  group_by(Category) %>%
  summarise(
    net_Early = mean(early_logFC, na.rm = TRUE),
    net_Late  = mean(late_logFC,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(net_Early, net_Late),
               names_to = "Stage", values_to = "Mean_logFC") %>%
  mutate(Stage = recode(Stage,
                        net_ER = "ER vs WT",
                        net_LR  = "LR vs Early"))

p_net <- ggplot(net, aes(x = Category, y = Mean_logFC, fill = Stage)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_classic(base_size = 13) +
  labs(title = "Net change by category (mean logFC)",
       x = NULL, y = "Mean logFC")

print(p_net)

## =========================
## 9) 输出一个汇总表（方便你写结果/补图注）
## =========================
df_out <- df %>%
  select(pathway, Category,
         ER_dir, early_logFC, early_FDR,
         LR_dir,  late_logFC,  late_FDR) %>%
  arrange(Category, pathway)

print(df_out)


