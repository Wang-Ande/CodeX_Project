library(openxlsx)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Data input ----
# Load data (replace with BeatAML2 "Raw" GitHub links from the webpage)
clinical <- read.xlsx("./01_data/Clinical Summary.xlsx")
auc      <- read.csv("01_data/Inhibitor AUC values.csv")
mut      <- read.csv("01_data/WES targeted Sequencing Mutation Calls.csv")
expr     <- read.csv("./01_data/Normalized Expression.csv")   # typically genes x samples or long format
expr     <- expr[,-c(1,3,4)]
# 1. 处理 NA 值
expr <- expr[!is.na(expr$display_label), ]  # remove NA rows
# 2. 处理重复基因
# 1️⃣ 先算所有行的均值（向量化，非常快）
row_mean <- rowMeans(expr[ , colnames(expr) != "display_label"])
# 2️⃣ 按 display_label + 均值排序
ord <- order(expr$display_label, -row_mean)
expr_sorted <- expr[ord, ]
# 3️⃣ 每个 display_label 只保留第一行（即均值最大）
expr_clean <- expr_sorted[!duplicated(expr_sorted$display_label), ]

write.csv(expr_clean, file = "./01_data/Venetocalx/Normalized Expression cleaned.csv", row.names = FALSE)

# clean version
expr <- read.csv("./01_data/Venetocalx/Normalized Expression cleaned.csv", row.names = 1)
expr <- as.data.frame(expr_clean)
rownames(expr) <- expr$display_label
# Group input ----
# 根据文献中的介绍，按照分别取AUC的前后20%作为resistance和sensitivity
# 病人样本仅取在sensitive or resistant中的样本
#  Get gilteritinib AUC
gilt_auc <- auc %>%
  filter(inhibitor == "Venetoclax") %>%   # adjust field name
  select(sample_id = dbgap_rnaseq_sample, auc)
flt3_itd_samples <- clinical %>%
  filter(`FLT3-ITD` == "positive") %>%   # adjust field name to actual column
  select(sample_id = dbgap_rnaseq_sample)
gilt_auc <- gilt_auc[gilt_auc$sample_id%in%flt3_itd_samples$sample_id,]

df <- gilt_auc
df <- df %>%
  mutate(group = case_when(
    auc <= quantile(auc, 0.20, na.rm = TRUE) ~ "Sensitivity",    
    auc >= quantile(auc, 0.80, na.rm = TRUE) ~ "Resistance",
    TRUE ~ "Moderate"
  ))                                    # quantile按升序排列c(1,2,3,4)

# 查看 n
table(df$group)
select_Response <- df
Gene_select <- expr[,colnames(expr)%in%select_Response$sample_id]

# 筛选只在sensitive和resistant中的转录组测序样本
filtered_id <- select_Response[select_Response$group!= "Moderate", "sample_id"]
Gene_select <- expr[,colnames(expr)%in%filtered_id]

select_Response <- select_Response[select_Response$sample_id %in% colnames(Gene_select),]
table(select_Response$group)

# res output
write.csv(Gene_select, file = "./01_data/Venetoclax/Venetoclax 20propsubgroup Norm data.csv")
write.xlsx(select_Response, file = "./01_data/Venetoclax/Venetoclax 20propsubgroup Response.xlsx")


# Filter low expr  
# approach 1
min_sam <- 3   # 一半样本中表达
filtered_data <- expr[rowSums(expr > 1) >= min_sam, ] # 1392
