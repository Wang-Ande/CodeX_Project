# install.packages("multiUS")
library(multiUS)
library(RColorBrewer)
library(dplyr)
library(data.table)
library(ggpubr)
library(ggrepel)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(matrixStats) 
Sys.setenv(LANGUAGE="en")
options(stringsAsFactors = F)
# 1. Data input ----
data_input <- read.csv("./01_data/AML_PBMC_pg.csv")
colnames(data_input)
data_input <- data_input[,-c(1,3)]
colnames(data_input)[c(1,2:ncol(data_input))] <- c("Gene",
                                                 paste0("AML_",1:170)
                                                 )

data_input$Gene <- unlist(lapply(data_input$Gene,function(a) unlist(strsplit(a,split = ";"))[1]))
##----------------------------去除重复的蛋白------------------------------------
# 基于缺失值删除，按照缺失数量排序、去重
SampleType <- c("AML")
index2 <- min(grep(paste0("^", SampleType, collapse="|"), colnames(data_input))) # 从样本列开始
data_input[,(index2):ncol(data_input)][data_input[,(index2):ncol(data_input)]=="NaN"] <- NA
missing <- data.frame(num=row.names(data_input),missing=apply(data_input[,2:ncol(data_input)],1,function(x) sum(is.na(x))))
missing <- missing[order(missing$missing,decreasing = F),]
data_input <- data_input[match(missing$num,row.names(data_input)),]
data_input <- data_input[which(!duplicated(data_input$Gene)),]
data_input <- data_input[!is.na(data_input$Gene), ] # 删除gene列为 NA的 
row.names(data_input) <- data_input$Gene
data_input <- data_input[,-1]
write.csv(data_input, file = "./01_data/Total_pg_Report.csv")

# 查看批次效应
pdf("./03_result/01_QC/Boxplot_before_norm.pdf", width = 10, height = 5)
cols <- brewer.pal(8, "Set1")
par(mar = c(10, 4, 4, 2))  # 下边界加大
boxplot(log2(data_input),
        col=cols,
        las=2,outline=F,
        ylab="log2(intensity)",
        main="Before Normalization",
        xaxt="n")   # 不画默认 x 轴

dev.off()

# 计算NA值比例
NA_ratio <- colSums(is.na(data_input))/dim(data_input)[1]
NA_ratio <- as.data.frame(NA_ratio)
NA_ratio$Samples <- rownames(NA_ratio)

# 设定NA值比例分类
NA_ratio[NA_ratio$NA_ratio < 0.2,"group"] <- "Good"
NA_ratio[NA_ratio$NA_ratio < 0.5&NA_ratio$NA_ratio > 0.2,"group"] <- "OK"
NA_ratio[NA_ratio$NA_ratio < 0.8&NA_ratio$NA_ratio > 0.5,"group"] <- "Bad"
NA_ratio[NA_ratio$NA_ratio > 0.8,"group"] <- "Remove"
NA_ratio$group <- factor(NA_ratio$group,levels = c("Good","OK","Bad","Remove"))

# 绘制NA值比例分布图
library(ggplot2)
p <- ggplot(NA_ratio,aes(x = NA_ratio, y = reorder(Samples,NA_ratio,decreasing = T), fill = group)) + 
  geom_bar(stat = "identity", color = "black") + 
  geom_text(aes(label = sprintf("%.2f", NA_ratio)),size = 2, hjust = -0.5) + 
  scale_fill_manual(
    name = "QC",
    values = c("Good" = "#62c882",
               "OK" = "#b5e281",
               "Bad" = "#febe80",
               "Remove" = "#bb2022")) + 
  theme_bw() + 
  labs( x = "Number of missing proteins",
        y = "Samples") +
  scale_x_continuous(
    limits = c(0, max(NA_ratio$NA_ratio) * 1.12)  # 右边多留 15%
  )

print(p)
ggsave(filename = "Sample_NA_Ratio.pdf", device = pdf, plot = p,
       path = "./03_result/01_QC/", dpi = 300)
dev.off()

# 去除缺失比例较高样本
bad_samples <- NA_ratio[NA_ratio$group%in%c("Bad", "Remove"),2]
data_input <- data_input[,!colnames(data_input)%in%bad_samples]
pdf("./03_result/01_QC/Boxplot_filter.pdf", width = 10, height = 5)
cols <- brewer.pal(8, "Set1")
par(mar = c(10, 4, 4, 2))  # 下边界加大
boxplot(log2(data_input),
        col=cols,
        las=2,outline=F,
        ylab="log2(intensity)",
        main="After filter",
        xaxt="n")   # 不画默认 x 轴

dev.off()

#---------------------2. 缺失值填补------------------
# 是否进行log2运算？
data_input <- log2(data_input)

# 是否进行median normalization运算？
# 计算每列的中位数
median_values <- apply(data_input, 2, median, na.rm = TRUE)
# 计算全局中位数（所有列中位数的中位数）
global_median <- median(median_values, na.rm = TRUE)
# 计算归一化因子：每列中位数 / 全局中位数
nor_facter <- median_values / global_median
# 对数据进行全局中位数归一化
nor_data_input <- data.frame(t(t(data_input) / nor_facter), check.names = FALSE)

pdf("./03_result/01_QC/Boxplot_after_norm.pdf", width = 10, height = 5)
cols <- brewer.pal(8, "Set1")
par(mar = c(10, 4, 4, 2))  # 下边界加大
boxplot(nor_data_input,
        col=cols,
        las=2,outline=F,
        ylab="log2(intensity)",
        main="After Normalization",
        xaxt="n")   # 不画默认 x 轴
# 手动画 45 度倾斜标签
# text(x = 1:ncol(nor_data_input),
#      y = par("usr")[3] - 0.5,   # 往下移动一些
#      labels = colnames(nor_data_input),
#      srt = 45,                  # 45°倾斜
#      adj = 1,
#      xpd = TRUE)
dev.off()

## ------------------ 0) 输入 ------------------
# nor_data_input：行=gene，列=样本（形如 EV_R1, OE_A175fs_R2, ...）
# SampleType：c("EV", "OE_A175fs", "OE_R402C", "OE_WT")
#在一组中缺失值60%（根据n数量选择）及以上的位点则使用填补
dat <- nor_data_input
SampleType <- c("AML")

#--------------------筛选可定量蛋白---------------
index <- list()
for (i in 1:nrow(data.frame(SampleType))) {
  num <- which(colnames(dat) %like% SampleType[i])
  index[[i]] <- apply(dat[,num], 1, function(x) sum(is.na(x)))*100/nrow(data.frame(num))
}

index1 <- data.frame(do.call(cbind,index),check.names = FALSE)
colnames(index1) <- c("AML")

#过滤至少有一组缺失值不高于25%的位点
dat_1 <- data.frame(dat[which(apply(index1, 1, function(x) sum(x<=30))>0),], check.names = FALSE)
nrow(dat);nrow(dat_1)

## ------------------ 1)（可选）删除指定离群样本 ------------------
# 若你已人工判断出异常样本：
#drop_samples <- c("M13_VR_sgK1_1", "M13_VR_sgSC_1")  # 例如 c("OE_R402C_R2", "EV_R3")
#dat_1 <- dat_1[, !(colnames(dat_1) %in% drop_samples)]


## ------------------- 2) 创建分组索引列表 ---------------------------------
SampleType <- c("AML")
group_cols <- lapply(SampleType, function(group) {
  grep(paste0("^", group, "_"), colnames(dat_1))
})
group_cols

## ------------------- 3) 先10% 分位数填充 ---------------------------------
# 遍历每个蛋白（行）
for (i in 1:nrow(dat_1)) {
  # 遍历每个处理组
  for (j in seq_along(SampleType)) {
    # 提取当前组的所有样本值
    group_values <- dat_1[i, group_cols[[j]]]
    
    # 计算缺失值比例
    na_ratio <- sum(is.na(group_values)) / length(group_values)
    
    # 当缺失比例 ≥ 75% 且存在非缺失值时 比例动态变化，根据组内平行重复个数
    if (na_ratio >= 0.60) {
      # 计算该组内所有非缺失值的 10% 分位数（即低丰度信号）
      fill_val <- quantile(
        unlist(dat_1[, group_cols[[j]]]), 
        probs = 0.1, 
        na.rm = TRUE
      )
      # or
      # 组内最小值填充
      # fill_val <- min(
      #   unlist(dat_1[, group_cols[[j]]]),
      #   na.rm = TRUE
      #   )
      
      dat_1[i, group_cols[[j]]][is.na(group_values)] <- fill_val
    }
  }
}

## ------------------- 4) 剩余缺失值seqKNN填补 ---------------------------------
# library(DMwR2)
# data <- knnImputation(dat)
library(multiUS)
data_fill <- multiUS::seqKNNimp(data = dat_1, k = 10)
min(data_fill)
data_fill <- 2^data_fill
pdf("./03_result/01_QC/Boxplot_after_fill.pdf", width = 10, height = 5)
cols <- brewer.pal(8, "Set1")
par(mar = c(10, 4, 4, 2))  # 下边界加大
boxplot(log2(data_fill),
        col=cols,
        las=2,outline=F,
        ylab="log2(intensity)",
        main="After fill",
        xaxt="n")   # 不画默认 x 轴
# 手动画 45 度倾斜标签
# text(x = 1:ncol(data_fill),
#      y = par("usr")[3] - 0.5,   # 往下移动一些
#      labels = colnames(data_fill),
#      srt = 45,                  # 45°倾斜
#      adj = 1,
#      xpd = TRUE)
dev.off()

write.csv(data_fill,file = "./01_data/Total_pg_norm_filled.csv")

# 根据箱线图结果判断是否需要再次归一化
global_median <- median(apply(data_fill, 2, median, na.rm = TRUE))
nor_facter <- apply(data_fill, 2, median, na.rm = TRUE) / global_median
data_fill_norm <- data.frame(t(t(data_fill) / nor_facter), check.names = FALSE)
pdf("./03_result/01_QC/Boxplot_after_fill_norm.pdf", width = 10, height = 5)
par(mar = c(10, 4, 4, 2))  # 下边界加大
cols <- brewer.pal(8, "Set1")
boxplot(log2(data_fill_norm),
        col=cols,las=2,
        outline=F,
        ylab="log2(intensity)",
        main="NA Filled Norm",
        xaxt="n")   # 不画默认 x 轴
# 手动画 45 度倾斜标签
# text(x = 1:ncol(data_fill_norm),
#      y = par("usr")[3] - 0.5,   # 往下移动一些
#      labels = colnames(data_fill_norm),
#      srt = 90,                  # 45°倾斜
#      adj = 1,
#      xpd = TRUE)
dev.off()
write.csv(data_fill_norm, file = "./01_data/Total_pg_norm_filled_normed.csv")

# PCA ----
pca_data <- log2(data_fill)
a1=t(pca_data)
groups=data.frame(group=rep(c("AML"),c(150)))
rownames(groups) = row.names(a1)

data1 <- cbind(groups,a1)
df_pca <- PCA(data1[,-1], graph = FALSE)
data1$group <- factor(data1$group, levels = c("AML"))

dir <- "./03_result/01_QC/"
pdf(file = paste0(dir,"PCA.pdf"), # _de-protein
    width = 8,
    height = 7.5)
fviz_pca_ind(df_pca,
             geom = c("point", "text"),
             col.ind = data1$group,
             palette = c("#E60012"),
             addEllipses = FALSE,
             mean.point = FALSE,           # 去除分组的中心点
             title="",
             repel = TRUE, 
             legend.title="Group")+
  theme_bw()+
  # 修改图例
  theme(legend.title = element_text(size = 14, face = "plain"), # 图例标题大小
        legend.text  = element_text(size = 12),                # 图例文字大小
        legend.position = "right")+           
  guides(color = guide_legend(override.aes = list(size=5))) + # 放大图例点
  theme(legend.text = element_text(size = 12))
dev.off()

# Pearson ----
cor_data <- data_fill
res <- round(cor(cor_data,method="pearson",use="pairwise.complete.obs"),digits = 2)
min(res);max(res[which(res<1)])
res <- as.matrix(res)
mycol <- c("#4DBBD57F","#3C54887F","#DC00007F")

dir_cor <- "./03_result/01_QC/"
pdf(file = paste0(dir_cor,"Pearson_corr.pdf"), # _de-protein
    width = 6,
    height = 5)
corrplot(res,
         method ="circle",#指定相关系数以圆的形式展示
         addgrid.col="black",#方框的颜色
         tl.col="black",#字体颜色
         tl.pos = "d",#字体放在对角线
         tl.cex=0.4,#字体的大小
         cl.cex=0.5,
         number.cex=0.8,
         is.corr = F,
         col=colorRampPalette(mycol)(110),
         col.lim=c(0,1),
         type ="lower",#只展示矩阵的下半部分
         addCoefasPercent = F)#以小数点的形式表示相关系数

corrplot(res,
         method ="number",#指定相关系数以圆的形式展示
         addgrid.col="black",#方框的颜色
         tl.col="black",#字体颜色
         tl.pos = "d",#字体放在对角线
         tl.cex=0.4,#字体的大小
         cl.cex=0.5,
         number.cex=0.8,
         add = TRUE,#上下两个图叠加在一起
         is.corr = F,
         col=colorRampPalette(mycol)(110),
         col.lim=c(0,1),
         type ="upper",#只展示矩阵的上半部分
         addCoefasPercent = F)#以小数点的形式表示相关系数
dev.off()


library(matrixStats)
library(dplyr)
library(tidyr)
library(ggplot2)

# CV ----
cv_data <- data_fill
## ----- 0) 可选：删除异常样本列 -----
## 若 data 是 log2 值，保留你这句：转回原尺度再算CV
# stastic_lin <- as.matrix(2^cv_data)           # 如果 data 是 log2 矩阵
# drop_samples <- c("EV_R3", "OE_WT_R3")     # 想删谁就写谁，可以多个
# stastic_lin <- stastic_lin[, !colnames(stastic_lin) %in% drop_samples]
stastic_lin <- cv_data
## ----- 1) 动态匹配各组列（不再用固定3列的索引）-----
SampleType <- c("AML")

group_cols <- setNames(
  lapply(SampleType, function(g){
    grep(paste0("^", g, "_\\d+$"), colnames(stastic_lin))
  }),
  SampleType
)
group_cols
## 可选：把没有任何样本列的组过滤掉
group_cols <- group_cols[sapply(group_cols, length) > 0]

## ----- 2) 计算每组 CV（行SD/行均值*100），<2个样本则返回NA -----
CV_list <- lapply(names(group_cols), function(g) {
  cols <- group_cols[[g]]
  if (length(cols) < 2) {
    data.frame(cv = rep(NA_real_, nrow(stastic_lin)))
  } else {
    m <- as.matrix(stastic_lin[, cols, drop = FALSE])  # 关键：转 matrix
    data.frame(
      cv = rowSds(m, na.rm = TRUE) * 100 / rowMeans(m, na.rm = TRUE)
    )
  }
})
names(CV_list) <- names(group_cols)

CV1 <- do.call(cbind, CV_list)
colnames(CV1) <- names(group_cols)
CV1$SYMBOL <- rownames(stastic_lin)

## ----- 3) 长表用于作图 -----
data1 <- CV1 |>
  pivot_longer(
    cols = all_of(names(group_cols)),
    names_to = "group",
    values_to = "CV"
  ) |>
  mutate(group = factor(group, levels = SampleType[SampleType %in% names(group_cols)]))

## ----- 4) 画图（x 轴和中位数标签都不再硬编码）-----
# mycol <- c("#82A1A7","#5C7A86","#D1BB81","#D0BDBD","#B16D6E","#F0C8B8")
mycol <- c("#DC00007F")
levels(data1$group)
# 控制样本排列顺序
data1$group <- factor(data1$group, levels = c("AML"))

dir_cv <- "./03_result/01_QC/"
if (!dir.exists(dir_cv)) dir.create(dir_cv, recursive = TRUE)

p <- ggplot(data1, aes(x = group, y = CV, fill = group)) +
  geom_violin(alpha = 0.5, show.legend = FALSE) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = colorRampPalette(mycol)(length(levels(data1$group)))) +
  labs(y = "coefficient of variation(%)", x = "", title = "") +
  theme_classic() +
  theme(legend.text = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1))

# 计算每组的中位数，用于标注
med_tbl <- data1 |>
  group_by(group) |>
  summarise(med = median(CV, na.rm = TRUE), .groups = "drop")

# 取一个合理的标注高度（例如上四分位数的 1.1 倍）
y_annot <- quantile(data1$CV, 0.9, na.rm = TRUE)

p <- p + geom_text(data = med_tbl,
                   aes(x = group, y = y_annot, label = sprintf("%.2f", med)),
                   inherit.aes = FALSE)
print(p)
ggsave(file.path(dir_cv, "CV.pdf"), plot = p, width = 2, height = 5)
dev.off()

