library(readr)
library(dplyr)
library(openxlsx)
library(limma)
library(org.Hs.eg.db) # 假设是人类数据

# 1. 处理蛋白组基因名
# alias2SymbolTable 会查找别名并返回官方 Symbol
lsc_df$official_gene <- alias2SymbolTable(lsc_df$lsc17_genes, species = "Hs")

# 2. 处理转录组基因名
prot_mat$official_gene <- alias2SymbolTable(rownames(prot_mat), species = "Hs")

# 3. 再次取交集（记得处理转换后可能产生的 NA）
common_genes0 <- intersect(lsc_df$lsc17_genes, rownames(prot_mat))
common_genes <- intersect(lsc_df$official_gene, prot_mat$official_gene)
diff_genes <- setdiff(official_gene1,official_gene2)

# 没有官方名的NA转为原用名
DE1$official_gene <- ifelse(is.na(DE1$official_gene), DE1$gene, DE1$official_gene)
DE2$official_gene <- ifelse(is.na(DE2$official_gene), DE2$gene, DE2$official_gene)

