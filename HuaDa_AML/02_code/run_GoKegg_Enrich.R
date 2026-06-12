run_GoKegg_Enrich <- function(Genelist,
                              Output_dir,
                              pvalue = 0.05,
                              logfc = 1,
                              go_db = "org.Hs.eg.db",
                              kegg_db = 'hsa') {
  library(clusterProfiler)
  library(ggplot2)
  library(openxlsx)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(KEGG.db)
  
  # 检查列名 ----
  required_cols <- c("gene", "logFC", "P.Value")
  if (!all(required_cols %in% colnames(Genelist))) {
    stop("Input Genelist must contain columns: gene, logFC, P.Value")
  }
  
  # 创建输出目录 ----
  dir.create(Output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 按上下调分组 ----
  message("Filtering gene sets...")
  up_genes   <- subset(Genelist, logFC >  logfc & P.Value < pvalue)
  down_genes <- subset(Genelist, logFC < -logfc & P.Value < pvalue)
  all_genes  <- subset(Genelist, abs(logFC) > logfc & P.Value < pvalue)
  
  # 内部函数：富集分析模块 ----
  do_enrichment <- function(df, tag) {
    if (nrow(df) == 0) {
      warning(paste0("No genes passed filter for ", tag))
      return(NULL)
    }
    
    gene_entrez <- tryCatch({
      bitr(df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = get(go_db))
    }, error = function(e) return(NULL))
    
    if (is.null(gene_entrez) || nrow(gene_entrez) == 0) {
      warning(paste0("No ENTREZID found for ", tag))
      return(NULL)
    }
    
    enrich_go <- enrichGO(gene = gene_entrez$ENTREZID,
                          OrgDb = get(go_db),
                          keyType = "ENTREZID",
                          ont = "ALL",
                          pvalueCutoff = 0.5,
                          qvalueCutoff = 1,
                          readable = TRUE)
    
    enrich_kegg <- enrichKEGG(gene = gene_entrez$ENTREZID,
                              organism = kegg_db,
                              keyType = "kegg",
                              pvalueCutoff = 1,
                              qvalueCutoff = 1,
                              use_internal_data = TRUE)
    enrich_kegg <- setReadable(enrich_kegg, OrgDb = get(go_db), keyType = "ENTREZID")
    
    # 写出 GO 结果
    write.xlsx(enrich_go@result,
               file = file.path(Output_dir, paste0("GO_", tag, ".xlsx")),
               rowNames = FALSE)
    pdf(file = file.path(Output_dir, paste0("GO_", tag, ".pdf")),
        width = 5.5, height = 7)
    print(dotplot(enrich_go, showCategory = 5, color = "pvalue", split = "ONTOLOGY") +
            facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free'))
    dev.off()
    
    # 写出 KEGG 结果
    write.xlsx(enrich_kegg@result,
               file = file.path(Output_dir, paste0("KEGG_", tag, ".xlsx")),
               rowNames = FALSE)
    pdf(file = file.path(Output_dir, paste0("KEGG_", tag, ".pdf")),
        width = 6, height = 5)
    print(dotplot(enrich_kegg, showCategory = 10, color = "pvalue"))
    dev.off()
    
    return(list(GO = enrich_go,
                KEGG = enrich_kegg,
                ENTREZ = gene_entrez$ENTREZID))
  }
  
  # 跑三组富集 ----
  message("Running GO/KEGG for downregulated genes...")
  result_down <- do_enrichment(down_genes, "down")
  
  message("Running GO/KEGG for upregulated genes...")
  result_up <- do_enrichment(up_genes, "up")
  
  message("Running GO/KEGG for all DE genes...")
  result_all <- do_enrichment(all_genes, "all_DE")
  
  # 报告汇总函数 ----
  summarize_res <- function(res) {
    if (is.null(res)) return(NULL)
    go_res   <- res$GO@result
    kegg_res <- res$KEGG@result
    
    list(
      n_gene = length(unique(res$ENTREZ)),
      go_count = nrow(go_res),
      kegg_count = nrow(kegg_res),
      go_sig = subset(go_res, pvalue < 0.05),
      kegg_sig = subset(kegg_res, pvalue < 0.05)
    )
  }
  
  s_down <- summarize_res(result_down)
  s_up   <- summarize_res(result_up)
  s_all  <- summarize_res(result_all)
  
  # 统一生成整合报告 ----
  report_file <- file.path(Output_dir, "GO_KEGG_report.txt")
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  sink(report_file)
  cat("============================================\n")
  cat("        GO/KEGG 富集分析报告（整合版）\n")
  cat("============================================\n")
  cat("分析时间: ", timestamp, "\n\n")
  
  cat("【logFC阈值】\n")
  cat("logFC阈值: ", logfc, "\n\n")
  cat("【pvalue阈值】\n")
  cat("pvalue阈值: ", pvalue, "\n\n")
  
  write_section <- function(s, tag_cn) {
    cat("============================================\n")
    cat("【", tag_cn, "】\n\n", sep = "")
    
    if (is.null(s)) {
      cat("无可用结果（无基因通过阈值或无法映射ENTREZID）。\n\n")
      return()
    }
    
    cat("【总基因数量】\n")
    cat("总基因数: ", s$n_gene, "\n\n")
    
    cat("【GO富集结果】\n")
    cat("总通路数: ", s$go_count, "\n")
    cat("显著通路 (p.value < 0.05): ", nrow(s$go_sig), "\n")
    if (nrow(s$go_sig) > 0) {
      cat("\n显著GO通路前5:\n")
      print(head(s$go_sig[, c("ONTOLOGY","ID","Description","pvalue")], 5))
    }
    cat("\n--------------------------------------------\n")
    
    cat("【KEGG富集结果】\n")
    cat("总通路数: ", s$kegg_count, "\n")
    cat("显著通路 (p.value < 0.05): ", nrow(s$kegg_sig), "\n")
    if (nrow(s$kegg_sig) > 0) {
      cat("\n显著KEGG通路前5:\n")
      print(head(s$kegg_sig[, c("ID","Description","pvalue")], 5))
    }
    cat("\n\n")
  }
  
  write_section(s_down, "下调基因 (down)")
  write_section(s_up,   "上调基因 (up)")
  write_section(s_all,  "全部差异基因 (all_DE)")
  
  cat("输出文件示例:\n")
  cat("- GO 富集表: GO_down.xlsx / GO_up.xlsx / GO_all_DE.xlsx\n")
  cat("- GO 富集图: GO_down.pdf  / GO_up.pdf  / GO_all_DE.pdf\n")
  cat("- KEGG 富集表: KEGG_down.xlsx / KEGG_up.xlsx / KEGG_all_DE.xlsx\n")
  cat("- KEGG 富集图: KEGG_down.pdf  / KEGG_up.pdf  / KEGG_all_DE.pdf\n\n")
  cat("报告文件: ", report_file, "\n")
  cat("============================================\n")
  sink()
  
  return(list(
    down = result_down,
    up   = result_up,
    all  = result_all
  ))
}
