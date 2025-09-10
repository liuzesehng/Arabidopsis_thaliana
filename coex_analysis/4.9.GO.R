library(clusterProfiler)
library(org.At.tair.db)  # 拟南芥物种数据库
library(ggplot2)
library(dplyr)
library(stringr)

# 定义一个函数来执行GO富集分析并保存结果
perform_enrichment <- function(genes, output_prefix) {
  ontologies <- c("BP", "MF", "CC")
  for (ont in ontologies) {
    enrich.go <- enrichGO(gene = genes,
                          OrgDb = org.At.tair.db,
                          keyType = "TAIR",  # 使用基因符号
                          ont = ont,  # BP, MF, 或 CC
                          pAdjustMethod = 'fdr',  # 指定 p 值校正方法
                          pvalueCutoff = 0.05,  # 指定 p 值阈值（可指定 1 以输出全部）
                          qvalueCutoff = 0.2,  # 指定 q 值阈值（可指定 1 以输出全部）
                          readable = TRUE)  # 是否将基因 ID 转换为基因名
        # 检查富集结果是否为空
    if (!is.null(enrich.go) && nrow(enrich.go@result) > 0) {
      # 1. 保存原始的、未简化的富集结果
      write.table(enrich.go, paste0(output_prefix, ".", ont, ".enrich.go.original.txt"), sep = '\t', row.names = FALSE, quote = FALSE)
      # 2. 使用simplify()函数去除冗余条目
      # cutoff=0.7: 语义相似性得分大于0.7的条目被认为是冗余的
      # by="p.adjust": 当多个条目冗余时，保留p.adjust值最小的那个
      enrich.go.simplified <- simplify(enrich.go, cutoff = 0.7, by = "p.adjust", select_fun = min)
      # 3. 保存简化后的富集结果
      write.table(enrich.go.simplified, paste0(output_prefix, ".", ont, ".enrich.go.simplified.txt"), sep = '\t', row.names = FALSE, quote = FALSE)
      
      # 4. 绘制气泡图
      if (nrow(enrich.go.simplified@result) > 0) {
        create_bubble_plot(enrich.go.simplified, output_prefix, ont)
      }
      
      cat(paste("分析完成:", ont, "，原始结果和简化结果已保存。\n"))
    } else {
      cat(paste("注意:", ont, "没有富集到任何结果。\n"))
    }
  }
}

# 定义气泡图绘制函数
create_bubble_plot <- function(enrich_result, output_prefix, ontology) {
  # 提取前20个最显著的结果
  plot_data <- enrich_result@result %>%
    arrange(p.adjust) %>%
    head(20) %>%
    mutate(
      # 截断长的Description以便显示
      Description = str_wrap(Description, width = 50),
      # 计算富集比例
      GeneRatio_numeric = sapply(GeneRatio, function(x) {
        nums <- as.numeric(unlist(strsplit(x, "/")))
        return(nums[1]/nums[2])
      }),
      # 负log10转换p.adjust值用于气泡颜色
      neg_log10_padj = -log10(p.adjust)
    )
  
  # 创建气泡图
  p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = reorder(Description, GeneRatio_numeric))) +
    geom_point(aes(size = Count, color = neg_log10_padj), alpha = 0.7) +
    scale_size_continuous(range = c(3, 12), name = "Gene Count") +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adjust)") +
    labs(
      x = "Gene Ratio",
      y = "GO Terms"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold")
    ) +
    guides(
      size = guide_legend(override.aes = list(alpha = 1), title = "Gene Count"),
      color = guide_colorbar(title = "-log10(p.adjust)")
    )
  
  # 保存图片 - 仅保存PDF格式
  output_file <- paste0(gsub("\\.txt$", "", output_prefix), ".", ontology, ".bubble_plot.pdf")
  ggsave(output_file, plot = p, width = 12, height = 9, dpi = 1200)
  
  cat(paste("气泡图已保存:", output_file, "\n"))
}


# 读取基因列表并执行富集分析
# all
genes <- read.delim('total/all.gene_id.txt', header = TRUE, stringsAsFactors = FALSE)[[1]]
perform_enrichment(genes, 'total/all.gene_id')

# a
genes <- read.delim('a/a.gene_id.txt', header = TRUE, stringsAsFactors = FALSE)[[1]]
perform_enrichment(genes, 'a/a.gene_id')

# b
genes <- read.delim('b/b.gene_id.txt', header = TRUE, stringsAsFactors = FALSE)[[1]]
perform_enrichment(genes, 'b/b.gene_id')

