library(clusterProfiler)
library(org.At.tair.db)  # 拟南芥物种数据库
library(ggplot2)
library(dplyr)
library(stringr)
library(grid)

# 定义一个函数来执行GO富集分析并保存结果
perform_enrichment <- function(genes, output_prefix) {
  ontologies <- c("BP", "MF", "CC")
  plot_results <- list()
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
      
      # 4. 收集绘图数据，后续统一绘制合并图
      if (nrow(enrich.go.simplified@result) > 0) {
        plot_results[[ont]] <- enrich.go.simplified@result
      }
      
      cat(paste("分析完成:", ont, "，原始结果和简化结果已保存。\n"))
    } else {
      cat(paste("注意:", ont, "没有富集到任何结果。\n"))
    }
  }

  create_combined_bubble_plot(plot_results, output_prefix)
}

# 定义合并气泡图绘制函数
create_combined_bubble_plot <- function(enrich_results, output_prefix) {
  top_n_terms <- 10

  if (length(enrich_results) == 0) {
    cat(paste("注意:", output_prefix, "没有可用于绘图的GO富集结果。\n"))
    return(invisible(NULL))
  }

  plot_data <- bind_rows(lapply(names(enrich_results), function(ontology) {
    enrich_results[[ontology]] %>%
      arrange(p.adjust) %>%
      head(top_n_terms) %>%
      mutate(Ontology = ontology)
  }))

  if (nrow(plot_data) == 0) {
    cat(paste("注意:", output_prefix, "没有可用于绘图的GO富集结果。\n"))
    return(invisible(NULL))
  }

  plot_data <- plot_data %>%
    mutate(
      # 截断长的Description以便显示
      Description = str_wrap(Description, width = 34),
      # 计算富集比例
      GeneRatio_numeric = sapply(GeneRatio, function(x) {
        nums <- as.numeric(unlist(strsplit(x, "/")))
        return(nums[1]/nums[2])
      }),
      # 负log10转换p.adjust值用于气泡颜色
      neg_log10_padj = -log10(p.adjust),
      Ontology = factor(Ontology, levels = c("BP", "MF", "CC"))
    )

  plot_data <- plot_data %>%
    group_by(Ontology) %>%
    arrange(GeneRatio_numeric, .by_group = TRUE) %>%
    mutate(
      Description_plot = factor(
        paste(Description, Ontology, sep = "__"),
        levels = paste(Description, Ontology, sep = "__")
      )
    ) %>%
    ungroup()
  
  label_map <- setNames(as.character(plot_data$Description), as.character(plot_data$Description_plot))

  p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description_plot)) +
    geom_point(aes(size = Count, color = neg_log10_padj), alpha = 0.7) +
    scale_size_continuous(range = c(3, 12), name = "Gene Count") +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adjust)") +
    scale_y_discrete(labels = label_map) +
    facet_grid(. ~ Ontology, scales = "free_y", space = "free_y") +
    labs(
      x = "Gene Ratio",
      y = paste0("Top ", top_n_terms, " GO Terms Per Ontology")
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 13, face = "bold", lineheight = 1.08, margin = margin(r = 12)),
      axis.text.x = element_text(size = 14, angle = 30, hjust = 1, vjust = 1, face = "bold", margin = margin(t = 10)),
      axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 12)),
      axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 14)),
      strip.text = element_text(size = 16, face = "bold"),
      panel.spacing.x = unit(2, "lines"),
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 13, face = "bold"),
      plot.margin = margin(t = 18, r = 24, b = 22, l = 28)
    ) +
    guides(
      size = guide_legend(override.aes = list(alpha = 1), title = "Gene Count"),
      color = guide_colorbar(title = "-log10(p.adjust)")
    )
  
  output_file <- paste0(gsub("\\.txt$", "", output_prefix), ".GO.bubble_plot.pdf")
  ggsave(output_file, plot = p, width = 22, height = 11, dpi = 1200)
  
  cat(paste("GO合并气泡图已保存:", output_file, "\n"))
}


# 读取基因列表并执行富集分析
# all
genes <- read.delim('total/yellow_module_genes.txt', header = TRUE, stringsAsFactors = FALSE)[[1]]
perform_enrichment(genes, 'total/total.gene_id')
genes_module <- read.delim('total/TF_in_module.txt', header = TRUE, stringsAsFactors = FALSE)[[2]]
perform_enrichment(genes_module, 'total/total.gene_id_TF')
genes_negative <- read.delim('total/negative_grey_module_genes.txt', header = TRUE, stringsAsFactors = FALSE)[[1]]
perform_enrichment(genes_negative, 'total/total.gene_id_negative')
genes_negative_module <- read.delim('total/TF_negative_in_module.txt', header = TRUE, stringsAsFactors = FALSE)[[2]]
perform_enrichment(genes_negative_module, 'total/total.gene_id_negative_TF')

# 10C
genes <- read.delim('10C/10C.greenyellow_module_genes.txt', header = TRUE, stringsAsFactors = FALSE)[[1]]
perform_enrichment(genes, '10C/10C.gene_id')
genes_module <- read.delim('10C/TF_in_module.txt', header = TRUE, stringsAsFactors = FALSE)[[2]]
perform_enrichment(genes_module, '10C/10C.gene_id_TF')
genes_negative <- read.delim('10C/10C.negative_green_module_genes.txt', header = TRUE, stringsAsFactors = FALSE)[[1]]
perform_enrichment(genes_negative, '10C/10C.gene_id_negative')
genes_negative_module <- read.delim('10C/TF_negative_in_module.txt', header = TRUE, stringsAsFactors = FALSE)[[2]]
perform_enrichment(genes_negative_module, '10C/10C.gene_id_negative_TF')

# 16C
genes <- read.delim('16C/16C.pink_module_genes.txt', header = TRUE, stringsAsFactors = FALSE)[[1]]
perform_enrichment(genes, '16C/16C.gene_id')
genes_module <- read.delim('16C/TF_in_module.txt', header = TRUE, stringsAsFactors = FALSE)[[2]]
perform_enrichment(genes_module, '16C/16C.gene_id_TF')
genes_negative <- read.delim('16C/16C.negative_salmon_module_genes.txt', header = TRUE, stringsAsFactors = FALSE)[[1]]
perform_enrichment(genes_negative, '16C/16C.gene_id_negative')
# genes_negative_module <- read.delim('16C/TF_negative_in_module.txt', header = TRUE, stringsAsFactors = FALSE)[[2]]
# perform_enrichment(genes_negative_module, '16C/16C.gene_id_negative_TF')

# 22C
genes <- read.delim('22C/22C.red_module_genes.txt', header = TRUE, stringsAsFactors = FALSE)[[1]]
perform_enrichment(genes, '22C/22C.gene_id')
genes_module <- read.delim('22C/TF_in_module.txt', header = TRUE, stringsAsFactors = FALSE)[[2]]
perform_enrichment(genes_module, '22C/22C.gene_id_TF')
genes_negative <- read.delim('22C/22C.negative_magenta_module_genes.txt', header = TRUE, stringsAsFactors = FALSE)[[1]]
perform_enrichment(genes_negative, '22C/22C.gene_id_negative')
genes_negative_module <- read.delim('22C/TF_negative_in_module.txt', header = TRUE, stringsAsFactors = FALSE)[[2]]
perform_enrichment(genes_negative_module, '22C/22C.gene_id_negative_TF')
