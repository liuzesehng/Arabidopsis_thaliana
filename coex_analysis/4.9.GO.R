library(clusterProfiler)
library(org.At.tair.db)  # 拟南芥物种数据库

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
      cat(paste("分析完成:", ont, "，原始结果和简化结果已保存。\n"))
    } else {
      cat(paste("注意:", ont, "没有富集到任何结果。\n"))
    }
  }
}


# 读取基因列表并执行富集分析
genes <- read.delim('total/all.gene_id.txt', header = TRUE, stringsAsFactors = FALSE)[[1]]
perform_enrichment(genes, 'total/all.gene_id.txt')

genes <- read.delim('total/gene.TF_ids.txt', header = TRUE, stringsAsFactors = FALSE)[[1]]
perform_enrichment(genes, 'total/gene.TF_ids')
library(ggplot2)
library(dplyr)

# 定义一个函数来读取和处理富集结果
read_and_process_results <- function(prefix, ont) {
  read.delim(paste0(prefix, ".", ont, ".enrich.go.txt"), header = TRUE, stringsAsFactors = FALSE)
}

# 读取和处理BP结果
enrich.go.10C.BP <- read_and_process_results('10C/10C.nonTF_ids', 'BP')
enrich.go.16C.BP <- read_and_process_results('16C/16C.nonTF_ids', 'BP')
enrich.go.22C.a.BP <- read_and_process_results('22C/22C.a.nonTF_ids', 'BP')
enrich.go.22C.b.BP <- read_and_process_results('22C/22C.b.nonTF_ids', 'BP')
enrich.go.mid.a.BP <- read_and_process_results('mid/mid.a.nonTF_ids', 'BP')
enrich.go.mid.b.BP <- read_and_process_results('mid/mid.b.nonTF_ids', 'BP')
enrich.go.high.a.BP <- read_and_process_results('high/high.a.nonTF_ids', 'BP')
enrich.go.high.b.BP <- read_and_process_results('high/high.b.nonTF_ids', 'BP')

enrich.go.10C.TF.BP <- read_and_process_results('10C/10C.TF_ids', 'BP')
enrich.go.16C.TF.BP <- read_and_process_results('16C/16C.TF_ids', 'BP')
enrich.go.22C.a.TF.BP <- read_and_process_results('22C/22C.a.TF_ids', 'BP')
enrich.go.22C.b.TF.BP <- read_and_process_results('22C/22C.b.TF_ids', 'BP')
enrich.go.mid.a.TF.BP <- read_and_process_results('mid/mid.a.TF_ids', 'BP')
enrich.go.mid.b.TF.BP <- read_and_process_results('mid/mid.b.TF_ids', 'BP')
enrich.go.high.a.TF.BP <- read_and_process_results('high/high.a.TF_ids', 'BP')
enrich.go.high.b.TF.BP <- read_and_process_results('high/high.b.TF_ids', 'BP')

# 提取前 10 个富集结果
top10_10C.BP.results <- enrich.go.10C.BP[1:10, ]
top10_16C.BP.results <- enrich.go.16C.BP[1:10, ]
top10_22C.a.BP.results <- enrich.go.22C.a.BP[1:10, ]
top10_22C.b.BP.results <- enrich.go.22C.b.BP[1:10, ]
top10_mid.a.BP.results <- enrich.go.mid.a.BP[1:10, ]
top10_mid.b.BP.results <- enrich.go.mid.b.BP[1:10, ]
top10_high.a.BP.results <- enrich.go.high.a.BP[1:10, ]
top10_high.b.BP.results <- enrich.go.high.b.BP[1:10, ]

top10_10C.TF.BP.results <- enrich.go.10C.TF.BP[1:10, ]
top10_16C.TF.BP.results <- enrich.go.16C.TF.BP[1:10, ]
top10_22C.a.TF.BP.results <- enrich.go.22C.a.TF.BP[1:10, ]
top10_22C.b.TF.BP.results <- enrich.go.22C.b.TF.BP[1:10, ]
top10_mid.a.TF.BP.results <- enrich.go.mid.a.TF.BP[1:10, ]
top10_mid.b.TF.BP.results <- enrich.go.mid.b.TF.BP[1:10, ]
top10_high.a.TF.BP.results <- enrich.go.high.a.TF.BP[1:10, ]
top10_high.b.TF.BP.results <- enrich.go.high.b.TF.BP[1:10, ]

# 假设有多个样本的富集结果 (list of enrichGO results)
sample_results_BP_nonTF <- list(
  "10" = top10_10C.BP.results,  # enrichGO 的结果
  "16" = top10_16C.BP.results,
  "22_a" = top10_22C.a.BP.results,
  "22_b" = top10_22C.b.BP.results,
  "mid_a" = top10_mid.a.BP.results,
  "mid_b" = top10_mid.b.BP.results,
  "high_a" = top10_high.a.BP.results,
  "high_b" = top10_high.b.BP.results
)

sample_results_BP_TF <- list(
  "10" = top10_10C.TF.BP.results,  # enrichGO 的结果
  "16" = top10_16C.TF.BP.results,
  "22_a" = top10_22C.a.TF.BP.results,
  "22_b" = top10_22C.b.TF.BP.results,
  "mid_a" = top10_mid.a.TF.BP.results,
  "mid_b" = top10_mid.b.TF.BP.results,
  "high_a" = top10_high.a.TF.BP.results,
  "high_b" = top10_high.b.TF.BP.results
)


# 将多个样本的富集结果合并到一个数据框
combined_results_BP_nonTF <- do.call(rbind, lapply(names(sample_results_BP_nonTF), function(sample) {
  result <- sample_results_BP_nonTF[[sample]]  # 提取 enrichGO 结果
  result$Sample <- sample  # 添加样本名称
  return(result)
}))

combined_results_BP_TF <- do.call(rbind, lapply(names(sample_results_BP_TF), function(sample) {
  result <- sample_results_BP_TF[[sample]]  # 提取 enrichGO 结果
  result$Sample <- sample  # 添加样本名称
  return(result)
}))

# 将 p.adjust 转换为 -log10(p.adjust)
combined_results_BP_nonTF <- combined_results_BP_nonTF %>%
  mutate(log_padj = -log10(p.adjust))  # 将 p.adjust 转换为 -log10(p.adjust)

combined_results_BP_TF <- combined_results_BP_TF %>%
  mutate(log_padj = -log10(p.adjust))  # 将 p.adjust 转换为 -log10(p.adjust)

# 将 GeneRatio 转换为数值类型
combined_results_BP_nonTF$GeneRatio <- sapply(strsplit(combined_results_BP_nonTF$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
combined_results_BP_TF$GeneRatio <- sapply(strsplit(combined_results_BP_TF$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))

# 删除含有 NA 的行
combined_results_BP_nonTF <- na.omit(combined_results_BP_nonTF)
combined_results_BP_TF <- na.omit(combined_results_BP_TF)

# 替换 NA 为默认值
combined_results_BP_nonTF <- combined_results_BP_nonTF %>%
  mutate(
    log_padj = ifelse(is.na(log_padj), 0, log_padj),  # 替换 log_padj 的 NA 为 0
    GeneRatio = ifelse(is.na(GeneRatio), 0, GeneRatio),  # 替换 GeneRatio 的 NA 为 0
    #Description = ifelse(is.na(Description), 0, Description)  # 替换 Description 的 NA 为 0
  )

combined_results_BP_TF <- combined_results_BP_TF %>%
  mutate(
    log_padj = ifelse(is.na(log_padj), 0, log_padj),  # 替换 log_padj 的 NA 为 0
    GeneRatio = ifelse(is.na(GeneRatio), 0, GeneRatio),  # 替换 GeneRatio 的 NA 为 0
    #Description = ifelse(is.na(Description), 0, Description)  # 替换 Description 的 NA 为 0
  )

# 绘制BP点图
ggplot(combined_results_BP_nonTF, aes(x = Sample, y = Description)) +
  geom_point(aes(size = GeneRatio, color = log_padj)) +  # 根据 count 调整点的大小，根据 log_padj 调整颜色
  scale_color_gradient(low = "blue", high = "red") +  # 自定义颜色渐变
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), #文本加粗
    axis.text.y = element_text(size = 10, face = "bold"),  # y 轴文本加粗
    legend.title = element_text(face = "bold"),  # 图例标题加粗
    legend.text = element_text(face = "bold")  # 图例文本加粗
  ) +
  labs(
    x = "",
    y = "",
    size = "GeneRatio",  # 表示点大小的注释
    color = "-log10(p.adjust)"  # 表示颜色的注释
  )


ggsave("GO_Enrichment_BP_nonTF_Across_Samples.pdf", width = 10, height = 8)


ggplot(combined_results_BP_TF, aes(x = Sample, y = Description)) +
  geom_point(aes(size = GeneRatio, color = log_padj)) +  # 根据 count 调整点的大小，根据 log_padj 调整颜色
  scale_color_gradient(low = "blue", high = "red") +  # 自定义颜色渐变
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), #文本加粗
    axis.text.y = element_text(size = 10, face = "bold"),  # y 轴文本加粗
    legend.title = element_text(face = "bold"),  # 图例标题加粗
    legend.text = element_text(face = "bold")  # 图例文本加粗
  ) +
  labs(
    x = "",
    y = "",
    size = "GeneRatio",  # 表示点大小的注释
    color = "-log10(p.adjust)"  # 表示颜色的注释
  )

ggsave("GO_Enrichment_BP_TF_Across_Samples.pdf", width = 10, height = 8)

# 读取和处理MF结果
enrich.go.10C.MF <- read_and_process_results('10C/10C.nonTF_ids', 'MF')
enrich.go.16C.MF <- read_and_process_results('16C/16C.nonTF_ids', 'MF')
enrich.go.22C.a.MF <- read_and_process_results('22C/22C.a.nonTF_ids', 'MF')
enrich.go.22C.b.MF <- read_and_process_results('22C/22C.b.nonTF_ids', 'MF')
enrich.go.mid.a.MF <- read_and_process_results('mid/mid.a.nonTF_ids', 'MF')
enrich.go.mid.b.MF <- read_and_process_results('mid/mid.b.nonTF_ids', 'MF')
enrich.go.high.a.MF <- read_and_process_results('high/high.a.nonTF_ids', 'MF')
enrich.go.high.b.MF <- read_and_process_results('high/high.b.nonTF_ids', 'MF')

enrich.go.10C.TF.MF <- read_and_process_results('10C/10C.TF_ids', 'MF')
enrich.go.16C.TF.MF <- read_and_process_results('16C/16C.TF_ids', 'MF')
enrich.go.22C.a.TF.MF <- read_and_process_results('22C/22C.a.TF_ids', 'MF')
enrich.go.22C.b.TF.MF <- read_and_process_results('22C/22C.b.TF_ids', 'MF')
enrich.go.mid.a.TF.MF <- read_and_process_results('mid/mid.a.TF_ids', 'MF')
enrich.go.mid.b.TF.MF <- read_and_process_results('mid/mid.b.TF_ids', 'MF')
enrich.go.high.a.TF.MF <- read_and_process_results('high/high.a.TF_ids', 'MF')
enrich.go.high.b.TF.MF <- read_and_process_results('high/high.b.TF_ids', 'MF')

# 提取前 10 个富集结果
top10_10C.MF.results <- enrich.go.10C.MF[1:20, ]
top10_16C.MF.results <- enrich.go.16C.MF[120, ]
top10_22C.a.MF.results <- enrich.go.22C.a.MF[1:20, ]
top10_22C.b.MF.results <- enrich.go.22C.b.MF[1:20, ]
top10_mid.a.MF.results <- enrich.go.mid.a.MF[1:20, ]
top10_mid.b.MF.results <- enrich.go.mid.b.MF[1:20, ]
top10_high.a.MF.results <- enrich.go.high.a.MF[1:20, ]
top10_high.b.MF.results <- enrich.go.high.b.MF[1:20, ]

top10_10C.TF.MF.results <- enrich.go.10C.TF.MF[1:20, ]
top10_16C.TF.MF.results <- enrich.go.16C.TF.MF[1:20, ]
top10_22C.a.TF.MF.results <- enrich.go.22C.a.TF.MF[1:20, ]
top10_22C.b.TF.MF.results <- enrich.go.22C.b.TF.MF[1:20, ]
top10_mid.a.TF.MF.results <- enrich.go.mid.a.TF.MF[1:20, ]
top10_mid.b.TF.MF.results <- enrich.go.mid.b.TF.MF[1:20, ]
top10_high.a.TF.MF.results <- enrich.go.high.a.TF.MF[1:20, ]
top10_high.b.TF.MF.results <- enrich.go.high.b.TF.MF[1:20, ]

# 假设有多个样本的富集结果 (list of enrichGO results)
sample_results_MF_nonTF <- list(
  "10" = top10_10C.MF.results,  # enrichGO 的结果
  "16" = top10_16C.MF.results,
  "22_a" = top10_22C.a.MF.results,
  "22_b" = top10_22C.b.MF.results,
  "mid_a" = top10_mid.a.MF.results,
  "mid_b" = top10_mid.b.MF.results,
  "high_a" = top10_high.a.MF.results,
  "high_b" = top10_high.b.MF.results
)

sample_results_MF_TF <- list(
  "10" = top10_10C.TF.MF.results,  # enrichGO 的结果
  "16" = top10_16C.TF.MF.results,
  "22_a" = top10_22C.a.TF.MF.results,
  "22_b" = top10_22C.b.TF.MF.results,
  "mid_a" = top10_mid.a.TF.MF.results,
  "mid_b" = top10_mid.b.TF.MF.results,
  "high_a" = top10_high.a.TF.MF.results,
  "high_b" = top10_high.b.TF.MF.results
)

# 将多个样本的富集结果合并到一个数据框
combined_results_MF_nonTF <- do.call(rbind, lapply(names(sample_results_MF_nonTF), function(sample) {
  result <- sample_results_MF_nonTF[[sample]]  # 提取 enrichGO 结果
  result$Sample <- sample  # 添加样本名称
  return(result)
}))

combined_results_MF_TF <- do.call(rbind, lapply(names(sample_results_MF_TF), function(sample) {
  result <- sample_results_MF_TF[[sample]]  # 提取 enrichGO 结果
  result$Sample <- sample  # 添加样本名称
  return(result)
}))

combined_results_MF_nonTF <- combined_results_MF_nonTF %>%
  mutate(log_padj = -log10(p.adjust))  # 将 p.adjust 转换为 -log10(p.adjust)

combined_results_MF_TF <- combined_results_MF_TF %>%
  mutate(log_padj = -log10(p.adjust))  # 将 p.adjust 转换为 -log10(p.adjust)

# 将 GeneRatio 转换为数值类型
combined_results_MF_nonTF$GeneRatio <- sapply(strsplit(combined_results_MF_nonTF$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
combined_results_MF_TF$GeneRatio <- sapply(strsplit(combined_results_MF_TF$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))

# 删除含有 NA 的行
combined_results_MF_nonTF <- na.omit(combined_results_MF_nonTF)
combined_results_MF_TF <- na.omit(combined_results_MF_TF)

# 绘制MF点图
ggplot(combined_results_MF_nonTF, aes(x = Sample, y = Description)) +
  geom_point(aes(size = GeneRatio, color = log_padj)) +  # 根据 count 调整点的大小，根据 log_padj 调整颜色
  scale_color_gradient(low = "blue", high = "red") +  # 自定义颜色渐变
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), #文本加粗
    axis.text.y = element_text(size = 10, face = "bold"),  # y 轴文本加粗
    legend.title = element_text(face = "bold"),  # 图例标题加粗
    legend.text = element_text(face = "bold")  # 图例文本加粗
  ) +
  labs(
    x = "",
    y = "",
    size = "GeneRatio",  # 表示点大小的注释
    color = "-log10(p.adjust)"  # 表示颜色的注释
  )

ggsave("GO_Enrichment_MF_nonTF_Across_Samples.pdf", width = 12, height = 10)

ggplot(combined_results_MF_TF, aes(x = Sample, y = Description)) +
  geom_point(aes(size = GeneRatio, color = log_padj)) +  # 根据 count 调整点的大小，根据 log_padj 调整颜色
  scale_color_gradient(low = "blue", high = "red") +  # 自定义颜色渐变
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), #文本加粗
    axis.text.y = element_text(size = 10, face = "bold"),  # y 轴文本加粗
    legend.title = element_text(face = "bold"),  # 图例标题加粗
    legend.text = element_text(face = "bold")  # 图例文本加粗
  ) +
  labs(
    x = "",
    y = "",
    size = "GeneRatio",  # 表示点大小的注释
    color = "-log10(p.adjust)"  # 表示颜色的注释
  )

ggsave("GO_Enrichment_MF_TF_Across_Samples.pdf", width = 10, height = 8)

# 读取和处理CC结果
enrich.go.10C.CC <- read_and_process_results('10C/10C.nonTF_ids', 'CC')
enrich.go.16C.CC <- read_and_process_results('16C/16C.nonTF_ids', 'CC')
enrich.go.22C.a.CC <- read_and_process_results('22C/22C.a.nonTF_ids', 'CC')
enrich.go.22C.b.CC <- read_and_process_results('22C/22C.b.nonTF_ids', 'CC')
enrich.go.mid.a.CC <- read_and_process_results('mid/mid.a.nonTF_ids', 'CC')
enrich.go.mid.b.CC <- read_and_process_results('mid/mid.b.nonTF_ids', 'CC')
enrich.go.high.a.CC <- read_and_process_results('high/high.a.nonTF_ids', 'CC')
enrich.go.high.b.CC <- read_and_process_results('high/high.b.nonTF_ids', 'CC')


# 提取前 10 个富集结果
top10_10C.CC.results <- enrich.go.10C.CC[1:10, ]
top10_16C.CC.results <- enrich.go.16C.CC[1:10, ]
top10_22C.a.CC.results <- enrich.go.22C.a.CC[1:10, ]
top10_22C.b.CC.results <- enrich.go.22C.b.CC[1:10, ]
top10_mid.a.CC.results <- enrich.go.mid.a.CC[1:10, ]
top10_mid.b.CC.results <- enrich.go.mid.b.CC[1:10, ]
top10_high.a.CC.results <- enrich.go.high.a.CC[1:10, ]
top10_high.b.CC.results <- enrich.go.high.b.CC[1:10, ]


# 假设有多个样本的富集结果 (list of enrichGO results)
sample_results_CC_nonTF <- list(
  "10" = top10_10C.CC.results,  # enrichGO 的结果
  "16" = top10_16C.CC.results,
  "22_a" = top10_22C.a.CC.results,
  "22_b" = top10_22C.b.CC.results,
  "mid_a" = top10_mid.a.CC.results,
  "mid_b" = top10_mid.b.CC.results,
  "high_a" = top10_high.a.CC.results,
  "high_b" = top10_high.b.CC.results
)


# 将多个样本的富集结果合并到一个数据框
combined_results_CC_nonTF <- do.call(rbind, lapply(names(sample_results_CC_nonTF), function(sample) {
  result <- sample_results_CC_nonTF[[sample]]  # 提取 enrichGO 结果
  result$Sample <- sample  # 添加样本名称
  return(result)
}))


combined_results_CC_nonTF <- combined_results_CC_nonTF %>%
  mutate(log_padj = -log10(p.adjust))  # 将 p.adjust 转换为 -log10(p.adjust)


# 将 GeneRatio 转换为数值类型
combined_results_CC_nonTF$GeneRatio <- sapply(strsplit(combined_results_CC_nonTF$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))

# 删除含有 NA 的行
combined_results_CC_nonTF <- na.omit(combined_results_CC_nonTF)

# 绘制CC点图
ggplot(combined_results_CC_nonTF, aes(x = Sample, y = Description)) +
  geom_point(aes(size = GeneRatio, color = log_padj)) +  # 根据 count 调整点的大小，根据 log_padj 调整颜色
  scale_color_gradient(low = "blue", high = "red") +  # 自定义颜色渐变
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), #文本加粗
    axis.text.y = element_text(size = 10, face = "bold"),  # y 轴文本加粗
    legend.title = element_text(face = "bold"),  # 图例标题加粗
    legend.text = element_text(face = "bold")  # 图例文本加粗
  ) +
  labs(
    x = "",
    y = "",
    size = "GeneRatio",  # 表示点大小的注释
    color = "-log10(p.adjust)"  # 表示颜色的注释
  ) 

ggsave("GO_Enrichment_CC_nonTF_Across_Samples.pdf", width = 7.5, height = 8)



