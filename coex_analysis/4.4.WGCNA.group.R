library(parallel)
library(WGCNA)
options(stringsAsFactors = FALSE) 
enableWGCNAThreads() ## 打开多线程解释
library(ggplot2)

# Rca_10C data
#1.数据预处理
dataExprtotal <- read.table("RCA.10C.tran.txt", header=T, comment.char = "", check.names=F)

# 基因TPM表达矩阵筛选：
# 行是基因，列是样品；保留 TPM > 1 的样本数占总样本数 10% 及以上的基因
expr_mat <- dataExprtotal[, -1]
threshold <- ceiling(ncol(expr_mat) * 0.1)
dataExprtotal <- dataExprtotal[
  rowSums(expr_mat > 1, na.rm = TRUE) >= threshold,
]

# 提取表达矩阵，并保留“行=基因，列=样品”的格式
rca_total <- as.data.frame(dataExprtotal[, -1])
rownames(rca_total) <- dataExprtotal$Name
# 检查维度
dim(rca_total)
# rca_total
#2.筛选方差前75%的基因##
# m.vars <- apply(rca_total,1,var)
# rca_total.upper <- rca_total[which(m.vars > quantile(m.vars, probs = seq(0, 1, 0.25))[4]),]
# dim(rca_total.upper)

#3.聚类前数据转置##
rca_total=as.data.frame(t(rca_total));
# 统计基因数和样品数
nGenes = ncol(rca_total)
nSamples = nrow(rca_total)
dim(rca_total)

#4.样本聚类检查离群值##
gsg = goodSamplesGenes(rca_total, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(rca_total)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(rca_total)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  rca_total = rca_total[gsg$goodSamples, gsg$goodGenes]
}

#查看是否有离群样品
sampleTree = hclust(dist(rca_total), method = "average")
# 使用 pdf 保存图像
pdf("10C/rca.10C.cluster.pdf", height = 20, width = 30)  # 打开 PDF 设备
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
dev.off()  # 关闭设备，保存文件
save(rca_total, file = "10C/rca.10C.Rdata")

head(rca_total)[,1:8]

# 5.软阈值的预设范围
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
# 自动计算推荐的软阈值
sft <- pickSoftThreshold(rca_total, powerVector = powers)
#Warning message:executing %dopar% sequentially: no parallel backend registered 
# 推荐值。如果是NA，就需要画图来自己挑选
sft$powerEstimate #最佳beta

# 创建PDF文件
pdf("10C/10C.soft_threshold_example.pdf", width = 10, height = 5)  # 设置宽度为10，高度为5

# 设置图形布局：1行2列
par(mfrow = c(1, 2))

# 文字大小
cex1 = 0.9

# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2", 
     type = "n", 
     main = "Scale independence")  # 设置坐标和标题
text(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     labels = powers, 
     cex = cex1, 
     col = "red")  # 添加文本标注
abline(h = 0.80, col = "red")  # 添加R^2=0.80的参考线

# Soft threshold与平均连通性
plot(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", 
     type = "n", 
     main = "Mean connectivity")  # 设置坐标和标题
text(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     labels = powers, 
     cex = cex1, 
     col = "red")  # 添加文本标注

# 关闭PDF文件
dev.off()

#6.构建网络
net <- blockwiseModules(rca_total, power = 16, maxBlockSize = 20000,,networkType = "signed",
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, saveTOMFileBase = "10C/femaleMouseTOM", verbose = 3)
table(net$colors)

#7.模块可视化
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# 打开 PDF 设备 
pdf("10C/10C.module_dendrogram.pdf", height = 15, width = 15) 
# 绘制图形 
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], 
                    "Module colors", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05) 
# 关闭 PDF 设备 
dev.off()

#8. 计算模块特征向量Mes
table(moduleColors)
MEs=net$MEs
geneTree =net$dendrograms[[1]]

#9.模块的导出
#主要模块里面的基因直接的相互作用关系信息可以导出到cytoscape等网络可视化软件。
TOM=TOMsimilarityFromExpr(rca_total, power=16)
modules= c("greenyellow")
probes = colnames(rca_total)
inModule =is.finite(match(moduleColors,modules));
modProbes=probes[inModule] #确定保留下来的。
# 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
## 模块对应的基因关系矩阵 

#10.导出网络到cytoscape
cyt = exportNetworkToCytoscape(
       modTOM,
      edgeFile = paste("10C/10C.CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
      nodeFile = paste("10C/10C.CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
      weighted = TRUE,
      threshold = 0.05,
      nodeNames = modProbes, 
      nodeAttr = moduleColors[inModule])
      
#11.绘制（合并后）模块之间相关性图
MEDiss <- 1-cor(MEs)#距离矩阵
METree <- hclust(as.dist(MEDiss), method = "average")
tiff(file="10C/10C.METree. moduleColors.tiff",width=17,height=15, units="cm", compression="lzw", res=1200)
p1 <- plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
print(p1)
dev.off()

#12.可视化-画模块之间的热图
dissTOM = 1-TOMsimilarityFromExpr(rca_total, power = 16);
tiff(file="10C/10C.moduleColors.module_heatmap.tiff",width=17,height=17, units="cm", compression="lzw", res=1200)
plotTOM <- dissTOM^7##为了更显著，用7次方
p1 <- TOMplot(plotTOM, geneTree, moduleColors, main="Network heatmap plot, all genes")
print(p1)
dev.off()

#13.可视化-画模块之间的热图
tiff(file="10C/10C.moduleColors.Eigengene adjacency heatmap.tiff",width=17,height=21, units="cm", compression="lzw", res=1200)
p2 <- plotEigengeneNetworks(MEs, 
                      "Eigengene adjacency heatmap", 
                      marHeatmap = c(4,4,4,4), 
                     plotDendrograms = T, 
                      xLabelsAngle = 90) 
print(p2)
dev.off()
# marHeatmap 设置下、左、上、右的边距
#根据基因间表达量进行聚类所得到的各模块间的相关性图

#14.输出模块基因!!!!
datME=moduleEigengenes(rca_total, moduleColors)[[1]]
#color1=rep("pink",dim(dataExpr)[[2]])#仅仅获得pink的gene
color1=as.character(moduleColors)#获得全部module对应的基因
datKME=signedKME(rca_total, datME)
dataExpr1=as.data.frame(t(rca_total));

datSummary=rownames(dataExpr1)
datout=data.frame(datSummary,colorNEW=color1,datKME )
#write.table(datout, "pink_gene_module.xls", sep="\t", row.names=F,quote=F)
write.table(datout, "10C/10C.module_gene_module.xls", sep="\t", row.names=F,quote=F)

# 输出greenyellow模块的基因
greenyellow_genes <- datout[datout$colorNEW == "greenyellow", ]
write.table(greenyellow_genes, "10C/10C.greenyellow_module_genes.txt", sep="\t", row.names=F, quote=F)
cat("Greenyellow模块包含", nrow(greenyellow_genes), "个基因\n")

#15.找与ME11模块负相关最高的模块
# 计算模块间的相关性
module_cor <- cor(MEs, use = "pairwise.complete.obs")
# 提取ME11与其他模块的相关性
ME11_cor <- module_cor["ME11", ]
# 找到与ME11负相关最高的模块（排除ME11自身）
ME11_cor_others <- ME11_cor[names(ME11_cor) != "ME11"]
most_neg_module <- names(which.min(ME11_cor_others))
cat("与ME11负相关最高的模块是:", most_neg_module, "\n")
cat("相关系数为:", ME11_cor_others[most_neg_module], "\n")

#主要模块里面的基因直接的相互作用关系信息可以导出到cytoscape等网络可视化软件。
TOM=TOMsimilarityFromExpr(rca_total, power=16)
modules= c("green") #与ME11负相关最高的模块是green
probes = colnames(rca_total)
inModule =is.finite(match(moduleColors,modules));
modProbes=probes[inModule] #确定保留下来的。
# 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
## 模块对应的基因关系矩阵 

# 导出网络到cytoscape
cyt = exportNetworkToCytoscape(
       modTOM,
      edgeFile = paste("10C/10C.negative.CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
      nodeFile = paste("10C/10C.negative.CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
      weighted = TRUE,
      threshold = 0.05,
      nodeNames = modProbes, 
      nodeAttr = moduleColors[inModule])

# 输出green模块的基因
green_genes <- datout[datout$colorNEW == "green", ]
write.table(green_genes, "10C/10C.negative_green_module_genes.txt", sep="\t", row.names=F, quote=F)
cat("Green模块包含", nrow(green_genes), "个基因\n")

# #16.提取该负相关模块的基因
# target_module_color <- "grey"
# target_genes <- colnames(rca_total)[moduleColors == target_module_color]
# cat("该模块包含", length(target_genes), "个基因\n")

# #17.提取RCA基因表达数据并计算基因与Rca的相关性
# # RCA是基因ID AT2G39730，从表达数据中提取其表达值
# RCA_gene <- "AT2G39730"

# # 检查RCA基因是否在数据中
# if (!RCA_gene %in% colnames(rca_total)) {
#   stop(paste("RCA基因", RCA_gene, "不在表达数据中！"))
# }

# # 提取RCA基因的表达值作为向量
# Rca <- rca_total[, RCA_gene]
# cat("成功提取RCA基因 (", RCA_gene, ") 的表达数据\n")

# # 计算目标模块基因与Rca的相关性
# cor_results <- corAndPvalue(rca_total[, target_genes], Rca, use = "pairwise.complete.obs")
# gene_rca_cor <- cor_results$cor
# gene_rca_pval <- cor_results$p

# # FDR校正
# gene_rca_fdr <- p.adjust(gene_rca_pval, method = "fdr")

# # 创建结果数据框
# gene_rca_results <- data.frame(
#   Gene = target_genes,
#   Correlation = as.vector(gene_rca_cor),
#   P_value = as.vector(gene_rca_pval),
#   FDR_q = gene_rca_fdr
# )

# #18.筛选符合条件的基因 (r ≤ -0.3, FDR q < 0.05)
# significant_genes <- gene_rca_results[
#   !is.na(gene_rca_results$Correlation) & 
#   !is.na(gene_rca_results$FDR_q) &
#   gene_rca_results$Correlation <= -0.3 & 
#   gene_rca_results$FDR_q < 0.05, 
# ]

# cat("筛选出", nrow(significant_genes), "个符合条件的基因 (r ≤ -0.3, FDR q < 0.05)\n")

# # 只有在有结果时才进行排序和输出
# if (nrow(significant_genes) > 0) {
#   # 按相关系数排序
#   significant_genes <- significant_genes[order(significant_genes$Correlation), ]
  
#   cat("\n前几个基因:\n")
#   print(head(significant_genes))
  
#   # 导出结果
#   write.table(significant_genes, 
#               "10C/ME2_negative_correlated_genes.txt", 
#               sep = "\t", 
#               row.names = FALSE, 
#               quote = FALSE)
  
#   cat("\n结果已导出到: 10C/ME2_negative_correlated_genes.txt\n")
# } else {
#   cat("警告：没有找到符合条件的基因！\n")
# }

# Rca_16C data
#1.数据预处理
dataExprtotal <- read.table("RCA.16C.tran.txt", header=T, comment.char = "", check.names=F)

# 基因TPM表达矩阵筛选：
# 行是基因，列是样品；保留 TPM > 1 的样本数占总样本数 10% 及以上的基因
expr_mat <- dataExprtotal[, -1]
threshold <- ceiling(ncol(expr_mat) * 0.1)
dataExprtotal <- dataExprtotal[
  rowSums(expr_mat > 1, na.rm = TRUE) >= threshold,
]

# 提取表达矩阵，并保留“行=基因，列=样品”的格式
rca_total <- as.data.frame(dataExprtotal[, -1])
rownames(rca_total) <- dataExprtotal$Name
# 检查维度
dim(rca_total)
# rca_total
#2.筛选方差前25%的基因##
# m.vars <- apply(rca_total,1,var)
# rca_total.upper <- rca_total[which(m.vars > quantile(m.vars, probs = seq(0, 1, 0.25))[4]),]
# dim(rca_total.upper)

#3.聚类前数据转置##
rca_total=as.data.frame(t(rca_total));
# 统计基因数和样品数
nGenes = ncol(rca_total)
nSamples = nrow(rca_total)
dim(rca_total)

#4.样本聚类检查离群值##
gsg = goodSamplesGenes(rca_total, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(rca_total)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(rca_total)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  rca_total = rca_total[gsg$goodSamples, gsg$goodGenes]
}

#查看是否有离群样品
sampleTree = hclust(dist(rca_total), method = "average")
# 使用 pdf 保存图像
pdf("16C/rca.16C.cluster.pdf", height = 20, width = 30)  # 打开 PDF 设备
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
dev.off()  # 关闭设备，保存文件
save(rca_total, file = "16C/rca.16C.Rdata")

head(rca_total)[,1:8]

# 5.软阈值的预设范围
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
# 自动计算推荐的软阈值
sft <- pickSoftThreshold(rca_total, powerVector = powers)
#Warning message:executing %dopar% sequentially: no parallel backend registered 
# 推荐值。如果是NA，就需要画图来自己挑选
sft$powerEstimate #最佳beta

# 创建PDF文件
pdf("16C/16C.soft_threshold_example.pdf", width = 10, height = 5)  # 设置宽度为10，高度为5

# 设置图形布局：1行2列
par(mfrow = c(1, 2))

# 文字大小
cex1 = 0.9

# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2", 
     type = "n", 
     main = "Scale independence")  # 设置坐标和标题
text(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     labels = powers, 
     cex = cex1, 
     col = "red")  # 添加文本标注
abline(h = 0.80, col = "red")  # 添加R^2=0.80的参考线

# Soft threshold与平均连通性
plot(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", 
     type = "n", 
     main = "Mean connectivity")  # 设置坐标和标题
text(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     labels = powers, 
     cex = cex1, 
     col = "red")  # 添加文本标注

# 关闭PDF文件
dev.off()

#6.构建网络
net <- blockwiseModules(rca_total, power = 20, maxBlockSize = 20000,networkType = "signed",
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, saveTOMFileBase = "16C/femaleMouseTOM", verbose = 3)
table(net$colors)

#7.模块可视化
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# 打开 PDF 设备 
pdf("16C/module_dendrogram.pdf", height = 15, width = 15) 
# 绘制图形 
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], 
                    "Module colors", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05) 
# 关闭 PDF 设备 
dev.off()

#8. 计算模块特征向量Mes
table(moduleColors)
MEs=net$MEs
geneTree =net$dendrograms[[1]]

#9.模块的导出
#主要模块里面的基因直接的相互作用关系信息可以导出到cytoscape等网络可视化软件。
TOM=TOMsimilarityFromExpr(rca_total, power=20)
modules= c("pink")
probes = colnames(rca_total)
inModule =is.finite(match(moduleColors,modules));
modProbes=probes[inModule] #确定保留下来的。
# 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
## 模块对应的基因关系矩阵 

#10.导出网络到cytoscape
cyt = exportNetworkToCytoscape(
       modTOM,
      edgeFile = paste("16C/16C.CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
      nodeFile = paste("16C/16C.CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
      weighted = TRUE,
      threshold = 0.05,
      nodeNames = modProbes, 
      nodeAttr = moduleColors[inModule])
      
#11.绘制（合并后）模块之间相关性图
MEDiss <- 1-cor(MEs)#距离矩阵
METree <- hclust(as.dist(MEDiss), method = "average")
tiff(file="16C/16C.METree. moduleColors.tiff",width=17,height=15, units="cm", compression="lzw", res=1200)
p1 <- plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
print(p1)
dev.off()

#12.可视化-画模块之间的热图
dissTOM = 1-TOMsimilarityFromExpr(rca_total, power = 20);
tiff(file="16C/16C.moduleColors.module_heatmap.tiff",width=17,height=17, units="cm", compression="lzw", res=1200)
plotTOM <- dissTOM^7##为了更显著，用7次方
p1 <- TOMplot(plotTOM, geneTree, moduleColors, main="Network heatmap plot, all genes")
print(p1)
dev.off()

#13.可视化-画模块之间的热图
tiff(file="16C/16C.moduleColors.Eigengene adjacency heatmap.tiff",width=17,height=21, units="cm", compression="lzw", res=1200)
p2 <- plotEigengeneNetworks(MEs, 
                      "Eigengene adjacency heatmap", 
                      marHeatmap = c(4,4,4,4), 
                     plotDendrograms = T, 
                      xLabelsAngle = 90) 
print(p2)
dev.off()
# marHeatmap 设置下、左、上、右的边距
#根据基因间表达量进行聚类所得到的各模块间的相关性图

#14.输出模块基因!!!!
datME=moduleEigengenes(rca_total, moduleColors)[[1]]
#color1=rep("pink",dim(dataExpr)[[2]])#仅仅获得pink的gene
color1=as.character(moduleColors)#获得全部module对应的基因
datKME=signedKME(rca_total, datME)
dataExpr1=as.data.frame(t(rca_total));

datSummary=rownames(dataExpr1)
datout=data.frame(datSummary,colorNEW=color1,datKME )
#write.table(datout, "pink_gene_module.xls", sep="\t", row.names=F,quote=F)
write.table(datout, "16C/16C.module_gene_module.xls", sep="\t", row.names=F,quote=F)

# 输出pink模块的基因
pink_genes <- datout[datout$colorNEW == "pink", ]
write.table(pink_genes, "16C/16C.pink_module_genes.txt", sep="\t", row.names=F, quote=F)
cat("Pink模块包含", nrow(pink_genes), "个基因\n")

#15.找与ME8模块负相关最高的模块
# 计算模块间的相关性
module_cor <- cor(MEs, use = "pairwise.complete.obs")
# 提取ME8与其他模块的相关性
ME8_cor <- module_cor["ME8", ]
# 找到与ME8负相关最高的模块（排除ME8自身）
ME8_cor_others <- ME8_cor[names(ME8_cor) != "ME8"]
most_neg_module <- names(which.min(ME8_cor_others))
cat("与ME8负相关最高的模块是:", most_neg_module, "\n")
cat("相关系数为:", ME8_cor_others[most_neg_module], "\n")

#主要模块里面的基因直接的相互作用关系信息可以导出到cytoscape等网络可视化软件。
TOM=TOMsimilarityFromExpr(rca_total, power=20)
modules= c("salmon") #与ME8负相关最高的模块是salmon
probes = colnames(rca_total)
inModule =is.finite(match(moduleColors,modules));
modProbes=probes[inModule] #确定保留下来的。
# 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
## 模块对应的基因关系矩阵 

# 导出网络到cytoscape
cyt = exportNetworkToCytoscape(
       modTOM,
      edgeFile = paste("16C/16C.negative.CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
      nodeFile = paste("16C/16C.negative.CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
      weighted = TRUE,
      threshold = 0.05,
      nodeNames = modProbes, 
      nodeAttr = moduleColors[inModule])

# 输出salmon模块的基因
salmon_genes <- datout[datout$colorNEW == "salmon", ]
write.table(salmon_genes, "16C/16C.negative_salmon_module_genes.txt", sep="\t", row.names=F, quote=F)
cat("Salmon模块包含", nrow(salmon_genes), "个基因\n")

# #16.提取该负相关模块的基因
# target_module_color <- "pink"
# target_genes <- colnames(rca_total)[moduleColors == target_module_color]
# cat("该模块包含", length(target_genes), "个基因\n")

# #17.提取RCA基因表达数据并计算基因与Rca的相关性
# # RCA是基因ID AT2G39730，从表达数据中提取其表达值
# RCA_gene <- "AT2G39730"

# # 检查RCA基因是否在数据中
# if (!RCA_gene %in% colnames(rca_total)) {
#   stop(paste("RCA基因", RCA_gene, "不在表达数据中！"))
# }

# # 提取RCA基因的表达值作为向量
# Rca <- rca_total[, RCA_gene]
# cat("成功提取RCA基因 (", RCA_gene, ") 的表达数据\n")

# # 计算目标模块基因与Rca的相关性
# cor_results <- corAndPvalue(rca_total[, target_genes], Rca, use = "pairwise.complete.obs")
# gene_rca_cor <- cor_results$cor
# gene_rca_pval <- cor_results$p

# # FDR校正
# gene_rca_fdr <- p.adjust(gene_rca_pval, method = "fdr")

# # 创建结果数据框
# gene_rca_results <- data.frame(
#   Gene = target_genes,
#   Correlation = as.vector(gene_rca_cor),
#   P_value = as.vector(gene_rca_pval),
#   FDR_q = gene_rca_fdr
# )

# #18.筛选符合条件的基因 (r ≤ -0.3, FDR q < 0.05)
# significant_genes <- gene_rca_results[
#   !is.na(gene_rca_results$Correlation) & 
#   !is.na(gene_rca_results$FDR_q) &
#   gene_rca_results$Correlation <= -0.3 & 
#   gene_rca_results$FDR_q < 0.05, 
# ]

# cat("筛选出", nrow(significant_genes), "个符合条件的基因 (r ≤ -0.3, FDR q < 0.05)\n")

# # 只有在有结果时才进行排序和输出
# if (nrow(significant_genes) > 0) {
#   # 按相关系数排序
#   significant_genes <- significant_genes[order(significant_genes$Correlation), ]
  
#   cat("\n前几个基因:\n")
#   print(head(significant_genes))
  
#   # 导出结果
#   write.table(significant_genes, 
#               "16C/ME2_negative_correlated_genes.txt", 
#               sep = "\t", 
#               row.names = FALSE, 
#               quote = FALSE)
  
#   cat("\n结果已导出到: 16C/ME2_negative_correlated_genes.txt\n")
# } else {
#   cat("警告：没有找到符合条件的基因！\n")
# }

# Rca_22C data
#1.数据预处理
dataExprtotal <- read.table("RCA.22C.tran.txt", header=T, comment.char = "", check.names=F)

# 基因TPM表达矩阵筛选：
# 行是基因，列是样品；保留 TPM > 1 的样本数占总样本数 10% 及以上的基因
expr_mat <- dataExprtotal[, -1]
threshold <- ceiling(ncol(expr_mat) * 0.1)
dataExprtotal <- dataExprtotal[
  rowSums(expr_mat > 1, na.rm = TRUE) >= threshold,
]

# 提取表达矩阵，并保留“行=基因，列=样品”的格式
rca_total <- as.data.frame(dataExprtotal[, -1])
rownames(rca_total) <- dataExprtotal$Name
# 检查维度
dim(rca_total)
# rca_total
#2.筛选方差前25%的基因##
# m.vars <- apply(rca_total,1,var)
# rca_total.upper <- rca_total[which(m.vars > quantile(m.vars, probs = seq(0, 1, 0.25))[4]),]
# dim(rca_total.upper)

#3.聚类前数据转置##
rca_total=as.data.frame(t(rca_total));
# 统计基因数和样品数
nGenes = ncol(rca_total)
nSamples = nrow(rca_total)
dim(rca_total)

#4.样本聚类检查离群值##
gsg = goodSamplesGenes(rca_total, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(rca_total)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(rca_total)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  rca_total = rca_total[gsg$goodSamples, gsg$goodGenes]
}

#查看是否有离群样品
sampleTree = hclust(dist(rca_total), method = "average")
# 使用 pdf 保存图像
pdf("22C/rca.22C.cluster.pdf", height = 20, width = 30)  # 打开 PDF 设备
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
dev.off()  # 关闭设备，保存文件
save(rca_total, file = "22C/rca.22C.Rdata")

head(rca_total)[,1:8]

# 5.软阈值的预设范围
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
# 自动计算推荐的软阈值
sft <- pickSoftThreshold(rca_total, powerVector = powers)
#Warning message:executing %dopar% sequentially: no parallel backend registered 
# 推荐值。如果是NA，就需要画图来自己挑选
sft$powerEstimate #最佳beta

# 创建PDF文件
pdf("22C/22C.soft_threshold_example.pdf", width = 10, height = 5)  # 设置宽度为10，高度为5

# 设置图形布局：1行2列
par(mfrow = c(1, 2))

# 文字大小
cex1 = 0.9

# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2", 
     type = "n", 
     main = "Scale independence")  # 设置坐标和标题
text(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     labels = powers, 
     cex = cex1, 
     col = "red")  # 添加文本标注
abline(h = 0.80, col = "red")  # 添加R^2=0.80的参考线

# Soft threshold与平均连通性
plot(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", 
     type = "n", 
     main = "Mean connectivity")  # 设置坐标和标题
text(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     labels = powers, 
     cex = cex1, 
     col = "red")  # 添加文本标注

# 关闭PDF文件
dev.off()

#6.构建网络
net <- blockwiseModules(rca_total, power = 9, maxBlockSize = 20000,networkType = "signed",
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, saveTOMFileBase = "22C/femaleMouseTOM", verbose = 3)
table(net$colors)

#7.模块可视化
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# 打开 PDF 设备 
pdf("22C/module_dendrogram.pdf", height = 15, width = 15) 
# 绘制图形 
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], 
                    "Module colors", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05) 
# 关闭 PDF 设备 
dev.off()

#8. 计算模块特征向量Mes
table(moduleColors)
MEs=net$MEs
geneTree =net$dendrograms[[1]]

#9.模块的导出
#主要模块里面的基因直接的相互作用关系信息可以导出到cytoscape等网络可视化软件。
TOM=TOMsimilarityFromExpr(rca_total, power=9)
modules= c("red")
probes = colnames(rca_total)
inModule =is.finite(match(moduleColors,modules));
modProbes=probes[inModule] #确定保留下来的。
# 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
## 模块对应的基因关系矩阵 

#10.导出网络到cytoscape
cyt = exportNetworkToCytoscape(
       modTOM,
      edgeFile = paste("22C/22C.CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
      nodeFile = paste("22C/22C.CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
      weighted = TRUE,
      threshold = 0.05,
      nodeNames = modProbes, 
      nodeAttr = moduleColors[inModule])
      
#11.绘制（合并后）模块之间相关性图
MEDiss <- 1-cor(MEs)#距离矩阵
METree <- hclust(as.dist(MEDiss), method = "average")
tiff(file="22C/22C.METree. moduleColors.tiff",width=17,height=15, units="cm", compression="lzw", res=1200)
p1 <- plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
print(p1)
dev.off()

#12.可视化-画模块之间的热图
dissTOM = 1-TOMsimilarityFromExpr(rca_total, power = 9);
tiff(file="22C/22C.moduleColors.module_heatmap.tiff",width=17,height=17, units="cm", compression="lzw", res=1200)
plotTOM <- dissTOM^7##为了更显著，用7次方
p1 <- TOMplot(plotTOM, geneTree, moduleColors, main="Network heatmap plot, all genes")
print(p1)
dev.off()

#13.可视化-画模块之间的热图
tiff(file="22C/22C.moduleColors.Eigengene adjacency heatmap.tiff",width=17,height=21, units="cm", compression="lzw", res=1200)
p2 <- plotEigengeneNetworks(MEs, 
                      "Eigengene adjacency heatmap", 
                      marHeatmap = c(4,4,4,4), 
                     plotDendrograms = T, 
                      xLabelsAngle = 90) 
print(p2)
dev.off()
# marHeatmap 设置下、左、上、右的边距
#根据基因间表达量进行聚类所得到的各模块间的相关性图

#14.输出模块基因!!!!
datME=moduleEigengenes(rca_total, moduleColors)[[1]]
#color1=rep("pink",dim(dataExpr)[[2]])#仅仅获得pink的gene
color1=as.character(moduleColors)#获得全部module对应的基因
datKME=signedKME(rca_total, datME)
dataExpr1=as.data.frame(t(rca_total));

datSummary=rownames(dataExpr1)
datout=data.frame(datSummary,colorNEW=color1,datKME )
#write.table(datout, "pink_gene_module.xls", sep="\t", row.names=F,quote=F)
write.table(datout, "22C/22C.module_gene_module.xls", sep="\t", row.names=F,quote=F)

# 输出red模块的基因
red_genes <- datout[datout$colorNEW == "red", ]
write.table(red_genes, "22C/22C.red_module_genes.txt", sep="\t", row.names=F, quote=F)
cat("Red模块包含", nrow(red_genes), "个基因\n")

#15.找与ME6模块负相关最高的模块
# 计算模块间的相关性
module_cor <- cor(MEs, use = "pairwise.complete.obs")
# 提取ME6与其他模块的相关性
ME6_cor <- module_cor["ME6", ]
# 找到与ME6负相关最高的模块（排除ME6自身）
ME6_cor_others <- ME6_cor[names(ME6_cor) != "ME6"]
most_neg_module <- names(which.min(ME6_cor_others))
cat("与ME6负相关最高的模块是:", most_neg_module, "\n")
cat("相关系数为:", ME6_cor_others[most_neg_module], "\n")

#主要模块里面的基因直接的相互作用关系信息可以导出到cytoscape等网络可视化软件。
TOM=TOMsimilarityFromExpr(rca_total, power=9)
modules= c("magenta") #与ME6负相关最高的模块是magenta
probes = colnames(rca_total)
inModule =is.finite(match(moduleColors,modules));
modProbes=probes[inModule] #确定保留下来的。
# 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
## 模块对应的基因关系矩阵 

# 导出网络到cytoscape
cyt = exportNetworkToCytoscape(
       modTOM,
      edgeFile = paste("22C/22C.negative.CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
      nodeFile = paste("22C/22C.negative.CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
      weighted = TRUE,
      threshold = 0.05,
      nodeNames = modProbes, 
      nodeAttr = moduleColors[inModule])

# 输出magenta模块的基因
magenta_genes <- datout[datout$colorNEW == "magenta", ]
write.table(magenta_genes, "22C/22C.negative_magenta_module_genes.txt", sep="\t", row.names=F, quote=F)
cat("Magenta模块包含", nrow(magenta_genes), "个基因\n")

# #16.提取该负相关模块的基因
# target_module_color <- "black"
# target_genes <- colnames(rca_total)[moduleColors == target_module_color]
# cat("该模块包含", length(target_genes), "个基因\n")

# #17.提取RCA基因表达数据并计算基因与Rca的相关性
# # RCA是基因ID AT2G39730，从表达数据中提取其表达值
# RCA_gene <- "AT2G39730"

# # 检查RCA基因是否在数据中
# if (!RCA_gene %in% colnames(rca_total)) {
#   stop(paste("RCA基因", RCA_gene, "不在表达数据中！"))
# }

# # 提取RCA基因的表达值作为向量
# Rca <- rca_total[, RCA_gene]
# cat("成功提取RCA基因 (", RCA_gene, ") 的表达数据\n")

# # 计算目标模块基因与Rca的相关性
# cor_results <- corAndPvalue(rca_total[, target_genes], Rca, use = "pairwise.complete.obs")
# gene_rca_cor <- cor_results$cor
# gene_rca_pval <- cor_results$p

# # FDR校正
# gene_rca_fdr <- p.adjust(gene_rca_pval, method = "fdr")

# # 创建结果数据框
# gene_rca_results <- data.frame(
#   Gene = target_genes,
#   Correlation = as.vector(gene_rca_cor),
#   P_value = as.vector(gene_rca_pval),
#   FDR_q = gene_rca_fdr
# )

# #18.筛选符合条件的基因 (r ≤ -0.3, FDR q < 0.05)
# significant_genes <- gene_rca_results[
#   !is.na(gene_rca_results$Correlation) & 
#   !is.na(gene_rca_results$FDR_q) &
#   gene_rca_results$Correlation <= -0.3 & 
#   gene_rca_results$FDR_q < 0.05, 
# ]

# cat("筛选出", nrow(significant_genes), "个符合条件的基因 (r ≤ -0.3, FDR q < 0.05)\n")

# # 只有在有结果时才进行排序和输出
# if (nrow(significant_genes) > 0) {
#   # 按相关系数排序
#   significant_genes <- significant_genes[order(significant_genes$Correlation), ]
  
#   cat("\n前几个基因:\n")
#   print(head(significant_genes))
  
#   # 导出结果
#   write.table(significant_genes, 
#               "22C/ME4_negative_correlated_genes.txt", 
#               sep = "\t", 
#               row.names = FALSE, 
#               quote = FALSE)
  
#   cat("\n结果已导出到: 22C/ME4_negative_correlated_genes.txt\n")
# } else {
#   cat("警告：没有找到符合条件的基因！\n")
# }

