library(parallel)
library(WGCNA)
options(stringsAsFactors = FALSE) 
enableWGCNAThreads() ## 打开多线程解释
library(ggplot2)

# total Rca data
#1.数据预处理
dataExprtotal <- read.table("RCA.tran.txt", header=T, comment.char = "", check.names=F)

rca_total <- as.data.frame(dataExprtotal[, -1])
rownames(rca_total) <- dataExprtotal$Name


# rca_total
#2.筛选方差前25%的基因##
m.vars <- apply(rca_total,1,var)
rca_total.upper <- rca_total[which(m.vars > quantile(m.vars, probs = seq(0, 1, 0.25))[4]),]
dim(rca_total.upper)

#3.聚类前数据转置##
rca_total=as.data.frame(t(rca_total.upper));

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
pdf("total/rca.cluster.pdf", height = 20, width = 30)  # 打开 PDF 设备
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
dev.off()  # 关闭设备，保存文件
save(rca_total, file = "total/rca.Rdata")

head(rca_total)[,1:8]

# 5.软阈值的预设范围
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
# 自动计算推荐的软阈值
sft <- pickSoftThreshold(rca_total, powerVector = powers)
#Warning message:executing %dopar% sequentially: no parallel backend registered 
# 推荐值。如果是NA，就需要画图来自己挑选
sft$powerEstimate #最佳beta

# 创建PDF文件
pdf("total/soft_threshold_example.pdf", width = 10, height = 5)  # 设置宽度为10，高度为5

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
net <- blockwiseModules(rca_total, power = 14, maxBlockSize = 20000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, saveTOMFileBase = "total/femaleMouseTOM", verbose = 3)
table(net$colors)

#7.模块可视化
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# 打开 PDF 设备 
pdf("total/module_dendrogram.pdf", height = 15, width = 15) 
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
TOM=TOMsimilarityFromExpr(rca_total, power=14)
modules= c("turquoise")
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
      edgeFile = paste("total/CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
      nodeFile = paste("total/CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
      weighted = TRUE,
      threshold = 0.05,
      nodeNames = modProbes, 
      nodeAttr = moduleColors[inModule])
      
#11.绘制（合并后）模块之间相关性图
MEDiss <- 1-cor(MEs)#距离矩阵
METree <- hclust(as.dist(MEDiss), method = "average")
tiff(file="total/METree. moduleColors.tiff",width=17,height=15, units="cm", compression="lzw", res=1200)
p1 <- plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
print(p1)
dev.off()

#12.可视化-画模块之间的热图
dissTOM = 1-TOMsimilarityFromExpr(rca_total, power = 14);
tiff(file="total/moduleColors.module_heatmap.tiff",width=17,height=17, units="cm", compression="lzw", res=1200)
plotTOM <- dissTOM^7##为了更显著，用7次方
p1 <- TOMplot(plotTOM, geneTree, moduleColors, main="Network heatmap plot, all genes")
print(p1)
dev.off()

#13.可视化-画模块之间的热图
tiff(file="total/moduleColors.Eigengene adjacency heatmap.tiff",width=17,height=21, units="cm", compression="lzw", res=1200)
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
write.table(datout, "total/all.module_gene_module.xls", sep="\t", row.names=F,quote=F)

