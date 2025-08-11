# 设置编码
Sys.setlocale("LC_ALL", "en_US.UTF-8")
options(encoding = "UTF-8")

# 加载必要的包
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(viridis)

# 读取数据
data <- read.table("Alt.RCA.me_snp.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

# 数据预处理
## 提取CG_1到CG_226列和温度信息（第1列）
cg_cols <- grep('^CG_[0-9]+$', colnames(data), value = TRUE)
cg_values <- as.matrix(data[, cg_cols])
temperature <- data[, 1]

# 执行PCA（如果数据维度允许）
## 基于其他列创建更多特征用于PCA
feature_cols <- c("Tem..", "Lat", "Long", "alt.m", cg_cols)
pca_data <- data[, feature_cols]
pca_data <- pca_data[complete.cases(pca_data), ]

# 去除方差为零的列
nzv <- apply(pca_data, 2, function(x) var(x, na.rm = TRUE) != 0)
pca_data_nzv <- pca_data[, nzv]

# 执行PCA
pca_result <- prcomp(pca_data_nzv, scale. = TRUE)

# 使用PCA结果对样本进行聚类
pca_coords <- pca_result$x[, 1:2]
dist_matrix <- dist(pca_coords)
hclust_result <- hclust(dist_matrix, method = "ward.D2")

# 根据聚类结果重新排序数据
cluster_order <- hclust_result$order
ordered_indices <- which(complete.cases(data[, feature_cols]))[cluster_order]

# 创建热图数据 - 转置矩阵使样本成为列
## 创建热图数据 - CG_1到CG_226，样本为列

heatmap_data <- t(cg_values[ordered_indices, ])
colnames(heatmap_data) <- ordered_indices
rownames(heatmap_data) <- cg_cols

# 自动调整热图和PDF尺寸
num_samples <- ncol(heatmap_data)
num_cg <- nrow(heatmap_data)
# 每个样本宽度0.25cm，最小宽度12cm，最大30cm
pdf_width <- min(max(num_samples * 0.25, 12), 30)
# 每个CG高度0.18cm，最小高度8cm，最大20cm
pdf_height <- min(max(num_cg * 0.18, 8), 20)

# 创建温度分组注释
temp_groups <- data$Tem..[ordered_indices]
temp_annotation <- data.frame(
  Temperature = factor(temp_groups, levels = c("10", "16", "22"))
)
rownames(temp_annotation) <- colnames(heatmap_data)

# 定义颜色
temp_colors <- c("10" = "blue", "16" = "yellow", "22" = "red")

# 使用ComplexHeatmap包创建热图
## 创建热图主体（无log10转换，直接用原始CG值）
ht <- Heatmap(
  heatmap_data,
  name = "CG",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  col = colorRamp2(
    seq(min(heatmap_data), max(heatmap_data), length.out = 100),
    viridis(100)
  ),
  height = unit(pdf_height, "cm"),
  heatmap_legend_param = list(
    title = "CG",
    title_position = "topcenter",
    legend_direction = "vertical",
    legend_width = unit(4, "cm"),
    labels_gp = gpar(fontsize = 8),
    at = pretty(range(heatmap_data), n = 5),
    labels = sprintf("%.1f", pretty(range(heatmap_data), n = 5))
  )
)

# 创建温度注释 - 底部注释
bottom_ha <- HeatmapAnnotation(
  Temperature = temp_annotation$Temperature,
  col = list(Temperature = temp_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  height = unit(0.8, "cm")
)

# 创建温度图例
temp_legend <- Legend(
  title = "Tem/°C",
  labels = c("10", "16", "22"),
  legend_gp = gpar(fill = c("blue", "yellow", "red")),
  title_gp = gpar(fontsize = 10),
  labels_gp = gpar(fontsize = 8)
)

# 组合热图和注释
final_heatmap <- ht
final_heatmap <- final_heatmap %v% bottom_ha

# 保存为PDF
pdf_file <- "rca_CG_pca_heatmap.pdf"
cairo_pdf(pdf_file, width = pdf_width, height = pdf_height, fallback_resolution = 1200)

# 创建图例列表
lgd_list <- packLegend(
  temp_legend,
  direction = "vertical"
)

# 绘制热图
draw(final_heatmap,
     heatmap_legend_list = list(lgd_list),
     heatmap_legend_side = "right",
     legend_gap = unit(1, "cm"),
     padding = unit(c(2, 6, 5, 2), "cm"))

dev.off()
cat("热图已成功保存到:", pdf_file, "\n")

# 输出PCA结果摘要
cat("PCA Summary:\n")
summary(pca_result)
cat("\nProportion of variance explained by first two components:\n")
cat("PC1:", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 2), "%\n")
cat("PC2:", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 2), "%\n")

# 输出样本数量统计
cat("\nSample counts by temperature:\n")
print(table(data$Tem..))


