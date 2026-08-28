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
data <- read.table("Alt.RCA.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 数据预处理
# 提取TPM值（第5列）和温度信息（第1列）
tpm_values <- data$TPM
temperature <- data[, 1]

# 执行PCA（如果数据维度允许）
# 基于其他列创建更多特征用于PCA
feature_cols <- c("Tem..", "Lat", "Long", "alt.m", "TPM")
pca_data <- data[, feature_cols]
pca_data <- pca_data[complete.cases(pca_data), ]

# 执行PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# 使用PCA结果对样本进行聚类
pca_coords <- pca_result$x[, 1:2]
dist_matrix <- dist(pca_coords)
hclust_result <- hclust(dist_matrix, method = "ward.D2")

# 根据聚类结果重新排序数据
cluster_order <- hclust_result$order
ordered_indices <- which(complete.cases(data[, feature_cols]))[cluster_order]

# 创建热图数据 - 转置矩阵使样本成为列
heatmap_data <- matrix(data$TPM[ordered_indices], nrow = 1)
colnames(heatmap_data) <- ordered_indices
rownames(heatmap_data) <- "total"

# 对TPM值进行log10转换
heatmap_data_log <- log10(heatmap_data)  # 加1避免log(0)

# 创建温度分组注释
temp_groups <- data$Tem..[ordered_indices]
temp_annotation <- data.frame(
  Temperature = factor(temp_groups, levels = c("10", "16", "22"))
)
rownames(temp_annotation) <- colnames(heatmap_data_log)

# 定义颜色
temp_colors <- c("10" = "blue", "16" = "yellow", "22" = "red")

# 使用ComplexHeatmap包创建热图
# 创建热图主体
ht <- Heatmap(
  heatmap_data_log,
  name = "log10(TPM)",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 12),  # Y轴标签字体大小
  row_names_side = "right",  # Y轴标签放在右侧
  col = colorRamp2(
    seq(min(heatmap_data_log, na.rm = TRUE), max(heatmap_data_log, na.rm = TRUE), length.out = 100),
    viridis(100)
  ),
  heatmap_legend_param = list(
    title = "log10(Tpm)",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 14),  # 增大标题字体
    labels_gp = gpar(fontsize = 12),  # 增大标签字体
    legend_direction = "vertical",
    legend_height = unit(8, "cm"),  # 图例条高度，与热图主体高度相匹配
    grid_width = unit(1, "cm"),  # 增加图例条的实际宽度
    at = pretty(range(heatmap_data_log, na.rm = TRUE), n = 5),
    labels = sprintf("%.1f", pretty(range(heatmap_data_log, na.rm = TRUE), n = 5))
  )
)

# 创建温度注释 - 底部注释
bottom_ha <- HeatmapAnnotation(
  Temperature = temp_annotation$Temperature,
  col = list(Temperature = temp_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  height = unit(1, "cm")
)

# 创建温度图例
temp_legend <- Legend(
  title = "Tem/°C",
  labels = c("10", "16", "22"),
  legend_gp = gpar(fill = c("blue", "yellow", "red")),
  title_gp = gpar(fontsize = 14),
  labels_gp = gpar(fontsize = 12),
  grid_height = unit(0.8, "cm"),
  grid_width = unit(0.8, "cm"),  # 增加图例条的实际宽度
)

# 组合热图和注释
final_heatmap <- ht
final_heatmap <- final_heatmap %v% bottom_ha

# 保存为PDF
pdf_file <- "rca_pca_heatmap.pdf"
cairo_pdf(pdf_file, width = 12, height = 9, fallback_resolution = 1200)

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
     padding = unit(c(2, 2, 2, 1), "cm"))

dev.off()
cat("热图已成功保存到:", pdf_file, "\n")
