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
# 提取CG、CHG、CHH所有列
cg_cols <- grep('^CG_[0-9]+$', colnames(data), value = TRUE)
chg_cols <- grep('^CHG_[0-9]+$', colnames(data), value = TRUE)
chh_cols <- grep('^CHH_[0-9]+$', colnames(data), value = TRUE)

# 合并所有甲基化位点
all_meth_cols <- c(cg_cols, chg_cols, chh_cols)
meth_values <- as.matrix(data[, all_meth_cols])

# 执行PCA（基于第一个文件的方法）
feature_cols <- c("Tem..", "Lat", "Long", "alt.m", all_meth_cols)
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

# 创建合并的热图数据 - 所有甲基化位点为行，样本为列
heatmap_data <- t(meth_values[ordered_indices, ])
colnames(heatmap_data) <- ordered_indices
rownames(heatmap_data) <- all_meth_cols

# 创建温度分组注释
temp_groups <- data$Tem..[ordered_indices]
temp_annotation <- data.frame(
  Temperature = factor(temp_groups, levels = c("10", "16", "22"))
)
rownames(temp_annotation) <- colnames(heatmap_data)

# 定义颜色
temp_colors <- c("10" = "blue", "16" = "yellow", "22" = "red")

# 创建行分组注释（用于Y轴分组）
# 创建甲基化类型分组
meth_types <- c(rep("CG", length(cg_cols)), 
                rep("CHG", length(chg_cols)), 
                rep("CHH", length(chh_cols)))

# 创建行注释
row_annotation <- rowAnnotation(
  Type = factor(meth_types, levels = c("CG", "CHG", "CHH")),
  col = list(Type = c("CG" = "#E31A1C", "CHG" = "#33A02C", "CHH" = "#1F78B4")),
  width = unit(0.5, "cm"),
  show_annotation_name = FALSE,  # 不显示"Type"标签
  annotation_name_gp = gpar(fontsize = 12),
  show_legend = FALSE  # 不显示甲基化类型图例
)

# 使用ComplexHeatmap包创建热图
ht <- Heatmap(
  heatmap_data,
  name = "Methylation",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,  # 由于位点太多，不显示行名
  show_column_names = FALSE,
  left_annotation = row_annotation,
  col = colorRamp2(
    seq(min(heatmap_data, na.rm = TRUE), max(heatmap_data, na.rm = TRUE), length.out = 100),
    viridis(100)
  ),
  heatmap_legend_param = list(
    title = "5mC(%)",
    legend_direction = "vertical",
    title_position = "topleft",
    legend_height = unit(8, "cm"),
    grid_width = unit(1, "cm"),
    title_gp = gpar(fontsize = 14),
    labels_gp = gpar(fontsize = 12),
    at = pretty(range(heatmap_data, na.rm = TRUE), n = 5),
    labels = sprintf("%.1f", pretty(range(heatmap_data, na.rm = TRUE), n = 5))
  ),
  # 添加行分割线以区分不同的甲基化类型
  row_split = factor(meth_types, levels = c("CG", "CHG", "CHH")),
  row_gap = unit(2, "mm"),
  border = TRUE
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
  grid_width = unit(0.8, "cm")
)

# 组合热图和注释
final_heatmap <- ht %v% bottom_ha

# 保存为PDF
pdf_file <- "rca_combined_methylation_heatmap.pdf"
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
cat("合并甲基化热图已成功保存到:", pdf_file, "\n")
