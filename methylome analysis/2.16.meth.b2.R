# 加载必要的包
library(ggplot2)
library(dplyr)
library(ggsignif)  # 添加这个包用于显著性标记
library(multcompView)  # 添加用于自动生成显著性字母的包

# 读取数据
data <- read.table("Alt.RCA.me_snp.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

# 查看数据结构
str(data)
head(data)


# CHH_terminator-β2TPM列
# 清理数据，去除缺失值
data_clean_chh_terminator <- data %>%
  filter(!is.na(CHH_terminator) & !is.na(β2_TPM..))

# 区间
start_val <- min(data_clean_chh_terminator$CHH_terminator, na.rm = TRUE)
end_val <- ceiling(max(data_clean_chh_terminator$CHH_terminator, na.rm = TRUE) * 10) / 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 0.05)
} else {
  breaks_seq <- seq(start_val, end_val, by = 0.05)
}
data_clean_chh_terminator$CHH_terminator_group <- cut(
  data_clean_chh_terminator$CHH_terminator,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_chh_terminator$CHH_terminator_group)
data_clean_chh_terminator$CHH_terminator_group_label <- break_labels[as.character(data_clean_chh_terminator$CHH_terminator_group)]

# 计算每个CHH_terminator区间的中点，用于x轴位置
data_clean_chh_terminator$CHH_terminator_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_chh_terminator$CHH_terminator_group)) + 0.025

# 进行方差分析和多重比较检验
aov_result <- aov(β2_TPM.. ~ CHH_terminator_group_label, data = data_clean_chh_terminator)
tukey_result <- TukeyHSD(aov_result)

# 提取p值并手动分配显著性字母
p_values <- tukey_result$CHH_terminator_group_label[,"p adj"]
group_names <- unique(data_clean_chh_terminator$CHH_terminator_group_label)

significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("a", "b", "b", "ab", "ab", "ab", "ab"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean_chh_terminator %>%
  group_by(CHH_terminator_group_label) %>%
  summarise(max_β2_TPM.. = max(β2_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CHH_terminator_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean_chh_terminator, aes(x = factor(CHH_terminator_group_label, levels = break_labels), y = β2_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CHH_terminator_group_label, levels = break_labels), 
                y = max_β2_TPM.. * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CHH Terminator",
    y = "b2TPM"
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    # 显示完整的边框
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    # 设置刻度线朝外
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    # 确保刻度线朝外
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5)
  ) +
  # 设置y轴范围，为字母标记留出空间
  scale_y_continuous(limits = c(0, max(data_clean_chh_terminator$β2_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_CHH_terminator_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean_chh_terminator %>%
  group_by(CHH_terminator_group, CHH_terminator_group_label) %>%
  summarise(
    mean_β2_TPM.. = mean(β2_TPM.., na.rm = TRUE),
    median_β2_TPM.. = median(β2_TPM.., na.rm = TRUE),
    sd_β2_TPM.. = sd(β2_TPM.., na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CHH_terminator_group_label) %>%
  left_join(significance_letters, by = c("CHH_terminator_group_label" = "group"))

print("各纬度区间的β2_TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== β2_TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== β2_TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的纬度范围
cat("β2_TPM原始数据纬度范围:", min(data_clean_chh_terminator$Lat, na.rm = TRUE), "到", max(data_clean_chh_terminator$Lat, na.rm = TRUE), "\n")
cat("β2_TPM总数据点数:", nrow(data_clean_chh_terminator), "\n")
cat("β2_TPM纬度区间数:", length(unique(data_clean_chh_terminator$Lat_group)), "\n")


# CHH-β2TPM列
# 清理数据，去除缺失值
data_clean_chh <- data %>%
  filter(!is.na(CHH) & !is.na(β2_TPM..))

# 区间
start_val <- min(data_clean_chh$CHH, na.rm = TRUE)
end_val <- ceiling(max(data_clean_chh$CHH, na.rm = TRUE) * 100) / 100
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 0.001)
} else {
  breaks_seq <- seq(start_val, end_val, by = 0.001)
}
data_clean_chh$CHH_group <- cut(
  data_clean_chh$CHH,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_chh$CHH_group)
data_clean_chh$CHH_group_label <- break_labels[as.character(data_clean_chh$CHH_group)]

# 计算每个CHH区间的中点，用于x轴位置
data_clean_chh$CHH_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_chh$CHH_group)) + 0.0005

# 进行方差分析和多重比较检验
aov_result <- aov(β2_TPM.. ~ CHH_group_label, data = data_clean_chh)
tukey_result <- TukeyHSD(aov_result)

# 提取p值并手动分配显著性字母
p_values <- tukey_result$CHH_group_label[,"p adj"]
group_names <- unique(data_clean_chh$CHH_group_label)

significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("ab", "ab", "a", "b", "b", "b", "b", "ab", "ab", "ab", "ab"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean_chh %>%
  group_by(CHH_group_label) %>%
  summarise(max_β2_TPM.. = max(β2_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CHH_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean_chh, aes(x = factor(CHH_group_label, levels = break_labels), y = β2_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CHH_group_label, levels = break_labels), 
                y = max_β2_TPM.. * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CHH",
    y = "b2TPM"
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    # 显示完整的边框
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    # 设置刻度线朝外
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    # 确保刻度线朝外
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5)
  ) +
  # 设置y轴范围，为字母标记留出空间
  scale_y_continuous(limits = c(0, max(data_clean_chh$β2_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_CHH_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean_chh %>%
  group_by(CHH_group, CHH_group_label) %>%
  summarise(
    mean_β2_TPM.. = mean(β2_TPM.., na.rm = TRUE),
    median_β2_TPM.. = median(β2_TPM.., na.rm = TRUE),
    sd_β2_TPM.. = sd(β2_TPM.., na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CHH_group_label) %>%
  left_join(significance_letters, by = c("CHH_group_label" = "group"))

print("各纬度区间的β2_TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== β2_TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== β2_TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的纬度范围
cat("β2_TPM原始数据纬度范围:", min(data_clean_chh$Lat, na.rm = TRUE), "到", max(data_clean_chh$Lat, na.rm = TRUE), "\n")
cat("β2_TPM总数据点数:", nrow(data_clean_chh), "\n")
cat("β2_TPM纬度区间数:", length(unique(data_clean_chh$Lat_group)), "\n")


# CG_promoter-β2TPM列
# 清理数据，去除缺失值
data_clean_cg_promoter <- data %>%
  filter(!is.na(CG_promoter) & !is.na(β2_TPM..))

# 区间
start_val <- min(data_clean_cg_promoter$CG_promoter, na.rm = TRUE)
end_val <- ceiling(max(data_clean_cg_promoter$CG_promoter, na.rm = TRUE) * 10) / 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 0.05)
} else {
  breaks_seq <- seq(start_val, end_val, by = 0.05)
}
data_clean_cg_promoter$CG_promoter_group <- cut(
  data_clean_cg_promoter$CG_promoter,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_cg_promoter$CG_promoter_group)
data_clean_cg_promoter$CG_promoter_group_label <- break_labels[as.character(data_clean_cg_promoter$CG_promoter_group)]

# 计算每个CG_promoter区间的中点，用于x轴位置
data_clean_cg_promoter$CG_promoter_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_cg_promoter$CG_promoter_group)) + 0.025

# 进行方差分析和多重比较检验
aov_result <- aov(β2_TPM.. ~ CG_promoter_group_label, data = data_clean_cg_promoter)
tukey_result <- TukeyHSD(aov_result)

# 提取p值并手动分配显著性字母
p_values <- tukey_result$CG_promoter_group_label[,"p adj"]
group_names <- unique(data_clean_cg_promoter$CG_promoter_group_label)

significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("c", "c", "b", "a", "a", "a", "abc", "a", "abc", "c"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean_cg_promoter %>%
  group_by(CG_promoter_group_label) %>%
  summarise(max_β2_TPM.. = max(β2_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CG_promoter_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean_cg_promoter, aes(x = factor(CG_promoter_group_label, levels = break_labels), y = β2_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CG_promoter_group_label, levels = break_labels), 
                y = max_β2_TPM.. * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CG Promoter",
    y = "b2TPM.."
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    # 显示完整的边框
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    # 设置刻度线朝外
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    # 确保刻度线朝外
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5)
  ) +
  # 设置y轴范围，为字母标记留出空间
  scale_y_continuous(limits = c(0, max(data_clean_cg_promoter$β2_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_CG_promoter_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean_cg_promoter %>%
  group_by(CG_promoter_group, CG_promoter_group_label) %>%
  summarise(
    mean_β2_TPM.. = mean(β2_TPM.., na.rm = TRUE),
    median_β2_TPM.. = median(β2_TPM.., na.rm = TRUE),
    sd_β2_TPM.. = sd(β2_TPM.., na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CG_promoter_group_label) %>%
  left_join(significance_letters, by = c("CG_promoter_group_label" = "group"))

print("各纬度区间的β2_TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== β2_TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== β2_TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的纬度范围
cat("β2_TPM原始数据纬度范围:", min(data_clean_cg_promoter$Lat, na.rm = TRUE), "到", max(data_clean_cg_promoter$Lat, na.rm = TRUE), "\n")
cat("β2_TPM总数据点数:", nrow(data_clean_cg_promoter), "\n")
cat("β2_TPM纬度区间数:", length(unique(data_clean_cg_promoter$Lat_group)), "\n")


# CG_terminator-β2TPM列
# 清理数据，去除缺失值
data_clean_cg_terminator <- data %>%
  filter(!is.na(CG_terminator) & !is.na(β2_TPM..))

# 区间
start_val <- min(data_clean_cg_terminator$CG_terminator, na.rm = TRUE)
end_val <- ceiling(max(data_clean_cg_terminator$CG_terminator, na.rm = TRUE) * 10) / 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 0.05)
} else {
  breaks_seq <- seq(start_val, end_val, by = 0.05)
}
data_clean_cg_terminator$CG_terminator_group <- cut(
  data_clean_cg_terminator$CG_terminator,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_cg_terminator$CG_terminator_group)
data_clean_cg_terminator$CG_terminator_group_label <- break_labels[as.character(data_clean_cg_terminator$CG_terminator_group)]

# 计算每个CG_terminator区间的中点，用于x轴位置
data_clean_cg_terminator$CG_terminator_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_cg_terminator$CG_terminator_group)) + 0.025

# 进行方差分析和多重比较检验
aov_result <- aov(β2_TPM.. ~ CG_terminator_group_label, data = data_clean_cg_terminator)
tukey_result <- TukeyHSD(aov_result)

# 提取p值并手动分配显著性字母
p_values <- tukey_result$CG_terminator_group_label[,"p adj"]
group_names <- unique(data_clean_cg_terminator$CG_terminator_group_label)

significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("a", "a", "b", "c", "c", "c", "a", "abc", "abc"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean_cg_terminator %>%
  group_by(CG_terminator_group_label) %>%
  summarise(max_β2_TPM.. = max(β2_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CG_terminator_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean_cg_terminator, aes(x = factor(CG_terminator_group_label, levels = break_labels), y = β2_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CG_terminator_group_label, levels = break_labels), 
                y = max_β2_TPM.. * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CG Terminator",
    y = "b2TPM"
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    # 显示完整的边框
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    # 设置刻度线朝外
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    # 确保刻度线朝外
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5)
  ) +
  # 设置y轴范围，为字母标记留出空间
  scale_y_continuous(limits = c(0, max(data_clean_cg_terminator$β2_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_CG_terminator_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean_cg_terminator %>%
  group_by(CG_terminator_group, CG_terminator_group_label) %>%
  summarise(
    mean_β2_TPM.. = mean(β2_TPM.., na.rm = TRUE),
    median_β2_TPM.. = median(β2_TPM.., na.rm = TRUE),
    sd_β2_TPM.. = sd(β2_TPM.., na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CG_terminator_group_label) %>%
  left_join(significance_letters, by = c("CG_terminator_group_label" = "group"))

print("各纬度区间的β2_TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== β2_TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== β2_TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的纬度范围
cat("β2_TPM原始数据纬度范围:", min(data_clean_cg_terminator$Lat, na.rm = TRUE), "到", max(data_clean_cg_terminator$Lat, na.rm = TRUE), "\n")
cat("β2_TPM总数据点数:", nrow(data_clean_cg_terminator), "\n")
cat("β2_TPM纬度区间数:", length(unique(data_clean_cg_terminator$Lat_group)), "\n")


# CG_116-β2TPM列
# 清理数据，去除缺失值
data_clean_cg_116 <- data %>%
  filter(!is.na(CG_116) & !is.na(β2_TPM..))

# 区间
start_val <- min(data_clean_cg_116$CG_116, na.rm = TRUE)
end_val <- ceiling(max(data_clean_cg_116$CG_116, na.rm = TRUE) / 10) * 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 10)
} else {
  breaks_seq <- seq(start_val, end_val, by = 10)
}
data_clean_cg_116$CG_116_group <- cut(
  data_clean_cg_116$CG_116,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_cg_116$CG_116_group)
data_clean_cg_116$CG_116_group_label <- break_labels[as.character(data_clean_cg_116$CG_116_group)]

# 计算每个CG_116区间的中点，用于x轴位置
data_clean_cg_116$CG_116_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_cg_116$CG_116_group)) + 5

# 进行方差分析和多重比较检验
aov_result <- aov(β2_TPM.. ~ CG_116_group_label, data = data_clean_cg_116)
tukey_result <- TukeyHSD(aov_result)

# 提取p值并手动分配显著性字母
p_values <- tukey_result$CG_116_group_label[,"p adj"]
group_names <- unique(data_clean_cg_116$CG_116_group_label)

significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("a", "a", "a", "a", "a", "a", "a", "a", "a"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean_cg_116 %>%
  group_by(CG_116_group_label) %>%
  summarise(max_β2_TPM.. = max(β2_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CG_116_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean_cg_116, aes(x = factor(CG_116_group_label, levels = break_labels), y = β2_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CG_116_group_label, levels = break_labels), 
                y = max_β2_TPM.. * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CG 116",
    y = "b2TPM"
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    # 显示完整的边框
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    # 设置刻度线朝外
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    # 确保刻度线朝外
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5)
  ) +
  # 设置y轴范围，为字母标记留出空间
  scale_y_continuous(limits = c(0, max(data_clean_cg_116$β2_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_CG_116_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean_cg_116 %>%
  group_by(CG_116_group, CG_116_group_label) %>%
  summarise(
    mean_β2_TPM.. = mean(β2_TPM.., na.rm = TRUE),
    median_β2_TPM.. = median(β2_TPM.., na.rm = TRUE),
    sd_β2_TPM.. = sd(β2_TPM.., na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CG_116_group_label) %>%
  left_join(significance_letters, by = c("CG_116_group_label" = "group"))

print("各纬度区间的β2_TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== β2_TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== β2_TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的纬度范围
cat("β2_TPM原始数据纬度范围:", min(data_clean_cg_116$Lat, na.rm = TRUE), "到", max(data_clean_cg_116$Lat, na.rm = TRUE), "\n")
cat("β2_TPM总数据点数:", nrow(data_clean_cg_116), "\n")
cat("β2_TPM纬度区间数:", length(unique(data_clean_cg_116$Lat_group)), "\n")


# CG_44-β2TPM列
# 清理数据，去除缺失值
data_clean_cg_44 <- data %>%
  filter(!is.na(CG_44) & !is.na(β2_TPM..))

# 区间
start_val <- min(data_clean_cg_44$CG_44, na.rm = TRUE)
end_val <- ceiling(max(data_clean_cg_44$CG_44, na.rm = TRUE) / 10) * 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 10)
} else {
  breaks_seq <- seq(start_val, end_val, by = 10)
}
data_clean_cg_44$CG_44_group <- cut(
  data_clean_cg_44$CG_44,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_cg_44$CG_44_group)
data_clean_cg_44$CG_44_group_label <- break_labels[as.character(data_clean_cg_44$CG_44_group)]

# 计算每个CG_44区间的中点，用于x轴位置
data_clean_cg_44$CG_44_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_cg_44$CG_44_group)) + 5

# 进行方差分析和多重比较检验
aov_result <- aov(β2_TPM.. ~ CG_44_group_label, data = data_clean_cg_44)
tukey_result <- TukeyHSD(aov_result)

# 提取p值并手动分配显著性字母
p_values <- tukey_result$CG_44_group_label[,"p adj"]
group_names <- unique(data_clean_cg_44$CG_44_group_label)

significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("b", "ab", "b", "b", "ab", "ab", "a", "ab", "a", "a"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean_cg_44 %>%
  group_by(CG_44_group_label) %>%
  summarise(max_β2_TPM.. = max(β2_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CG_44_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean_cg_44, aes(x = factor(CG_44_group_label, levels = break_labels), y = β2_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CG_44_group_label, levels = break_labels), 
                y = max_β2_TPM.. * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CG 44",
    y = "b2TPM"
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    # 显示完整的边框
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    # 设置刻度线朝外
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    # 确保刻度线朝外
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5)
  ) +
  # 设置y轴范围，为字母标记留出空间
  scale_y_continuous(limits = c(0, max(data_clean_cg_44$β2_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_CG_44_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean_cg_44 %>%
  group_by(CG_44_group, CG_44_group_label) %>%
  summarise(
    mean_β2_TPM.. = mean(β2_TPM.., na.rm = TRUE),
    median_β2_TPM.. = median(β2_TPM.., na.rm = TRUE),
    sd_β2_TPM.. = sd(β2_TPM.., na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CG_44_group_label) %>%
  left_join(significance_letters, by = c("CG_44_group_label" = "group"))

print("各纬度区间的β2_TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== β2_TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== β2_TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的纬度范围
cat("β2_TPM原始数据纬度范围:", min(data_clean_cg_44$Lat, na.rm = TRUE), "到", max(data_clean_cg_44$Lat, na.rm = TRUE), "\n")
cat("β2_TPM总数据点数:", nrow(data_clean_cg_44), "\n")
cat("β2_TPM纬度区间数:", length(unique(data_clean_cg_44$Lat_group)), "\n")

