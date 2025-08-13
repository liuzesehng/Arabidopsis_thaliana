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

# CHG_terminator-TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(CHG_terminator) & !is.na(TPM))

# 区间
start_val <- min(data_clean$CHG_terminator, na.rm = TRUE)
end_val <- ceiling(max(data_clean$CHG_terminator, na.rm = TRUE) * 10) / 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 0.05)
} else {
  breaks_seq <- seq(start_val, end_val, by = 0.05)
}
data_clean$CHG_terminator_group <- cut(
  data_clean$CHG_terminator,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$CHG_terminator_group)
data_clean$CHG_terminator_group_label <- break_labels[as.character(data_clean$CHG_terminator_group)]

# 计算每个CHG_terminator区间的中点，用于x轴位置
data_clean$CHG_terminator_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean$CHG_terminator_group)) + 0.025

# 进行方差分析和多重比较检验
aov_result <- aov(TPM ~ CHG_terminator_group_label, data = data_clean)
tukey_result <- TukeyHSD(aov_result)

# 提取p值并手动分配显著性字母
p_values <- tukey_result$CHG_terminator_group_label[,"p adj"]
group_names <- unique(data_clean$CHG_terminator_group_label)

significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("a", "a", "a", "a", "a", "a", "a", "a"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(CHG_terminator_group_label) %>%
  summarise(max_TPM = max(TPM, na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CHG_terminator_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(CHG_terminator_group_label, levels = break_labels), y = TPM)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CHG_terminator_group_label, levels = break_labels), 
                y = max_TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CHG Terminator",
    y = "TPM"
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
  scale_y_continuous(limits = c(0, max(data_clean$TPM, na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("TPM_vs_CHG_terminator_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean %>%
  group_by(CHG_terminator_group, CHG_terminator_group_label) %>%
  summarise(
    mean_TPM = mean(TPM, na.rm = TRUE),
    median_TPM = median(TPM, na.rm = TRUE),
    sd_TPM = sd(TPM, na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CHG_terminator_group_label) %>%
  left_join(significance_letters, by = c("CHG_terminator_group_label" = "group"))

print("各纬度区间的TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的纬度范围
cat("TPM原始数据纬度范围:", min(data_clean$Lat, na.rm = TRUE), "到", max(data_clean$Lat, na.rm = TRUE), "\n")
cat("TPM总数据点数:", nrow(data_clean), "\n")
cat("TPM纬度区间数:", length(unique(data_clean$Lat_group)), "\n")


# CG-TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(CG) & !is.na(TPM))

# 区间
start_val <- min(data_clean$CG, na.rm = TRUE)
end_val <- ceiling(max(data_clean$CG, na.rm = TRUE) * 10) / 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 0.05)
} else {
  breaks_seq <- seq(start_val, end_val, by = 0.05)
}
data_clean$CG_group <- cut(
  data_clean$CG,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$CG_group)
data_clean$CG_group_label <- break_labels[as.character(data_clean$CG_group)]

# 计算每个纬度区间的中点，用于x轴位置
data_clean$CG_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean$CG_group)) + 0.025

# 进行方差分析和多重比较检验
aov_result <- aov(TPM ~ CG_group_label, data = data_clean)
tukey_result <- TukeyHSD(aov_result)

# 提取p值并手动分配显著性字母
p_values <- tukey_result$CG_group_label[,"p adj"]
group_names <- unique(data_clean$CG_group_label)

significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("a", "b", "b", "ab", "ab", "ab"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(CG_group_label) %>%
  summarise(max_TPM = max(TPM, na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CG_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(CG_group_label, levels = break_labels), y = TPM)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CG_group_label, levels = break_labels), 
                y = max_TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CG",
    y = "Tpm"
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
  scale_y_continuous(limits = c(0, max(data_clean$TPM, na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("TPM_vs_CG.boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean %>%
  group_by(CG_group, CG_group_label) %>%
  summarise(
    mean_TPM = mean(TPM, na.rm = TRUE),
    median_TPM = median(TPM, na.rm = TRUE),
    sd_TPM = sd(TPM, na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CG_group_label) %>%
  left_join(significance_letters, by = c("CG_group_label" = "group"))

print("各CG区间的TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的纬度范围
cat("TPM原始数据纬度范围:", min(data_clean$CG, na.rm = TRUE), "到", max(data_clean$CG, na.rm = TRUE), "\n")
cat("TPM总数据点数:", nrow(data_clean), "\n")
cat("TPM纬度区间数:", length(unique(data_clean$CG_group)), "\n")


# CG_95_TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(CG_95) & !is.na(TPM))

# 区间
start_val <- min(data_clean$CG_95, na.rm = TRUE)
end_val <- ceiling(max(data_clean$CG_95, na.rm = TRUE) / 10) * 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 10)
} else {
  breaks_seq <- seq(start_val, end_val, by = 10)
}
data_clean$CG_95_group <- cut(
  data_clean$CG_95,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$CG_95_group)
data_clean$CG_95_group_label <- break_labels[as.character(data_clean$CG_95_group)]

# 计算每个CG_95区间的中点，用于x轴位置
data_clean$CG_95_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean$CG_95_group)) + 5

# 进行方差分析和多重比较检验
aov_result <- aov(TPM ~ CG_95_group_label, data = data_clean)
tukey_result <- TukeyHSD(aov_result)

# 提取p值并手动分配显著性字母
p_values <- tukey_result$CG_95_group_label[,"p adj"]
group_names <- unique(data_clean$CG_95_group_label)

significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("a", "a", "ab", "ab", "ab", "b", "ab", "ab", "ab", "ab"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(CG_95_group_label) %>%
  summarise(max_TPM = max(TPM, na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CG_95_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(CG_95_group_label, levels = break_labels), y = TPM)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CG_95_group_label, levels = break_labels), 
                y = max_TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CG_95",
    y = "TPM"
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
  scale_y_continuous(limits = c(0, max(data_clean$TPM, na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("TPM_vs_CG_95.box.plot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean %>%
  group_by(CG_95_group, CG_95_group_label) %>%
  summarise(
    mean_TPM = mean(TPM, na.rm = TRUE),
    median_TPM = median(TPM, na.rm = TRUE),
    sd_TPM = sd(TPM, na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CG_95_group_label) %>%
  left_join(significance_letters, by = c("CG_95_group_label" = "group"))

print("各CG_95区间的TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== CG_95_TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== CG_95_TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的CG_95范围
cat("CG_95_TPM原始数据范围:", min(data_clean$CG_95, na.rm = TRUE), "到", max(data_clean$CG_95, na.rm = TRUE), "\n")
cat("CG_95_TPM总数据点数:", nrow(data_clean), "\n")
cat("CG_95_TPM区间数:", length(unique(data_clean$CG_95_group)), "\n")


# CHH_terminator_TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(CHH_terminator) & !is.na(TPM))

# 区间
start_val <- min(data_clean$CHH_terminator, na.rm = TRUE)
end_val <- ceiling(max(data_clean$CHH_terminator, na.rm = TRUE) * 10) / 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 0.05)
} else {
  breaks_seq <- seq(start_val, end_val, by = 0.05)
}
data_clean$CHH_terminator_group <- cut(
  data_clean$CHH_terminator,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$CHH_terminator_group)
data_clean$CHH_terminator_group_label <- break_labels[as.character(data_clean$CHH_terminator_group)]

# 计算每个CHH_terminator区间的中点，用于x轴位置
data_clean$CHH_terminator_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean$CHH_terminator_group)) + 0.025

# 进行方差分析和多重比较检验
aov_result <- aov(TPM ~ CHH_terminator_group_label, data = data_clean)
tukey_result <- TukeyHSD(aov_result)

# 提取p值并手动分配显著性字母
p_values <- tukey_result$CHH_terminator_group_label[,"p adj"]
group_names <- unique(data_clean$CHH_terminator_group_label)

significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("b", "a", "ab", "ab", "ab", "ab", "ab"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(CHH_terminator_group_label) %>%
  summarise(max_TPM = max(TPM, na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CHH_terminator_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(CHH_terminator_group_label, levels = break_labels), y = TPM)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CHH_terminator_group_label, levels = break_labels), 
                y = max_TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CHH Terminator",
    y = "TPM"
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
  scale_y_continuous(limits = c(0, max(data_clean$TPM, na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("TPM_vs_CHH_terminator.boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean %>%
  group_by(CHH_terminator_group, CHH_terminator_group_label) %>%
  summarise(
    mean_TPM = mean(TPM, na.rm = TRUE),
    median_TPM = median(TPM, na.rm = TRUE),
    sd_TPM = sd(TPM, na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CHH_terminator_group_label) %>%
  left_join(significance_letters, by = c("CHH_terminator_group_label" = "group"))

print("各CHH_terminator区间的TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== CHH_terminator_TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== CHH_terminator_TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的CHH_terminator范围
cat("CHH_terminator_TPM原始数据范围:", min(data_clean$CHH_terminator, na.rm = TRUE), "到", max(data_clean$CHH_terminator, na.rm = TRUE), "\n")
cat("CHH_terminator_TPM总数据点数:", nrow(data_clean), "\n")
cat("CHH_terminator_TPM区间数:", length(unique(data_clean$CHH_terminator_group)), "\n")


# CHG_promoter_TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(CHG_promoter) & !is.na(TPM))

# 区间
start_val <- min(data_clean$CHG_promoter, na.rm = TRUE)
end_val <- ceiling(max(data_clean$CHG_promoter, na.rm = TRUE) * 100) / 100
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 0.005)
} else {
  breaks_seq <- seq(start_val, end_val, by = 0.005)
}
data_clean$CHG_promoter_group <- cut(
  data_clean$CHG_promoter,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$CHG_promoter_group)
data_clean$CHG_promoter_group_label <- break_labels[as.character(data_clean$CHG_promoter_group)]

# 计算每个CHG_promoter区间的中点，用于x轴位置
data_clean$CHG_promoter_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean$CHG_promoter_group)) + 0.0025

# 进行方差分析和多重比较检验
aov_result <- aov(TPM ~ CHG_promoter_group_label, data = data_clean)
tukey_result <- TukeyHSD(aov_result)

# 提取p值并手动分配显著性字母
p_values <- tukey_result$CHG_promoter_group_label[,"p adj"]
group_names <- unique(data_clean$CHG_promoter_group_label)

significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("a", "a", "a", "a"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(CHG_promoter_group_label) %>%
  summarise(max_TPM = max(TPM, na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CHG_promoter_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(CHG_promoter_group_label, levels = break_labels), y = TPM)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CHG_promoter_group_label, levels = break_labels), 
                y = max_TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CHG Promoter",
    y = "TPM"
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
  scale_y_continuous(limits = c(0, max(data_clean$TPM, na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("TPM_vs_CHG_promoter.boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean %>%
  group_by(CHG_promoter_group, CHG_promoter_group_label) %>%
  summarise(
    mean_TPM = mean(TPM, na.rm = TRUE),
    median_TPM = median(TPM, na.rm = TRUE),
    sd_TPM = sd(TPM, na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CHG_promoter_group_label) %>%
  left_join(significance_letters, by = c("CHG_promoter_group_label" = "group"))

print("各CHG_promoter区间的TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== CHG_promoter_TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== CHG_promoter_TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的CHG_promoter范围
cat("CHG_promoter_TPM原始数据范围:", min(data_clean$CHG_promoter, na.rm = TRUE), "到", max(data_clean$CHG_promoter, na.rm = TRUE), "\n")
cat("CHG_promoter_TPM总数据点数:", nrow(data_clean), "\n")
cat("CHG_promoter_TPM区间数:", length(unique(data_clean$CHG_promoter_group)), "\n")


# CHH_132_TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(CHH_132) & !is.na(TPM))

# 区间
start_val <- min(data_clean$CHH_132, na.rm = TRUE)
end_val <- ceiling(max(data_clean$CHH_132, na.rm = TRUE) / 10) * 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 5)
} else {
  breaks_seq <- seq(start_val, end_val, by = 5)
}
data_clean$CHH_132_group <- cut(
  data_clean$CHH_132,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$CHH_132_group)
data_clean$CHH_132_group_label <- break_labels[as.character(data_clean$CHH_132_group)]

# 计算每个CHH_132区间的中点，用于x轴位置
data_clean$CHH_132_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean$CHH_132_group)) + 5

# 进行方差分析和多重比较检验
aov_result <- aov(TPM ~ CHH_132_group_label, data = data_clean)
tukey_result <- TukeyHSD(aov_result)

# 提取p值并手动分配显著性字母
p_values <- tukey_result$CHH_132_group_label[,"p adj"]
group_names <- unique(data_clean$CHH_132_group_label)

significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("b", "a", "ab", "ab", "ab", "ab", "ab"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(CHH_132_group_label) %>%
  summarise(max_TPM = max(TPM, na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CHH_132_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(CHH_132_group_label, levels = break_labels), y = TPM)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CHH_132_group_label, levels = break_labels), 
                y = max_TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CHH_132",
    y = "TPM"
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
  scale_y_continuous(limits = c(0, max(data_clean$TPM, na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("TPM_vs_CHH_132.boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean %>%
  group_by(CHH_132_group, CHH_132_group_label) %>%
  summarise(
    mean_TPM = mean(TPM, na.rm = TRUE),
    median_TPM = median(TPM, na.rm = TRUE),
    sd_TPM = sd(TPM, na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CHH_132_group_label) %>%
  left_join(significance_letters, by = c("CHH_132_group_label" = "group"))

print("各CHH_132区间的TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== CHH_132_TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== CHH_132_TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的CHH_132范围
cat("CHH_132_TPM原始数据范围:", min(data_clean$CHH_132, na.rm = TRUE), "到", max(data_clean$CHH_132, na.rm = TRUE), "\n")
cat("CHH_132_TPM总数据点数:", nrow(data_clean), "\n")
cat("CHH_132_TPM区间数:", length(unique(data_clean$CHH_132_group)), "\n")


