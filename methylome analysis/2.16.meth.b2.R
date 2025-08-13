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


# ===== CHH_terminator_β2TPM 分析 =====
# 数据预处理和分组
data_clean <- data %>%
  filter(!is.na(CHH_terminator) & !is.na(β2TPM) & is.finite(CHH_terminator) & is.finite(β2TPM))

# 创建 CHH_terminator 的分组
breaks <- quantile(data_clean$CHH_terminator, probs = seq(0, 1, by = 1/8), na.rm = TRUE)
data_clean$CHH_terminator_group <- cut(data_clean$CHH_terminator, breaks = breaks, include.lowest = TRUE)

# 创建分组标签
break_labels <- paste0("[", round(breaks[-length(breaks)], 3), ", ", round(breaks[-1], 3), "]")
break_labels[length(break_labels)] <- paste0("[", round(breaks[length(breaks)-1], 3), ", ", round(breaks[length(breaks)], 3), "]")
data_clean$CHH_terminator_group_label <- cut(data_clean$CHH_terminator, breaks = breaks, labels = break_labels, include.lowest = TRUE)

# 进行方差分析
aov_result <- aov(β2TPM ~ CHH_terminator_group_label, data = data_clean)

# 进行Tukey多重比较检验
tukey_result <- TukeyHSD(aov_result)

# 获取分组名称
group_names <- levels(data_clean$CHH_terminator_group_label)

# 为每个分组分配显著性字母（这里需要根据Tukey检验结果手动调整）
significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("a", "a", "a", "a", "a", "a", "a", "a"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(CHH_terminator_group_label) %>%
  summarise(max_b2TPM = max(β2TPM, na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CHH_terminator_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(CHH_terminator_group_label, levels = break_labels), y = β2TPM)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CHH_terminator_group_label, levels = break_labels), 
                y = max_b2TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CHH Terminator",
    y = "β2TPM"
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
    # 设置坐标轴刻度线
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5)
  ) +
  # 设置y轴范围，为字母标记留出空间
  scale_y_continuous(limits = c(0, max(data_clean$β2TPM, na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_CHH_terminator_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean %>%
  group_by(CHH_terminator_group, CHH_terminator_group_label) %>%
  summarise(
    mean_b2TPM = mean(β2TPM, na.rm = TRUE),
    median_b2TPM = median(β2TPM, na.rm = TRUE),
    sd_b2TPM = sd(β2TPM, na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CHH_terminator_group_label) %>%
  left_join(significance_letters, by = c("CHH_terminator_group_label" = "group"))

print("各纬度区间的β2TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== β2TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== β2TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的纬度范围
cat("β2TPM原始数据纬度范围:", min(data_clean$Lat, na.rm = TRUE), "到", max(data_clean$Lat, na.rm = TRUE), "\n")
cat("β2TPM总数据点数:", nrow(data_clean), "\n")
cat("β2TPM纬度区间数:", length(unique(data_clean$Lat_group)), "\n")


# ===== CHH_β2TPM 分析 =====
# 数据预处理和分组
data_clean <- data %>%
  filter(!is.na(CHH) & !is.na(β2TPM) & is.finite(CHH) & is.finite(β2TPM))

# 创建 CHH 的分组
breaks <- quantile(data_clean$CHH, probs = seq(0, 1, by = 1/8), na.rm = TRUE)
data_clean$CHH_group <- cut(data_clean$CHH, breaks = breaks, include.lowest = TRUE)

# 创建分组标签
break_labels <- paste0("[", round(breaks[-length(breaks)], 3), ", ", round(breaks[-1], 3), "]")
break_labels[length(break_labels)] <- paste0("[", round(breaks[length(breaks)-1], 3), ", ", round(breaks[length(breaks)], 3), "]")
data_clean$CHH_group_label <- cut(data_clean$CHH, breaks = breaks, labels = break_labels, include.lowest = TRUE)

# 进行方差分析
aov_result <- aov(β2TPM ~ CHH_group_label, data = data_clean)

# 进行Tukey多重比较检验
tukey_result <- TukeyHSD(aov_result)

# 获取分组名称
group_names <- levels(data_clean$CHH_group_label)

# 为每个分组分配显著性字母（这里需要根据Tukey检验结果手动调整）
significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("a", "a", "a", "a", "a", "a", "a", "a"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(CHH_group_label) %>%
  summarise(max_b2TPM = max(β2TPM, na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CHH_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(CHH_group_label, levels = break_labels), y = β2TPM)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CHH_group_label, levels = break_labels), 
                y = max_b2TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CHH",
    y = "β2TPM"
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
    # 设置坐标轴刻度线
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5)
  ) +
  # 设置y轴范围，为字母标记留出空间
  scale_y_continuous(limits = c(0, max(data_clean$β2TPM, na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_CHH_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean %>%
  group_by(CHH_group, CHH_group_label) %>%
  summarise(
    mean_b2TPM = mean(β2TPM, na.rm = TRUE),
    median_b2TPM = median(β2TPM, na.rm = TRUE),
    sd_b2TPM = sd(β2TPM, na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CHH_group_label) %>%
  left_join(significance_letters, by = c("CHH_group_label" = "group"))

print("各纬度区间的β2TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== β2TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== β2TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的纬度范围
cat("β2TPM原始数据纬度范围:", min(data_clean$Lat, na.rm = TRUE), "到", max(data_clean$Lat, na.rm = TRUE), "\n")
cat("β2TPM总数据点数:", nrow(data_clean), "\n")
cat("β2TPM纬度区间数:", length(unique(data_clean$Lat_group)), "\n")


# ===== CG_promoter_β2TPM 分析 =====
# 数据预处理和分组
data_clean <- data %>%
  filter(!is.na(CG_promoter) & !is.na(β2TPM) & is.finite(CG_promoter) & is.finite(β2TPM))

# 创建 CG_promoter 的分组
breaks <- quantile(data_clean$CG_promoter, probs = seq(0, 1, by = 1/8), na.rm = TRUE)
data_clean$CG_promoter_group <- cut(data_clean$CG_promoter, breaks = breaks, include.lowest = TRUE)

# 创建分组标签
break_labels <- paste0("[", round(breaks[-length(breaks)], 3), ", ", round(breaks[-1], 3), "]")
break_labels[length(break_labels)] <- paste0("[", round(breaks[length(breaks)-1], 3), ", ", round(breaks[length(breaks)], 3), "]")
data_clean$CG_promoter_group_label <- cut(data_clean$CG_promoter, breaks = breaks, labels = break_labels, include.lowest = TRUE)

# 进行方差分析
aov_result <- aov(β2TPM ~ CG_promoter_group_label, data = data_clean)

# 进行Tukey多重比较检验
tukey_result <- TukeyHSD(aov_result)

# 获取分组名称
group_names <- levels(data_clean$CG_promoter_group_label)

# 为每个分组分配显著性字母（这里需要根据Tukey检验结果手动调整）
significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("a", "a", "a", "a", "a", "a", "a", "a"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(CG_promoter_group_label) %>%
  summarise(max_b2TPM = max(β2TPM, na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CG_promoter_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(CG_promoter_group_label, levels = break_labels), y = β2TPM)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CG_promoter_group_label, levels = break_labels), 
                y = max_b2TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CG Promoter",
    y = "β2TPM"
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
    # 设置坐标轴刻度线
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5)
  ) +
  # 设置y轴范围，为字母标记留出空间
  scale_y_continuous(limits = c(0, max(data_clean$β2TPM, na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_CG_promoter_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean %>%
  group_by(CG_promoter_group, CG_promoter_group_label) %>%
  summarise(
    mean_b2TPM = mean(β2TPM, na.rm = TRUE),
    median_b2TPM = median(β2TPM, na.rm = TRUE),
    sd_b2TPM = sd(β2TPM, na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CG_promoter_group_label) %>%
  left_join(significance_letters, by = c("CG_promoter_group_label" = "group"))

print("各纬度区间的β2TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== β2TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== β2TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的纬度范围
cat("β2TPM原始数据纬度范围:", min(data_clean$Lat, na.rm = TRUE), "到", max(data_clean$Lat, na.rm = TRUE), "\n")
cat("β2TPM总数据点数:", nrow(data_clean), "\n")
cat("β2TPM纬度区间数:", length(unique(data_clean$Lat_group)), "\n")


# ===== CG_terminator_β2TPM 分析 =====
# 数据预处理和分组
data_clean <- data %>%
  filter(!is.na(CG_terminator) & !is.na(β2TPM) & is.finite(CG_terminator) & is.finite(β2TPM))

# 创建 CG_terminator 的分组
breaks <- quantile(data_clean$CG_terminator, probs = seq(0, 1, by = 1/8), na.rm = TRUE)
data_clean$CG_terminator_group <- cut(data_clean$CG_terminator, breaks = breaks, include.lowest = TRUE)

# 创建分组标签
break_labels <- paste0("[", round(breaks[-length(breaks)], 3), ", ", round(breaks[-1], 3), "]")
break_labels[length(break_labels)] <- paste0("[", round(breaks[length(breaks)-1], 3), ", ", round(breaks[length(breaks)], 3), "]")
data_clean$CG_terminator_group_label <- cut(data_clean$CG_terminator, breaks = breaks, labels = break_labels, include.lowest = TRUE)

# 进行方差分析
aov_result <- aov(β2TPM ~ CG_terminator_group_label, data = data_clean)

# 进行Tukey多重比较检验
tukey_result <- TukeyHSD(aov_result)

# 获取分组名称
group_names <- levels(data_clean$CG_terminator_group_label)

# 为每个分组分配显著性字母（这里需要根据Tukey检验结果手动调整）
significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("a", "a", "a", "a", "a", "a", "a", "a"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(CG_terminator_group_label) %>%
  summarise(max_b2TPM = max(β2TPM, na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CG_terminator_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(CG_terminator_group_label, levels = break_labels), y = β2TPM)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CG_terminator_group_label, levels = break_labels), 
                y = max_b2TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CG Terminator",
    y = "β2TPM"
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
    # 设置坐标轴刻度线
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5)
  ) +
  # 设置y轴范围，为字母标记留出空间
  scale_y_continuous(limits = c(0, max(data_clean$β2TPM, na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_CG_terminator_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean %>%
  group_by(CG_terminator_group, CG_terminator_group_label) %>%
  summarise(
    mean_b2TPM = mean(β2TPM, na.rm = TRUE),
    median_b2TPM = median(β2TPM, na.rm = TRUE),
    sd_b2TPM = sd(β2TPM, na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CG_terminator_group_label) %>%
  left_join(significance_letters, by = c("CG_terminator_group_label" = "group"))

print("各纬度区间的β2TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== β2TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== β2TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的纬度范围
cat("β2TPM原始数据纬度范围:", min(data_clean$Lat, na.rm = TRUE), "到", max(data_clean$Lat, na.rm = TRUE), "\n")
cat("β2TPM总数据点数:", nrow(data_clean), "\n")
cat("β2TPM纬度区间数:", length(unique(data_clean$Lat_group)), "\n")


# ===== CG_116_β2TPM 分析 =====
# 数据预处理和分组
data_clean <- data %>%
  filter(!is.na(CG_116) & !is.na(β2TPM) & is.finite(CG_116) & is.finite(β2TPM))

# 创建 CG_116 的分组
breaks <- quantile(data_clean$CG_116, probs = seq(0, 1, by = 1/8), na.rm = TRUE)
data_clean$CG_116_group <- cut(data_clean$CG_116, breaks = breaks, include.lowest = TRUE)

# 创建分组标签
break_labels <- paste0("[", round(breaks[-length(breaks)], 3), ", ", round(breaks[-1], 3), "]")
break_labels[length(break_labels)] <- paste0("[", round(breaks[length(breaks)-1], 3), ", ", round(breaks[length(breaks)], 3), "]")
data_clean$CG_116_group_label <- cut(data_clean$CG_116, breaks = breaks, labels = break_labels, include.lowest = TRUE)

# 进行方差分析
aov_result <- aov(β2TPM ~ CG_116_group_label, data = data_clean)

# 进行Tukey多重比较检验
tukey_result <- TukeyHSD(aov_result)

# 获取分组名称
group_names <- levels(data_clean$CG_116_group_label)

# 为每个分组分配显著性字母（这里需要根据Tukey检验结果手动调整）
significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("a", "a", "a", "a", "a", "a", "a", "a"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(CG_116_group_label) %>%
  summarise(max_b2TPM = max(β2TPM, na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CG_116_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(CG_116_group_label, levels = break_labels), y = β2TPM)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CG_116_group_label, levels = break_labels), 
                y = max_b2TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CG 116",
    y = "β2TPM"
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
    # 设置坐标轴刻度线
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5)
  ) +
  # 设置y轴范围，为字母标记留出空间
  scale_y_continuous(limits = c(0, max(data_clean$β2TPM, na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_CG_116_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean %>%
  group_by(CG_116_group, CG_116_group_label) %>%
  summarise(
    mean_b2TPM = mean(β2TPM, na.rm = TRUE),
    median_b2TPM = median(β2TPM, na.rm = TRUE),
    sd_b2TPM = sd(β2TPM, na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CG_116_group_label) %>%
  left_join(significance_letters, by = c("CG_116_group_label" = "group"))

print("各纬度区间的β2TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== β2TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== β2TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的纬度范围
cat("β2TPM原始数据纬度范围:", min(data_clean$Lat, na.rm = TRUE), "到", max(data_clean$Lat, na.rm = TRUE), "\n")
cat("β2TPM总数据点数:", nrow(data_clean), "\n")
cat("β2TPM纬度区间数:", length(unique(data_clean$Lat_group)), "\n")


# ===== CG_44_β2TPM 分析 =====
# 数据预处理和分组
data_clean <- data %>%
  filter(!is.na(CG_44) & !is.na(β2TPM) & is.finite(CG_44) & is.finite(β2TPM))

# 创建 CG_44 的分组
breaks <- quantile(data_clean$CG_44, probs = seq(0, 1, by = 1/8), na.rm = TRUE)
data_clean$CG_44_group <- cut(data_clean$CG_44, breaks = breaks, include.lowest = TRUE)

# 创建分组标签
break_labels <- paste0("[", round(breaks[-length(breaks)], 3), ", ", round(breaks[-1], 3), "]")
break_labels[length(break_labels)] <- paste0("[", round(breaks[length(breaks)-1], 3), ", ", round(breaks[length(breaks)], 3), "]")
data_clean$CG_44_group_label <- cut(data_clean$CG_44, breaks = breaks, labels = break_labels, include.lowest = TRUE)

# 进行方差分析
aov_result <- aov(β2TPM ~ CG_44_group_label, data = data_clean)

# 进行Tukey多重比较检验
tukey_result <- TukeyHSD(aov_result)

# 获取分组名称
group_names <- levels(data_clean$CG_44_group_label)

# 为每个分组分配显著性字母（这里需要根据Tukey检验结果手动调整）
significance_letters <- data.frame(
  group = sort(group_names),
  letter = c("a", "a", "a", "a", "a", "a", "a", "a"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(CG_44_group_label) %>%
  summarise(max_b2TPM = max(β2TPM, na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("CG_44_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(CG_44_group_label, levels = break_labels), y = β2TPM)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(CG_44_group_label, levels = break_labels), 
                y = max_b2TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CG 44",
    y = "β2TPM"
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
    # 设置坐标轴刻度线
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5)
  ) +
  # 设置y轴范围，为字母标记留出空间
  scale_y_continuous(limits = c(0, max(data_clean$β2TPM, na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_CG_44_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary <- data_clean %>%
  group_by(CG_44_group, CG_44_group_label) %>%
  summarise(
    mean_b2TPM = mean(β2TPM, na.rm = TRUE),
    median_b2TPM = median(β2TPM, na.rm = TRUE),
    sd_b2TPM = sd(β2TPM, na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CG_44_group_label) %>%
  left_join(significance_letters, by = c("CG_44_group_label" = "group"))

print("各纬度区间的β2TPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== β2TPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== β2TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的纬度范围
cat("β2TPM原始数据纬度范围:", min(data_clean$Lat, na.rm = TRUE), "到", max(data_clean$Lat, na.rm = TRUE), "\n")
cat("β2TPM总数据点数:", nrow(data_clean), "\n")
cat("β2TPM纬度区间数:", length(unique(data_clean$Lat_group)), "\n")
