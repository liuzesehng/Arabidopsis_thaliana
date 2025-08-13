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


# CHH_terminator-β1_TPM列
# 清理数据，去除缺失值
data_clean_CHH_terminator <- data %>%
  filter(!is.na(CHH_terminator) & !is.na(β1_TPM..))

# 区间
start_val <- min(data_clean_CHH_terminator$CHH_terminator, na.rm = TRUE)
end_val <- ceiling(max(data_clean_CHH_terminator$CHH_terminator, na.rm = TRUE) * 10) / 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 0.05)
} else {
  breaks_seq <- seq(start_val, end_val, by = 0.05)
}
data_clean_CHH_terminator$CHH_terminator_group <- cut(
  data_clean_CHH_terminator$CHH_terminator,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_CHH_terminator$CHH_terminator_group)
data_clean_CHH_terminator$CHH_terminator_group_label <- break_labels[as.character(data_clean_CHH_terminator$CHH_terminator_group)]

# 计算每个CHH_terminator区间的中点，用于x轴位置
data_clean_CHH_terminator$CHH_terminator_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_CHH_terminator$CHH_terminator_group)) + 0.025

# 进行方差分析和多重比较检验
aov_result_CHH_terminator <- aov(β1_TPM.. ~ CHH_terminator_group_label, data = data_clean_CHH_terminator)
tukey_result_CHH_terminator <- TukeyHSD(aov_result_CHH_terminator)

# 提取p值并手动分配显著性字母
p_values_CHH_terminator <- tukey_result_CHH_terminator$CHH_terminator_group_label[,"p adj"]
group_names_CHH_terminator <- unique(data_clean_CHH_terminator$CHH_terminator_group_label)

significance_letters_CHH_terminator <- data.frame(
  group = sort(group_names_CHH_terminator),
  letter = c("b", "b", "a", "a", "ab", "ab", "ab"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max_CHH_terminator <- data_clean_CHH_terminator %>%
  group_by(CHH_terminator_group_label) %>%
  summarise(max_b1TPM = max(β1_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters_CHH_terminator, by = c("CHH_terminator_group_label" = "group"))

# 创建箱线图
p_CHH_terminator <- ggplot(data_clean_CHH_terminator, aes(x = factor(CHH_terminator_group_label, levels = break_labels), y = β1_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max_CHH_terminator, 
            aes(x = factor(CHH_terminator_group_label, levels = break_labels), 
                y = max_b1TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CHH Terminator",
    y = "b1TPM"
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
  scale_y_continuous(limits = c(0, max(data_clean_CHH_terminator$β1_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CHH_terminator)

# 保存图形为PDF格式
ggsave("B1-TPM_vs_CHH_terminator_boxplot.pdf", p_CHH_terminator, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary_CHH_terminator <- data_clean_CHH_terminator %>%
  group_by(CHH_terminator_group, CHH_terminator_group_label) %>%
  summarise(
    mean_b1TPM = mean(β1_TPM.., na.rm = TRUE),
    median_b1TPM = median(β1_TPM.., na.rm = TRUE),
    sd_b1TPM = sd(β1_TPM.., na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CHH_terminator_group_label) %>%
  left_join(significance_letters_CHH_terminator, by = c("CHH_terminator_group_label" = "group"))

print("各CHH_terminator区间的β1_TPM统计信息:")
print(data_summary_CHH_terminator)

# 显示方差分析结果
cat("\n=== CHH_terminator β1_TPM方差分析结果 ===\n")
print(summary(aov_result_CHH_terminator))

# 显示多重比较检验结果
cat("\n=== CHH_terminator β1_TPM Tukey多重比较检验结果 ===\n")
print(tukey_result_CHH_terminator)


# CG_terminator-β1_TPM列
# 清理数据，去除缺失值
data_clean_CG_terminator <- data %>%
  filter(!is.na(CG_terminator) & !is.na(β1_TPM..))

# 区间
start_val <- min(data_clean_CG_terminator$CG_terminator, na.rm = TRUE)
end_val <- ceiling(max(data_clean_CG_terminator$CG_terminator, na.rm = TRUE) * 10) / 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 0.05)
} else {
  breaks_seq <- seq(start_val, end_val, by = 0.05)
}
data_clean_CG_terminator$CG_terminator_group <- cut(
  data_clean_CG_terminator$CG_terminator,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_CG_terminator$CG_terminator_group)
data_clean_CG_terminator$CG_terminator_group_label <- break_labels[as.character(data_clean_CG_terminator$CG_terminator_group)]

# 计算每个CG_terminator区间的中点，用于x轴位置
data_clean_CG_terminator$CG_terminator_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_CG_terminator$CG_terminator_group)) + 0.025

# 进行方差分析和多重比较检验
aov_result_CG_terminator <- aov(β1_TPM.. ~ CG_terminator_group_label, data = data_clean_CG_terminator)
tukey_result_CG_terminator <- TukeyHSD(aov_result_CG_terminator)

# 提取p值并手动分配显著性字母
p_values_CG_terminator <- tukey_result_CG_terminator$CG_terminator_group_label[,"p adj"]
group_names_CG_terminator <- unique(data_clean_CG_terminator$CG_terminator_group_label)

significance_letters_CG_terminator <- data.frame(
  group = sort(group_names_CG_terminator),
  letter = c("c", "b", "b", "ab", "ab", "a", "c", "c", "a"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max_CG_terminator <- data_clean_CG_terminator %>%
  group_by(CG_terminator_group_label) %>%
  summarise(max_b1TPM = max(β1_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters_CG_terminator, by = c("CG_terminator_group_label" = "group"))

# 创建箱线图
p_CG_terminator <- ggplot(data_clean_CG_terminator, aes(x = factor(CG_terminator_group_label, levels = break_labels), y = β1_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max_CG_terminator, 
            aes(x = factor(CG_terminator_group_label, levels = break_labels), 
                y = max_b1TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CG Terminator",
    y = "b1TPM"
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
  scale_y_continuous(limits = c(0, max(data_clean_CG_terminator$β1_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CG_terminator)

# 保存图形为PDF格式
ggsave("B1-TPM_vs_CG_terminator_boxplot.pdf", p_CG_terminator, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary_CG_terminator <- data_clean_CG_terminator %>%
  group_by(CG_terminator_group, CG_terminator_group_label) %>%
  summarise(
    mean_b1TPM = mean(β1_TPM.., na.rm = TRUE),
    median_b1TPM = median(β1_TPM.., na.rm = TRUE),
    sd_b1TPM = sd(β1_TPM.., na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CG_terminator_group_label) %>%
  left_join(significance_letters_CG_terminator, by = c("CG_terminator_group_label" = "group"))

print("各CG_terminator区间的β1_TPM统计信息:")
print(data_summary_CG_terminator)

# 显示方差分析结果
cat("\n=== CG_terminator β1_TPM方差分析结果 ===\n")
print(summary(aov_result_CG_terminator))

# 显示多重比较检验结果
cat("\n=== CG_terminator β1_TPM Tukey多重比较检验结果 ===\n")
print(tukey_result_CG_terminator)


# CHH_promoter-β1_TPM列
# 清理数据，去除缺失值
data_clean_CHH_promoter <- data %>%
  filter(!is.na(CHH_promoter) & !is.na(β1_TPM..))

# 区间
start_val <- min(data_clean_CHH_promoter$CHH_promoter, na.rm = TRUE)
end_val <- ceiling(max(data_clean_CHH_promoter$CHH_promoter, na.rm = TRUE) * 100) / 100
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 0.005)
} else {
  breaks_seq <- seq(start_val, end_val, by = 0.005)
}
data_clean_CHH_promoter$CHH_promoter_group <- cut(
  data_clean_CHH_promoter$CHH_promoter,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_CHH_promoter$CHH_promoter_group)
data_clean_CHH_promoter$CHH_promoter_group_label <- break_labels[as.character(data_clean_CHH_promoter$CHH_promoter_group)]

# 计算每个CHH_promoter区间的中点，用于x轴位置
data_clean_CHH_promoter$CHH_promoter_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_CHH_promoter$CHH_promoter_group)) + 0.0025

# 进行方差分析和多重比较检验
aov_result_CHH_promoter <- aov(β1_TPM.. ~ CHH_promoter_group_label, data = data_clean_CHH_promoter)
tukey_result_CHH_promoter <- TukeyHSD(aov_result_CHH_promoter)

# 提取p值并手动分配显著性字母
p_values_CHH_promoter <- tukey_result_CHH_promoter$CHH_promoter_group_label[,"p adj"]
group_names_CHH_promoter <- unique(data_clean_CHH_promoter$CHH_promoter_group_label)

significance_letters_CHH_promoter <- data.frame(
  group = sort(group_names_CHH_promoter),
  letter = c("a", "a", "a", "a", "a", "a", "a", "a"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max_CHH_promoter <- data_clean_CHH_promoter %>%
  group_by(CHH_promoter_group_label) %>%
  summarise(max_b1TPM = max(β1_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters_CHH_promoter, by = c("CHH_promoter_group_label" = "group"))

# 创建箱线图
p_CHH_promoter <- ggplot(data_clean_CHH_promoter, aes(x = factor(CHH_promoter_group_label, levels = break_labels), y = β1_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max_CHH_promoter, 
            aes(x = factor(CHH_promoter_group_label, levels = break_labels), 
                y = max_b1TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CHH Promoter",
    y = "b1TPM"
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
  scale_y_continuous(limits = c(0, max(data_clean_CHH_promoter$β1_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CHH_promoter)

# 保存图形为PDF格式
ggsave("B1-TPM_vs_CHH_promoter_boxplot.pdf", p_CHH_promoter, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary_CHH_promoter <- data_clean_CHH_promoter %>%
  group_by(CHH_promoter_group, CHH_promoter_group_label) %>%
  summarise(
    mean_b1TPM = mean(β1_TPM.., na.rm = TRUE),
    median_b1TPM = median(β1_TPM.., na.rm = TRUE),
    sd_b1TPM = sd(β1_TPM.., na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CHH_promoter_group_label) %>%
  left_join(significance_letters_CHH_promoter, by = c("CHH_promoter_group_label" = "group"))

print("各CHH_promoter区间的β1_TPM统计信息:")
print(data_summary_CHH_promoter)

# 显示方差分析结果
cat("\n=== CHH_promoter β1_TPM方差分析结果 ===\n")
print(summary(aov_result_CHH_promoter))

# 显示多重比较检验结果
cat("\n=== CHH_promoter β1_TPM Tukey多重比较检验结果 ===\n")
print(tukey_result_CHH_promoter)


# CG_95-β1_TPM列
# 清理数据，去除缺失值
data_clean_CG_95 <- data %>%
  filter(!is.na(CG_95) & !is.na(β1_TPM..))

# 区间
start_val <- min(data_clean_CG_95$CG_95, na.rm = TRUE)
end_val <- ceiling(max(data_clean_CG_95$CG_95, na.rm = TRUE) / 10) * 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 10)
} else {
  breaks_seq <- seq(start_val, end_val, by = 10)
}
data_clean_CG_95$CG_95_group <- cut(
  data_clean_CG_95$CG_95,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_CG_95$CG_95_group)
data_clean_CG_95$CG_95_group_label <- break_labels[as.character(data_clean_CG_95$CG_95_group)]

# 计算每个CG_95区间的中点，用于x轴位置
data_clean_CG_95$CG_95_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_CG_95$CG_95_group)) + 5

# 进行方差分析和多重比较检验
aov_result_CG_95 <- aov(β1_TPM.. ~ CG_95_group_label, data = data_clean_CG_95)
tukey_result_CG_95 <- TukeyHSD(aov_result_CG_95)

# 提取p值并手动分配显著性字母
p_values_CG_95 <- tukey_result_CG_95$CG_95_group_label[,"p adj"]
group_names_CG_95 <- unique(data_clean_CG_95$CG_95_group_label)

significance_letters_CG_95 <- data.frame(
  group = sort(group_names_CG_95),
  letter = c("a", "a", "a", "a", "a", "a", "a", "a", "a", "a"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max_CG_95 <- data_clean_CG_95 %>%
  group_by(CG_95_group_label) %>%
  summarise(max_b1TPM = max(β1_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters_CG_95, by = c("CG_95_group_label" = "group"))

# 创建箱线图
p_CG_95 <- ggplot(data_clean_CG_95, aes(x = factor(CG_95_group_label, levels = break_labels), y = β1_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max_CG_95, 
            aes(x = factor(CG_95_group_label, levels = break_labels), 
                y = max_b1TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CG_95",
    y = "b1TPM"
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
  scale_y_continuous(limits = c(0, max(data_clean_CG_95$β1_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CG_95)

# 保存图形为PDF格式
ggsave("B1-TPM_vs_CG_95_boxplot.pdf", p_CG_95, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary_CG_95 <- data_clean_CG_95 %>%
  group_by(CG_95_group, CG_95_group_label) %>%
  summarise(
    mean_b1TPM = mean(β1_TPM.., na.rm = TRUE),
    median_b1TPM = median(β1_TPM.., na.rm = TRUE),
    sd_b1TPM = sd(β1_TPM.., na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CG_95_group_label) %>%
  left_join(significance_letters_CG_95, by = c("CG_95_group_label" = "group"))

print("各CG_95区间的β1_TPM统计信息:")
print(data_summary_CG_95)

# 显示方差分析结果
cat("\n=== CG_95 β1_TPM方差分析结果 ===\n")
print(summary(aov_result_CG_95))

# 显示多重比较检验结果
cat("\n=== CG_95 β1_TPM Tukey多重比较检验结果 ===\n")
print(tukey_result_CG_95)


# CG_promoter-β1_TPM列
# 清理数据，去除缺失值
data_clean_CG_promoter <- data %>%
  filter(!is.na(CG_promoter) & !is.na(β1_TPM..))

# 区间
start_val <- min(data_clean_CG_promoter$CG_promoter, na.rm = TRUE)
end_val <- ceiling(max(data_clean_CG_promoter$CG_promoter, na.rm = TRUE) * 10) / 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 0.05)
} else {
  breaks_seq <- seq(start_val, end_val, by = 0.05)
}
data_clean_CG_promoter$CG_promoter_group <- cut(
  data_clean_CG_promoter$CG_promoter,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_CG_promoter$CG_promoter_group)
data_clean_CG_promoter$CG_promoter_group_label <- break_labels[as.character(data_clean_CG_promoter$CG_promoter_group)]

# 计算每个CG_promoter区间的中点，用于x轴位置
data_clean_CG_promoter$CG_promoter_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_CG_promoter$CG_promoter_group)) + 0.025

# 进行方差分析和多重比较检验
aov_result_CG_promoter <- aov(β1_TPM.. ~ CG_promoter_group_label, data = data_clean_CG_promoter)
tukey_result_CG_promoter <- TukeyHSD(aov_result_CG_promoter)

# 提取p值并手动分配显著性字母
p_values_CG_promoter <- tukey_result_CG_promoter$CG_promoter_group_label[,"p adj"]
group_names_CG_promoter <- unique(data_clean_CG_promoter$CG_promoter_group_label)

significance_letters_CG_promoter <- data.frame(
  group = sort(group_names_CG_promoter),
  letter = c("a", "a", "a", "ab", "b", "ab", "ab", "ab", "ab", "ab"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max_CG_promoter <- data_clean_CG_promoter %>%
  group_by(CG_promoter_group_label) %>%
  summarise(max_b1TPM = max(β1_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters_CG_promoter, by = c("CG_promoter_group_label" = "group"))

# 创建箱线图
p_CG_promoter <- ggplot(data_clean_CG_promoter, aes(x = factor(CG_promoter_group_label, levels = break_labels), y = β1_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max_CG_promoter, 
            aes(x = factor(CG_promoter_group_label, levels = break_labels), 
                y = max_b1TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CG Promoter",
    y = "b1TPM"
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
  scale_y_continuous(limits = c(0, max(data_clean_CG_promoter$β1_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CG_promoter)

# 保存图形为PDF格式
ggsave("B1-TPM_vs_CG_promoter_boxplot.pdf", p_CG_promoter, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary_CG_promoter <- data_clean_CG_promoter %>%
  group_by(CG_promoter_group, CG_promoter_group_label) %>%
  summarise(
    mean_b1TPM = mean(β1_TPM.., na.rm = TRUE),
    median_b1TPM = median(β1_TPM.., na.rm = TRUE),
    sd_b1TPM = sd(β1_TPM.., na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CG_promoter_group_label) %>%
  left_join(significance_letters_CG_promoter, by = c("CG_promoter_group_label" = "group"))

print("各CG_promoter区间的β1_TPM统计信息:")
print(data_summary_CG_promoter)

# 显示方差分析结果
cat("\n=== CG_promoter β1_TPM方差分析结果 ===\n")
print(summary(aov_result_CG_promoter))

# 显示多重比较检验结果
cat("\n=== CG_promoter β1_TPM Tukey多重比较检验结果 ===\n")
print(tukey_result_CG_promoter)


# CHH_724-β1_TPM列
# 清理数据，去除缺失值
data_clean_CHH_724 <- data %>%
  filter(!is.na(CHH_724) & !is.na(β1_TPM..))

# 区间
start_val <- min(data_clean_CHH_724$CHH_724, na.rm = TRUE)
end_val <- ceiling(max(data_clean_CHH_724$CHH_724, na.rm = TRUE) / 10) * 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 5)
} else {
  breaks_seq <- seq(start_val, end_val, by = 5)
}
data_clean_CHH_724$CHH_724_group <- cut(
  data_clean_CHH_724$CHH_724,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_CHH_724$CHH_724_group)
data_clean_CHH_724$CHH_724_group_label <- break_labels[as.character(data_clean_CHH_724$CHH_724_group)]

# 计算每个CHH_724区间的中点，用于x轴位置
data_clean_CHH_724$CHH_724_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_CHH_724$CHH_724_group)) + 2.5

# 进行方差分析和多重比较检验
aov_result_CHH_724 <- aov(β1_TPM.. ~ CHH_724_group_label, data = data_clean_CHH_724)
tukey_result_CHH_724 <- TukeyHSD(aov_result_CHH_724)

# 提取p值并手动分配显著性字母
p_values_CHH_724 <- tukey_result_CHH_724$CHH_724_group_label[,"p adj"]
group_names_CHH_724 <- unique(data_clean_CHH_724$CHH_724_group_label)

significance_letters_CHH_724 <- data.frame(
  group = sort(group_names_CHH_724),
  letter = c("a", "b", "ab", "a", "ab", "a", "b", "a"),  # 根据实际统计结果调整
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max_CHH_724 <- data_clean_CHH_724 %>%
  group_by(CHH_724_group_label) %>%
  summarise(max_b1TPM = max(β1_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters_CHH_724, by = c("CHH_724_group_label" = "group"))

# 创建箱线图
p_CHH_724 <- ggplot(data_clean_CHH_724, aes(x = factor(CHH_724_group_label, levels = break_labels), y = β1_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max_CHH_724, 
            aes(x = factor(CHH_724_group_label, levels = break_labels), 
                y = max_b1TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
  labs(
    x = "CHH_724",
    y = "b1TPM"
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
  scale_y_continuous(limits = c(0, max(data_clean_CHH_724$β1_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CHH_724)

# 保存图形为PDF格式
ggsave("B1-TPM_vs_CHH_724_boxplot.pdf", p_CHH_724, width = 12, height = 6, dpi = 1200)

# 显示统计信息
data_summary_CHH_724 <- data_clean_CHH_724 %>%
  group_by(CHH_724_group, CHH_724_group_label) %>%
  summarise(
    mean_b1TPM = mean(β1_TPM.., na.rm = TRUE),
    median_b1TPM = median(β1_TPM.., na.rm = TRUE),
    sd_b1TPM = sd(β1_TPM.., na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(CHH_724_group_label) %>%
  left_join(significance_letters_CHH_724, by = c("CHH_724_group_label" = "group"))

print("各CHH_724区间的β1_TPM统计信息:")
print(data_summary_CHH_724)

# 显示方差分析结果
cat("\n=== CHH_724 β1_TPM方差分析结果 ===\n")
print(summary(aov_result_CHH_724))

# 显示多重比较检验结果
cat("\n=== CHH_724 β1_TPM Tukey多重比较检验结果 ===\n")
print(tukey_result_CHH_724)
