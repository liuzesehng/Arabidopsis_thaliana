# 加载必要的包
library(ggplot2)
library(dplyr)

# 读取数据
data <- read.table("Alt.RCA.me_snp.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

# 查看数据结构
str(data)
head(data)


# CHH_terminator-βTPM列
# 清理数据，去除缺失值
data_clean_CHH_terminator <- data %>%
  filter(!is.na(CHH_terminator) & !is.na(β_TPM..))

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

# 创建箱线图
p_CHH_terminator <- ggplot(data_clean_CHH_terminator, aes(x = factor(CHH_terminator_group_label, levels = break_labels), y = β_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "CHH Terminator",
    y = "bTPM"
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
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
  scale_y_continuous(limits = c(0, max(data_clean_CHH_terminator$β_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CHH_terminator)

# 保存图形为PDF格式
ggsave("B-TPM_vs_CHH_terminator_boxplot.pdf", p_CHH_terminator, width = 12, height = 9, dpi = 1200)


# CHG-β_TPM列
# 清理数据，去除缺失值
data_clean_CHG <- data %>%
  filter(!is.na(CHG) & !is.na(β_TPM..))

# 区间
start_val <- min(data_clean_CHG$CHG, na.rm = TRUE)
end_val <- ceiling(max(data_clean_CHG$CHG, na.rm = TRUE) * 10) / 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 0.005)
} else {
  breaks_seq <- seq(start_val, end_val, by = 0.005)
}
data_clean_CHG$CHG_group <- cut(
  data_clean_CHG$CHG,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_CHG$CHG_group)
data_clean_CHG$CHG_group_label <- break_labels[as.character(data_clean_CHG$CHG_group)]

# 计算每个CHG区间的中点，用于x轴位置
data_clean_CHG$CHG_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_CHG$CHG_group)) + 0.0025

# 创建箱线图
p_CHG <- ggplot(data_clean_CHG, aes(x = factor(CHG_group_label, levels = break_labels), y = β_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "CHG",
    y = "bTPM"
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
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
  scale_y_continuous(limits = c(0, max(data_clean_CHG$β_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CHG)

# 保存图形为PDF格式
ggsave("B-TPM_vs_CHG_boxplot.pdf", p_CHG, width = 12, height = 9, dpi = 1200)


# CHG_162-β_TPM列
# 清理数据，去除缺失值
data_clean_CHG_162 <- data %>%
  filter(!is.na(CHG_162) & !is.na(β_TPM..))

# 区间
start_val <- min(data_clean_CHG_162$CHG_162, na.rm = TRUE)
end_val <- ceiling(max(data_clean_CHG_162$CHG_162, na.rm = TRUE) / 10) * 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 5)
} else {
  breaks_seq <- seq(start_val, end_val, by = 5)
}
data_clean_CHG_162$CHG_162_group <- cut(
  data_clean_CHG_162$CHG_162,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_CHG_162$CHG_162_group)
data_clean_CHG_162$CHG_162_group_label <- break_labels[as.character(data_clean_CHG_162$CHG_162_group)]

# 计算每个CHG_162区间的中点，用于x轴位置
data_clean_CHG_162$CHG_162_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_CHG_162$CHG_162_group)) + 2.5

# 创建箱线图
p_CHG_162 <- ggplot(data_clean_CHG_162, aes(x = factor(CHG_162_group_label, levels = break_labels), y = β_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "CHG_162",
    y = "bTPM"
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
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
  scale_y_continuous(limits = c(0, max(data_clean_CHG_162$β_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CHG_162)

# 保存图形为PDF格式
ggsave("B-TPM_vs_CHG_162_boxplot.pdf", p_CHG_162, width = 12, height = 9, dpi = 1200)


# CHH-β_TPM列
# 清理数据，去除缺失值
data_clean_CHH <- data %>%
  filter(!is.na(CHH) & !is.na(β_TPM..))

# 区间
start_val <- min(data_clean_CHH$CHH, na.rm = TRUE)
end_val <- ceiling(max(data_clean_CHH$CHH, na.rm = TRUE) * 100) / 100
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 0.005)
} else {
  breaks_seq <- seq(start_val, end_val, by = 0.005)
}
data_clean_CHH$CHH_group <- cut(
  data_clean_CHH$CHH,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_CHH$CHH_group)
data_clean_CHH$CHH_group_label <- break_labels[as.character(data_clean_CHH$CHH_group)]

# 计算每个CHH区间的中点，用于x轴位置
data_clean_CHH$CHH_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_CHH$CHH_group)) + 0.0025

# 创建箱线图
p_CHH <- ggplot(data_clean_CHH, aes(x = factor(CHH_group_label, levels = break_labels), y = β_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "CHH",
    y = "bTPM"
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
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
  scale_y_continuous(limits = c(0, max(data_clean_CHH$β_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CHH)

# 保存图形为PDF格式
ggsave("B-TPM_vs_CHH_boxplot.pdf", p_CHH, width = 12, height = 9, dpi = 1200)


# CHH_1259-β_TPM列
# 清理数据，去除缺失值
data_clean_CHH_1259 <- data %>%
  filter(!is.na(CHH_1259) & !is.na(β_TPM..))

# 区间
start_val <- min(data_clean_CHH_1259$CHH_1259, na.rm = TRUE)
end_val <- ceiling(max(data_clean_CHH_1259$CHH_1259, na.rm = TRUE) / 10) * 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 5)
} else {
  breaks_seq <- seq(start_val, end_val, by = 5)
}
data_clean_CHH_1259$CHH_1259_group <- cut(
  data_clean_CHH_1259$CHH_1259,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_CHH_1259$CHH_1259_group)
data_clean_CHH_1259$CHH_1259_group_label <- break_labels[as.character(data_clean_CHH_1259$CHH_1259_group)]

# 计算每个CHH_1259区间的中点，用于x轴位置
data_clean_CHH_1259$CHH_1259_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_CHH_1259$CHH_1259_group)) + 2.5

# 创建箱线图
p_CHH_1259 <- ggplot(data_clean_CHH_1259, aes(x = factor(CHH_1259_group_label, levels = break_labels), y = β_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "CHH_1259",
    y = "bTPM"
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
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
  scale_y_continuous(limits = c(0, max(data_clean_CHH_1259$β_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CHH_1259)

# 保存图形为PDF格式
ggsave("B-TPM_vs_CHH_1259_boxplot.pdf", p_CHH_1259, width = 12, height = 9, dpi = 1200)


# CHH_306-β_TPM列
# 清理数据，去除缺失值
data_clean_CHH_306 <- data %>%
  filter(!is.na(CHH_306) & !is.na(β_TPM..))

# 区间
start_val <- min(data_clean_CHH_306$CHH_306, na.rm = TRUE)
end_val <- ceiling(max(data_clean_CHH_306$CHH_306, na.rm = TRUE) / 10) * 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 5)
} else {
  breaks_seq <- seq(start_val, end_val, by = 5)
}
data_clean_CHH_306$CHH_306_group <- cut(
  data_clean_CHH_306$CHH_306,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_CHH_306$CHH_306_group)
data_clean_CHH_306$CHH_306_group_label <- break_labels[as.character(data_clean_CHH_306$CHH_306_group)]

# 计算每个CHH_306区间的中点，用于x轴位置
data_clean_CHH_306$CHH_306_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_CHH_306$CHH_306_group)) + 2.5

# 创建箱线图
p_CHH_306 <- ggplot(data_clean_CHH_306, aes(x = factor(CHH_306_group_label, levels = break_labels), y = β_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "CHH_306",
    y = "bTPM"
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
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
  scale_y_continuous(limits = c(0, max(data_clean_CHH_306$β_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CHH_306)

# 保存图形为PDF格式
ggsave("B-TPM_vs_CHH_306_boxplot.pdf", p_CHH_306, width = 12, height = 9, dpi = 1200)

