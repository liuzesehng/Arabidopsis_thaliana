# 加载必要的包
library(ggplot2)
library(dplyr)

# 读取数据
data <- read.table("Alt.RCA.me_snp.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

# 查看数据结构
str(data)
head(data)

# CHH_terminator-α_TPM列
# 清理数据，去除缺失值
data_clean_CHH_terminator <- data %>%
  filter(!is.na(CHH_terminator) & !is.na(α_TPM..))

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
p_CHH_terminator <- ggplot(data_clean_CHH_terminator, aes(x = factor(CHH_terminator_group_label, levels = break_labels), y = α_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "CHH Terminator",
    y = "aTPM"
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
  scale_y_continuous(limits = c(0, max(data_clean_CHH_terminator$α_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CHH_terminator)

# 保存图形为PDF格式
ggsave("A-TPM_vs_CHH_terminator_boxplot.pdf", p_CHH_terminator, width = 12, height = 9, dpi = 1200)


# CHH-α_TPM列
# 清理数据，去除缺失值
data_clean_CHH <- data %>%
  filter(!is.na(CHH) & !is.na(α_TPM..))

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
p_CHH <- ggplot(data_clean_CHH, aes(x = factor(CHH_group_label, levels = break_labels), y = α_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "CHH",
    y = "aTPM"
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
  scale_y_continuous(limits = c(0, max(data_clean_CHH$α_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CHH)

# 保存图形为PDF格式
ggsave("A-TPM_vs_CHH_boxplot.pdf", p_CHH, width = 12, height = 9, dpi = 1200)


# CG_terminator-α_TPM列
# 清理数据，去除缺失值
data_clean_CG_terminator <- data %>%
  filter(!is.na(CG_terminator) & !is.na(α_TPM..))

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

# 创建箱线图
p_CG_terminator <- ggplot(data_clean_CG_terminator, aes(x = factor(CG_terminator_group_label, levels = break_labels), y = α_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "CG Terminator",
    y = "aTPM"
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
  scale_y_continuous(limits = c(0, max(data_clean_CG_terminator$α_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CG_terminator)

# 保存图形为PDF格式
ggsave("A-TPM_vs_CG_terminator_boxplot.pdf", p_CG_terminator, width = 12, height = 9, dpi = 1200)


# CG_95-α_TPM列
# 清理数据，去除缺失值
data_clean_CG_95 <- data %>%
  filter(!is.na(CG_95) & !is.na(α_TPM..))

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

# 创建箱线图
p_CG_95 <- ggplot(data_clean_CG_95, aes(x = factor(CG_95_group_label, levels = break_labels), y = α_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "CG_95",
    y = "aTPM"
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
  scale_y_continuous(limits = c(0, max(data_clean_CG_95$α_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CG_95)

# 保存图形为PDF格式
ggsave("A-TPM_vs_CG_95_boxplot.pdf", p_CG_95, width = 12, height = 9, dpi = 1200)


# CHH_promoter-α_TPM列
# 清理数据，去除缺失值
data_clean_CHH_promoter <- data %>%
  filter(!is.na(CHH_promoter) & !is.na(α_TPM..))

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

# 创建箱线图
p_CHH_promoter <- ggplot(data_clean_CHH_promoter, aes(x = factor(CHH_promoter_group_label, levels = break_labels), y = α_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "CHH Promoter",
    y = "aTPM"
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
  scale_y_continuous(limits = c(0, max(data_clean_CHH_promoter$α_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CHH_promoter)

# 保存图形为PDF格式
ggsave("A-TPM_vs_CHH_promoter_boxplot.pdf", p_CHH_promoter, width = 12, height = 9, dpi = 1200)


# CHH_327-α_TPM列
# 清理数据，去除缺失值
data_clean_CHH_327 <- data %>%
  filter(!is.na(CHH_327) & !is.na(α_TPM..))

# 区间
start_val <- min(data_clean_CHH_327$CHH_327, na.rm = TRUE)
end_val <- ceiling(max(data_clean_CHH_327$CHH_327, na.rm = TRUE) / 10) * 10
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 5)
} else {
  breaks_seq <- seq(start_val, end_val, by = 5)
}
data_clean_CHH_327$CHH_327_group <- cut(
  data_clean_CHH_327$CHH_327,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean_CHH_327$CHH_327_group)
data_clean_CHH_327$CHH_327_group_label <- break_labels[as.character(data_clean_CHH_327$CHH_327_group)]

# 计算每个CHH_327区间的中点，用于x轴位置
data_clean_CHH_327$CHH_327_midpoint <- as.numeric(gsub("\\[([0-9.]+),.*", "\\1", data_clean_CHH_327$CHH_327_group)) + 2.5

# 创建箱线图
p_CHH_327 <- ggplot(data_clean_CHH_327, aes(x = factor(CHH_327_group_label, levels = break_labels), y = α_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "CHH_327",
    y = "aTPM"
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
  scale_y_continuous(limits = c(0, max(data_clean_CHH_327$α_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p_CHH_327)

# 保存图形为PDF格式
ggsave("A-TPM_vs_CHH_327_boxplot.pdf", p_CHH_327, width = 12, height = 9, dpi = 1200)

