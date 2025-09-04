# 加载必要的包
library(ggplot2)
library(dplyr)

# 读取数据
data <- read.table("Alt.RCA.tsv", sep="\t", header=TRUE)

# 查看数据结构
str(data)
head(data)

# Long-TPM列
# 清理数据，去除缺失值并限制经度范围
data_clean <- data %>%
  filter(!is.na(Long) & !is.na(TPM)) %>%
  filter(Long >= -120 & Long < 180)

# 创建经度区间分组（以20度为区间，从-120到180）
data_clean$Long_group <- cut(data_clean$Long, 
                           breaks = seq(-120, 180, by = 20),
                           include.lowest = TRUE,
                           right = FALSE)

# 创建更简洁的标签格式
breaks_seq <- seq(-120, 180, by = 20)
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$Long_group)
data_clean$Long_group_label <- break_labels[as.character(data_clean$Long_group)]

# 计算每个经度区间的中点，用于x轴位置
data_clean$Long_midpoint <- as.numeric(gsub("\\[([0-9.-]+),.*", "\\1", data_clean$Long_group)) + 10

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(Long_group_label, levels = break_labels), y = TPM)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "Longitude",
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
ggsave("TPM_vs_Longitude_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)


# Long-αTPM列
# 清理数据，去除缺失值并限制经度范围
data_clean <- data %>%
  filter(!is.na(Long) & !is.na(α_TPM..)) %>%
  filter(Long >= -120 & Long < 180)

# 创建经度区间分组（以20度为区间，从-120到180）
data_clean$Long_group <- cut(data_clean$Long, 
                           breaks = seq(-120, 180, by = 20),
                           include.lowest = TRUE,
                           right = FALSE)

# 创建更简洁的标签格式
breaks_seq <- seq(-120, 180, by = 20)
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$Long_group)
data_clean$Long_group_label <- break_labels[as.character(data_clean$Long_group)]

# 计算每个经度区间的中点，用于x轴位置
data_clean$Long_midpoint <- as.numeric(gsub("\\[([0-9.-]+),.*", "\\1", data_clean$Long_group)) + 10

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(Long_group_label, levels = break_labels), y = α_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "Longitude",
    y = "aTpm"
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
  scale_y_continuous(limits = c(0, max(data_clean$α_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("A-TPM_vs_Longitude_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)


# Long-βTPM列
# 清理数据，去除缺失值并限制经度范围
data_clean <- data %>%
  filter(!is.na(Long) & !is.na(β_TPM..)) %>%
  filter(Long >= -120 & Long < 180)

# 创建经度区间分组（以20度为区间，从-120到180）
data_clean$Long_group <- cut(data_clean$Long,
                            breaks = seq(-120, 180, by = 20),
                            include.lowest = TRUE,
                            right = FALSE)

# 创建更简洁的标签格式
breaks_seq <- seq(-120, 180, by = 20)
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$Long_group)
data_clean$Long_group_label <- break_labels[as.character(data_clean$Long_group)]

# 计算每个经度区间的中点，用于x轴位置
data_clean$Long_midpoint <- as.numeric(gsub("\\[([0-9.-]+),.*", "\\1", data_clean$Long_group)) + 10

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(Long_group_label, levels = break_labels), y = β_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "Longitude",
    y = "bTpm"
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
  scale_y_continuous(limits = c(0, max(data_clean$β_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B-TPM_vs_Longitude_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)


# Long-B1TPM列
# 清理数据，去除缺失值并限制经度范围
data_clean <- data %>%
  filter(!is.na(Long) & !is.na(β1_TPM..)) %>%
  filter(Long >= -120 & Long < 180)

# 创建经度区间分组（以20度为区间，从-120到180）
data_clean$Long_group <- cut(data_clean$Long, 
                           breaks = seq(-120, 180, by = 20),
                           include.lowest = TRUE,
                           right = FALSE)

# 创建更简洁的标签格式
breaks_seq <- seq(-120, 180, by = 20)
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$Long_group)
data_clean$Long_group_label <- break_labels[as.character(data_clean$Long_group)]

# 计算每个经度区间的中点，用于x轴位置
data_clean$Long_midpoint <- as.numeric(gsub("\\[([0-9.-]+),.*", "\\1", data_clean$Long_group)) + 10

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(Long_group_label, levels = break_labels), y = β1_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "Longitude",
    y = "b1Tpm"
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
  scale_y_continuous(limits = c(0, max(data_clean$β1_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B1-TPM_vs_Longitude_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)


# Lat-B2TPM列
# 清理数据，去除缺失值并限制经度范围
data_clean <- data %>%
  filter(!is.na(Long) & !is.na(β2_TPM..)) %>%
  filter(Long >= -120 & Long < 180)

# 创建经度区间分组（以20度为区间，从-120到180）
data_clean$Long_group <- cut(data_clean$Long, 
                           breaks = seq(-120, 180, by = 20),
                           include.lowest = TRUE,
                           right = FALSE)

# 创建更简洁的标签格式
breaks_seq <- seq(-120, 180, by = 20)
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$Long_group)
data_clean$Long_group_label <- break_labels[as.character(data_clean$Long_group)]

# 计算每个经度区间的中点，用于x轴位置
data_clean$Long_midpoint <- as.numeric(gsub("\\[([0-9.-]+),.*", "\\1", data_clean$Long_group)) + 10

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(Long_group_label, levels = break_labels), y = β2_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "Longitude",
    y = "b2Tpm"
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
  scale_y_continuous(limits = c(0, max(data_clean$β2_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_Longitude_boxplot.pdf", p, width = 12, height = 6, dpi = 1200)

