# 加载必要的包
library(ggplot2)
library(dplyr)


# 读取数据
data <- read.table("Alt.RCA.tsv", sep="\t", header=TRUE)

# 查看数据结构
str(data)
head(data)

# Alt-TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(alt.m) & !is.na(TPM))

# 区间
start_val <- min(data_clean$alt.m, na.rm = TRUE)
end_val <- ceiling(max(data_clean$alt.m, na.rm = TRUE) / 100) * 100
if (end_val < start_val) {
  # 处理特殊情况，强制区间为10
  breaks_seq <- c(start_val, start_val + 400)
} else {
  breaks_seq <- seq(start_val, end_val, by = 400)
}
data_clean$Alt_group <- cut(
  data_clean$alt.m,
  breaks = breaks_seq,
  include.lowest = TRUE,
  right = FALSE
)

# 创建更简洁的标签格式（30-35）
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$Alt_group)
data_clean$Alt_group_label <- break_labels[as.character(data_clean$Alt_group)]

# 创建海拔区间分组（以400米为区间）
data_clean$Alt_group <- cut(data_clean$alt.m, 
                           breaks = seq(-20, ceiling(max(data_clean$alt.m, na.rm = TRUE)/400)*400, by = 400),
                           include.lowest = TRUE,
                           right = FALSE)

# 创建更简洁的标签格式
breaks_seq <- seq(-20, ceiling(max(data_clean$alt.m, na.rm = TRUE)/400)*400, by = 400)
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$Alt_group)
data_clean$Alt_group_label <- break_labels[as.character(data_clean$Alt_group)]

# 计算每个海拔区间的中点，用于x轴位置
data_clean$Alt_midpoint <- as.numeric(gsub("\\[([0-9.-]+),.*", "\\1", data_clean$Alt_group)) + 200

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(Alt_group_label, levels = break_labels), y = TPM)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "Altitude",
    y = "Tpm"
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
  scale_y_continuous(limits = c(0, max(data_clean$TPM, na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("TPM_vs_Altitude_boxplot.pdf", p, width = 12, height = 9, dpi = 1200)


# Alt-αTPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(alt.m) & !is.na(α_TPM..))

# 创建海拔区间分组（以400米为区间）
data_clean$Alt_group <- cut(data_clean$alt.m, 
                           breaks = seq(-20, ceiling(max(data_clean$alt.m, na.rm = TRUE)/400)*400, by = 400),
                           include.lowest = TRUE,
                           right = FALSE)

# 创建更简洁的标签格式
breaks_seq <- seq(-20, ceiling(max(data_clean$alt.m, na.rm = TRUE)/400)*400, by = 400)
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$Alt_group)
data_clean$Alt_group_label <- break_labels[as.character(data_clean$Alt_group)]

# 计算每个海拔区间的中点，用于x轴位置
data_clean$Alt_midpoint <- as.numeric(gsub("\\[([0-9.-]+),.*", "\\1", data_clean$Alt_group)) + 200

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(Alt_group_label, levels = break_labels), y = α_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "Altitude",
    y = "aTpm"
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
  scale_y_continuous(limits = c(0, max(data_clean$α_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("A-TPM_vs_Altitude_boxplot.pdf", p, width = 12, height = 9, dpi = 1200)


# Alt-βTPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(alt.m) & !is.na(β_TPM..))

# 创建海拔区间分组（以400米为区间）
data_clean$Alt_group <- cut(data_clean$alt.m, 
                            breaks = seq(-20, ceiling(max(data_clean$alt.m, na.rm = TRUE)/400)*400, by = 400),
                            include.lowest = TRUE,
                            right = FALSE)

# 创建更简洁的标签格式
breaks_seq <- seq(-20, ceiling(max(data_clean$alt.m, na.rm = TRUE)/400)*400, by = 400)
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$Alt_group)
data_clean$Alt_group_label <- break_labels[as.character(data_clean$Alt_group)]

# 计算每个海拔区间的中点，用于x轴位置
data_clean$Alt_midpoint <- as.numeric(gsub("\\[([0-9.-]+),.*", "\\1", data_clean$Alt_group)) + 200

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(Alt_group_label, levels = break_labels), y = β_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "Altitude",
    y = "bTpm"
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
  scale_y_continuous(limits = c(0, max(data_clean$β_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B-TPM_vs_Altitude_boxplot.pdf", p, width = 12, height = 9, dpi = 1200)


# Alt-β1TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(alt.m) & !is.na(β1_TPM..))

# 创建海拔区间分组（以400米为区间）
data_clean$Alt_group <- cut(data_clean$alt.m, 
                           breaks = seq(-20, ceiling(max(data_clean$alt.m, na.rm = TRUE)/400)*400, by = 400),
                           include.lowest = TRUE,
                           right = FALSE)

# 创建更简洁的标签格式
breaks_seq <- seq(-20, ceiling(max(data_clean$alt.m, na.rm = TRUE)/400)*400, by = 400)
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$Alt_group)
data_clean$Alt_group_label <- break_labels[as.character(data_clean$Alt_group)]

# 计算每个海拔区间的中点，用于x轴位置
data_clean$Alt_midpoint <- as.numeric(gsub("\\[([0-9.-]+),.*", "\\1", data_clean$Alt_group)) + 200

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(Alt_group_label, levels = break_labels), y = β1_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "Altitude",
    y = "b1Tpm"
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
  scale_y_continuous(limits = c(0, max(data_clean$β1_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B1-TPM_vs_Altitude_boxplot.pdf", p, width = 12, height = 9, dpi = 1200)


# Alt-β2TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(alt.m) & !is.na(β2_TPM..))

# 创建海拔区间分组（以400米为区间）
data_clean$Alt_group <- cut(data_clean$alt.m, 
                           breaks = seq(-20, ceiling(max(data_clean$alt.m, na.rm = TRUE)/400)*400, by = 400),
                           include.lowest = TRUE,
                           right = FALSE)

# 创建更简洁的标签格式
breaks_seq <- seq(-20, ceiling(max(data_clean$alt.m, na.rm = TRUE)/400)*400, by = 400)
break_labels <- paste0(breaks_seq[-length(breaks_seq)], "-", breaks_seq[-1])
names(break_labels) <- levels(data_clean$Alt_group)
data_clean$Alt_group_label <- break_labels[as.character(data_clean$Alt_group)]

# 计算每个海拔区间的中点，用于x轴位置
data_clean$Alt_midpoint <- as.numeric(gsub("\\[([0-9.-]+),.*", "\\1", data_clean$Alt_group)) + 200

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(Alt_group_label, levels = break_labels), y = β2_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  labs(
    x = "Altitude",
    y = "b2Tpm"
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
  scale_y_continuous(limits = c(0, max(data_clean$β2_TPM.., na.rm = TRUE) * 1.3))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_Altitude_boxplot.pdf", p, width = 12, height = 9, dpi = 1200)

