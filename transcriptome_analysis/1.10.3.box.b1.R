# 加载必要的包
library(ggplot2)
library(dplyr)
library(readr)
library(ggsignif)
library(car)
library(dunn.test)

# 读取数据
data <- read.table("RCA.climate.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 查看数据结构
str(data)
head(data)

# 检查温度列的唯一值
unique(data$`tem`)

# 将温度列转换为因子，确保正确排序（保留所有温度）
data$Temperature <- factor(data$`tem`, levels = c("10", "16", "22"))

# 创建用于统计分析的数据集（只包含10和16）
data_stats <- data %>% filter(`tem` %in% c("10", "16"))
data_stats$Temperature <- factor(data_stats$`tem`, levels = c("10", "16"))

# 检查数据
table(data$Temperature)
table(data_stats$Temperature)
summary(data$β1..)

# 统计摘要（只对10和16进行）
cat("\n=== 各温度组的β1/%TPM统计摘要 ===\n")
summary_stats <- data_stats %>%
  group_by(Temperature) %>%
  summarise(
    count = n(),
    mean = mean(β1.., na.rm = TRUE),
    median = median(β1.., na.rm = TRUE),
    sd = sd(β1.., na.rm = TRUE),
    min = min(β1.., na.rm = TRUE),
    max = max(β1.., na.rm = TRUE),
    q25 = quantile(β1.., 0.25, na.rm = TRUE),
    q75 = quantile(β1.., 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_stats)

# 进行正态性检验和方差齐性检验（只对10和16）
cat("\n=== 正态性检验 (Shapiro-Wilk test) ===\n")
shapiro_results <- data_stats %>%
  group_by(Temperature) %>%
  summarise(
    shapiro_p = shapiro.test(β1..)$p.value,
    .groups = 'drop'
  )
print(shapiro_results)

# 对残差进行正态性检验
temp_aov <- aov(β1.. ~ Temperature, data = data_stats)
cat("\n=== 残差正态性检验 ===\n")
residuals_shapiro <- shapiro.test(residuals(temp_aov))
cat("残差Shapiro-Wilk检验 p值:", residuals_shapiro$p.value, "\n")

# 方差齐性检验 (Levene's test)
cat("\n=== 方差齐性检验 (Levene's test) ===\n")
levene_result <- leveneTest(β1.. ~ Temperature, data = data_stats)
print(levene_result)

# Bartlett检验 (对正态分布数据)
cat("\n=== 方差齐性检验 (Bartlett's test) ===\n")
bartlett_result <- bartlett.test(β1.. ~ Temperature, data = data_stats)
print(bartlett_result)

# 使用 Welch's t-test (适用于方差不齐的情况)
cat("\n=== Welch's t-test (韦尔奇 t 检验) ===\n")
welch_result <- t.test(β1.. ~ Temperature, data = data_stats, var.equal = FALSE)
print(welch_result)
cat("t 统计量:", welch_result$statistic, "\n")
cat("自由度:", welch_result$parameter, "\n")
cat("p 值:", welch_result$p.value, "\n")
cat("95% 置信区间:", welch_result$conf.int, "\n")

# 根据 Welch's t-test 结果确定显著性标记
welch_annotation <- ifelse(welch_result$p.value < 0.001, "***",
                           ifelse(welch_result$p.value < 0.01, "**",
                                 ifelse(welch_result$p.value < 0.05, "*", "ns")))

# 创建箱线图（包含所有温度10、16、22）
p <- ggplot(data, aes(x = Temperature, y = β1.., fill = Temperature)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(width = 0.5)) +
  # 标准的 geom_boxplot，它会画出箱体和简洁的箱须
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, position = position_dodge(width = 0.5), 
               width = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  scale_fill_manual(values = c("10" = "blue", "16" = "yellow", "22" = "red")) +
  # 添加显著性标记（只标注10和16之间）
  geom_signif(comparisons = list(c("10", "16")),
              annotations = welch_annotation,
              y_position = c(max(data$β1.., na.rm = TRUE) * 1.1),
              tip_length = 0.02,
              textsize = 6) +
  labs(
    y = "b1/%",
    fill = "Tem/°C"
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    # 图例设置
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    # 去除x轴标题
    axis.title.x = element_blank(),
    # 显示完整的边框
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    # 设置刻度线朝外
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    # 确保刻度线朝外
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5)
  ) +
  # 设置y轴范围
  scale_y_continuous(limits = c(0, max(data$β1.., na.rm = TRUE) * 1.4))


# 显示图形
print(p)

# 保存图形
ggsave("temperature_b1_per_tpm_boxplot.pdf", plot = p, width = 12, height = 9, dpi = 1200)


# 统计摘要（只对10和16进行）
cat("\n=== 各温度组的b1TPM统计摘要 ===\n")
summary_stats <- data_stats %>%
  group_by(Temperature) %>%
  summarise(
    count = n(),
    mean = mean(β1, na.rm = TRUE),
    median = median(β1, na.rm = TRUE),
    sd = sd(β1, na.rm = TRUE),
    min = min(β1, na.rm = TRUE),
    max = max(β1, na.rm = TRUE),
    q25 = quantile(β1, 0.25, na.rm = TRUE),
    q75 = quantile(β1, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_stats)

# 进行正态性检验和方差齐性检验（只对10和16）
cat("\n=== 正态性检验 (Shapiro-Wilk test) ===\n")
shapiro_results <- data_stats %>%
  group_by(Temperature) %>%
  summarise(
    shapiro_p = shapiro.test(β1)$p.value,
    .groups = 'drop'
  )
print(shapiro_results)

# 对残差进行正态性检验
temp_aov <- aov(β1 ~ Temperature, data = data_stats)
cat("\n=== 残差正态性检验 ===\n")
residuals_shapiro <- shapiro.test(residuals(temp_aov))
cat("残差Shapiro-Wilk检验 p值:", residuals_shapiro$p.value, "\n")

# 方差齐性检验 (Levene's test)
cat("\n=== 方差齐性检验 (Levene's test) ===\n")
levene_result <- leveneTest(β1 ~ Temperature, data = data_stats)
print(levene_result)

# Bartlett检验 (对正态分布数据)
cat("\n=== 方差齐性检验 (Bartlett's test) ===\n")
bartlett_result <- bartlett.test(β1 ~ Temperature, data = data_stats)
print(bartlett_result)

# 使用 Wilcoxon 秩和检验 (Mann-Whitney U test) - 适用于非正态分布数据
cat("\n=== Wilcoxon 秩和检验 (Mann-Whitney U test) ===\n")
wilcox_result <- wilcox.test(β1 ~ Temperature, data = data_stats, exact = FALSE)
print(wilcox_result)
cat("W 统计量:", wilcox_result$statistic, "\n")
cat("p 值:", wilcox_result$p.value, "\n")

# 根据 Wilcoxon 检验结果确定显著性标记
wilcox_annotation <- ifelse(wilcox_result$p.value < 0.001, "***",
                            ifelse(wilcox_result$p.value < 0.01, "**",
                                  ifelse(wilcox_result$p.value < 0.05, "*", "ns")))

# 创建箱线图（包含所有温度10、16、22）
p <- ggplot(data, aes(x = Temperature, y = β1, fill = Temperature)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(width = 0.5)) +
  # 标准的 geom_boxplot，它会画出箱体和简洁的箱须
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, position = position_dodge(width = 0.5), 
               width = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  scale_fill_manual(values = c("10" = "blue", "16" = "yellow", "22" = "red")) +
  # 添加显著性标记（只标注10和16之间）
  geom_signif(comparisons = list(c("10", "16")),
              annotations = wilcox_annotation,
              y_position = c(max(data$β1) * 1.1),
              tip_length = 0.02,
              textsize = 6) +
  labs(
    y = "b1/TPM",
    fill = "Tem/°C"
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    # 图例设置
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    # 去除x轴标题
    axis.title.x = element_blank(),
    # 显示完整的边框
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    # 设置刻度线朝外
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    # 确保刻度线朝外
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5)
  ) +
  # 设置y轴范围
  scale_y_continuous(limits = c(0, max(data$β1, na.rm = TRUE) * 1.4))


# 显示图形
print(p)

# 保存图形
ggsave("temperature_b1_tpm_boxplot.pdf", plot = p, width = 12, height = 9, dpi = 1200)

# 进行方差分析（仅供参考，因为假设不满足）
# cat("\n=== 方差分析结果（仅供参考，假设不满足） ===\n")
# aov_result <- aov(β1_TPM.. ~ Temperature, data = data)
# print(summary(aov_result))

# 进行多重比较检验
# cat("\n=== 多重比较检验 ===\n")
# tukey_result <- TukeyHSD(aov_result)
# print(tukey_result)
