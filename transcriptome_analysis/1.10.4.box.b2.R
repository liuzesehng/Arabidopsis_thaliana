# 设置编码
Sys.setlocale("LC_ALL", "en_US.UTF-8")
options(encoding = "UTF-8")

# 加载必要的包
library(ggplot2)
library(dplyr)
library(readr)
library(ggsignif)
library(car)
library(dunn.test)

# 读取数据
data <- read.table("Alt.RCA.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 查看数据结构
str(data)
head(data)

# 检查温度列的唯一值
unique(data$`Tem..`)

# 将温度列转换为因子，确保正确排序
data$Temperature <- factor(data$`Tem..`, levels = c("10", "16", "22"))

# 检查数据
table(data$Temperature)
summary(data$β2_TPM..)

# 创建箱线图
p <- ggplot(data, aes(x = Temperature, y = β2_TPM.., fill = Temperature)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(width = 0.5)) +
  # 标准的 geom_boxplot，它会画出箱体和简洁的箱须
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, position = position_dodge(width = 0.5), 
               width = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  scale_fill_manual(values = c("10" = "blue", "16" = "yellow", "22" = "red")) +
  # 添加显著性标记
  geom_signif(comparisons = list(c("10", "16"), c("16", "22"), c("10", "22")),
              annotations = c("ns", "***", "***"),
              y_position = c(max(data$β2_TPM..) * 1.1, 
                           max(data$β2_TPM..) * 1.2, 
                           max(data$β2_TPM..) * 1.3),
              tip_length = 0.02,
              textsize = 6) +
  labs(
    y = "b2Tpm",
    fill = "Tem/°C"
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    # 图例设置
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
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
  scale_y_continuous(limits = c(0, max(data$β2_TPM.., na.rm = TRUE) * 1.4))


# 显示图形
print(p)

# 保存图形
ggsave("temperature_b2_tpm_boxplot.pdf", plot = p, width = 12, height = 9, dpi = 1200)

# 统计摘要
cat("\n=== 各温度组的β2TPM统计摘要 ===\n")
summary_stats <- data %>%
  group_by(Temperature) %>%
  summarise(
    count = n(),
    mean = mean(β2_TPM.., na.rm = TRUE),
    median = median(β2_TPM.., na.rm = TRUE),
    sd = sd(β2_TPM.., na.rm = TRUE),
    min = min(β2_TPM.., na.rm = TRUE),
    max = max(β2_TPM.., na.rm = TRUE),
    q25 = quantile(β2_TPM.., 0.25, na.rm = TRUE),
    q75 = quantile(β2_TPM.., 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_stats)

# 进行正态性检验和方差齐性检验
cat("\n=== 正态性检验 (Shapiro-Wilk test) ===\n")
shapiro_results <- data %>%
  group_by(Temperature) %>%
  summarise(
    shapiro_p = shapiro.test(β2_TPM..)$p.value,
    .groups = 'drop'
  )
print(shapiro_results)

# 对残差进行正态性检验
temp_aov <- aov(β2_TPM.. ~ Temperature, data = data)
cat("\n=== 残差正态性检验 ===\n")
residuals_shapiro <- shapiro.test(residuals(temp_aov))
cat("残差Shapiro-Wilk检验 p值:", residuals_shapiro$p.value, "\n")

# 方差齐性检验 (Levene's test)
cat("\n=== 方差齐性检验 (Levene's test) ===\n")
levene_result <- leveneTest(β2_TPM.. ~ Temperature, data = data)
print(levene_result)

# Bartlett检验 (对正态分布数据)
cat("\n=== 方差齐性检验 (Bartlett's test) ===\n")
bartlett_result <- bartlett.test(β2_TPM.. ~ Temperature, data = data)
print(bartlett_result)

# 由于数据不符合正态性和方差齐性假设，使用非参数检验
cat("\n=== 非参数检验：Kruskal-Wallis检验 ===\n")
kruskal_result <- kruskal.test(β2_TPM.. ~ Temperature, data = data)
print(kruskal_result)

# 如果Kruskal-Wallis检验显著，进行事后检验
if(kruskal_result$p.value < 0.05) {
  cat("\n=== 非参数多重比较：Dunn检验 ===\n")
  dunn_result <- dunn.test(data$β2_TPM.., data$Temperature, method = "bonferroni")
  print(dunn_result)
  
} else {
  cat("Kruskal-Wallis检验不显著，无需进行事后检验\n")
}

# 进行方差分析（仅供参考，因为假设不满足）
# cat("\n=== 方差分析结果（仅供参考，假设不满足） ===\n")
# aov_result <- aov(β2_TPM.. ~ Temperature, data = data)
# print(summary(aov_result))

# 进行多重比较检验
# cat("\n=== 多重比较检验 ===\n")
# tukey_result <- TukeyHSD(aov_result)
# print(tukey_result)

