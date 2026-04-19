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
unique(data$`Temperature`)

# 将温度列转换为因子，确保正确排序（包含10、16、22）
data$Temperature <- factor(data$`Temperature`, levels = c("10", "16", "22"))

# 对total进行log1p转换以避免log(0)
data$total <- log1p(data$total)

# 检查数据
table(data$Temperature)
summary(data$total)

# 统计摘要
cat("\n=== 各温度组的total TPM统计摘要 ===\n")
summary_stats <- data %>%
  group_by(Temperature) %>%
  summarise(
    count = n(),
    mean = mean(total, na.rm = TRUE),
    median = median(total, na.rm = TRUE),
    sd = sd(total, na.rm = TRUE),
    min = min(total, na.rm = TRUE),
    max = max(total, na.rm = TRUE),
    q25 = quantile(total, 0.25, na.rm = TRUE),
    q75 = quantile(total, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_stats)

median_output <- summary_stats %>%
  select(Temperature, median) %>%
  rename(median_value = median)
write_tsv(median_output, "temperature_total_tpm_boxplot_median.tsv")

cat("\n=== 正态性检验 (Shapiro-Wilk test) ===\n")
shapiro_results <- data %>%
  group_by(Temperature) %>%
  summarise(
    shapiro_p = shapiro.test(total)$p.value,
    .groups = 'drop'
  )
print(shapiro_results)

temp_aov <- aov(total ~ Temperature, data = data)
cat("\n=== 残差正态性检验 ===\n")
residuals_shapiro <- shapiro.test(residuals(temp_aov))
cat("残差Shapiro-Wilk检验 p值:", residuals_shapiro$p.value, "\n")

cat("\n=== 方差齐性检验 (Levene's test) ===\n")
levene_result <- leveneTest(total ~ Temperature, data = data)
print(levene_result)

cat("\n=== 方差齐性检验 (Bartlett's test) ===\n")
bartlett_result <- bartlett.test(total ~ Temperature, data = data)
print(bartlett_result)

if (all(shapiro_results$shapiro_p > 0.05) && residuals_shapiro$p.value > 0.05 && levene_result$`Pr(>F)`[1] > 0.05) {
  cat("\n=== 单因素方差分析 (One-way ANOVA) ===\n")
  anova_result <- aov(total ~ Temperature, data = data)
  anova_summary <- summary(anova_result)
  print(anova_summary)
  overall_p <- anova_summary[[1]]$`Pr(>F)`[1]
  cat("\n=== Tukey HSD 多重比较检验 ===\n")
  tukey_result <- TukeyHSD(anova_result)
  print(tukey_result)
  comparison_pvals <- tukey_result$Temperature[, "p adj"]
  stat_method <- "ANOVA + Tukey HSD"
} else {
  cat("\n=== Kruskal-Wallis 秩和检验 ===\n")
  kruskal_result <- kruskal.test(total ~ Temperature, data = data)
  print(kruskal_result)
  overall_p <- kruskal_result$p.value
  cat("\n=== Dunn's 多重比较检验 ===\n")
  dunn_result <- dunn.test(data$total, data$Temperature, method="bonferroni")
  comparison_pvals <- dunn_result$P.adjusted
  stat_method <- "Kruskal-Wallis + Dunn's test"
}
cat("\n使用的统计方法:", stat_method, "\n")
cat("整体检验 overall P 值:", format.pval(overall_p, digits = 3), "\n")

get_significance <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "n.s.")))
}

comparisons <- list(c("10", "16"), c("16", "22"), c("10", "22"))
comparison_labels <- c("10_vs_16", "16_vs_22", "10_vs_22")

pairwise_pvals <- tibble(
  comparison = comparison_labels,
  adjusted_p_value = as.numeric(comparison_pvals)
)
write_tsv(pairwise_pvals, "temperature_total_tpm_boxplot_pairwise_pvalues.tsv")

significance_labels <- get_significance(comparison_pvals)
y_positions <- c(max(data$total, na.rm = TRUE) * 1.1,
                 max(data$total, na.rm = TRUE) * 1.2,
                 max(data$total, na.rm = TRUE) * 1.3)

# 创建箱线图（包含温度10、16、22）
p <- ggplot(data, aes(x = Temperature, y = total)) +
  stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(width = 0.5)) +
  geom_boxplot(fill = "#4472C4", color = "black", alpha = 1, outlier.shape = 16, outlier.size = 1, 
               position = position_dodge(width = 0.5), width = 0.5) +
  labs(
    y = expression("total"~"(ln(TPM+1))"),
    x = expression(Temperature~"("*degree*C*")")
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 26),
    axis.text.x = element_text(size = 26),
    # 图例设置
    legend.title = element_text(size = 28),
    legend.text = element_text(size = 26),
    # 去除右侧颜色图例
    legend.position = "none",
    # 显示完整的边框
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    # 设置刻度线朝外
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    # 确保刻度线朝外
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5)
  ) +
  # 设置y轴范围
  scale_y_continuous(
    limits = c(5, max(data$total, na.rm = TRUE) * 1.4),
    breaks = scales::pretty_breaks(n = 8)
  )

p <- p + geom_signif(
  comparisons = comparisons,
  annotations = significance_labels,
  y_position = y_positions,
  tip_length = 0.02,
  textsize = 9
)


# 显示图形
print(p)

# 保存图形
ggsave("temperature_total_tpm_boxplot.pdf", plot = p, width = 12, height = 9, dpi = 1200)

# 进行方差分析（仅供参考，因为假设不满足）
# cat("\n=== 方差分析结果（仅供参考，假设不满足） ===\n")
# aov_result <- aov(TPM ~ Temperature, data = data)
# print(summary(aov_result))

# 进行多重比较检验
# cat("\n=== 多重比较检验 ===\n")
# tukey_result <- TukeyHSD(aov_result)
# print(tukey_result)
