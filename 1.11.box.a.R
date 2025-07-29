# 设置编码
Sys.setlocale("LC_ALL", "en_US.UTF-8")
options(encoding = "UTF-8")

# 加载必要的包
library(ggplot2)
library(dplyr)
library(readr)
library(ggsignif)

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
summary(data$α_TPM..)

# 创建箱线图
p <- ggplot(data, aes(x = Temperature, y = α_TPM.., fill = Temperature)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(width = 0.5)) +
  # 标准的 geom_boxplot，它会画出箱体和简洁的箱须
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, position = position_dodge(width = 0.5), 
               width = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  scale_fill_manual(values = c("10" = "blue", "16" = "yellow", "22" = "red")) +
  # 添加显著性标记
  geom_signif(comparisons = list(c("10", "16"), c("16", "22"), c("10", "22")),
              annotations = c("***", "***", "***"),
              y_position = c(max(data$α_TPM..) * 1.1, 
                           max(data$α_TPM..) * 1.2, 
                           max(data$α_TPM..) * 1.3),
              tip_length = 0.02,
              textsize = 4) +
  labs(
    y = "aTpm",
    fill = "Tem/°C"
  ) +
  theme_minimal() +
  theme(
    # 去除所有网格线
    panel.grid = element_blank(),
    # 设置坐标轴标题和文本
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(size = 8),
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
  scale_y_continuous(limits = c(0, max(data$α_TPM.., na.rm = TRUE) * 1.4))


# 显示图形
print(p)

# 保存图形
ggsave("temperature_a_tpm_boxplot.pdf", plot = p, width = 10, height = 6, dpi = 1200)

# 统计摘要
cat("\n=== 各温度组的αTPM统计摘要 ===\n")
summary_stats <- data %>%
  group_by(Temperature) %>%
  summarise(
    count = n(),
    mean = mean(α_TPM.., na.rm = TRUE),
    median = median(α_TPM.., na.rm = TRUE),
    sd = sd(α_TPM.., na.rm = TRUE),
    min = min(α_TPM.., na.rm = TRUE),
    max = max(α_TPM.., na.rm = TRUE),
    q25 = quantile(α_TPM.., 0.25, na.rm = TRUE),
    q75 = quantile(α_TPM.., 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_stats)

# 进行方差分析
cat("\n=== 方差分析结果 ===\n")
aov_result <- aov(α_TPM.. ~ Temperature, data = data)
print(summary(aov_result))

# 进行多重比较检验
cat("\n=== 多重比较检验 ===\n")
tukey_result <- TukeyHSD(aov_result)
print(tukey_result)
