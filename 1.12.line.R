# 加载必要的包
library(ggplot2)
library(dplyr)

# 读取数据
data <- read.table("Alt.RCA.tsv", sep="\t", header=TRUE)

# 查看数据结构
str(data)
head(data)

# Lat-TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(Lat) & !is.na(TPM)) %>%
  filter(Lat >= 30)  # 只保留纬度≥30的数据

# 创建纬度区间分组（以2度为区间）
data_clean$Lat_group <- cut(data_clean$Lat, 
                           breaks = seq(30, ceiling(max(data_clean$Lat, na.rm = TRUE)/2)*2, by = 2),
                           include.lowest = TRUE,
                           right = FALSE)

# 计算每个纬度区间的TPM平均值
data_summary <- data_clean %>%
  group_by(Lat_group) %>%
  summarise(
    mean_TPM = mean(TPM, na.rm = TRUE),
    mean_Lat = mean(Lat, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(mean_Lat)

# 绘制折线图
p <- ggplot(data_summary, aes(x = mean_Lat, y = mean_TPM)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    x = "Latitude",
    y = "Tpm"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_blank(),  # 移除默认的轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = seq(30, 65, by = 2)) +
  scale_y_continuous(expand = c(0.02, 0))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("TPM_vs_Latitude.pdf", p, width = 10, height = 6, dpi = 1200)

# 显示统计信息
print("各纬度区间的统计信息:")
print(data_summary)

# 显示原始数据的纬度范围
cat("原始数据纬度范围:", min(data_clean$Lat, na.rm = TRUE), "到", max(data_clean$Lat, na.rm = TRUE), "\n")
cat("总数据点数:", nrow(data_clean), "\n")


# Lat-ATPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(Lat) & !is.na(α_TPM..)) %>%
  filter(Lat >= 30)  # 只保留纬度≥30的数据

# 创建纬度区间分组（以2度为区间）
data_clean$Lat_group <- cut(data_clean$Lat, 
                           breaks = seq(30, ceiling(max(data_clean$Lat, na.rm = TRUE)/2)*2, by = 2),
                           include.lowest = TRUE,
                           right = FALSE)

# 计算每个纬度区间的TPM平均值
data_summary <- data_clean %>%
  group_by(Lat_group) %>%
  summarise(
    mean_TPM = mean(α_TPM.., na.rm = TRUE),
    mean_Lat = mean(Lat, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(mean_Lat)

# 绘制折线图
p <- ggplot(data_summary, aes(x = mean_Lat, y = mean_TPM)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    x = "Latitude",
    y = "aTpm"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_blank(),  # 移除默认的轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = seq(30, 65, by = 2)) +
  scale_y_continuous(expand = c(0.02, 0))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("A-TPM_vs_Latitude.pdf", p, width = 10, height = 6, dpi = 1200)

# 显示统计信息
print("各纬度区间的统计信息:")
print(data_summary)

# 显示原始数据的纬度范围
cat("原始数据纬度范围:", min(data_clean$Lat, na.rm = TRUE), "到", max(data_clean$Lat, na.rm = TRUE), "\n")
cat("总数据点数:", nrow(data_clean), "\n")

# Lat-BTPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(Lat) & !is.na(β_TPM..)) %>%
  filter(Lat >= 30)  # 只保留纬度≥30的数据

# 创建纬度区间分组（以2度为区间）
data_clean$Lat_group <- cut(data_clean$Lat,
                            breaks = seq(30, ceiling(max(data_clean$Lat, na.rm = TRUE)/2)*2, by = 2),
                            include.lowest = TRUE,
                            right = FALSE)
                  
# 计算每个纬度区间的TPM平均值
data_summary <- data_clean %>%
  group_by(Lat_group) %>%
  summarise(
    mean_TPM = mean(β_TPM.., na.rm = TRUE),
    mean_Lat = mean(Lat, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(mean_Lat)

# 绘制折线图
p <- ggplot(data_summary, aes(x = mean_Lat, y = mean_TPM)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    x = "Latitude",
    y = "bTpm"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_blank(),  # 移除默认的轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = seq(30, 65, by = 2)) +
  scale_y_continuous(expand = c(0.02, 0))

# 显示图形  
print(p)

# 保存图形为PDF格式
ggsave("B-TPM_vs_Latitude.pdf", p, width = 10, height = 6, dpi = 1200)

# 显示统计信息
print("各纬度区间的统计信息:")
print(data_summary)

# 显示原始数据的纬度范围
cat("原始数据纬度范围:", min(data_clean$Lat, na.rm = TRUE), "到", max(data_clean$Lat, na.rm = TRUE), "\n")
cat("总数据点数:", nrow(data_clean), "\n")

# Lat-B1TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(Lat) & !is.na(β1_TPM..)) %>%
  filter(Lat >= 30)  # 只保留纬度≥30的数据

# 创建纬度区间分组（以2度为区间）
data_clean$Lat_group <- cut(data_clean$Lat, 
                           breaks = seq(30, ceiling(max(data_clean$Lat, na.rm = TRUE)/2)*2, by = 2),
                           include.lowest = TRUE,
                           right = FALSE)

# 计算每个纬度区间的TPM平均值
data_summary <- data_clean %>%
  group_by(Lat_group) %>%
  summarise(
    mean_TPM = mean(β1_TPM.., na.rm = TRUE),
    mean_Lat = mean(Lat, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(mean_Lat)

# 绘制折线图
p <- ggplot(data_summary, aes(x = mean_Lat, y = mean_TPM)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    x = "Latitude",
    y = "b1Tpm"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_blank(),  # 移除默认的轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = seq(30, 65, by = 2)) +
  scale_y_continuous(expand = c(0.02, 0))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B1-TPM_vs_Latitude.pdf", p, width = 10, height = 6, dpi = 1200)

# 显示统计信息
print("各纬度区间的统计信息:")
print(data_summary)

# 显示原始数据的纬度范围
cat("原始数据纬度范围:", min(data_clean$Lat, na.rm = TRUE), "到", max(data_clean$Lat, na.rm = TRUE), "\n")
cat("总数据点数:", nrow(data_clean), "\n")


# Lat-B2TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(Lat) & !is.na(β2_TPM..)) %>%
  filter(Lat >= 30)  # 只保留纬度≥30的数据

# 创建纬度区间分组（以2度为区间）
data_clean$Lat_group <- cut(data_clean$Lat, 
                           breaks = seq(30, ceiling(max(data_clean$Lat, na.rm = TRUE)/2)*2, by = 2),
                           include.lowest = TRUE,
                           right = FALSE)

# 计算每个纬度区间的TPM平均值
data_summary <- data_clean %>%
  group_by(Lat_group) %>%
  summarise(
    mean_TPM = mean(β2_TPM.., na.rm = TRUE),
    mean_Lat = mean(Lat, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(mean_Lat)

# 绘制折线图
p <- ggplot(data_summary, aes(x = mean_Lat, y = mean_TPM)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    x = "Latitude",
    y = "b2Tpm"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_blank(),  # 移除默认的轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = seq(30, 65, by = 2)) +
  scale_y_continuous(expand = c(0.02, 0))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_Latitude.pdf", p, width = 10, height = 6, dpi = 1200)

# 显示统计信息
print("各纬度区间的统计信息:")
print(data_summary)

# 显示原始数据的纬度范围
cat("原始数据纬度范围:", min(data_clean$Lat, na.rm = TRUE), "到", max(data_clean$Lat, na.rm = TRUE), "\n")
cat("总数据点数:", nrow(data_clean), "\n")

# Long-TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(Long) & !is.na(TPM))

# 创建经度区间分组（以2度为区间）
data_clean$Long_group <- cut(data_clean$Long, 
                           breaks = seq(-120, ceiling(max(data_clean$Long, na.rm = TRUE)/10)*10, by = 10),
                           include.lowest = TRUE,
                           right = FALSE)

# 计算每个经度区间的TPM平均值
data_summary <- data_clean %>%
  group_by(Long_group) %>%
  summarise(
    mean_TPM = mean(TPM, na.rm = TRUE),
    mean_Long = mean(Long, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(mean_Long)

# 绘制折线图
p <- ggplot(data_summary, aes(x = mean_Long, y = mean_TPM)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    x = "Longitude",
    y = "Tpm"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_blank(),  # 移除默认的轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = seq(-120, 180, by = 10)) +
  scale_y_continuous(expand = c(0.02, 0))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("TPM_vs_Longitude.pdf", p, width = 10, height = 6, dpi = 1200)

# 显示统计信息
print("各经度区间的统计信息:")
print(data_summary)

# 显示原始数据的经度范围
cat("原始数据经度范围:", min(data_clean$Long, na.rm = TRUE), "到", max(data_clean$Long, na.rm = TRUE), "\n")
cat("总数据点数:", nrow(data_clean), "\n")


# Long-ATPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(Long) & !is.na(α_TPM..))

# 创建经度区间分组（以10度为区间）
data_clean$Long_group <- cut(data_clean$Long, 
                           breaks = seq(-120, ceiling(max(data_clean$Long, na.rm = TRUE)/10)*10, by = 10),
                           include.lowest = TRUE,
                           right = FALSE)

# 计算每个经度区间的TPM平均值
data_summary <- data_clean %>%
  group_by(Long_group) %>%
  summarise(
    mean_TPM = mean(α_TPM.., na.rm = TRUE),
    mean_Long = mean(Long, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(mean_Long)

# 绘制折线图
p <- ggplot(data_summary, aes(x = mean_Long, y = mean_TPM)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    x = "Longitude",
    y = "aTpm"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_blank(),  # 移除默认的轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = seq(-120, 180, by = 10)) +
  scale_y_continuous(expand = c(0.02, 0))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("A-TPM_vs_Longitude.pdf", p, width = 10, height = 6, dpi = 1200)

# 显示统计信息
print("各经度区间的统计信息:")
print(data_summary)

# 显示原始数据的经度范围
cat("原始数据经度范围:", min(data_clean$Long, na.rm = TRUE), "到", max(data_clean$Long, na.rm = TRUE), "\n")
cat("总数据点数:", nrow(data_clean), "\n")

# Long-BTPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(Long) & !is.na(β_TPM..))

# 创建经度区间分组（以10度为区间）
data_clean$Long_group <- cut(data_clean$Long,
                            breaks = seq(-120, ceiling(max(data_clean$Long, na.rm = TRUE)/10)*10, by = 10),
                            include.lowest = TRUE,
                            right = FALSE)

# 计算每个经度区间的TPM平均值
data_summary <- data_clean %>%
  group_by(Long_group) %>%
  summarise(
    mean_TPM = mean(β_TPM.., na.rm = TRUE),
    mean_Long = mean(Long, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(mean_Long)

# 绘制折线图
p <- ggplot(data_summary, aes(x = mean_Long, y = mean_TPM)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    x = "Longitude",
    y = "bTpm"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_blank(),  # 移除默认的轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = seq(-120, 180, by = 10)) +
  scale_y_continuous(expand = c(0.02, 0))

# 显示图形  
print(p)

# 保存图形为PDF格式
ggsave("B-TPM_vs_Longitude.pdf", p, width = 10, height = 6, dpi = 1200)

# 显示统计信息
print("各经度区间的统计信息:")
print(data_summary)

# 显示原始数据的经度范围
cat("原始数据经度范围:", min(data_clean$Long, na.rm = TRUE), "到", max(data_clean$Long, na.rm = TRUE), "\n")
cat("总数据点数:", nrow(data_clean), "\n")

# Long-B1TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(Long) & !is.na(β1_TPM..))

# 创建经度区间分组（以10度为区间）
data_clean$Long_group <- cut(data_clean$Long, 
                           breaks = seq(-120, ceiling(max(data_clean$Long, na.rm = TRUE)/10)*10, by = 10),
                           include.lowest = TRUE,
                           right = FALSE)

# 计算每个经度区间的TPM平均值
data_summary <- data_clean %>%
  group_by(Long_group) %>%
  summarise(
    mean_TPM = mean(β1_TPM.., na.rm = TRUE),
    mean_Long = mean(Long, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(mean_Long)

# 绘制折线图
p <- ggplot(data_summary, aes(x = mean_Long, y = mean_TPM)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    x = "Longitude",
    y = "b1Tpm"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_blank(),  # 移除默认的轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = seq(-120, 180, by = 10)) +
  scale_y_continuous(expand = c(0.02, 0))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B1-TPM_vs_Longitude.pdf", p, width = 10, height = 6, dpi = 1200)

# 显示统计信息
print("各经度区间的统计信息:")
print(data_summary)

# 显示原始数据的经度范围
cat("原始数据经度范围:", min(data_clean$Long, na.rm = TRUE), "到", max(data_clean$Long, na.rm = TRUE), "\n")
cat("总数据点数:", nrow(data_clean), "\n")


# Lat-B2TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(Long) & !is.na(β2_TPM..))

# 创建纬度区间分组（以10度为区间）
data_clean$Long_group <- cut(data_clean$Long, 
                           breaks = seq(-120, ceiling(max(data_clean$Long, na.rm = TRUE)/10)*10, by = 10),
                           include.lowest = TRUE,
                           right = FALSE)

# 计算每个纬度区间的TPM平均值
data_summary <- data_clean %>%
  group_by(Long_group) %>%
  summarise(
    mean_TPM = mean(β2_TPM.., na.rm = TRUE),
    mean_Long = mean(Long, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(mean_Long)

# 绘制折线图
p <- ggplot(data_summary, aes(x = mean_Long, y = mean_TPM)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    x = "Longitude",
    y = "b2Tpm"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_blank(),  # 移除默认的轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = seq(-120, 180, by = 10)) +
  scale_y_continuous(expand = c(0.02, 0))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_Longitude.pdf", p, width = 10, height = 6, dpi = 1200)

# 显示统计信息
print("各经度区间的统计信息:")
print(data_summary)

# 显示原始数据的经度范围
cat("原始数据经度范围:", min(data_clean$Long, na.rm = TRUE), "到", max(data_clean$Long, na.rm = TRUE), "\n")
cat("总数据点数:", nrow(data_clean), "\n")


# Alt-TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(alt.m) & !is.na(TPM))

# 创建海拔区间分组（以200米为区间）
data_clean$Alt_group <- cut(data_clean$alt.m, 
                           breaks = seq(-20, ceiling(max(data_clean$alt.m, na.rm = TRUE)/200)*200, by = 200),
                           include.lowest = TRUE,
                           right = FALSE)

# 计算每个海拔区间的TPM平均值
data_summary <- data_clean %>%
  group_by(Alt_group) %>%
  summarise(
    mean_TPM = mean(TPM, na.rm = TRUE),
    mean_Alt = mean(alt.m, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(mean_Alt)

# 绘制折线图
p <- ggplot(data_summary, aes(x = mean_Alt, y = mean_TPM)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    x = "Altitude",
    y = "Tpm"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_blank(),  # 移除默认的轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = seq(-20, 4200, by = 200)) +
  scale_y_continuous(expand = c(0.02, 0))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("TPM_vs_Altitude.pdf", p, width = 10, height = 6, dpi = 1200)

# 显示统计信息
print("各海拔区间的统计信息:")
print(data_summary)

# 显示原始数据的海拔范围
cat("原始数据海拔范围:", min(data_clean$alt.m, na.rm = TRUE), "到", max(data_clean$alt.m, na.rm = TRUE), "\n")
cat("总数据点数:", nrow(data_clean), "\n")


# Alt-ATPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(alt.m) & !is.na(α_TPM..))

# 创建海拔区间分组（以200米为区间）
data_clean$Alt_group <- cut(data_clean$alt.m, 
                           breaks = seq(-20, ceiling(max(data_clean$alt.m, na.rm = TRUE)/200)*200, by = 200),
                           include.lowest = TRUE,
                           right = FALSE)

# 计算每个海拔区间的TPM平均值
data_summary <- data_clean %>%
  group_by(Alt_group) %>%
  summarise(
    mean_TPM = mean(α_TPM.., na.rm = TRUE),
    mean_Alt = mean(alt.m, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(mean_Alt)

# 绘制折线图
p <- ggplot(data_summary, aes(x = mean_Alt, y = mean_TPM)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    x = "Altitude",
    y = "aTpm"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_blank(),  # 移除默认的轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = seq(-20, 4200, by = 200)) +
  scale_y_continuous(expand = c(0.02, 0))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("A-TPM_vs_Altitude.pdf", p, width = 10, height = 6, dpi = 1200)

# 显示统计信息
print("各海拔区间的统计信息:")
print(data_summary)

# 显示原始数据的海拔范围
cat("原始数据海拔范围:", min(data_clean$alt.m, na.rm = TRUE), "到", max(data_clean$alt.m, na.rm = TRUE), "\n")
cat("总数据点数:", nrow(data_clean), "\n")

# Alt-BTPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(alt.m) & !is.na(β_TPM..))

# 创建海拔区间分组（以200米为区间）
data_clean$Alt_group <- cut(data_clean$alt.m,
                            breaks = seq(-20, ceiling(max(data_clean$alt.m, na.rm = TRUE)/200)*200, by = 200),
                            include.lowest = TRUE,
                            right = FALSE)

# 计算每个海拔区间的TPM平均值
data_summary <- data_clean %>%
  group_by(Alt_group) %>%
  summarise(
    mean_TPM = mean(β_TPM.., na.rm = TRUE),
    mean_Alt = mean(alt.m, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(mean_Alt)

# 绘制折线图
p <- ggplot(data_summary, aes(x = mean_Alt, y = mean_TPM)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    x = "Altitude",
    y = "bTpm"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_blank(),  # 移除默认的轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = seq(-20, 4200, by = 200)) +
  scale_y_continuous(expand = c(0.02, 0))

# 显示图形  
print(p)

# 保存图形为PDF格式
ggsave("B-TPM_vs_Altitude.pdf", p, width = 10, height = 6, dpi = 1200)

# 显示统计信息
print("各海拔区间的统计信息:")
print(data_summary)

# 显示原始数据的海拔范围
cat("原始数据海拔范围:", min(data_clean$alt.m, na.rm = TRUE), "到", max(data_clean$alt.m, na.rm = TRUE), "\n")
cat("总数据点数:", nrow(data_clean), "\n")

# Alt-B1TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(alt.m) & !is.na(β1_TPM..))

# 创建海拔区间分组（以200米为区间）
data_clean$Alt_group <- cut(data_clean$alt.m, 
                           breaks = seq(-20, ceiling(max(data_clean$alt.m, na.rm = TRUE)/200)*200, by = 200),
                           include.lowest = TRUE,
                           right = FALSE)

# 计算每个海拔区间的TPM平均值
data_summary <- data_clean %>%
  group_by(Alt_group) %>%
  summarise(
    mean_TPM = mean(β1_TPM.., na.rm = TRUE),
    mean_Alt = mean(alt.m, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(mean_Alt)

# 绘制折线图
p <- ggplot(data_summary, aes(x = mean_Alt, y = mean_TPM)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    x = "Altitude",
    y = "b1Tpm"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_blank(),  # 移除默认的轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = seq(-20, 4200, by = 200)) +
  scale_y_continuous(expand = c(0.02, 0))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B1-TPM_vs_Altitude.pdf", p, width = 10, height = 6, dpi = 1200)

# 显示统计信息
print("各海拔区间的统计信息:")
print(data_summary)

# 显示原始数据的海拔范围
cat("原始数据海拔范围:", min(data_clean$alt.m, na.rm = TRUE), "到", max(data_clean$alt.m, na.rm = TRUE), "\n")
cat("总数据点数:", nrow(data_clean), "\n")


# Alt-B2TPM列
# 清理数据，去除缺失值
data_clean <- data %>%
  filter(!is.na(alt.m) & !is.na(β2_TPM..))

# 创建海拔区间分组（以200m为区间）
data_clean$Alt_group <- cut(data_clean$alt.m, 
                           breaks = seq(-20, ceiling(max(data_clean$alt.m, na.rm = TRUE)/200)*200, by = 200),
                           include.lowest = TRUE,
                           right = FALSE)

# 计算每个海拔区间的TPM平均值
data_summary <- data_clean %>%
  group_by(Alt_group) %>%
  summarise(
    mean_TPM = mean(β2_TPM.., na.rm = TRUE),
    mean_Alt = mean(alt.m, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(mean_Alt)

# 绘制折线图
p <- ggplot(data_summary, aes(x = mean_Alt, y = mean_TPM)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    x = "Altitude",
    y = "b2Tpm"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_blank(),  # 移除默认的轴线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = seq(-20, 4200, by = 200)) +
  scale_y_continuous(expand = c(0.02, 0))

# 显示图形
print(p)

# 保存图形为PDF格式
ggsave("B2-TPM_vs_Altitude.pdf", p, width = 10, height = 6, dpi = 1200)

# 显示统计信息
print("各海拔区间的统计信息:")
print(data_summary)

# 显示原始数据的海拔范围
cat("原始数据海拔范围:", min(data_clean$alt.m, na.rm = TRUE), "到", max(data_clean$alt.m, na.rm = TRUE), "\n")
cat("总数据点数:", nrow(data_clean), "\n")

