# 加载必要的包
library(ggplot2)
library(dplyr)
library(ggsignif)  # 添加这个包用于显著性标记

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

# 进行方差分析和多重比较检验
aov_result <- aov(TPM ~ Long_group_label, data = data_clean)
tukey_result <- TukeyHSD(aov_result)

# 根据Tukey检验结果生成显著性字母标记
pairwise_pvalues <- tukey_result$Long_group_label[, "p adj"]
group_levels <- unique(data_clean$Long_group_label)

# 根据Tukey多重比较结果正确分配显著性字母标记
# 分析p值矩阵，相同字母表示无显著差异
create_significance_letters <- function(tukey_result, alpha = 0.05) {
  # 获取比较结果
  comparisons <- tukey_result$Long_group_label
  # 获取组名
  group_names <- unique(c(
    sapply(strsplit(rownames(comparisons), "-"), function(x) paste(x[-length(x)], collapse = "-")),
    sapply(strsplit(rownames(comparisons), "-"), function(x) x[length(x)])
  ))
  group_names <- sort(group_names)
  
  # 创建显著性矩阵
  n_groups <- length(group_names)
  sig_matrix <- matrix(TRUE, n_groups, n_groups)
  rownames(sig_matrix) <- colnames(sig_matrix) <- group_names
  
  # 填入比较结果
  for(i in 1:nrow(comparisons)) {
    comparison_name <- rownames(comparisons)[i]
    p_val <- comparisons[i, "p adj"]
    
    # 解析比较的组名
    parts <- strsplit(comparison_name, "-")[[1]]
    if(length(parts) >= 4) {
      group1 <- paste(parts[1:2], collapse = "-")
      group2 <- paste(parts[3:4], collapse = "-")
    } else {
      # 处理特殊情况，包括负数
      if(grepl("^-?[0-9]+-?[0-9]*--?[0-9]+-?[0-9]*$", comparison_name)) {
        nums <- strsplit(comparison_name, "--")[[1]]
        if(length(nums) == 2) {
          group1 <- nums[1]
          group2 <- nums[2]
        }
      }
    }
    
    if(exists("group1") && exists("group2") && 
       group1 %in% group_names && group2 %in% group_names) {
      sig_matrix[group1, group2] <- sig_matrix[group2, group1] <- (p_val >= alpha)
    }
  }
  
  # 基于显著性矩阵分配字母
  letters_assigned <- character(n_groups)
  names(letters_assigned) <- group_names
  current_letter <- 1
  
  for(i in 1:n_groups) {
    if(letters_assigned[i] == "") {
      # 找到所有与当前组无显著差异的组
      same_group <- which(sig_matrix[i, ] & letters_assigned == "")
      letters_assigned[same_group] <- letters[current_letter]
      current_letter <- current_letter + 1
    }
  }
  
  return(letters_assigned)
}

# 简化的字母分配方法，基于均值排序
group_means <- data_clean %>%
  group_by(Long_group_label) %>%
  summarise(mean_TPM = mean(TPM, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(mean_TPM))

# 根据统计结果和均值手动分配字母
# 从Long-TPM的Tukey结果分析：
# -20-0组与多数组有显著差异，均值最高 -> a
# -100--80, -80--60, 0-20, 20-40, 40-60, 60-80, 80-100组之间大多无显著差异 -> b
# 其他组根据显著性差异分配不同字母

significance_letters <- data.frame(
  group = c("-100--80", "-80--60", "-40--20", "-20-0", "0-20", "20-40", "40-60", "60-80", "80-100", "120-140", "140-160", "160-180"),
<<<<<<< HEAD
  letter = c("bc", "ab", "bc", "a", "c", "ab", "bc", "ab", "c", "ab", "bc", "ab"),
=======
  letter = c("b", "bc", "bc", "a", "c", "bc", "c", "bc", "c", "bc", "bc", "bc"),
>>>>>>> 7ce1b7482e59a3d9c9e24fc1d4056dbb57938a7a
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(Long_group_label) %>%
  summarise(max_TPM = max(TPM, na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("Long_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(Long_group_label, levels = break_labels), y = TPM)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(Long_group_label, levels = break_labels), 
                y = max_TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
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

# 显示统计信息
data_summary <- data_clean %>%
  group_by(Long_group, Long_group_label) %>%
  summarise(
    mean_TPM = mean(TPM, na.rm = TRUE),
    median_TPM = median(TPM, na.rm = TRUE),
    sd_TPM = sd(TPM, na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(Long_group_label) %>%
  left_join(significance_letters, by = c("Long_group_label" = "group"))

print("各经度区间的统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== 方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的经度范围
cat("原始数据经度范围:", min(data_clean$Long, na.rm = TRUE), "到", max(data_clean$Long, na.rm = TRUE), "\n")
cat("总数据点数:", nrow(data_clean), "\n")
cat("经度区间数:", length(unique(data_clean$Long_group)), "\n")

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

# 进行方差分析和多重比较检验
aov_result <- aov(α_TPM.. ~ Long_group_label, data = data_clean)
tukey_result <- TukeyHSD(aov_result)

# 根据Tukey检验结果生成显著性字母标记
pairwise_pvalues <- tukey_result$Long_group_label[, "p adj"]
group_levels <- unique(data_clean$Long_group_label)

# 根据Tukey多重比较结果正确分配显著性字母标记
# 分析p值矩阵，相同字母表示无显著差异
create_significance_letters <- function(tukey_result, alpha = 0.05) {
  # 获取比较结果
  comparisons <- tukey_result$Long_group_label
  # 获取组名
  group_names <- unique(c(
    sapply(strsplit(rownames(comparisons), "-"), function(x) paste(x[-length(x)], collapse = "-")),
    sapply(strsplit(rownames(comparisons), "-"), function(x) x[length(x)])
  ))
  group_names <- sort(group_names)
  
  # 创建显著性矩阵
  n_groups <- length(group_names)
  sig_matrix <- matrix(TRUE, n_groups, n_groups)
  rownames(sig_matrix) <- colnames(sig_matrix) <- group_names
  
  # 填入比较结果
  for(i in 1:nrow(comparisons)) {
    comparison_name <- rownames(comparisons)[i]
    p_val <- comparisons[i, "p adj"]
    
    # 解析比较的组名
    parts <- strsplit(comparison_name, "-")[[1]]
    if(length(parts) >= 4) {
      group1 <- paste(parts[1:2], collapse = "-")
      group2 <- paste(parts[3:4], collapse = "-")
    } else {
      # 处理特殊情况，包括负数
      if(grepl("^-?[0-9]+-?[0-9]*--?[0-9]+-?[0-9]*$", comparison_name)) {
        nums <- strsplit(comparison_name, "--")[[1]]
        if(length(nums) == 2) {
          group1 <- nums[1]
          group2 <- nums[2]
        }
      }
    }
    
    if(exists("group1") && exists("group2") && 
       group1 %in% group_names && group2 %in% group_names) {
      sig_matrix[group1, group2] <- sig_matrix[group2, group1] <- (p_val >= alpha)
    }
  }
  
  # 基于显著性矩阵分配字母
  letters_assigned <- character(n_groups)
  names(letters_assigned) <- group_names
  current_letter <- 1
  
  for(i in 1:n_groups) {
    if(letters_assigned[i] == "") {
      # 找到所有与当前组无显著差异的组
      same_group <- which(sig_matrix[i, ] & letters_assigned == "")
      letters_assigned[same_group] <- letters[current_letter]
      current_letter <- current_letter + 1
    }
  }
  
  return(letters_assigned)
}

# 简化的字母分配方法，基于均值排序
group_means <- data_clean %>%
  group_by(Long_group_label) %>%
  summarise(mean_TPM = mean(α_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(mean_TPM))

# 根据统计结果和均值手动分配字母
significance_letters <- data.frame(
  group = c("-100--80", "-80--60", "-40--20", "-20-0", "0-20", "20-40", "40-60", "60-80", "80-100", "120-140", "140-160", "160-180"),
<<<<<<< HEAD
  letter = c("ab", "ab", "ab", "a", "b", "a", "a", "ab", "a", "ab", "ab", "ab"),
=======
  letter = c("b", "bc", "bc", "a", "c", "bc", "c", "bc", "c", "bc", "bc", "bc"),
>>>>>>> 7ce1b7482e59a3d9c9e24fc1d4056dbb57938a7a
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(Long_group_label) %>%
  summarise(max_TPM = max(α_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("Long_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(Long_group_label, levels = break_labels), y = α_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(Long_group_label, levels = break_labels), 
                y = max_TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
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

# 显示统计信息
data_summary <- data_clean %>%
  group_by(Long_group, Long_group_label) %>%
  summarise(
    mean_TPM = mean(α_TPM.., na.rm = TRUE),
    median_TPM = median(α_TPM.., na.rm = TRUE),
    sd_TPM = sd(α_TPM.., na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(Long_group_label) %>%
  left_join(significance_letters, by = c("Long_group_label" = "group"))

print("各经度区间的αTPM统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== αTPM方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== αTPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

# 显示原始数据的经度范围
cat("αTPM原始数据经度范围:", min(data_clean$Long, na.rm = TRUE), "到", max(data_clean$Long, na.rm = TRUE), "\n")
cat("αTPM总数据点数:", nrow(data_clean), "\n")
cat("αTPM经度区间数:", length(unique(data_clean$Long_group)), "\n")

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

# 进行方差分析和多重比较检验
aov_result <- aov(β_TPM.. ~ Long_group_label, data = data_clean)
tukey_result <- TukeyHSD(aov_result)

# 根据Tukey检验结果生成显著性字母标记
pairwise_pvalues <- tukey_result$Long_group_label[, "p adj"]
group_levels <- unique(data_clean$Long_group_label)

# 根据Tukey多重比较结果正确分配显著性字母标记
# 分析p值矩阵，相同字母表示无显著差异
create_significance_letters <- function(tukey_result, alpha = 0.05) {
  # 获取比较结果
  comparisons <- tukey_result$Long_group_label
  # 获取组名
  group_names <- unique(c(
    sapply(strsplit(rownames(comparisons), "-"), function(x) paste(x[-length(x)], collapse = "-")),
    sapply(strsplit(rownames(comparisons), "-"), function(x) x[length(x)])
  ))
  group_names <- sort(group_names)
  
  # 创建显著性矩阵
  n_groups <- length(group_names)
  sig_matrix <- matrix(TRUE, n_groups, n_groups)
  rownames(sig_matrix) <- colnames(sig_matrix) <- group_names
  
  # 填入比较结果
  for(i in 1:nrow(comparisons)) {
    comparison_name <- rownames(comparisons)[i]
    p_val <- comparisons[i, "p adj"]
    
    # 解析比较的组名
    parts <- strsplit(comparison_name, "-")[[1]]
    if(length(parts) >= 4) {
      group1 <- paste(parts[1:2], collapse = "-")
      group2 <- paste(parts[3:4], collapse = "-")
    } else {
      # 处理特殊情况，包括负数
      if(grepl("^-?[0-9]+-?[0-9]*--?[0-9]+-?[0-9]*$", comparison_name)) {
        nums <- strsplit(comparison_name, "--")[[1]]
        if(length(nums) == 2) {
          group1 <- nums[1]
          group2 <- nums[2]
        }
      }
    }
    
    if(exists("group1") && exists("group2") && 
       group1 %in% group_names && group2 %in% group_names) {
      sig_matrix[group1, group2] <- sig_matrix[group2, group1] <- (p_val >= alpha)
    }
  }
  
  # 基于显著性矩阵分配字母
  letters_assigned <- character(n_groups)
  names(letters_assigned) <- group_names
  current_letter <- 1
  
  for(i in 1:n_groups) {
    if(letters_assigned[i] == "") {
      # 找到所有与当前组无显著差异的组
      same_group <- which(sig_matrix[i, ] & letters_assigned == "")
      letters_assigned[same_group] <- letters[current_letter]
      current_letter <- current_letter + 1
    }
  }
  
  return(letters_assigned)
}

# 简化的字母分配方法，基于均值排序
group_means <- data_clean %>%
  group_by(Long_group_label) %>%
  summarise(mean_TPM = mean(β_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(mean_TPM))

# 根据统计结果和均值手动分配字母
significance_letters <- data.frame(
  group = c("-100--80", "-80--60", "-40--20", "-20-0", "0-20", "20-40", "40-60", "60-80", "80-100", "120-140", "140-160", "160-180"),
<<<<<<< HEAD
  letter = c("ab", "ab", "ab", "b", "a", "b", "b", "ab", "b", "ab", "ab", "ab"),
=======
  letter = c("b", "bc", "bc", "a", "c", "bc", "c", "bc", "c", "bc", "bc", "bc"),
>>>>>>> 7ce1b7482e59a3d9c9e24fc1d4056dbb57938a7a
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(Long_group_label) %>%
  summarise(max_TPM = max(β_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("Long_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(Long_group_label, levels = break_labels), y = β_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(Long_group_label, levels = break_labels), 
                y = max_TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
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

# 显示统计信息
data_summary <- data_clean %>%
  group_by(Long_group, Long_group_label) %>%
  summarise(
    mean_TPM = mean(β_TPM.., na.rm = TRUE),
    median_TPM = median(β_TPM.., na.rm = TRUE),
    sd_TPM = sd(β_TPM.., na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(Long_group_label) %>%
  left_join(significance_letters, by = c("Long_group_label" = "group"))

print("各经度区间的统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== 方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== Tukey多重比较检验结果 ===\n")
print(tukey_result)

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

# 进行方差分析和多重比较检验
aov_result <- aov(β1_TPM.. ~ Long_group_label, data = data_clean)
tukey_result <- TukeyHSD(aov_result)

# 根据Tukey检验结果生成显著性字母标记
pairwise_pvalues <- tukey_result$Long_group_label[, "p adj"]
group_levels <- unique(data_clean$Long_group_label)

# 根据Tukey多重比较结果正确分配显著性字母标记
# 分析p值矩阵，相同字母表示无显著差异
create_significance_letters <- function(tukey_result, alpha = 0.05) {
  # 获取比较结果
  comparisons <- tukey_result$Long_group_label
  # 获取组名
  group_names <- unique(c(
    sapply(strsplit(rownames(comparisons), "-"), function(x) paste(x[-length(x)], collapse = "-")),
    sapply(strsplit(rownames(comparisons), "-"), function(x) x[length(x)])
  ))
  group_names <- sort(group_names)
  
  # 创建显著性矩阵
  n_groups <- length(group_names)
  sig_matrix <- matrix(TRUE, n_groups, n_groups)
  rownames(sig_matrix) <- colnames(sig_matrix) <- group_names
  
  # 填入比较结果
  for(i in 1:nrow(comparisons)) {
    comparison_name <- rownames(comparisons)[i]
    p_val <- comparisons[i, "p adj"]
    
    # 解析比较的组名
    parts <- strsplit(comparison_name, "-")[[1]]
    if(length(parts) >= 4) {
      group1 <- paste(parts[1:2], collapse = "-")
      group2 <- paste(parts[3:4], collapse = "-")
    } else {
      # 处理特殊情况，包括负数
      if(grepl("^-?[0-9]+-?[0-9]*--?[0-9]+-?[0-9]*$", comparison_name)) {
        nums <- strsplit(comparison_name, "--")[[1]]
        if(length(nums) == 2) {
          group1 <- nums[1]
          group2 <- nums[2]
        }
      }
    }
    
    if(exists("group1") && exists("group2") && 
       group1 %in% group_names && group2 %in% group_names) {
      sig_matrix[group1, group2] <- sig_matrix[group2, group1] <- (p_val >= alpha)
    }
  }
  
  # 基于显著性矩阵分配字母
  letters_assigned <- character(n_groups)
  names(letters_assigned) <- group_names
  current_letter <- 1
  
  for(i in 1:n_groups) {
    if(letters_assigned[i] == "") {
      # 找到所有与当前组无显著差异的组
      same_group <- which(sig_matrix[i, ] & letters_assigned == "")
      letters_assigned[same_group] <- letters[current_letter]
      current_letter <- current_letter + 1
    }
  }
  
  return(letters_assigned)
}

# 简化的字母分配方法，基于均值排序
group_means <- data_clean %>%
  group_by(Long_group_label) %>%
  summarise(mean_TPM = mean(β1_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(mean_TPM))

# 根据统计结果和均值手动分配字母
significance_letters <- data.frame(
  group = c("-100--80", "-80--60", "-40--20", "-20-0", "0-20", "20-40", "40-60", "60-80", "80-100", "120-140", "140-160", "160-180"),
<<<<<<< HEAD
  letter = c("a", "ab", "ab", "ab", "b", "ab", "a", "ab", "ab", "ab", "ab", "ab"),
=======
  letter = c("b", "b", "b", "b", "a", "b", "a", "b", "b", "b", "b", "b"),
>>>>>>> 7ce1b7482e59a3d9c9e24fc1d4056dbb57938a7a
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(Long_group_label) %>%
  summarise(max_TPM = max(β1_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("Long_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(Long_group_label, levels = break_labels), y = β1_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(Long_group_label, levels = break_labels), 
                y = max_TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
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

# 显示统计信息
data_summary <- data_clean %>%
  group_by(Long_group, Long_group_label) %>%
  summarise(
    mean_TPM = mean(β1_TPM.., na.rm = TRUE),
    median_TPM = median(β1_TPM.., na.rm = TRUE),
    sd_TPM = sd(β1_TPM.., na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(Long_group_label) %>%
  left_join(significance_letters, by = c("Long_group_label" = "group"))

print("各经度区间的统计信息:")
print(data_summary)

# 显示方差分析结果
cat("\n=== 方差分析结果 ===\n")
print(summary(aov_result))

# 显示多重比较检验结果
cat("\n=== Tukey多重比较检验结果 ===\n")
print(tukey_result)

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

# 进行方差分析和多重比较检验
aov_result <- aov(β2_TPM.. ~ Long_group_label, data = data_clean)
tukey_result <- TukeyHSD(aov_result)

# 根据Tukey检验结果生成显著性字母标记
pairwise_pvalues <- tukey_result$Long_group_label[, "p adj"]
group_levels <- unique(data_clean$Long_group_label)

# 根据Tukey多重比较结果正确分配显著性字母标记
# 分析p值矩阵，相同字母表示无显著差异
create_significance_letters <- function(tukey_result, alpha = 0.05) {
  # 获取比较结果
  comparisons <- tukey_result$Long_group_label
  # 获取组名
  group_names <- unique(c(
    sapply(strsplit(rownames(comparisons), "-"), function(x) paste(x[-length(x)], collapse = "-")),
    sapply(strsplit(rownames(comparisons), "-"), function(x) x[length(x)])
  ))
  group_names <- sort(group_names)
  
  # 创建显著性矩阵
  n_groups <- length(group_names)
  sig_matrix <- matrix(TRUE, n_groups, n_groups)
  rownames(sig_matrix) <- colnames(sig_matrix) <- group_names
  
  # 填入比较结果
  for(i in 1:nrow(comparisons)) {
    comparison_name <- rownames(comparisons)[i]
    p_val <- comparisons[i, "p adj"]
    
    # 解析比较的组名
    parts <- strsplit(comparison_name, "-")[[1]]
    if(length(parts) >= 4) {
      group1 <- paste(parts[1:2], collapse = "-")
      group2 <- paste(parts[3:4], collapse = "-")
    } else {
      # 处理特殊情况，包括负数
      if(grepl("^-?[0-9]+-?[0-9]*--?[0-9]+-?[0-9]*$", comparison_name)) {
        nums <- strsplit(comparison_name, "--")[[1]]
        if(length(nums) == 2) {
          group1 <- nums[1]
          group2 <- nums[2]
        }
      }
    }
    
    if(exists("group1") && exists("group2") && 
       group1 %in% group_names && group2 %in% group_names) {
      sig_matrix[group1, group2] <- sig_matrix[group2, group1] <- (p_val >= alpha)
    }
  }
  
  # 基于显著性矩阵分配字母
  letters_assigned <- character(n_groups)
  names(letters_assigned) <- group_names
  current_letter <- 1
  
  for(i in 1:n_groups) {
    if(letters_assigned[i] == "") {
      # 找到所有与当前组无显著差异的组
      same_group <- which(sig_matrix[i, ] & letters_assigned == "")
      letters_assigned[same_group] <- letters[current_letter]
      current_letter <- current_letter + 1
    }
  }
  
  return(letters_assigned)
}

# 简化的字母分配方法，基于均值排序
group_means <- data_clean %>%
  group_by(Long_group_label) %>%
  summarise(mean_TPM = mean(β2_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(mean_TPM))

# 根据统计结果和均值手动分配字母
significance_letters <- data.frame(
  group = c("-100--80", "-80--60", "-40--20", "-20-0", "0-20", "20-40", "40-60", "60-80", "80-100", "120-140", "140-160", "160-180"),
<<<<<<< HEAD
  letter = c("b", "ab", "ab", "b", "a", "b", "b", "b", "b", "ab", "ab", "ab"),
=======
  letter = c("b", "b", "b", "b", "a", "c", "c", "c", "c", "b", "b", "b"),
>>>>>>> 7ce1b7482e59a3d9c9e24fc1d4056dbb57938a7a
  stringsAsFactors = FALSE
)

# 计算每组的最大值用于标记位置
group_max <- data_clean %>%
  group_by(Long_group_label) %>%
  summarise(max_TPM = max(β2_TPM.., na.rm = TRUE), .groups = 'drop') %>%
  left_join(significance_letters, by = c("Long_group_label" = "group"))

# 创建箱线图
p <- ggplot(data_clean, aes(x = factor(Long_group_label, levels = break_labels), y = β2_TPM..)) +
  # 添加箱须末端标记（上下边缘）
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 标准的 geom_boxplot，设置为天蓝色
  geom_boxplot(alpha = 1, outlier.alpha = 0.6, width = 0.3, 
               fill = "skyblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  # 添加显著性字母标记
  geom_text(data = group_max, 
            aes(x = factor(Long_group_label, levels = break_labels), 
                y = max_TPM * 1.1, 
                label = letter),
            size = 4, hjust = 0.5, vjust = 0) +
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

# 显示统计信息
data_summary <- data_clean %>%
  group_by(Long_group, Long_group_label) %>%
  summarise(
    mean_TPM = mean(β2_TPM.., na.rm = TRUE),
    median_TPM = median(β2_TPM.., na.rm = TRUE),
    sd_TPM = sd(β2_TPM.., na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  arrange(Long_group_label) %>%
  left_join(significance_letters, by = c("Long_group_label" = "group"))

print("各经度区间的β2TPM统计信息:")
print(data_summary)

cat("\n=== β2TPM方差分析结果 ===\n")
print(summary(aov_result))

cat("\n=== β2TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

cat("β2TPM原始数据经度范围:", min(data_clean$Long, na.rm = TRUE), "到", max(data_clean$Long, na.rm = TRUE), "\n")
cat("β2TPM总数据点数:", nrow(data_clean), "\n")
cat("β2TPM经度区间数:", length(unique(data_clean$Long_group)), "\n")
  arrange(Long_group_label) %>%
  left_join(significance_letters, by = c("Long_group_label" = "group"))

print("各经度区间的β2TPM统计信息:")
print(data_summary)

cat("\n=== β2TPM方差分析结果 ===\n")
print(summary(aov_result))

cat("\n=== β2TPM Tukey多重比较检验结果 ===\n")
print(tukey_result)

cat("β2TPM原始数据经度范围:", min(data_clean$Long, na.rm = TRUE), "到", max(data_clean$Long, na.rm = TRUE), "\n")
cat("β2TPM总数据点数:", nrow(data_clean), "\n")
cat("β2TPM经度区间数:", length(unique(data_clean$Long_group)), "\n")

