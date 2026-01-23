# ==============================================================================
# 脚本功能: 分析甲基化和SNP位点与温度的关系
# 分析流程:
# 1. 读取特征文件(甲基化CG/CHG/CHH和SNP)
# 2. 从Alt.snp_meth.tsv中提取对应位点和温度信息
# 3. 对每个位点进行正态性检验
# 4. 根据正态性选择统计检验方法(ANOVA或Kruskal-Wallis)
# 5. 筛选显著位点(P < 0.05或FDR校正后Q < 0.05)
# 6. 对显著位点作图
# ==============================================================================

# 加载必要的包
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(car)
library(pheatmap)
library(RColorBrewer)
library(viridis)

# 设置工作目录
setwd("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana")

# ==============================================================================
# 定义函数
# ==============================================================================

# 读取特征文件并返回特征名称列表
read_feature_files <- function(folder_path, pattern) {
  files <- list.files(folder_path, pattern = pattern, full.names = TRUE)
  all_features <- c()
  
  for (file in files) {
    data <- read.csv(file, stringsAsFactors = FALSE)
    if ("Feature" %in% colnames(data)) {
      all_features <- c(all_features, data$Feature)
    }
  }
  
  return(unique(all_features))
}

# 对单个位点进行统计检验
test_site_temperature <- function(site_name, site_data, temperature_data) {
  # 合并数据
  df <- data.frame(
    Temperature = temperature_data,
    Value = site_data,
    stringsAsFactors = FALSE
  )
  
  # 去除NA值
  df <- df[!is.na(df$Value), ]
  
  # 确保至少有两个温度组
  if (length(unique(df$Temperature)) < 2) {
    return(list(
      site = site_name,
      test_method = "insufficient_data",
      p_value = NA,
      is_normal = NA
    ))
  }
  
  # 检查每个温度组的样本量
  group_counts <- table(df$Temperature)
  if (any(group_counts < 3)) {
    return(list(
      site = site_name,
      test_method = "insufficient_samples",
      p_value = NA,
      is_normal = NA
    ))
  }
  
  # 正态性检验 - 检查每个温度组的正态性
  is_normal <- TRUE
  temp_groups <- unique(df$Temperature)
  
  for (temp in temp_groups) {
    group_data <- df$Value[df$Temperature == temp]
    if (length(group_data) >= 3) {
      shapiro_result <- tryCatch(
        shapiro.test(group_data),
        error = function(e) list(p.value = 0)
      )
      if (shapiro_result$p.value < 0.05) {
        is_normal <- FALSE
        break
      }
    }
  }
  
  # 根据正态性选择检验方法
  if (is_normal) {
    # 使用ANOVA
    test_result <- tryCatch(
      {
        aov_result <- aov(Value ~ Temperature, data = df)
        summary_result <- summary(aov_result)
        list(
          site = site_name,
          test_method = "ANOVA",
          p_value = summary_result[[1]][["Pr(>F)"]][1],
          is_normal = TRUE
        )
      },
      error = function(e) {
        list(
          site = site_name,
          test_method = "ANOVA_error",
          p_value = NA,
          is_normal = TRUE
        )
      }
    )
  } else {
    # 使用Kruskal-Wallis检验
    test_result <- tryCatch(
      {
        kw_result <- kruskal.test(Value ~ Temperature, data = df)
        list(
          site = site_name,
          test_method = "Kruskal-Wallis",
          p_value = kw_result$p.value,
          is_normal = FALSE
        )
      },
      error = function(e) {
        list(
          site = site_name,
          test_method = "KW_error",
          p_value = NA,
          is_normal = FALSE
        )
      }
    )
  }
  
  return(test_result)
}

# 创建汇总表格，包含位点信息、P值、各温度均值和趋势
create_summary_table <- function(sites_data, alt_data, test_results_df) {
  summary_list <- list()
  
  for (i in 1:nrow(test_results_df)) {
    site <- as.character(test_results_df$site[i])
    p_value <- test_results_df$p_value[i]
    
    # 提取特征类型和数字
    feature_type <- sub("_.*", "", site)  # 提取前缀（CG, CHG, CHH, snp）
    feature_num <- as.numeric(sub(".*_(\\d+)$", "\\1", site))  # 提取数字
    
    # 获取位点数据
    site_values <- alt_data[[site]]
    temp_values <- alt_data$tem
    
    # 准备数据框
    df <- data.frame(
      Temperature = temp_values,
      Value = site_values,
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$Value), ]
    
    # 计算各温度的均值
    temp_means <- df %>%
      group_by(Temperature) %>%
      summarise(Mean = mean(Value, na.rm = TRUE), .groups = 'drop')
    
    # 获取10, 16, 22度的均值
    mean_10 <- ifelse(10 %in% temp_means$Temperature, 
                      temp_means$Mean[temp_means$Temperature == 10], NA)
    mean_16 <- ifelse(16 %in% temp_means$Temperature, 
                      temp_means$Mean[temp_means$Temperature == 16], NA)
    mean_22 <- ifelse(22 %in% temp_means$Temperature, 
                      temp_means$Mean[temp_means$Temperature == 22], NA)
    
    # 判断趋势
    # 如果有三个温度的数据，判断整体趋势
    if (!is.na(mean_10) && !is.na(mean_16) && !is.na(mean_22)) {
      if (mean_10 < mean_16 && mean_16 < mean_22) {
        trend <- "Up"
      } else if (mean_10 > mean_16 && mean_16 > mean_22) {
        trend <- "Down"
      } else if (mean_16 > mean_10 && mean_16 > mean_22) {
        trend <- "Peak_at_16C"
      } else if (mean_16 < mean_10 && mean_16 < mean_22) {
        trend <- "Valley_at_16C"
      } else {
        trend <- "Complex"
      }
    } else if (!is.na(mean_10) && !is.na(mean_16)) {
      if (mean_16 > mean_10) {
        trend <- "Up"
      } else if (mean_16 < mean_10) {
        trend <- "Down"
      } else {
        trend <- "Stable"
      }
    } else if (!is.na(mean_16) && !is.na(mean_22)) {
      if (mean_22 > mean_16) {
        trend <- "Up"
      } else if (mean_22 < mean_16) {
        trend <- "Down"
      } else {
        trend <- "Stable"
      }
    } else {
      trend <- "Unknown"
    }
    
    # 添加到列表
    summary_list[[i]] <- data.frame(
      Site_ID = site,
      Feature_Type = feature_type,
      Feature_Num = ifelse(is.na(feature_num), 0, feature_num),
      P_value = p_value,
      Mean_10C = ifelse(!is.na(mean_10), mean_10, "NA"),
      Mean_16C = ifelse(!is.na(mean_16), mean_16, "NA"),
      Mean_22C = ifelse(!is.na(mean_22), mean_22, "NA"),
      Trend = trend,
      stringsAsFactors = FALSE
    )
  }
  
  # 合并所有结果
  summary_df <- do.call(rbind, summary_list)
  
  # 按特征类型（SNP、CG、CHG、CHH）和数字排序
  summary_df$Feature_Type <- factor(summary_df$Feature_Type, 
                                    levels = c("snp", "CG", "CHG", "CHH"))
  summary_df <- summary_df %>%
    arrange(Feature_Type, Feature_Num)
  
  # 移除辅助列
  summary_df <- summary_df %>%
    select(-Feature_Type, -Feature_Num)
  
  return(summary_df)
}

# 为显著位点绘制热图
plot_heatmap <- function(sites_data, alt_data, output_file, 
                         site_names = NULL, max_sites = 50, feature_label = "") {
  # 如果没有指定位点名称，使用所有位点
  if (is.null(site_names)) {
    site_names <- sites_data$site
  }
  
  # 限制位点数量
  if (length(site_names) > max_sites) {
    cat("位点数量过多(", length(site_names), "), 仅显示前", max_sites, "个最显著的位点\n")
    site_names <- site_names[1:max_sites]
  }
  
  # 准备热图矩阵
  heatmap_matrix <- matrix(NA, nrow = length(site_names), ncol = 0)
  rownames(heatmap_matrix) <- site_names
  
  # 获取所有温度
  all_temps <- sort(unique(alt_data$tem))
  
  # 对每个温度组，计算每个位点的均值
  for (temp in all_temps) {
    temp_indices <- which(alt_data$tem == temp)
    temp_means <- sapply(site_names, function(site) {
      values <- alt_data[[site]][temp_indices]
      mean(values, na.rm = TRUE)
    })
    heatmap_matrix <- cbind(heatmap_matrix, temp_means)
  }
  
  # 设置列名
  colnames(heatmap_matrix) <- paste0(all_temps, "°C")
  
  # 标准化数据（按行进行Z-score标准化）
  heatmap_matrix_scaled <- t(apply(heatmap_matrix, 1, function(x) {
    if (all(is.na(x))) return(x)
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
  
  # 创建热图
  tryCatch({
    pdf(output_file, width = 8, height = max(6, length(site_names) * 0.15))
    
    pheatmap(
      heatmap_matrix_scaled,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      color = viridis(100, option = "D"),  # 使用viridis配色
      main = paste0("Temperature Effect on ", feature_label, " (Z-score)"),
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 10,
      border_color = "grey60",
      na_col = "grey90",
      angle_col = 0  # 列名水平显示
    )
    
    dev.off()
    cat("热图已保存至:", output_file, "\n")
  }, error = function(e) {
    cat("绘制热图失败:", e$message, "\n")
  })
  
  return(heatmap_matrix)
}

# 为多种特征类型绘制综合热图
plot_combined_heatmap <- function(all_sig_sites, alt_data, output_file, max_sites_per_type = 50) {
  # 收集所有位点数据
  all_site_names <- c()
  feature_annotations <- c()
  
  for (feature_name in names(all_sig_sites)) {
    sites <- all_sig_sites[[feature_name]]
    if (!is.null(sites) && nrow(sites) > 0) {
      # 按P值排序并限制数量
      sites <- sites[order(sites$p_value), ]
      n_sites <- min(max_sites_per_type, nrow(sites))
      selected_sites <- as.character(sites$site[1:n_sites])
      
      all_site_names <- c(all_site_names, selected_sites)
      feature_annotations <- c(feature_annotations, rep(feature_name, n_sites))
    }
  }
  
  if (length(all_site_names) == 0) {
    cat("没有可用的位点进行绘图\n")
    return(NULL)
  }
  
  cat("总共", length(all_site_names), "个位点将被绘制\n")
  
  # 准备热图矩阵
  heatmap_matrix <- matrix(NA, nrow = length(all_site_names), ncol = 0)
  rownames(heatmap_matrix) <- all_site_names
  
  # 获取所有温度
  all_temps <- sort(unique(alt_data$tem))
  
  # 对每个温度组，计算每个位点的均值
  for (temp in all_temps) {
    temp_indices <- which(alt_data$tem == temp)
    temp_means <- sapply(all_site_names, function(site) {
      values <- alt_data[[site]][temp_indices]
      mean(values, na.rm = TRUE)
    })
    heatmap_matrix <- cbind(heatmap_matrix, temp_means)
  }
  
  # 设置列名
  colnames(heatmap_matrix) <- paste0(all_temps, "°C")
  
  # 标准化数据（按行进行Z-score标准化）
  heatmap_matrix_scaled <- t(apply(heatmap_matrix, 1, function(x) {
    if (all(is.na(x))) return(x)
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
  
  # 准备注释数据框
  annotation_row <- data.frame(
    " " = factor(feature_annotations, levels = c("CG", "CHG", "CHH", "SNP")),
    check.names = FALSE
  )
  rownames(annotation_row) <- all_site_names
  
  # 设置注释颜色
  ann_colors <- list(
    " " = c(
      CG = "#E64B35",   # 红色
      CHG = "#4DBBD5",  # 蓝色
      CHH = "#00A087",  # 绿色
      SNP = "#F39B7F"   # 橙色
    )
  )
  
  # 创建热图
  tryCatch({
    pdf(output_file, width = 10, height = max(8, length(all_site_names) * 0.12))
    
    pheatmap(
      heatmap_matrix_scaled,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),  # 科研配色：红-黄-蓝
      main = "Temperature Effect on Methylation and SNP Sites",
      fontsize = 10,
      fontsize_row = 6,
      fontsize_col = 12,
      border_color = NA,
      na_col = "grey90",
      annotation_row = annotation_row,
      annotation_colors = ann_colors,
      angle_col = 0,  # 列名水平显示
      gaps_col = NULL,
      clustering_distance_rows = "euclidean",
      clustering_method = "complete"
    )
    
    dev.off()
    cat("综合热图已保存至:", output_file, "\n")
  }, error = function(e) {
    cat("绘制综合热图失败:", e$message, "\n")
  })
  
  return(heatmap_matrix)
}

# ==============================================================================
# 主分析流程
# ==============================================================================

cat("\n========================================\n")
cat("开始分析甲基化和SNP位点与温度的关系\n")
cat("========================================\n\n")

# 读取温度和甲基化/SNP数据
cat("读取Alt.snp_meth.tsv文件...\n")
alt_data <- read.table("list/RCA/Alt.snp_meth.tsv", 
                       header = TRUE, 
                       sep = "\t", 
                       stringsAsFactors = FALSE)

cat("数据维度:", nrow(alt_data), "行,", ncol(alt_data), "列\n")
cat("温度分组:", paste(unique(alt_data$tem), collapse = ", "), "\n\n")

# 定义要分析的文件夹
folders <- c(
  "list/xgboot/TPM_4.5/feature_data_extraction/Scoupled_specific",
  "list/xgboot/TPM_4.5/feature_data_extraction/Sexpr_only",
  "list/xgboot/TPM_4.5/feature_data_extraction/Ssplice_only"
)

# 定义特征类型
feature_types <- list(
  CG = list(pattern = "_CG_data.csv", y_label = "CG Methylation Level (%)"),
  CHG = list(pattern = "_CHG_data.csv", y_label = "CHG Methylation Level (%)"),
  CHH = list(pattern = "_CHH_data.csv", y_label = "CHH Methylation Level (%)"),
  SNP = list(pattern = "_snp_data.csv", y_label = "SNP Allele Frequency")
)

# 对每个文件夹和特征类型进行分析
for (folder in folders) {
  folder_name <- basename(folder)
  cat("\n========================================\n")
  cat("分析文件夹:", folder_name, "\n")
  cat("========================================\n\n")
  
  # 创建输出目录
  output_dir <- file.path("list/RCA", folder_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 存储每个特征类型的显著位点，用于后续综合绘图
  all_sig_sites <- list()
  
  for (feature_name in names(feature_types)) {
    cat("\n--- 分析特征类型:", feature_name, "---\n")
    
    feature_info <- feature_types[[feature_name]]
    
    # 读取特征文件
    features <- read_feature_files(folder, feature_info$pattern)
    
    if (length(features) == 0) {
      cat("未找到", feature_name, "特征文件,跳过...\n")
      all_sig_sites[[feature_name]] <- NULL
      next
    }
    
    cat("找到", length(features), "个特征位点\n")
    
    # 检查这些特征在alt_data中是否存在
    available_features <- intersect(features, colnames(alt_data))
    cat("在Alt.snp_meth.tsv中找到", length(available_features), "个位点\n")
    
    if (length(available_features) == 0) {
      cat("没有可用的位点进行分析,跳过...\n")
      all_sig_sites[[feature_name]] <- NULL
      next
    }
    
    # 对每个位点进行统计检验
    cat("开始统计检验...\n")
    test_results <- list()
    
    for (i in seq_along(available_features)) {
      site <- available_features[i]
      
      if (i %% 50 == 0) {
        cat("  处理进度:", i, "/", length(available_features), "\n")
      }
      
      site_values <- alt_data[[site]]
      temp_values <- alt_data$tem
      
      result <- test_site_temperature(site, site_values, temp_values)
      test_results[[i]] <- result
    }
    
    # 整理结果
    results_df <- do.call(rbind, lapply(test_results, as.data.frame))
    
    # 去除无效结果
    valid_results <- results_df[!is.na(results_df$p_value), ]
    
    cat("\n有效检验结果数:", nrow(valid_results), "\n")
    
    if (nrow(valid_results) == 0) {
      cat("没有有效的检验结果,跳过...\n")
      all_sig_sites[[feature_name]] <- NULL
      next
    }
    
    # FDR校正
    valid_results$q_value <- p.adjust(valid_results$p_value, method = "fdr")
    
    # 筛选显著位点 (P < 0.05)
    sig_sites_p <- valid_results[valid_results$p_value < 0.05, ]
    cat("P < 0.05 的显著位点数:", nrow(sig_sites_p), "\n")
    
    # 筛选显著位点 (Q < 0.05)
    sig_sites_q <- valid_results[valid_results$q_value < 0.05, ]
    cat("Q < 0.05 的显著位点数:", nrow(sig_sites_q), "\n")
    
    # 保存结果
    result_file <- file.path(output_dir, paste0(feature_name, "_test_results.tsv"))
    write.table(valid_results, result_file, 
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat("结果已保存至:", result_file, "\n")
    
    # 为显著位点创建汇总表格 (使用P < 0.05标准)
    if (nrow(sig_sites_p) > 0) {
      cat("\n创建显著位点汇总表格...\n")
      
      # 按P值排序
      sig_sites_p <- sig_sites_p[order(sig_sites_p$p_value), ]
      
      # 存储显著位点供后续综合绘图使用
      all_sig_sites[[feature_name]] <- sig_sites_p
      
      # 创建汇总表格
      summary_table <- create_summary_table(sig_sites_p, alt_data, sig_sites_p)
      
      # 保存汇总表格
      summary_file <- file.path(output_dir, paste0(feature_name, "_summary_table.tsv"))
      write.table(summary_table, summary_file, 
                  sep = "\t", row.names = FALSE, quote = FALSE)
      cat("汇总表格已保存至:", summary_file, "\n")
      
      # 输出趋势统计
      cat("\n趋势统计:\n")
      print(table(summary_table$Trend))
      
      # 绘制单独的热图（如果位点数量大于等于2）
      if (nrow(sig_sites_p) >= 2) {
        cat("\n开始绘制", feature_name, "热图...\n")
        heatmap_file <- file.path(output_dir, paste0(feature_name, "_heatmap.pdf"))
        
        plot_heatmap(
          sites_data = sig_sites_p,
          alt_data = alt_data,
          output_file = heatmap_file,
          site_names = as.character(sig_sites_p$site),
          max_sites = 50,
          feature_label = feature_name
        )
      } else {
        cat("\n位点数量不足(仅", nrow(sig_sites_p), "个)，跳过单独热图绘制\n")
      }
    } else {
      all_sig_sites[[feature_name]] <- NULL
    }
    
    # 输出统计摘要
    cat("\n统计检验方法分布:\n")
    print(table(valid_results$test_method))
    
    cat("\n正态性分布:\n")
    print(table(valid_results$is_normal))
  }
  
  # 绘制综合热图（包含所有特征类型）
  cat("\n========================================\n")
  cat("绘制综合热图（CG、CHG、CHH、SNP）\n")
  cat("========================================\n\n")
  
  combined_heatmap_file <- file.path(output_dir, "combined_all_features_heatmap.pdf")
  plot_combined_heatmap(
    all_sig_sites = all_sig_sites,
    alt_data = alt_data,
    output_file = combined_heatmap_file,
    max_sites_per_type = 30  # 每种特征类型最多显示30个位点
  )
}

cat("\n========================================\n")
cat("分析完成!\n")
cat("========================================\n")
