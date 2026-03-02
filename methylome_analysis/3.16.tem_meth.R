# ==============================================================================
# 脚本功能: 分析甲基化和SNP位点与温度的关系（按基因区域分区）
# 分析流程:
# 1. 读取特征文件(.csv格式)
# 2. 根据end列位置分为三个区域：上游、body区、下游
# 3. 从Alt.snp_meth.tsv中提取对应位点和温度信息
# 4. 对甲基化位点：Kruskal-Wallis检验 + Dunn两两比较
# 5. 对SNP位点：卡方检验 + 两两比较
# 6. 根据新分类规则进行分类
# 7. 三个区域分别绘制热图
# ==============================================================================

# 加载必要的包
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(FSA)  # for Dunn test
library(ggsignif)  # for significance markers

# 设置工作目录
setwd("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana")

# ==============================================================================
# 定义函数
# ==============================================================================

# 读取CSV文件并根据end列位置分区
read_and_classify_features <- function(folder_path) {
  files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  all_features_upstream <- list()
  all_features_body <- list()
  all_features_downstream <- list()
  
  for (file in files) {
    data <- read.csv(file, stringsAsFactors = FALSE)
    
    # 确保有Feature列
    if (!"Feature" %in% colnames(data)) {
      cat("警告: 文件", basename(file), "没有Feature列，跳过\n")
      next
    }
    
    # 对于SNP文件，使用Position列；对于甲基化文件，使用end列
    position_col <- NULL
    if ("end" %in% colnames(data)) {
      position_col <- "end"
    } else if ("Position" %in% colnames(data)) {
      position_col <- "Position"
    } else {
      cat("警告: 文件", basename(file), "没有end或Position列，跳过\n")
      next
    }
    
    # 根据位置列分区
    for (i in 1:nrow(data)) {
      pos <- data[[position_col]][i]
      feature_name <- data$Feature[i]
      
      # 分类到不同区域
      if (pos > 16573692 && pos <= 16575692) {
        # 上游
        all_features_upstream[[feature_name]] <- TRUE
      } else if (pos > 16570745 && pos <= 16573692) {
        # body区
        all_features_body[[feature_name]] <- TRUE
      } else if (pos > 16568745 && pos <= 16570745) {
        # 下游
        all_features_downstream[[feature_name]] <- TRUE
      }
    }
  }
  
  return(list(
    upstream = names(all_features_upstream),
    body = names(all_features_body),
    downstream = names(all_features_downstream)
  ))
}

# 对甲基化位点进行Kruskal-Wallis检验和Dunn两两比较
test_methylation_site <- function(site_name, site_data, temperature_data) {
  # 合并数据
  df <- data.frame(
    Temperature = as.factor(temperature_data),
    Value = site_data,
    stringsAsFactors = FALSE
  )
  
  # 去除NA值
  df <- df[!is.na(df$Value), ]
  
  # 确保至少有两个温度组
  if (length(unique(df$Temperature)) < 2) {
    return(NULL)
  }
  
  # 检查每个温度组的样本量
  group_counts <- table(df$Temperature)
  if (any(group_counts < 3)) {
    return(NULL)
  }
  
  # Kruskal-Wallis检验
  kw_result <- tryCatch(
    kruskal.test(Value ~ Temperature, data = df),
    error = function(e) NULL
  )
  
  if (is.null(kw_result)) {
    return(NULL)
  }
  
  # 计算各温度的均值
  temp_means <- df %>%
    group_by(Temperature) %>%
    summarise(Mean = mean(Value, na.rm = TRUE), .groups = 'drop')
  
  mean_10 <- ifelse(10 %in% temp_means$Temperature, 
                    temp_means$Mean[temp_means$Temperature == 10], NA)
  mean_16 <- ifelse(16 %in% temp_means$Temperature, 
                    temp_means$Mean[temp_means$Temperature == 16], NA)
  mean_22 <- ifelse(22 %in% temp_means$Temperature, 
                    temp_means$Mean[temp_means$Temperature == 22], NA)
  
  # Dunn两两比较
  dunn_result <- tryCatch(
    dunnTest(Value ~ Temperature, data = df, method = "bonferroni"),
    error = function(e) NULL
  )
  
  p_10_16 <- NA
  p_16_22 <- NA
  p_10_22 <- NA
  
  if (!is.null(dunn_result)) {
    comp_table <- dunn_result$res
    
    # 提取两两比较的P值
    for (i in 1:nrow(comp_table)) {
      comparison <- comp_table$Comparison[i]
      if (grepl("10.*16|16.*10", comparison)) {
        p_10_16 <- comp_table$P.adj[i]
      } else if (grepl("16.*22|22.*16", comparison)) {
        p_16_22 <- comp_table$P.adj[i]
      } else if (grepl("10.*22|22.*10", comparison)) {
        p_10_22 <- comp_table$P.adj[i]
      }
    }
  }
  
  # 根据规则分类
  trend <- classify_trend(p_10_16, p_16_22, mean_10, mean_16, mean_22)
  
  return(list(
    site = site_name,
    test_method = "Kruskal-Wallis + Dunn",
    overall_p_value = kw_result$p.value,
    p_10_16 = p_10_16,
    p_16_22 = p_16_22,
    p_10_22 = p_10_22,
    mean_10 = mean_10,
    mean_16 = mean_16,
    mean_22 = mean_22,
    trend = trend
  ))
}

# 对SNP位点进行卡方检验和两两比较
test_snp_site <- function(site_name, site_data, temperature_data) {
  # 合并数据
  df <- data.frame(
    Temperature = as.factor(temperature_data),
    Value = site_data,
    stringsAsFactors = FALSE
  )
  
  # 去除NA值
  df <- df[!is.na(df$Value), ]
  
  # 确保至少有两个温度组
  if (length(unique(df$Temperature)) < 2) {
    return(NULL)
  }
  
  # 检查每个温度组的样本量
  group_counts <- table(df$Temperature)
  if (any(group_counts < 3)) {
    return(NULL)
  }
  
  # 创建列联表（用于检验：基于计数）
  contingency_table <- table(df$Temperature, df$Value)
  
  # Fisher精确检验（使用计数数据）
  fisher_result <- tryCatch(
    fisher.test(contingency_table),
    error = function(e) NULL
  )
  
  if (is.null(fisher_result)) {
    return(NULL)
  }
  
  # 计算每个温度下SNP=1的频率（用于趋势分类）
  prop_df <- df %>%
    group_by(Temperature) %>%
    summarise(
      Prop = sum(Value == 1, na.rm = TRUE) / n(),
      .groups = 'drop'
    )
  
  prop_10 <- ifelse(10 %in% prop_df$Temperature, 
                    prop_df$Prop[prop_df$Temperature == 10], NA)
  prop_16 <- ifelse(16 %in% prop_df$Temperature, 
                    prop_df$Prop[prop_df$Temperature == 16], NA)
  prop_22 <- ifelse(22 %in% prop_df$Temperature, 
                    prop_df$Prop[prop_df$Temperature == 22], NA)
  
  # 两两Fisher精确检验（使用计数数据）
  p_10_16 <- NA
  p_16_22 <- NA
  p_10_22 <- NA
  
  temps <- unique(df$Temperature)
  if (length(temps) >= 2) {
    # 10 vs 16
    if ("10" %in% temps && "16" %in% temps) {
      df_10_16 <- df[df$Temperature %in% c("10", "16"), ]
      ct_10_16 <- table(df_10_16$Temperature, df_10_16$Value)  # 计数数据
      fisher_10_16 <- tryCatch(
        fisher.test(ct_10_16),
        error = function(e) NULL
      )
      if (!is.null(fisher_10_16)) p_10_16 <- fisher_10_16$p.value
    }
    
    # 16 vs 22
    if ("16" %in% temps && "22" %in% temps) {
      df_16_22 <- df[df$Temperature %in% c("16", "22"), ]
      ct_16_22 <- table(df_16_22$Temperature, df_16_22$Value)  # 计数数据
      fisher_16_22 <- tryCatch(
        fisher.test(ct_16_22),
        error = function(e) NULL
      )
      if (!is.null(fisher_16_22)) p_16_22 <- fisher_16_22$p.value
    }
    
    # 10 vs 22
    if ("10" %in% temps && "22" %in% temps) {
      df_10_22 <- df[df$Temperature %in% c("10", "22"), ]
      ct_10_22 <- table(df_10_22$Temperature, df_10_22$Value)  # 计数数据
      fisher_10_22 <- tryCatch(
        fisher.test(ct_10_22),
        error = function(e) NULL
      )
      if (!is.null(fisher_10_22)) p_10_22 <- fisher_10_22$p.value
    }
  }
  
  # 根据规则分类（使用频率数据）
  trend <- classify_trend(p_10_16, p_16_22, prop_10, prop_16, prop_22)
  
  return(list(
    site = site_name,
    test_method = "Fisher's Exact Test",
    overall_p_value = fisher_result$p.value,
    p_10_16 = p_10_16,
    p_16_22 = p_16_22,
    p_10_22 = p_10_22,
    mean_10 = prop_10,
    mean_16 = prop_16,
    mean_22 = prop_22,
    trend = trend
  ))
}

# 趋势分类函数
classify_trend <- function(p_10_16, p_16_22, mean_10, mean_16, mean_22) {
  # 检查是否有足够的数据
  # if (is.na(p_10_16) || is.na(p_16_22) || is.na(mean_10) || is.na(mean_16) || is.na(mean_22)) {
  #   return("Insufficient_Data")
  # }
  
  # 应用分类规则
  if (p_10_16 < 0.05 && p_16_22 < 0.05 && mean_10 < mean_16 && mean_16 < mean_22) {
    return("Up")
  } else if (p_10_16 < 0.05 && p_16_22 < 0.05 && mean_10 > mean_16 && mean_16 > mean_22) {
    return("Down")
  } else if (p_10_16 < 0.05 && p_16_22 < 0.05 && mean_10 < mean_16 && mean_16 > mean_22) {
    return("Peak_at_16C")
  } else if (p_10_16 < 0.05 && p_16_22 < 0.05 && mean_10 > mean_16 && mean_16 < mean_22) {
    return("Valley_at_16C")
  } else if ((p_10_16 > 0.05 || is.na(p_10_16)) && p_16_22 < 0.05 && mean_16 < mean_22) {
    return("Up-like_HighTemp")
  } else if (p_10_16 < 0.05 && p_16_22 > 0.05 && mean_10 < mean_16) {
    return("Up-like_LowMid")
  } else if ((p_10_16 > 0.05 || is.na(p_10_16)) && p_16_22 < 0.05 && mean_16 > mean_22) {
    return("Down-like_HighTemp")
  } else if (p_10_16 < 0.05 && p_16_22 > 0.05 && mean_10 > mean_16) {
    return("Down-like_LowTemp")
  } else {
    return("Other")
  }
}

# 为显著位点绘制热图（支持区域标注）
plot_heatmap <- function(alt_data, site_names, output_file, region_label = "", trend_label = "", region_annotation = NULL) {
  if (length(site_names) == 0) {
    cat("没有位点需要绘制\n")
    return(NULL)
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
      if (site %in% colnames(alt_data)) {
        values <- alt_data[[site]][temp_indices]
        mean(values, na.rm = TRUE)
      } else {
        NA
      }
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
  
  # 准备注释：区分CG、CHG、CHH、SNP
  feature_types <- sapply(site_names, function(site) {
    if (grepl("^CG_", site)) {
      return("CG")
    } else if (grepl("^CHG_", site)) {
      return("CHG")
    } else if (grepl("^CHH_", site)) {
      return("CHH")
    } else if (grepl("^snp", site)) {
      return("SNP")
    } else {
      return("Other")
    }
  })
  
  # 如果提供了区域注释，添加区域列
  if (!is.null(region_annotation)) {
    annotation_row <- data.frame(
      Region = factor(region_annotation, levels = c("Upstream", "Body", "Downstream")),
      Type = factor(feature_types, levels = c("CG", "CHG", "CHH", "SNP", "Other")),
      row.names = site_names,
      check.names = FALSE
    )
    
    ann_colors <- list(
      Region = c(
        Upstream = "#FF6B6B",
        Body = "#4ECDC4",
        Downstream = "#95E1D3"
      ),
      Type = c(
        CG = "#E64B35",
        CHG = "#4DBBD5",
        CHH = "#00A087",
        SNP = "#F39B7F",
        Other = "#8491B4"
      )
    )
  } else {
    annotation_row <- data.frame(
      " " = factor(feature_types, levels = c("CG", "CHG", "CHH", "SNP", "Other")),
      row.names = site_names,
      check.names = FALSE
    )
    
    ann_colors <- list(
      " " = c(
        CG = "#E64B35",
        CHG = "#4DBBD5",
        CHH = "#00A087",
        SNP = "#F39B7F",
        Other = "#8491B4"
      )
    )
  }
  
  # 创建热图
  tryCatch({
    pdf(output_file, width = 10, height = max(8, length(site_names) * 0.12))
    
    main_title <- paste0(region_label, " - ", trend_label, " (Z-score)")
    
    pheatmap(
      heatmap_matrix_scaled,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
      main = main_title,
      fontsize = 10,
      fontsize_row = 6,
      fontsize_col = 12,
      border_color = NA,
      na_col = "grey90",
      annotation_row = annotation_row,
      annotation_colors = ann_colors,
      annotation_names_row = FALSE,
      angle_col = 0,
      clustering_distance_rows = "euclidean",
      clustering_method = "complete"
    )
    
    dev.off()
    cat("热图已保存至:", output_file, "\n")
  }, error = function(e) {
    cat("绘制热图失败:", e$message, "\n")
  })
  
  return(heatmap_matrix)
}

# 为单个位点绘制箱线图
plot_boxplot <- function(alt_data, site_name, output_file, region_label = "", trend_label = "", 
                         p_10_16 = NA, p_16_22 = NA, p_10_22 = NA) {
  if (length(site_name) == 0) {
    cat("没有位点需要绘制\n")
    return(NULL)
  }
  
  # 提取位点数据
  if (!site_name %in% colnames(alt_data)) {
    cat("位点", site_name, "在数据中不存在\n")
    return(NULL)
  }
  
  # 准备数据
  plot_data <- data.frame(
    Temperature = factor(alt_data$tem, levels = c(10, 16, 22)),
    Value = alt_data[[site_name]],
    stringsAsFactors = FALSE
  )
  
  # 去除NA值
  plot_data <- plot_data[!is.na(plot_data$Value), ]
  
  if (nrow(plot_data) == 0) {
    cat("位点", site_name, "没有有效数据\n")
    return(NULL)
  }
  
  # 确定特征类型
  feature_type <- "Other"
  if (grepl("^CG_", site_name)) {
    feature_type <- "CG"
  } else if (grepl("^CHG_", site_name)) {
    feature_type <- "CHG"
  } else if (grepl("^CHH_", site_name)) {
    feature_type <- "CHH"
  } else if (grepl("^snp", site_name)) {
    feature_type <- "SNP"
  }
  
  # 准备显著性标记
  get_significance <- function(p_val) {
    if (is.na(p_val)) {
      return("ns")
    } else if (p_val < 0.001) {
      return("***")
    } else if (p_val < 0.01) {
      return("**")
    } else if (p_val < 0.05) {
      return("*")
    } else {
      return("ns")
    }
  }
  
  sig_10_16 <- get_significance(p_10_16)
  sig_16_22 <- get_significance(p_16_22)
  sig_10_22 <- get_significance(p_10_22)
  
  # 创建箱线图
  tryCatch({
    pdf(output_file, width = 8, height = 6)
    
    y_label <- ifelse(feature_type == "SNP", "SNP Allele Frequency(%)", "Methylation Level(%)")
    
    p <- ggplot(plot_data, aes(x = Temperature, y = Value, fill = Temperature)) +
      # 添加箱须末端标记
      stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(width = 0.5)) +
      # 标准的箱线图
      geom_boxplot(alpha = 1, outlier.shape = 16, outlier.size = 1, 
                   position = position_dodge(width = 0.5), width = 0.5) +
      scale_fill_manual(values = c("10" = "#4472C4", "16" = "#70AD47", "22" = "#FFC000")) +
      # 添加显著性标记
      geom_signif(comparisons = list(c("10", "16"), c("16", "22"), c("10", "22")),
                  annotations = c(sig_10_16, sig_16_22, sig_10_22),
                  y_position = c(max(plot_data$Value, na.rm = TRUE) * 1.1, 
                               max(plot_data$Value, na.rm = TRUE) * 1.2,
                               max(plot_data$Value, na.rm = TRUE) * 1.3),
                  tip_length = 0.02,
                  textsize = 5) +
      labs(
        title = paste0(region_label, " - ", trend_label, "\n", site_name, " (", feature_type, ")"),
        y = y_label,
        fill = "Tem(°C)"
      ) +
      theme_minimal() +
      theme(
        # 图标题
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
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
      scale_y_continuous(limits = c(0, max(plot_data$Value, na.rm = TRUE) * 1.4))
    
    print(p)
    
    dev.off()
    cat("箱线图已保存至:", output_file, "\n")
  }, error = function(e) {
    cat("绘制箱线图失败:", e$message, "\n")
  })
  
  return(invisible(NULL))
}

# ==============================================================================
# 主分析流程
# ==============================================================================

cat("\n========================================\n")
cat("开始分析甲基化和SNP位点与温度的关系（按基因区域分区）\n")
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
  "list/xgboot/TPM_4.5/feature_data_extraction/A_per_TPM",
  "list/xgboot/TPM_4.5/feature_data_extraction/A_TPM",
  "list/xgboot/TPM_4.5/feature_data_extraction/B_per_TPM",
  "list/xgboot/TPM_4.5/feature_data_extraction/B_TPM",
  "list/xgboot/TPM_4.5/feature_data_extraction/total_TPM"
)

# 对每个文件夹进行分析
for (folder in folders) {
  folder_name <- basename(folder)
  cat("\n========================================\n")
  cat("分析文件夹:", folder_name, "\n")
  cat("========================================\n\n")
  
  # 创建输出目录
  output_dir <- file.path("list/RCA", folder_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 读取并分类特征
  cat("读取CSV文件并根据end列位置分区...\n")
  classified_features <- read_and_classify_features(folder)
  
  cat("上游区域特征数:", length(classified_features$upstream), "\n")
  cat("Body区域特征数:", length(classified_features$body), "\n")
  cat("下游区域特征数:", length(classified_features$downstream), "\n\n")
  
  # 对三个区域分别进行分析
  regions <- list(
    upstream = list(features = classified_features$upstream, label = "Upstream"),
    body = list(features = classified_features$body, label = "Body"),
    downstream = list(features = classified_features$downstream, label = "Downstream")
  )
  
  # 存储所有区域的结果，用于汇总
  all_regions_results <- list()
  
  for (region_name in names(regions)) {
    region_info <- regions[[region_name]]
    features <- region_info$features
    region_label <- region_info$label
    
    cat("\n--- 分析区域:", region_label, "---\n")
    
    if (length(features) == 0) {
      cat("该区域没有特征，跳过...\n")
      next
    }
    
    # 检查特征在alt_data中是否存在
    available_features <- intersect(features, colnames(alt_data))
    cat("在Alt.snp_meth.tsv中找到", length(available_features), "个位点\n")
    
    if (length(available_features) == 0) {
      cat("没有可用的位点进行分析，跳过...\n")
      next
    }
    
    # 分离甲基化位点和SNP位点
    meth_features <- grep("^C", available_features, value = TRUE)
    snp_features <- grep("^snp", available_features, value = TRUE)
    
    cat("甲基化位点数:", length(meth_features), "\n")
    cat("SNP位点数:", length(snp_features), "\n")
    
    # 存储所有测试结果
    all_test_results <- list()
    
    # 对甲基化位点进行检验
    if (length(meth_features) > 0) {
      cat("\n处理甲基化位点...\n")
      
      for (i in seq_along(meth_features)) {
        site <- meth_features[i]
        
        if (i %% 50 == 0) {
          cat("  进度:", i, "/", length(meth_features), "\n")
        }
        
        site_values <- alt_data[[site]]
        temp_values <- alt_data$tem
        
        result <- test_methylation_site(site, site_values, temp_values)
        
        if (!is.null(result)) {
          all_test_results[[length(all_test_results) + 1]] <- result
        }
      }
    }
    
    # 对SNP位点进行检验
    if (length(snp_features) > 0) {
      cat("\n处理SNP位点...\n")
      
      for (i in seq_along(snp_features)) {
        site <- snp_features[i]
        
        if (i %% 50 == 0) {
          cat("  进度:", i, "/", length(snp_features), "\n")
        }
        
        site_values <- alt_data[[site]]
        temp_values <- alt_data$tem
        
        result <- test_snp_site(site, site_values, temp_values)
        
        if (!is.null(result)) {
          all_test_results[[length(all_test_results) + 1]] <- result
        }
      }
    }
    
    # 整理结果
    if (length(all_test_results) == 0) {
      cat("没有有效的检验结果，跳过...\n")
      next
    }
    
    results_df <- do.call(rbind, lapply(all_test_results, as.data.frame))
    
    cat("\n有效检验结果数:", nrow(results_df), "\n")
    
    # 保存完整结果
    result_file <- file.path(output_dir, paste0(region_name, "_test_results.tsv"))
    write.table(results_df, result_file, 
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat("结果已保存至:", result_file, "\n")
    
    # 筛选有显著趋势的位点（排除"Other"和"Insufficient_Data"）
    sig_results <- results_df[!results_df$trend %in% c("Other", "Insufficient_Data"), ]
    
    cat("有明确趋势的位点数:", nrow(sig_results), "\n")
    
    if (nrow(sig_results) == 0) {
      cat("没有显著位点，跳过绘图...\n")
      next
    }
    
    # 输出趋势统计
    cat("\n趋势统计:\n")
    print(table(sig_results$trend))
    
    # 保存显著结果摘要
    summary_file <- file.path(output_dir, paste0(region_name, "_significant_sites.tsv"))
    write.table(sig_results, summary_file, 
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat("显著位点已保存至:", summary_file, "\n")
    
    # 为每种趋势分别绘制单个区域热图
    trends <- unique(sig_results$trend)
    
    for (trend in trends) {
      trend_sites <- sig_results[sig_results$trend == trend, ]
      
      if (nrow(trend_sites) >= 2) {
        cat("\n开始绘制", region_label, "-", trend, "热图...\n")
        cat("  该趋势包含", nrow(trend_sites), "个位点\n")
        
        heatmap_file <- file.path(output_dir, paste0(region_name, "_", trend, "_heatmap.pdf"))
        
        plot_heatmap(
          alt_data = alt_data,
          site_names = as.character(trend_sites$site),
          output_file = heatmap_file,
          region_label = region_label,
          trend_label = trend
        )
      } else if (nrow(trend_sites) == 1) {
        cat("\n", trend, "趋势仅有1个位点，绘制箱线图\n")
        
        boxplot_file <- file.path(output_dir, paste0(region_name, "_", trend, "_boxplot.pdf"))
        
        # 获取该位点的p值
        site_result <- trend_sites[1, ]
        
        plot_boxplot(
          alt_data = alt_data,
          site_name = as.character(trend_sites$site[1]),
          output_file = boxplot_file,
          region_label = region_label,
          trend_label = trend,
          p_10_16 = site_result$p_10_16,
          p_16_22 = site_result$p_16_22,
          p_10_22 = site_result$p_10_22
        )
      } else {
        cat("\n", trend, "趋势位点数量为0，跳过绘图\n")
      }
    }
    
    # 存储该区域结果用于汇总
    all_regions_results[[region_name]] <- sig_results
  }
  
  # 绘制汇总热图（所有区域合并）
  cat("\n\n========================================\n")
  cat("绘制汇总热图（合并所有区域）\n")
  cat("========================================\n")
  
  if (length(all_regions_results) > 0) {
    # 收集所有趋势类型
    all_trends <- unique(unlist(lapply(all_regions_results, function(x) unique(x$trend))))
    
    for (trend in all_trends) {
      # 收集该趋势在所有区域的位点
      combined_sites <- c()
      combined_regions <- c()
      
      for (region_name in names(all_regions_results)) {
        region_results <- all_regions_results[[region_name]]
        trend_sites <- region_results[region_results$trend == trend, ]
        
        if (nrow(trend_sites) > 0) {
          combined_sites <- c(combined_sites, as.character(trend_sites$site))
          combined_regions <- c(combined_regions, rep(regions[[region_name]]$label, nrow(trend_sites)))
        }
      }
      
      if (length(combined_sites) >= 2) {
        cat("\n绘制", trend, "趋势汇总热图（包含所有区域）...\n")
        cat("  总共包含", length(combined_sites), "个位点\n")
        cat("  来自", length(unique(combined_regions)), "个区域\n")
        
        heatmap_file <- file.path(output_dir, paste0("combined_all_regions_", trend, "_heatmap.pdf"))
        
        plot_heatmap(
          alt_data = alt_data,
          site_names = combined_sites,
          output_file = heatmap_file,
          region_label = "All Regions",
          trend_label = trend,
          region_annotation = combined_regions
        )
      } else if (length(combined_sites) == 1) {
        cat("\n", trend, "趋势总位点数量为1个，绘制箱线图\n")
        
        boxplot_file <- file.path(output_dir, paste0("combined_all_regions_", trend, "_boxplot.pdf"))
        
        # 获取该位点的p值（从原始结果中查找）
        p_vals <- c(NA, NA, NA)
        for (region_name in names(all_regions_results)) {
          region_results <- all_regions_results[[region_name]]
          site_result <- region_results[region_results$site == combined_sites[1], ]
          if (nrow(site_result) > 0) {
            p_vals <- c(site_result$p_10_16, site_result$p_16_22, site_result$p_10_22)
            break
          }
        }
        
        plot_boxplot(
          alt_data = alt_data,
          site_name = combined_sites[1],
          output_file = boxplot_file,
          region_label = paste0("All Regions (", combined_regions[1], ")"),
          trend_label = trend,
          p_10_16 = p_vals[1],
          p_16_22 = p_vals[2],
          p_10_22 = p_vals[3]
        )
      } else {
        cat("\n", trend, "趋势总位点数量为0，跳过绘图\n")
      }
    }
  } else {
    cat("没有任何区域有显著结果，跳过汇总热图绘制\n")
  }
}

cat("\n========================================\n")
cat("分析完成!\n")
cat("========================================\n")
