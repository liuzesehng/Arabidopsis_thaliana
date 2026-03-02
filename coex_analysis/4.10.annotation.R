# 基因注释脚本
# 使用 org.At.tair.db 包为拟南芥基因添加功能注释
# ========================================

# 加载必要的包
library(org.At.tair.db)
library(GO.db)
library(dplyr)

# 定义批量注释函数（不包含家族信息）
annotate_genes <- function(gene_ids) {
  # 去重和去除NA值
  unique_genes <- unique(gene_ids[!is.na(gene_ids)])
  
  # 过滤掉空字符串和无效字符
  unique_genes <- unique_genes[unique_genes != "" & !is.na(unique_genes)]
  unique_genes <- trimws(unique_genes)  # 去除前后空格
  
  # 如果没有有效基因ID，返回空数据框
  if (length(unique_genes) == 0) {
    return(data.frame(
      Gene_ID = character(0),
      Symbol = character(0),
      Gene_Name = character(0),
      GO_Term = character(0),
      GO_Name = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # 获取基因符号（Symbol）
  symbols <- suppressMessages(mapIds(org.At.tair.db, 
                   keys = unique_genes,
                   column = "SYMBOL",
                   keytype = "TAIR",
                   multiVals = "first"))
  
  # 获取基因全称（Gene Name）
  genenames <- suppressMessages(mapIds(org.At.tair.db,
                     keys = unique_genes,
                     column = "GENENAME",
                     keytype = "TAIR",
                     multiVals = "first"))
  
  # 获取GO功能描述
  go_terms <- suppressMessages(mapIds(org.At.tair.db,
                    keys = unique_genes,
                    column = "GO",
                    keytype = "TAIR",
                    multiVals = "first"))
  
  # 获取GO term名称
  go_names <- sapply(go_terms, function(go_id) {
    if (is.na(go_id)) {
      return(NA)
    } else {
      tryCatch({
        Term(go_id)
      }, error = function(e) {
        return(NA)
      })
    }
  })
  
  # 创建注释数据框
  annotation_df <- data.frame(
    Gene_ID = unique_genes,
    Symbol = symbols,
    Gene_Name = genenames,
    GO_Term = go_terms,
    GO_Name = go_names,
    stringsAsFactors = FALSE
  )
  
  return(annotation_df)
}

# 定义从 TF 列表中获取家族信息的函数
annotate_genes_with_tf_family <- function(gene_ids, tf_list_file) {
  # 先调用基本注释函数
  annotation_df <- annotate_genes(gene_ids)
  
  # 读取 TF 家族列表
  if (file.exists(tf_list_file)) {
    tf_list <- read.table(tf_list_file, header = TRUE, sep = "\t", 
                         stringsAsFactors = FALSE)
    
    # 对每个基因ID匹配其家族
    families <- sapply(annotation_df$Gene_ID, function(gene_id) {
      matches <- tf_list[tf_list$Gene_ID == gene_id, "Family"]
      if (length(matches) > 0) {
        return(matches[1])
      } else {
        return(NA)
      }
    })
    
    # 在 Symbol 后添加 Family 列
    annotation_df <- data.frame(
      Gene_ID = annotation_df$Gene_ID,
      Symbol = annotation_df$Symbol,
      Family = families,
      Gene_Name = annotation_df$Gene_Name,
      GO_Term = annotation_df$GO_Term,
      GO_Name = annotation_df$GO_Name,
      stringsAsFactors = FALSE
    )
  } else {
    warning(sprintf("TF 列表文件不存在: %s", tf_list_file))
  }
  
  return(annotation_df)
}

# ========================================
# 处理 MCC.csv即gene_id.txt 文件
# ========================================
cat("正在处理 MCC.csvgene_id.txt 文件...\n")

# 定义需要处理的温度条件
temperatures <- c("10C", "16C", "22C", "total")

for (temp in temperatures) {
  mcc_file <- file.path(temp, paste0(temp, ".gene_id.txt"))
  
  # 检查文件是否存在
  if (file.exists(mcc_file)) {
    cat(sprintf("  处理文件: %s\n", mcc_file))
    
    # 读取文件
    mcc_data <- read.table(mcc_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # 如果超过20行则只保留前20行，否则全部保留
    if (nrow(mcc_data) > 20) {
      mcc_data <- mcc_data[1:20, ]
      cat(sprintf("    提取前20个基因进行注释\n"))
    } else {
      cat(sprintf("    共 %d 个基因，全部保留进行注释\n", nrow(mcc_data)))
    }
    
    # 检查第一列是否存在
    if (ncol(mcc_data) > 0) {
      # 提取第一列的基因ID
      gene_ids <- mcc_data[, 1]
      
      # 过滤和清理基因ID
      gene_ids <- as.character(gene_ids)
      gene_ids <- trimws(gene_ids)
      gene_ids <- gene_ids[gene_ids != "" & !is.na(gene_ids)]
      
      if (length(gene_ids) == 0) {
        cat(sprintf("    警告: 没有有效的基因ID\n"))
        next
      }
      
      cat(sprintf("    有效基因ID数量: %d\n", length(gene_ids)))
      
      # 进行注释
      tryCatch({
        annotations <- annotate_genes(gene_ids)
        
        # 将注释结果与原始数据合并
        # 通过第一列的基因ID进行匹配
        colnames(mcc_data)[1] <- "Gene_ID"  # 重命名第一列以便合并
        
        # 清理mcc_data中的基因ID
        mcc_data$Gene_ID <- trimws(as.character(mcc_data$Gene_ID))
        
        annotated_data <- merge(mcc_data, annotations, 
                               by = "Gene_ID", 
                               all.x = TRUE, 
                               sort = FALSE)
        
        # 调整列的顺序：Gene_ID, Symbol, Gene_Name, GO_Term, GO_Name
        final_cols <- c("Gene_ID", "Symbol", "Gene_Name", "GO_Term", "GO_Name")
        annotated_data <- annotated_data[, final_cols]
        
        # 保存注释后的文件
        output_file <- file.path(temp, paste0(temp, ".MCC.annotated.txt"))
        write.table(annotated_data, output_file, row.names = FALSE, sep = "\t", quote = FALSE)
        
        cat(sprintf("    -> 保存到: %s\n", output_file))
        cat(sprintf("    -> 注释了 %d 个基因\n", nrow(annotated_data)))
      }, error = function(e) {
        cat(sprintf("    错误: %s\n", e$message))
        cat(sprintf("    跳过该文件\n"))
      })
    } else {
      cat(sprintf("    警告: 文件为空或没有列\n"))
    }
  } else {
    cat(sprintf("  文件不存在，跳过: %s\n", mcc_file))
  }
}

cat("\n")

temperatures <- c("10C", "16C", "22C", "total")

for (temp in temperatures) {
  mcc_file <- file.path(temp, paste0(temp, ".negative.gene_id.txt"))
  
  # 检查文件是否存在
  if (file.exists(mcc_file)) {
    cat(sprintf("  处理文件: %s\n", mcc_file))
    
    # 读取文件
    mcc_data <- read.table(mcc_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # 如果超过20行则只保留前20行，否则全部保留
    if (nrow(mcc_data) > 20) {
      mcc_data <- mcc_data[1:20, ]
      cat(sprintf("    提取前20个基因进行注释\n"))
    } else {
      cat(sprintf("    共 %d 个基因，全部保留进行注释\n", nrow(mcc_data)))
    }
    
    # 检查第一列是否存在
    if (ncol(mcc_data) > 0) {
      # 提取第一列的基因ID
      gene_ids <- mcc_data[, 1]
      
      # 过滤和清理基因ID
      gene_ids <- as.character(gene_ids)
      gene_ids <- trimws(gene_ids)
      gene_ids <- gene_ids[gene_ids != "" & !is.na(gene_ids)]
      
      if (length(gene_ids) == 0) {
        cat(sprintf("    警告: 没有有效的基因ID\n"))
        next
      }
      
      cat(sprintf("    有效基因ID数量: %d\n", length(gene_ids)))
      
      # 进行注释
      tryCatch({
        annotations <- annotate_genes(gene_ids)
        
        # 将注释结果与原始数据合并
        # 通过第一列的基因ID进行匹配
        colnames(mcc_data)[1] <- "Gene_ID"  # 重命名第一列以便合并
        
        # 清理mcc_data中的基因ID
        mcc_data$Gene_ID <- trimws(as.character(mcc_data$Gene_ID))
        
        annotated_data <- merge(mcc_data, annotations, 
                               by = "Gene_ID", 
                               all.x = TRUE, 
                               sort = FALSE)
        
        # 调整列的顺序：Gene_ID, Symbol, Gene_Name, GO_Term, GO_Name
        final_cols <- c("Gene_ID", "Symbol", "Gene_Name", "GO_Term", "GO_Name")
        annotated_data <- annotated_data[, final_cols]
        
        # 保存注释后的文件
        output_file <- file.path(temp, paste0(temp, ".negative.MCC.annotated.txt"))
        write.table(annotated_data, output_file, row.names = FALSE, sep = "\t", quote = FALSE)
        
        cat(sprintf("    -> 保存到: %s\n", output_file))
        cat(sprintf("    -> 注释了 %d 个基因\n", nrow(annotated_data)))
      }, error = function(e) {
        cat(sprintf("    错误: %s\n", e$message))
        cat(sprintf("    跳过该文件\n"))
      })
    } else {
      cat(sprintf("    警告: 文件为空或没有列\n"))
    }
  } else {
    cat(sprintf("  文件不存在，跳过: %s\n", mcc_file))
  }
}

cat("\n")

# ========================================
# 处理 gene.TF_ids.txt 文件
# ========================================
cat("正在处理 gene.TF_ids.txt 文件...\n")

# TF家族列表文件路径
tf_list_file <- "Ath_TF_list.txt"

# 定义需要处理的文件
tf_files <- list(
  "10C" = file.path("10C", "gene.TF_ids.txt"),
  "10C.negative" = file.path("10C", "gene.TF_negative_ids.txt"),
  "16C" = file.path("16C", "gene.TF_ids.txt"),
  "22C.negative" = file.path("22C", "gene.TF_negative_ids.txt"),
  "total" = file.path("total", "gene.TF_ids.txt"),
  "total.negative" = file.path("total", "gene.TF_negative_ids.txt")
)

for (name in names(tf_files)) {
  tf_file <- tf_files[[name]]
  
  # 检查文件是否存在
  if (file.exists(tf_file)) {
    cat(sprintf("  处理文件: %s\n", tf_file))
    
    # 读取文件（制表符分隔）
    tf_data <- read.table(tf_file, header = TRUE, sep = "\t", 
                          stringsAsFactors = FALSE, comment.char = "",
                          quote = "", fill = TRUE)
    
    # 检查第一列是否存在
    if (ncol(tf_data) > 0) {
      # 提取第一列的基因ID
      gene_ids <- tf_data[, 1]
      
      # 过滤和清理基因ID
      gene_ids <- as.character(gene_ids)
      gene_ids <- trimws(gene_ids)
      gene_ids <- gene_ids[gene_ids != "" & !is.na(gene_ids)]
      
      if (length(gene_ids) == 0) {
        cat(sprintf("    警告: 没有有效的基因ID\n"))
        next
      }
      
      cat(sprintf("    有效基因ID数量: %d\n", length(gene_ids)))
      
      # 进行注释（包含TF家族信息）
      tryCatch({
        annotations <- annotate_genes_with_tf_family(gene_ids, tf_list_file)
        
        # 将注释结果与原始数据合并
        colnames(tf_data)[1] <- "Gene_ID"  # 重命名第一列以便合并
        
        # 清理tf_data中的基因ID
        tf_data$Gene_ID <- trimws(as.character(tf_data$Gene_ID))
        
        annotated_data <- merge(tf_data, annotations, 
                               by = "Gene_ID", 
                               all.x = TRUE, 
                               sort = FALSE)
        
        # 调整列的顺序：Gene_ID, Symbol, Family, Gene_Name, GO_Term, GO_Name
        final_cols <- c("Gene_ID", "Symbol", "Family", "Gene_Name", "GO_Term", "GO_Name")
        annotated_data <- annotated_data[, final_cols]
        
        # 保存注释后的文件
        output_file <- file.path(dirname(tf_file), 
                                paste0(name, ".gene.TF_ids.annotated.txt"))
        write.table(annotated_data, output_file, 
                   row.names = FALSE, sep = "\t", quote = FALSE)
        
        cat(sprintf("    -> 保存到: %s\n", output_file))
        cat(sprintf("    -> 注释了 %d 个基因\n", nrow(annotated_data)))
      }, error = function(e) {
        cat(sprintf("    错误: %s\n", e$message))
        cat(sprintf("    跳过该文件\n"))
      })
    } else {
      cat(sprintf("    警告: 文件为空或没有列\n"))
    }
  } else {
    cat(sprintf("  文件不存在，跳过: %s\n", tf_file))
  }
}

cat("\n")

# ========================================
# 处理 TF_in_module.txt 文件
# ========================================
cat("正在处理 TF_in_module.txt 文件...\n")

# 定义需要处理的文件
tf_module_files <- list(
  "total" = file.path("total", "TF_in_module.txt"),
  "total.negative" = file.path("total", "TF_negative_in_module.txt"),
  "22C" = file.path("22C", "TF_in_module.txt"),
  "22C.negative" = file.path("22C", "TF_negative_in_module.txt"),
  "16C" = file.path("16C", "TF_in_module.txt"),
  "10C" = file.path("10C", "TF_in_module.txt"),
  "10C.negative" = file.path("10C", "TF_negative_in_module.txt")
)

for (name in names(tf_module_files)) {
  tf_module_file <- tf_module_files[[name]]
  
  # 检查文件是否存在
  if (file.exists(tf_module_file)) {
    cat(sprintf("  处理文件: %s\n", tf_module_file))
    
    # 读取文件（制表符分隔）
    tf_module_data <- read.table(tf_module_file, header = TRUE, sep = "\t", 
                                 stringsAsFactors = FALSE, comment.char = "",
                                 quote = "", fill = TRUE)
    
    # 检查第二列是否存在
    if (ncol(tf_module_data) > 1) {
      # 提取第二列的基因ID
      gene_ids <- tf_module_data[, 2]
      
      # 过滤和清理基因ID
      gene_ids <- as.character(gene_ids)
      gene_ids <- trimws(gene_ids)
      gene_ids <- gene_ids[gene_ids != "" & !is.na(gene_ids)]
      
      if (length(gene_ids) == 0) {
        cat(sprintf("    警告: 没有有效的基因ID\n"))
        next
      }
      
      cat(sprintf("    有效基因ID数量: %d\n", length(gene_ids)))
      
      # 进行注释（包含TF家族信息）
      tryCatch({
        annotations <- annotate_genes_with_tf_family(gene_ids, tf_list_file)
        
        # 将注释结果与原始数据合并
        colnames(tf_module_data)[2] <- "Gene_ID"  # 重命名第二列以便合并
        
        # 清理tf_module_data中的基因ID
        tf_module_data$Gene_ID <- trimws(as.character(tf_module_data$Gene_ID))
        
        annotated_data <- merge(tf_module_data, annotations, 
                               by = "Gene_ID", 
                               all.x = TRUE, 
                               sort = FALSE)
        
        # 调整列的顺序：Gene_ID, Symbol, Family, Gene_Name, GO_Term, GO_Name
        final_cols <- c("Gene_ID", "Symbol", "Family", "Gene_Name", "GO_Term", "GO_Name")
        annotated_data <- annotated_data[, final_cols]
        
        # 保存注释后的文件
        output_file <- file.path(dirname(tf_module_file), 
                                paste0(name, ".TF_in_module.annotated.txt"))
        write.table(annotated_data, output_file, 
                   row.names = FALSE, sep = "\t", quote = FALSE)
        
        cat(sprintf("    -> 保存到: %s\n", output_file))
        cat(sprintf("    -> 注释了 %d 个基因\n", nrow(annotated_data)))
      }, error = function(e) {
        cat(sprintf("    错误: %s\n", e$message))
        cat(sprintf("    跳过该文件\n"))
      })
    } else {
      cat(sprintf("    警告: 文件为空或没有列\n"))
    }
  } else {
    cat(sprintf("  文件不存在，跳过: %s\n", tf_module_file))
  }
}

cat("\n所有文件处理完成！\n")


