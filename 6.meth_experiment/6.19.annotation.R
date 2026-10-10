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
    annotation_df$Family <- NA_character_
  }
  
  return(annotation_df)
}

# ========================================
# 注释 TF.non_coexpressed.txt 中的全部基因
# ========================================
input_file <- "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/meth_experiment/CHH_TF_loci/TF.non_coexpressed.txt"
output_file <- file.path(dirname(input_file), "TF.non_coexpressed.annotated.txt")
tf_list_file <- "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/coex2/Ath_TF_list.txt"

if (!file.exists(input_file)) {
  stop(sprintf("输入文件不存在: %s", input_file))
}

tf_data <- read.table(input_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, colClasses = "character",
                      comment.char = "", quote = "", check.names = FALSE)
if (!"tf_id" %in% colnames(tf_data)) {
  stop("输入文件缺少 tf_id 列")
}

gene_ids <- trimws(tf_data$tf_id)
gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]
if (length(gene_ids) == 0) {
  stop("输入文件没有有效的基因ID")
}

cat(sprintf("正在注释 %d 行、%d 个唯一基因...\n",
            length(gene_ids), length(unique(gene_ids))))
annotations <- annotate_genes_with_tf_family(gene_ids, tf_list_file)

# 按输入顺序取回注释，保留重复基因行，不限制基因数量。
annotated_data <- annotations[match(gene_ids, annotations$Gene_ID),
                              c("Gene_ID", "Symbol", "Family", "Gene_Name",
                                "GO_Term", "GO_Name"), drop = FALSE]
write.table(annotated_data, output_file, row.names = FALSE,
            sep = "\t", quote = FALSE)
cat(sprintf("注释完成，保存到: %s\n", output_file))
