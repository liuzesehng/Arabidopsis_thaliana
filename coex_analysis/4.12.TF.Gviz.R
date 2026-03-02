# mCHH 和 TF Motif 可视化脚本
# 基于 CHH_summary_table.tsv 文件绘制 Rca 基因上的 mCHH 位点和对应的 TF motif

library(Gviz)
library(GenomicRanges)
library(biomaRt)
library(GenomicFeatures)
library(txdbmaker)
library(tools)

# 设置工作目录
setwd("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/RCA")

# 1. 定义目标基因信息 (Rca = AT2G39730)
options(ucscChromosomeNames = FALSE)
target_chr <- "2"
target_gene_id <- "AT2G39730" # Rca

# 2. 加载 TxDb (从本地 GFF3 文件)
message("Building TxDb from local GFF3 file...")
gff3_file <- "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/refgen/Arabidopsis_thaliana.TAIR10.62.gff3"
txdb <- makeTxDbFromGFF(gff3_file, 
                        format = "gff3",
                        organism = "Arabidopsis thaliana")

# 3. 提取特定基因的转录本
tx_gr <- transcriptsBy(txdb, by = "gene")[[target_gene_id]]
if (!is.null(tx_gr)) {
    # 重新定义范围以聚焦于该基因
    target_chr <- as.character(seqnames(tx_gr)[1])
    target_start <- min(start(tx_gr)) - 2000 # 增加缓冲
    target_end <- max(end(tx_gr)) + 2000
    
    print(paste("Plotting Rca (AT2G39730) region:", target_chr, target_start, target_end))
}

# 4. 读取多个 intersection.csv 文件
base_dir <- "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana"
analysis_dir <- file.path(base_dir, "list/xgboot/TPM_4.5/feature_analysis")

# 查找所有intersection.csv文件
all_files <- list.files(analysis_dir, pattern = "_(SA|SB|SPα|SPβ|ST)_intersection\\.csv$", full.names = TRUE)
message(paste("Found", length(all_files), "intersection.csv files"))

# 读取并合并所有数据
all_data_list <- list()
all_columns <- character(0)

# 第一遍：收集所有可能的列名
for (file_path in all_files) {
    data <- read.csv(file_path, stringsAsFactors = FALSE, nrows = 1)
    all_columns <- union(all_columns, colnames(data))
}
message(paste("Total unique columns found:", length(all_columns)))

# 第二遍：读取数据并统一列结构
for (file_path in all_files) {
    # 提取文件类型信息
    filename <- basename(file_path)
    # 匹配模式: {type}_{SA|SB|SPα|SPβ|ST}_intersection.csv
    pattern_match <- regmatches(filename, regexec("^(.+)_(SA|SB|SPα|SPβ|ST)_intersection\\.csv$", filename))
    if (length(pattern_match[[1]]) == 3) {
        condition_type <- pattern_match[[1]][2]  # 如 "Down-like_HighTemp"
        methylation_type <- pattern_match[[1]][3]  # 如 "SA", "SB" 等
        
        # 读取文件
        data <- read.csv(file_path, stringsAsFactors = FALSE)
        if (nrow(data) > 0) {
            # 为缺失的列添加 NA
            for (col in all_columns) {
                if (!col %in% colnames(data)) {
                    data[[col]] <- NA
                }
            }
            # 确保列的顺序一致
            data <- data[, all_columns]
            
            data$condition_type <- condition_type
            data$methylation_type <- methylation_type
            data$file_source <- filename
            all_data_list[[length(all_data_list) + 1]] <- data
            message(paste("Read", nrow(data), "sites from:", filename))
        }
    }
}

# 合并所有数据
if (length(all_data_list) > 0) {
    summary_data <- do.call(rbind, all_data_list)
    message(paste("Total sites:", nrow(summary_data)))
} else {
    stop("No valid data found in intersection.csv files")
}

# 创建位点数据（添加 end 列）
mchh_data <- summary_data[, c("Site_ID", "chr", "position", "condition_type", "methylation_type")]
mchh_data$start <- mchh_data$position
mchh_data$end <- mchh_data$position

# 提取并统计位点类型
extract_site_type <- function(site_id) {
    # 使用不区分大小写的正则匹配
    type_match <- regmatches(site_id, regexec("^(CG|CHG|CHH|SNP|cg|chg|chh|snp)_", site_id))
    if (length(type_match[[1]]) >= 2) {
        # 统一转换为大写
        return(toupper(type_match[[1]][2]))
    }
    return("Unknown")
}

mchh_data$site_type <- sapply(mchh_data$Site_ID, extract_site_type)

# 统计各类型位点数量
site_type_counts <- table(mchh_data$site_type, mchh_data$condition_type)
message("\n=== Site Type Statistics ===")
print(site_type_counts)
message("===========================\n")

# 过滤至少有一个 motif_id 不为 NA 的行
motif_sites <- summary_data[!is.na(summary_data$motif_id_1) | 
                            !is.na(summary_data$motif_id_2) | 
                            !is.na(summary_data$motif_id_3) | 
                            !is.na(summary_data$motif_id_4), ]
message(paste("Found", nrow(motif_sites), "sites with at least one valid motif"))

# 将 motif 数据重构为长格式（每个 motif 一行）
motif_list <- list()
for (i in 1:nrow(motif_sites)) {
    site <- motif_sites[i, ]
    site_type <- extract_site_type(site$Site_ID)
    
    # 处理 motif_1
    if (!is.na(site$motif_id_1) && site$motif_id_1 != "NA") {
        motif_list[[length(motif_list) + 1]] <- data.frame(
            Site_ID = site$Site_ID,
            chr = site$chr,
            position = site$position,
            motif_id = site$motif_id_1,
            start = as.numeric(site$motif_start_1),
            end = as.numeric(site$motif_stop_1),
            motif_number = 1,
            condition_type = site$condition_type,
            methylation_type = site$methylation_type,
            site_type = site_type,
            stringsAsFactors = FALSE
        )
    }
    
    # 处理 motif_2
    if (!is.na(site$motif_id_2) && site$motif_id_2 != "NA") {
        motif_list[[length(motif_list) + 1]] <- data.frame(
            Site_ID = site$Site_ID,
            chr = site$chr,
            position = site$position,
            motif_id = site$motif_id_2,
            start = as.numeric(site$motif_start_2),
            end = as.numeric(site$motif_stop_2),
            motif_number = 2,
            condition_type = site$condition_type,
            methylation_type = site$methylation_type,
            site_type = site_type,
            stringsAsFactors = FALSE
        )
    }
    
    # 处理 motif_3
    if (!is.na(site$motif_id_3) && site$motif_id_3 != "NA") {
        motif_list[[length(motif_list) + 1]] <- data.frame(
            Site_ID = site$Site_ID,
            chr = site$chr,
            position = site$position,
            motif_id = site$motif_id_3,
            start = as.numeric(site$motif_start_3),
            end = as.numeric(site$motif_stop_3),
            motif_number = 3,
            condition_type = site$condition_type,
            methylation_type = site$methylation_type,
            site_type = site_type,
            stringsAsFactors = FALSE
        )
    }
    
    # 处理 motif_4
    if (!is.na(site$motif_id_4) && site$motif_id_4 != "NA") {
        motif_list[[length(motif_list) + 1]] <- data.frame(
            Site_ID = site$Site_ID,
            chr = site$chr,
            position = site$position,
            motif_id = site$motif_id_4,
            start = as.numeric(site$motif_start_4),
            end = as.numeric(site$motif_stop_4),
            motif_number = 4,
            condition_type = site$condition_type,
            methylation_type = site$methylation_type,
            site_type = site_type,
            stringsAsFactors = FALSE
        )
    }
}

# 合并所有 motif
if (length(motif_list) > 0) {
    motif_expanded <- do.call(rbind, motif_list)
    message(paste("Expanded to", nrow(motif_expanded), "individual motifs"))
} else {
    motif_expanded <- data.frame()
    message("No valid motifs found")
}

# 5. 提取基因结构信息用于位点注释
exons_by_tx <- exonsBy(txdb, by="tx")
cds_by_tx <- cdsBy(txdb, by="tx")
five_utrs <- fiveUTRsByTranscript(txdb)
three_utrs <- threeUTRsByTranscript(txdb)
rca_tx_names <- names(tx_gr)

# 位点注释函数
annotate_position <- function(pos, tx_names, exons_by_tx, cds_by_tx, five_utrs, three_utrs) {
    annotation <- "Intergenic"
    
    # 检查是否在任何转录本的外显子中
    for (tx_name in tx_names) {
        if (tx_name %in% names(exons_by_tx)) {
            exons <- exons_by_tx[[tx_name]]
            if (any(pos >= start(exons) & pos <= end(exons))) {
                # 在外显子中，进一步判断是否在CDS或UTR
                if (tx_name %in% names(cds_by_tx)) {
                    cds <- cds_by_tx[[tx_name]]
                    if (any(pos >= start(cds) & pos <= end(cds))) {
                        annotation <- "CDS"
                        break
                    }
                }
                if (tx_name %in% names(five_utrs)) {
                    utr5 <- five_utrs[[tx_name]]
                    if (any(pos >= start(utr5) & pos <= end(utr5))) {
                        annotation <- "5'UTR"
                        break
                    }
                }
                if (tx_name %in% names(three_utrs)) {
                    utr3 <- three_utrs[[tx_name]]
                    if (any(pos >= start(utr3) & pos <= end(utr3))) {
                        annotation <- "3'UTR"
                        break
                    }
                }
                if (annotation == "Intergenic") {
                    annotation <- "Exon"
                }
                break
            }
        }
    }
    
    # 如果不在外显子中，检查是否在基因范围内（则为内含子）
    if (annotation == "Intergenic") {
        gene_start <- min(start(tx_gr))
        gene_end <- max(end(tx_gr))
        gene_strand <- as.character(strand(tx_gr)[1])
        
        if (pos >= gene_start & pos <= gene_end) {
            annotation <- "Intron"
        } else if (gene_strand == "-") {
            # 负链基因：end右侧是上游，start左侧是下游
            if (pos > gene_end) {
                annotation <- "Upstream"
            } else if (pos < gene_start) {
                annotation <- "Downstream"
            }
        } else {
            # 正链基因：start左侧是上游，end右侧是下游
            if (pos < gene_start) {
                annotation <- "Upstream"
            } else if (pos > gene_end) {
                annotation <- "Downstream"
            }
        }
    }
    
    return(annotation)
}

# 6. 为位点添加位置注释（如果需要）
if (nrow(mchh_data) > 0) {
    mchh_data$Location <- sapply(mchh_data$start, function(pos) {
        annotate_position(pos, rca_tx_names, exons_by_tx, cds_by_tx, five_utrs, three_utrs)
    })
}

# 为每个 motif 添加位置注释
if (nrow(motif_expanded) > 0) {
    motif_expanded$Location <- sapply(motif_expanded$start, function(pos) {
        annotate_position(pos, rca_tx_names, exons_by_tx, cds_by_tx, five_utrs, three_utrs)
    })
}

# 获取所有的条件类型（用于生成多个图）
condition_types <- unique(summary_data$condition_type)
message(paste("Found", length(condition_types), "condition types:", paste(condition_types, collapse = ", ")))

# 为每种条件类型生成一个图
for (cond_type in condition_types) {
    message(paste("\n=== Processing condition:", cond_type, "==="))
    
    # 过滤当前条件的数据
    cond_mchh <- mchh_data[mchh_data$condition_type == cond_type, ]
    cond_motif <- motif_expanded[motif_expanded$condition_type == cond_type, ]
    
    if (nrow(cond_mchh) == 0) {
        message(paste("No data for condition:", cond_type))
        next
    }
    
    # 8. 创建基因结构轨道
    grtrack <- GeneRegionTrack(txdb, chromosome = target_chr, start = target_start, end = target_end, 
                               name = "Rca Gene", transcriptAnnotation = "transcript", 
                               showId = TRUE, geneSymbol = TRUE, 
                               collapseTranscripts = FALSE, stacking = "squish",
                               fill = "blue", col = "blue",
                               stackHeight = 0.3,
                               rotation.title = 0)
    
    # 过滤：仅保留 Rca 基因
    if (length(ranges(grtrack)) > 0) {
        mcols_df <- mcols(ranges(grtrack))
        if ("gene" %in% colnames(mcols_df)) {
            keep_idx <- which(toupper(unlist(mcols_df$gene)) == toupper(target_gene_id))
            if(length(keep_idx) > 0){
                ranges(grtrack) <- ranges(grtrack)[keep_idx]
            }
        }
    }
    
    # 坐标轴轨道
    axisTrack <- GenomeAxisTrack()
    
    # 9. 为每种甲基化类型创建轨道
    methylation_types <- c("SA", "SB", "SPα", "SPβ", "ST")
    meth_tracks <- list()
    
    # 定义颜色方案
    meth_colors <- c(
        "CG" = "red",
        "CHG" = "blue", 
        "CHH" = "orange",
        "SNP" = "purple"
    )
    
    for (meth_type in methylation_types) {
        # 过滤当前甲基化类型的位点
        meth_sites <- cond_mchh[cond_mchh$methylation_type == meth_type, ]
        
        if (nrow(meth_sites) > 0) {
            # 过滤到目标区域（染色体名称清理）
            meth_chr_clean <- gsub("Chr", "", meth_sites$chr)
            meth_filtered <- meth_sites[meth_chr_clean == target_chr & 
                                        meth_sites$start >= target_start & 
                                        meth_sites$end <= target_end, ]
            
            if (nrow(meth_filtered) > 0) {
                # 统计当前甲基化类型下的位点类型分布
                type_counts <- table(meth_filtered$site_type)
                message(paste(meth_type, "site types in region:", paste(names(type_counts), "=", type_counts, collapse=", ")))
                
                # 按位点类型分别创建 GRanges
                for (site_type in c("CG", "CHG", "CHH", "SNP")) {
                    type_sites <- meth_filtered[meth_filtered$site_type == site_type, ]
                    if (nrow(type_sites) > 0) {
                        meth_gr <- GRanges(
                            seqnames = target_chr,
                            ranges = IRanges(start = type_sites$start, end = type_sites$end),
                            location = type_sites$Location
                        )
                        
                        track_name <- paste0(meth_type, "_", site_type)
                        meth_tracks[[track_name]] <- AnnotationTrack(
                            meth_gr,
                            name = track_name,
                            shape = "box",
                            fill = meth_colors[site_type],
                            col = meth_colors[site_type],
                            stacking = "dense",
                            rotation.title = 0,
                            min.width = 3,
                            size = 1
                        )
                        message(paste("  -> Created track:", track_name, "with", nrow(type_sites), "sites"))
                    }
                }
            } else {
                message(paste(meth_type, "- No sites in target region"))
            }
        }
    }
    
    # 10. 创建 TF Motif 轨道
    motif_tracks <- list()
    if (nrow(cond_motif) > 0) {
        # 转换染色体名称并过滤
        cond_motif$chr_clean <- gsub("Chr", "", cond_motif$chr)
        motif_filtered <- cond_motif[cond_motif$chr_clean == target_chr & 
                                     cond_motif$start >= target_start & 
                                     cond_motif$end <= target_end, ]
        
        message(paste("Found", nrow(motif_filtered), "motifs in the target region"))
        
        # 获取所有唯一的TF ID，并为每个创建轨道
        unique_tfs <- unique(motif_filtered$motif_id)
        
        for (tf_id in unique_tfs) {
            tf_motifs <- motif_filtered[motif_filtered$motif_id == tf_id, ]
            
            if (nrow(tf_motifs) > 0) {
                motif_gr <- GRanges(
                    seqnames = target_chr,
                    ranges = IRanges(start = tf_motifs$start, end = tf_motifs$end),
                    motif_id = tf_motifs$motif_id,
                    location = tf_motifs$Location,
                    site_type = tf_motifs$site_type
                )
                
                motif_tracks[[tf_id]] <- AnnotationTrack(
                    motif_gr,
                    name = tf_id,  # 使用 TF ID 作为轨道名称
                    shape = "box",
                    fill = "darkgreen",
                    col = "darkgreen",
                    stacking = "dense",
                    rotation.title = 0,
                    min.width = 3,
                    size = 1
                )
                message(paste("Created TF track:", tf_id, "with", nrow(tf_motifs), "motifs"))
            }
        }
    }
    
    # 11. 绘图
    output_dir <- file.path(base_dir, "list/xgboot/TPM_4.5/feature_analysis/plots")
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    output_pdf <- file.path(output_dir, paste0(cond_type, "_Rca_Methylation_TF.pdf"))
    pdf(output_pdf, width = 15, height = max(8, 2 + length(meth_tracks) + length(motif_tracks)))
    
    # 组装轨道
    tracks_to_plot <- list(axisTrack, grtrack)
    track_sizes <- c(1, 5)
    
    # 添加甲基化轨道
    if (length(meth_tracks) > 0) {
        tracks_to_plot <- c(tracks_to_plot, meth_tracks)
        track_sizes <- c(track_sizes, rep(1, length(meth_tracks)))
    }
    
    # 添加 TF motif 轨道
    if (length(motif_tracks) > 0) {
        tracks_to_plot <- c(tracks_to_plot, motif_tracks)
        track_sizes <- c(track_sizes, rep(1, length(motif_tracks)))
    }
    
    plotTracks(tracks_to_plot, 
               from = target_start, 
               to = target_end, 
               chromosome = target_chr,
               sizes = track_sizes,
               background.title = "transparent", 
               col.title = "black",
               main = paste("Methylation Sites and TF Motifs in Rca (AT2G39730) -", cond_type))
    
    dev.off()
    message(paste("Saved plot to:", output_pdf))
}

message("\n=================================")
message("All plots generated successfully!")
message("=================================")
