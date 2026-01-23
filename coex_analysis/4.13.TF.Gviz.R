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

# 2. 加载 TxDb
message("Connecting to Ensembl Plants BioMart to build latest TxDb...")
txdb <- makeTxDbFromBiomart(biomart = "plants_mart", 
                            dataset = "athaliana_eg_gene",
                            host = "https://plants.ensembl.org")

# 3. 提取特定基因的转录本
tx_gr <- transcriptsBy(txdb, by = "gene")[[target_gene_id]]
if (!is.null(tx_gr)) {
    # 重新定义范围以聚焦于该基因
    target_chr <- as.character(seqnames(tx_gr)[1])
    target_start <- min(start(tx_gr)) - 2000 # 增加缓冲
    target_end <- max(end(tx_gr)) + 2000
    
    print(paste("Plotting Rca (AT2G39730) region:", target_chr, target_start, target_end))
}

# 4. 读取 CHH_summary_table.tsv（包含 mCHH 位点和 motif 信息）
input_file <- "Scoupled_specific/CHH_summary_table.tsv"
message(paste("Reading CHH and motif data from:", input_file))

summary_data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message(paste("Read", nrow(summary_data), "CHH sites from input file"))

# mCHH 位点数据直接从 summary_data 获取
mchh_data <- summary_data[, c("Site_ID", "chr", "start", "end")]

# 过滤至少有一个 motif_id 不为 NA 的行
motif_sites <- summary_data[!is.na(summary_data$motif_id_1) | 
                            !is.na(summary_data$motif_id_2) | 
                            !is.na(summary_data$motif_id_3) | 
                            !is.na(summary_data$motif_id_4), ]
message(paste("Found", nrow(motif_sites), "CHH sites with at least one valid motif"))

# 将 motif 数据重构为长格式（每个 motif 一行）
motif_list <- list()
for (i in 1:nrow(motif_sites)) {
    site <- motif_sites[i, ]
    
    # 处理 motif_1
    if (!is.na(site$motif_id_1)) {
        motif_list[[length(motif_list) + 1]] <- data.frame(
            Site_ID = site$Site_ID,
            chr = site$chr,
            motif_id = site$motif_id_1,
            start = site$motif_start_1,
            end = site$motif_stop_1,
            motif_number = 1,
            stringsAsFactors = FALSE
        )
    }
    
    # 处理 motif_2
    if (!is.na(site$motif_id_2)) {
        motif_list[[length(motif_list) + 1]] <- data.frame(
            Site_ID = site$Site_ID,
            chr = site$chr,
            motif_id = site$motif_id_2,
            start = site$motif_start_2,
            end = site$motif_stop_2,
            motif_number = 2,
            stringsAsFactors = FALSE
        )
    }
    
    # 处理 motif_3
    if (!is.na(site$motif_id_3)) {
        motif_list[[length(motif_list) + 1]] <- data.frame(
            Site_ID = site$Site_ID,
            chr = site$chr,
            motif_id = site$motif_id_3,
            start = site$motif_start_3,
            end = site$motif_stop_3,
            motif_number = 3,
            stringsAsFactors = FALSE
        )
    }
    
    # 处理 motif_4
    if (!is.na(site$motif_id_4)) {
        motif_list[[length(motif_list) + 1]] <- data.frame(
            Site_ID = site$Site_ID,
            chr = site$chr,
            motif_id = site$motif_id_4,
            start = site$motif_start_4,
            end = site$motif_stop_4,
            motif_number = 4,
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

# 6. 为 mCHH 位点添加位置注释
if (nrow(mchh_data) > 0) {
    mchh_data$Location <- sapply(mchh_data$start, function(pos) {  # CHH_summary_table.tsv 已是1-based
        annotate_position(pos, rca_tx_names, exons_by_tx, cds_by_tx, five_utrs, three_utrs)
    })
}

# 为每个 motif 添加位置注释
if (nrow(motif_expanded) > 0) {
    motif_expanded$Location <- sapply(motif_expanded$start, function(pos) {
        annotate_position(pos, rca_tx_names, exons_by_tx, cds_by_tx, five_utrs, three_utrs)
    })
    
    # 保存注释结果
    output_csv <- "Scoupled_specific/CHH_motif_annotations.csv"
    write.csv(motif_expanded, output_csv, row.names = FALSE)
    message(paste("Saved motif annotations to:", output_csv))
}

# 7. 创建基因结构轨道
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

# 8. 创建 mCHH 轨道
mchh_chr_clean <- gsub("Chr", "", mchh_data$chr)
mchh_filtered <- mchh_data[mchh_chr_clean == target_chr & 
                           mchh_data$start >= target_start & 
                           mchh_data$end <= target_end, ]
message(paste("Found", nrow(mchh_filtered), "mCHH sites in the target region"))

if (nrow(mchh_filtered) > 0) {
    mchh_gr <- GRanges(seqnames = target_chr,
                      ranges = IRanges(start = mchh_filtered$start,  # 已是1-based
                                     end = mchh_filtered$end),
                      location = mchh_filtered$Location)
    
    mchh_track <- AnnotationTrack(mchh_gr, 
                                 name = "mCHH",
                                 shape = "box",
                                 fill = "orange",
                                 col = "orange",
                                 stacking = "dense",
                                 rotation.title = 0,
                                 min.width = 3,
                                 size = 1)
}

# 9. 创建 Motif 轨道（按 motif_number 分组）
motif_tracks <- list()
if (nrow(motif_expanded) > 0) {
    # 转换染色体名称并过滤
    motif_expanded$chr_clean <- gsub("Chr", "", motif_expanded$chr)
    motif_filtered <- motif_expanded[motif_expanded$chr_clean == target_chr & 
                                     motif_expanded$start >= target_start & 
                                     motif_expanded$end <= target_end, ]
    
    message(paste("Found", nrow(motif_filtered), "motifs in the target region"))
    
    # 为每个 motif_number 创建独立轨道
    for (i in 1:4) {
        motif_subset <- motif_filtered[motif_filtered$motif_number == i, ]
        
        if (nrow(motif_subset) > 0) {
            motif_gr <- GRanges(seqnames = target_chr,
                               ranges = IRanges(start = motif_subset$start, 
                                              end = motif_subset$end),
                               motif_id = motif_subset$motif_id,
                               location = motif_subset$Location)
            
            motif_tracks[[i]] <- AnnotationTrack(motif_gr, 
                                                name = paste0("Motif_", i),
                                                shape = "box",
                                                fill = "darkgreen",
                                                col = "darkgreen",
                                                stacking = "dense",
                                                rotation.title = 0,
                                                min.width = 3,
                                                size = 1)
            message(paste("Created track for Motif", i, "with", nrow(motif_subset), "elements"))
        }
    }
}

# 10. 绘图
output_pdf <- "Scoupled_specific/Rca_CHH_and_Motifs.pdf"
pdf(output_pdf, width = 15, height = 8)

# 组装轨道
tracks_to_plot <- list(axisTrack, grtrack)
track_sizes <- c(1, 5)

# 添加 mCHH 轨道
if (exists("mchh_track")) {
    tracks_to_plot <- c(tracks_to_plot, list(mchh_track))
    track_sizes <- c(track_sizes, 1)
}

# 添加 motif 轨道
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
           main = "mCHH and TF Motifs in Rca (AT2G39730)")

dev.off()
message(paste("Saved plot to:", output_pdf))

# 11. 生成统计摘要
if (nrow(motif_filtered) > 0) {
    summary_stats <- as.data.frame(table(motif_filtered$Location, motif_filtered$motif_number))
    colnames(summary_stats) <- c("Location", "Motif_Number", "Count")
    summary_stats <- summary_stats[summary_stats$Count > 0, ]
    
    summary_csv <- "Scoupled_specific/CHH_motif_summary.csv"
    write.csv(summary_stats, summary_csv, row.names = FALSE)
    message(paste("Saved summary statistics to:", summary_csv))
    
    print(summary_stats)
}

print("Done.")
