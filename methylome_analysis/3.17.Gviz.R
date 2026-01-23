# 1. 安装与加载包
# BiocManager::install(c("Gviz", "GenomicRanges", "biomaRt", "GenomicFeatures"))
library(Gviz)
library(GenomicRanges)
library(biomaRt)
library(GenomicFeatures)
library(txdbmaker)
library(tools)

# 设置工作目录为脚本所在目录 (可选, 如果在命令行运行可能不需要)
setwd("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/RCA")

# 2. 定义基础路径
base_dir <- "../../list/xgboot/TPM_4.5/feature_data_extraction"
categories <- c("Scoupled_specific", "Sexpr_only", "Ssplice_only")

# 3. 定义目标基因信息 (Rca = AT2G39730)
# 根据TAIR10, Rca (AT2G39730) 位于 Chr2:16570746-16573692 (反义链)
# 这里定义一个可视化窗口，稍微扩大一点范围
options(ucscChromosomeNames = FALSE)
target_chr <- "2"
target_start <- 16570746
target_end <- 16573692
target_gene_id <- "AT2G39730" # Rca

# 加载 TxDb
# 使用 biomaRt 在线构建最新的 TxDb 对象 (替代旧的 TxDb 包)
message("Connecting to Ensembl Plants BioMart to build latest TxDb...")
txdb <- makeTxDbFromBiomart(biomart = "plants_mart", 
                            dataset = "athaliana_eg_gene",
                            host = "https://plants.ensembl.org")

# 创建基因结构轨道
# 自动从TxDb获取Rca结构
# method 1: by range (might include other overlapping genes)
# grtrack <- GeneRegionTrack(txdb, chromosome = target_chr, start = target_start, end = target_end, 
#                            name = "Gene Model", transcriptAnnotation = "gene")

# method 2: by gene id (more precise for the specific gene)
# 提取特定基因的转录本
tx_gr <- transcriptsBy(txdb, by = "gene")[[target_gene_id]]
if (!is.null(tx_gr)) {
    # 重新定义范围以聚焦于该基因
    target_chr <- as.character(seqnames(tx_gr)[1])
    target_start <- min(start(tx_gr)) - 2000 # 恢复 2000 缓冲
    target_end <- max(end(tx_gr)) + 2000
    
    print(paste("Plotting Rca (AT2G39730) region:", target_chr, target_start, target_end))
}

grtrack <- GeneRegionTrack(txdb, chromosome = target_chr, start = target_start, end = target_end, 
                           name = "Rca Gene", transcriptAnnotation = "transcript", 
                           showId = TRUE, geneSymbol = TRUE, 
                           collapseTranscripts = FALSE, stacking = "squish",
                           fill = "blue", col = "blue", # 统一蓝色
                           stackHeight = 0.3, # 减小转录本矩形高度(变细)
                           rotation.title = 0) # 轨道标题水平放置

# 过滤：仅保留 Rca 基因 (严格匹配 gene_id)
# 排除该区域的其他基因
# 注意：TxDb创建的轨道，mcols中通常包含 gene 列
if (length(ranges(grtrack)) > 0) {
    # 过滤掉非目标基因的转录本
    # 有时候 gene列可能是列表，需要处理
    mcols_df <- mcols(ranges(grtrack))
    if ("gene" %in% colnames(mcols_df)) {
        # 能够匹配到目标ID的保留（忽略大小写以防万一）
        keep_idx <- which(toupper(unlist(mcols_df$gene)) == toupper(target_gene_id))
        if(length(keep_idx) > 0){
             ranges(grtrack) <- ranges(grtrack)[keep_idx]
        }
    }
}

# 坐标轴轨道
axisTrack <- GenomeAxisTrack()

# 5. 提取基因结构信息用于位点注释
# 获取 Rca 基因的外显子、内含子、UTR 等信息
exons_by_tx <- exonsBy(txdb, by="tx")
cds_by_tx <- cdsBy(txdb, by="tx")
five_utrs <- fiveUTRsByTranscript(txdb)
three_utrs <- threeUTRsByTranscript(txdb)

# 获取 Rca 基因的转录本ID
rca_tx_names <- names(tx_gr)

# 位点注释函数
annotate_feature_position <- function(feature_gr, tx_names, exons_by_tx, cds_by_tx, five_utrs, three_utrs) {
    annotations <- character(length(feature_gr))
    
    for (i in seq_along(feature_gr)) {
        pos <- start(feature_gr[i])
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
                    # 如果在外显子但不在CDS或UTR，可能是未翻译区
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
        
        annotations[i] <- annotation
    }
    
    return(annotations)
}

# 4. 循环处理每个类别并绘图
for (cat in categories) {
    message(paste("Processing:", cat))
    
    input_dir <- file.path(base_dir, cat)
    
    # ---------------------------
    # 读取 SNP 数据
    # ---------------------------
    snp_file <- file.path(input_dir, paste0(cat, "_snp_data.csv"))
    snp_gr <- GRanges()
    
    if (file.exists(snp_file)) {
        snp_df <- read.csv(snp_file, check.names = FALSE)
        # 假设列名: Feature, #Chromosome, Position, ...
        # 需要处理 #Chromosome
        colnames(snp_df)[colnames(snp_df) == "#Chromosome"] <- "chr"
        colnames(snp_df)[colnames(snp_df) == "Position"] <- "start"
        
        # 匹配逻辑 (兼容带Chr和不带Chr的情况)
        input_chr_clean <- gsub("Chr", "", snp_df$chr)
        target_chr_clean <- gsub("Chr", "", target_chr)
        
        # 过滤只在目标染色体的点
        snp_df <- snp_df[input_chr_clean == target_chr_clean & snp_df$start >= target_start & snp_df$start <= target_end, ]
        
        if (nrow(snp_df) > 0) {
            # 强制 seqnames 与 target_chr 一致以确保 Gviz 绘图
            snp_gr <- GRanges(seqnames = target_chr, 
                              ranges = IRanges(start = snp_df$start, end = snp_df$start),
                              type = "SNP",
                              id = snp_df$Feature)
        }
    }
    
    # ---------------------------
    # 读取 甲基化 数据 (CG, CHG, CHH)
    # ---------------------------
    meth_types <- c("CG", "CHG", "CHH")
    meth_gr_list <- list()
    
    for (mtype in meth_types) {
        mfile <- file.path(input_dir, paste0(cat, "_", mtype, "_data.csv"))
        if (file.exists(mfile)) {
            mdf <- read.csv(mfile)
            # 假设列名: Feature, chr, start, end, CV
            
            # 统一染色体名称匹配逻辑
            input_chr_clean <- gsub("Chr", "", mdf$chr)
            target_chr_clean <- gsub("Chr", "", target_chr)
            
            # 过滤
            mdf <- mdf[input_chr_clean == target_chr_clean & mdf$start >= target_start & mdf$end <= target_end, ]
            
            if (nrow(mdf) > 0) {
                # Methylation is 0-based bed (start+1 for 1-based coordinate)
                # 强制 seqnames 与 target_chr 一致
                gr <- GRanges(seqnames = target_chr,
                              ranges = IRanges(start = mdf$start + 1, end = mdf$end),
                              type = paste0("m", mtype), # mCG, mCHG...
                              id = mdf$Feature,
                              score = mdf$CV)
                meth_gr_list[[mtype]] <- gr
            }
        }
    }
    
    # ---------------------------
    # 位点注释和统计
    # ---------------------------
    annotation_results <- data.frame(
        Feature = character(),
        Type = character(),
        Position = integer(),
        Annotation = character(),
        stringsAsFactors = FALSE
    )
    
    # 注释 SNP
    if (length(snp_gr) > 0) {
        snp_annotations <- annotate_feature_position(snp_gr, rca_tx_names, 
                                                      exons_by_tx, cds_by_tx, 
                                                      five_utrs, three_utrs)
        snp_results <- data.frame(
            Feature = snp_gr$id,
            Type = "SNP",
            Position = start(snp_gr),
            Annotation = snp_annotations,
            stringsAsFactors = FALSE
        )
        annotation_results <- rbind(annotation_results, snp_results)
    }
    
    # 注释甲基化位点
    for (mtype in meth_types) {
        if (!is.null(meth_gr_list[[mtype]])) {
            gr <- meth_gr_list[[mtype]]
            meth_annotations <- annotate_feature_position(gr, rca_tx_names,
                                                          exons_by_tx, cds_by_tx,
                                                          five_utrs, three_utrs)
            meth_results <- data.frame(
                Feature = gr$id,
                Type = paste0("m", mtype),
                Position = start(gr),
                Annotation = meth_annotations,
                stringsAsFactors = FALSE
            )
            annotation_results <- rbind(annotation_results, meth_results)
        }
    }
    
    # 保存注释结果
    if (nrow(annotation_results) > 0) {
        annotation_csv <- paste0("Rca_Feature_Annotations_", cat, ".csv")
        write.csv(annotation_results, annotation_csv, row.names = FALSE)
        message(paste("  Saved feature annotations to:", annotation_csv))
        
        # 生成统计摘要
        summary_stats <- as.data.frame(table(annotation_results$Type, annotation_results$Annotation))
        colnames(summary_stats) <- c("Feature_Type", "Location", "Count")
        summary_stats <- summary_stats[summary_stats$Count > 0, ]
        
        summary_csv <- paste0("Rca_Feature_Summary_", cat, ".csv")
        write.csv(summary_stats, summary_csv, row.names = FALSE)
        message(paste("  Saved summary statistics to:", summary_csv))
    }
    
    # ---------------------------
    # 创建轨道
    # ---------------------------
    tracks_to_plot <- list(axisTrack, grtrack)
    
    # SNP Track
    if (length(snp_gr) > 0) {
        atrack_snp <- AnnotationTrack(snp_gr, name = "SNPs", 
                                      shape = "box", fill = "red", stacking = "dense",
                                      rotation.title = 0) # 标题水平放置
        tracks_to_plot[[length(tracks_to_plot) + 1]] <- atrack_snp
    } else {
        # 空轨道占位或跳过
        message("  No SNPs found in region")
    }
    
    # Methylation Tracks (Separated)
    # Define colors for each type
    meth_colors <- c(CG = "darkgreen", CHG = "blue", CHH = "orange")
    
    for (mtype in meth_types) {
        if (!is.null(meth_gr_list[[mtype]])) {
            gr <- meth_gr_list[[mtype]]
            
            # Create track for this methylation type
            atrack_meth <- AnnotationTrack(gr, name = paste0("m", mtype),
                                           shape = "box", # 改回 box 或 ellipse 以确保显示，circle有时需要指定宽高
                                           stacking = "dense", 
                                           fill = meth_colors[mtype], col = meth_colors[mtype],
                                           min.width = 3, # 强制最小显示宽度，防止点太小看不见
                                           size = 1, # 增加轨道高度
                                           rotation.title = 0) # 标题水平放置
            tracks_to_plot[[length(tracks_to_plot) + 1]] <- atrack_meth
        } else {
             message(paste("  No", mtype, "methylation sites found in region"))
        }
    }
    
    # ---------------------------
    # 绘图并保存
    # ---------------------------
    output_pdf <- paste0("Rca_Gene_Map_", cat, ".pdf")
    # 增加PDF宽度(width=15)，以抵消增加缓冲带来的视觉压缩，保持基因的"长条形"美感
    pdf(output_pdf, width = 15, height = 5) # 稍微降低高度，拉长视觉比例
    
    # 计算轨道比例
    # Axis(1) + Gene(5) + Others(1 each) -> 增加基因轨道的相对高度
    track_sizes <- c(1, 5, rep(1, length(tracks_to_plot) - 2))
    
    # 绘制轨道
    plotTracks(tracks_to_plot, 
               from = target_start, 
               to = target_end, 
               chromosome = target_chr,
               sizes = track_sizes, # 设置比例
               background.title = "transparent", 
               col.title = "black",
               main = paste("Features in Rca (AT2G39730) -", cat))
    
    dev.off()
    message(paste("Saved plot to:", output_pdf))
}

print("Done.")


