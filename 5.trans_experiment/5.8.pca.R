#!/usr/bin/env Rscript
# Arabidopsis Col-0, 3 d temperature treatments: sample QC and exploratory PCA.
# Run: Rscript 5.8.pca.R [--input-dir DIR] [--output-dir DIR]
#                       [--top-n 1000] [--validate-only]
# Required: tximport, DESeq2, ggplot2, ggrepel, jsonlite, systemfonts; Cairo.
# Raw inputs are read only. No sample exclusion, jitter, ellipse or batch correction.

options(stringsAsFactors = FALSE, warn = 1)
main <- function() {
  base <- "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng"
  project <- file.path(base, "Arabidopsis_thaliana")
  opt <- list(input_dir = file.path(project, "list/rnadata"),
              output_dir = file.path(project, "list/RCA/PCA"),
              top_n = 1000L, validate_only = FALSE)
  args <- commandArgs(trailingOnly = TRUE)
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--help") {
      cat("Usage: Rscript 5.8.pca.R [--input-dir DIR] [--output-dir DIR]",
          "[--top-n 1000] [--validate-only]\n")
      return(invisible(NULL))
    }
    if (key == "--validate-only") opt$validate_only <- TRUE else {
      if (!key %in% c("--input-dir", "--output-dir", "--top-n") ||
          i == length(args)) stop("Unknown option or missing value: ", key)
      i <- i + 1L
      opt[[gsub("-", "_", substring(key, 3))]] <- args[i]
    }
    i <- i + 1L
  }
  opt$top_n <- suppressWarnings(as.numeric(opt$top_n))
  stopifnot(length(opt$top_n) == 1L, is.finite(opt$top_n), opt$top_n >= 2,
            opt$top_n == as.integer(opt$top_n))
  packages <- c("tximport", "DESeq2", "ggplot2", "ggrepel", "jsonlite", "systemfonts")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  if (!capabilities("cairo")) stop("R requires Cairo support for PDF and PNG export")
  suppressPackageStartupMessages(library(ggplot2))
  check <- function(ok, why) if (!isTRUE(ok)) stop(why, call. = FALSE)
  read_tsv <- function(path) read.delim(path, check.names = FALSE, quote = "",
                                       comment.char = "", stringsAsFactors = FALSE)
  norm <- function(path) normalizePath(path, mustWork = TRUE)
  json <- function(path) jsonlite::fromJSON(path)
  temps <- seq(4L, 32L, 4L)
  expected <- c("T4_2", "T4_4", "T4_5", "T8_1", "T8_2", "T8_4",
                "T12_2", "T12_3", "T12_6", "T16_3", "T16_5", "T16_6",
                "T20_1", "T20_2", "T20_4", "T24_1", "T24_2", "T24_4",
                "T28_1", "T28_4", "T28_6", "T32_1", "T32_3", "T32_4")
  files <- list.files(opt$input_dir, pattern = "^quant\\.genes\\.sf$", recursive = TRUE,
                      full.names = TRUE)
  check(length(files) == 24L, paste("Expected 24 quant.genes.sf files; found", length(files)))
  sample_ids <- basename(dirname(dirname(files)))
  check(!anyDuplicated(sample_ids) && setequal(sample_ids, expected),
        "Sample IDs do not match the verified 24-sample manifest")
  files <- files[match(expected, sample_ids)]
  names(files) <- expected
  samples <- data.frame(sample_id = expected,
                         temperature_C = as.integer(sub("^T([0-9]+)_.*", "\\1", expected)),
                         replicate_id = as.integer(sub("^T[0-9]+_", "", expected)),
                         quant_file = vapply(files, norm, character(1)))
  samples$is_control <- samples$temperature_C == 24L
  check(all(table(factor(samples$temperature_C, levels = temps)) == 3L),
        "Each of the eight temperatures must have exactly three samples")
  dirs <- dirname(files)
  cmd_files <- file.path(dirs, "cmd_info.json")
  meta_files <- file.path(dirs, "aux_info/meta_info.json")
  sample_source_files <- c(unname(files), cmd_files, meta_files)
  sample_source_md5 <- tools::md5sum(sample_source_files)
  cmds <- lapply(cmd_files, json)
  metas <- lapply(meta_files, json)
  for (j in seq_along(files)) {
    check(all(vapply(c("mates1", "mates2", "output"), function(k)
      grepl(paste0("(^|/)", expected[j], "(/|$)"), cmds[[j]][[k]]), logical(1))),
      paste("Salmon input/output sample mismatch:", expected[j]))
    check(length(metas[[j]]$quant_errors) == 0L,
          paste("Salmon quantification errors:", expected[j]))
  }
  index_paths <- unique(vapply(cmds, function(x) norm(x$index), character(1)))
  recorded_gene_maps <- unique(vapply(cmds, function(x) x$geneMap, character(1)))
  check(length(index_paths) == 1L && length(recorded_gene_maps) == 1L,
        "Samples used inconsistent Salmon indexes or recorded gene maps")
  index_file <- file.path(index_paths, "info.json")
  idx <- json(index_file)
  for (m in metas) {
    check(identical(m$index_seq_hash, idx$SeqHash) &&
            identical(m$index_name_hash, idx$NameHash), "Salmon reference hash mismatch")
  }
  ref <- dirname(index_paths)
  gff_path <- file.path(ref, "refgen/GCF_000001735.4_TAIR10.1_genomic.gff")
  gff <- read.delim(gff_path, header = FALSE, comment.char = "#", quote = "",
                    colClasses = "character")
  gene_records <- gff[[9]][grepl("(^|;)ID=gene-", gff[[9]])]
  annotated_genes <- unique(sub(".*(?:^|;)ID=gene-([^;]+).*", "\\1",
                                gene_records, perl = TRUE))
  rm(gff)
  quant <- lapply(files, read_tsv)
  ids <- sort(quant[[1]]$Name, method = "radix")
  check(!anyDuplicated(ids) && all(grepl("^AT[1-5CM]G[0-9]+$", ids)),
        "quant.genes.sf must contain unique Arabidopsis gene IDs, not transcript IDs")
  check(all(ids %in% annotated_genes), "Gene IDs absent from the reference GFF")
  for (j in seq_along(quant)) {
    q <- quant[[j]]
    check(identical(names(q), c("Name", "Length", "EffectiveLength", "TPM", "NumReads")),
          paste("Unexpected quant.genes.sf schema:", expected[j]))
    check(!anyDuplicated(q$Name) && setequal(q$Name, ids),
          paste("Gene set mismatch or duplicate gene ID:", expected[j]))
    check(all(vapply(q[-1], is.numeric, logical(1))) &&
            all(is.finite(as.matrix(q[-1]))) && all(as.matrix(q[-1]) >= 0) &&
            all(q$Length > 0) && all(q$EffectiveLength > 0) && sum(q$NumReads) > 0,
          paste("Invalid gene quantification values:", expected[j]))
    quant[[j]] <- q[match(ids, q$Name), , drop = FALSE]
  }
  samples$percent_mapped <- vapply(metas, function(x) x$percent_mapped, numeric(1))
  samples$fragments_processed <- vapply(metas, function(x) x$num_processed, numeric(1))
  samples$estimated_fragments <- vapply(quant, function(q) sum(q$NumReads), numeric(1))
  message("Validation passed: 24 quant.genes.sf files; ", length(ids),
          " unique gene IDs per sample; gene IDs verified against reference GFF.")
  if (opt$validate_only) return(invisible(samples))
  # txIn=FALSE imports each row as an already quantified gene. No tx2gene
  # aggregation is performed. Explicitly align rows: this tximport branch does
  # not itself verify matching gene order across input files.
  import_genes <- function(f) {
    q <- read_tsv(f)
    check(!anyDuplicated(q$Name) && setequal(q$Name, ids),
          paste("Gene IDs changed during import:", f))
    q[match(ids, q$Name), , drop = FALSE]
  }
  txi <- tximport::tximport(files, type = "salmon", txIn = FALSE, txOut = FALSE,
                            geneIdCol = "Name", countsFromAbundance = "no",
                            importer = import_genes)
  check(identical(rownames(txi$counts), ids) && identical(colnames(txi$counts), expected),
        "Gene or sample ordering changed during import")
  check(all(txi$counts == do.call(cbind, lapply(quant, function(q) q$NumReads))) &&
          all(txi$length == do.call(cbind, lapply(quant, function(q) q$EffectiveLength))),
        "Imported counts/lengths do not match quant.genes.sf")
  rm(quant)
  coldata <- data.frame(row.names = expected, temperature = factor(samples$temperature_C))
  dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = coldata, design = ~1)
  counts <- DESeq2::counts(dds)
  keep <- rowSums(counts >= 10) >= 3L
  check(sum(keep) >= 2L, "Too few genes pass the count filter")
  vst_messages <- character()
  vst <- withCallingHandlers(
    DESeq2::varianceStabilizingTransformation(dds[keep, ], blind = TRUE,
                                              fitType = "parametric"),
    message = function(m) vst_messages <<- c(vst_messages, conditionMessage(m)))
  # DESeqTransform does not retain dispersionFunction(). DESeq2 explicitly reports
  # its sole automatic parametric-to-local fallback; preserve that message log.
  fit_type <- if (any(grepl("local regression fit was automatically substituted",
                            vst_messages, fixed = TRUE))) "local" else "parametric"
  mat <- SummarizedExperiment::assay(vst)
  check(all(is.finite(mat)), "VST contains non-finite values")
  gene_var <- apply(mat, 1, var)
  eligible <- names(gene_var)[is.finite(gene_var) & gene_var > 0]
  ranked <- eligible[order(-gene_var[eligible], eligible, method = "radix")]
  selected <- head(ranked, opt$top_n)
  check(length(selected) >= 2L, "Too few nonconstant genes for PCA")
  run_pca <- function(gene_ids) {
    p <- prcomp(t(mat[gene_ids, , drop = FALSE]), center = TRUE, scale. = FALSE)
    # Make the otherwise arbitrary PC sign reproducible using the largest loading.
    signs <- vapply(seq_len(ncol(p$rotation)), function(k) {
      sign(p$rotation[which.max(abs(p$rotation[, k])), k])
    }, numeric(1))
    signs[signs == 0] <- 1
    p$x <- sweep(p$x, 2, signs, "*")
    p$rotation <- sweep(p$rotation, 2, signs, "*")
    p
  }
  pcs <- list(top = run_pca(selected), all = run_pca(sort(eligible)))
  variance <- function(p) {
    proportion <- p$sdev^2 / sum(p$sdev^2)
    data.frame(PC = colnames(p$x), eigenvalue = p$sdev^2,
               proportion = proportion, percent = 100 * proportion,
               cumulative_percent = 100 * cumsum(proportion))
  }
  variances <- lapply(pcs, variance)
  stopifnot(all(vapply(variances, function(v) abs(sum(v$proportion) - 1) < 1e-12,
                       logical(1))), identical(rownames(pcs$top$x), expected))
  score_table <- function(p) cbind(samples[, c("sample_id", "temperature_C", "replicate_id",
                                               "is_control")], as.data.frame(p$x))
  scores <- lapply(pcs, score_table)
  gene_table <- data.frame(gene_id = rownames(counts),
                           samples_count_ge10 = rowSums(counts >= 10),
                           passes_count_filter = keep,
                           vst_variance = gene_var[match(rownames(counts), names(gene_var))],
                           variance_rank = match(rownames(counts), ranked),
                           in_top_pca = rownames(counts) %in% selected,
                           in_all_pca = rownames(counts) %in% eligible)
  qc_metrics <- function(gene_ids, label) {
    m <- mat[gene_ids, , drop = FALSE]
    distance <- as.matrix(dist(t(m)))
    correlations <- cor(m, method = "pearson")
    rows <- lapply(seq_len(nrow(samples)), function(j) {
      others <- setdiff(seq_len(nrow(samples)), j)
      within <- others[samples$temperature_C[others] == samples$temperature_C[j]]
      nearest <- others[which.min(distance[j, others])]
      data.frame(analysis = label, sample_id = expected[j],
                 temperature_C = samples$temperature_C[j],
                 mean_within_group_distance = mean(distance[j, within]),
                 min_within_group_correlation = min(correlations[j, within]),
                 nearest_sample = expected[nearest],
                 nearest_sample_same_temperature = samples$temperature_C[nearest] ==
                   samples$temperature_C[j])
    })
    list(summary = do.call(rbind, rows), distances = distance, correlations = correlations)
  }
  qc <- list(top = qc_metrics(selected, "top_variable"),
              all = qc_metrics(eligible, "all_eligible"))
  available_fonts <- systemfonts::system_fonts()$family
  font <- if ("Arial" %in% available_fonts) "Arial" else if ("DejaVu Sans" %in%
    available_fonts) "DejaVu Sans" else "sans"
  palette <- setNames(c("#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677",
                         "#AA3377", "#BBBBBB", "#000000"), temps)
  labels <- paste0(temps, " °C")
  make_plot <- function(which, labeled = FALSE) {
    d <- scores[[which]]
    d$temperature <- factor(d$temperature_C, levels = temps)
    v <- variances[[which]]$percent
    p <- ggplot(d, aes(PC1, PC2, colour = temperature)) +
      geom_point(size = 3, alpha = 0.85) +
      scale_colour_manual(values = palette, breaks = as.character(temps), labels = labels,
                           drop = FALSE) +
      scale_x_continuous(expand = expansion(mult = if (labeled) 0.18 else 0.09)) +
      scale_y_continuous(expand = expansion(mult = if (labeled) 0.18 else 0.09)) +
      coord_cartesian() +
      labs(x = sprintf("PC1 (%.2f%%)", v[1]), y = sprintf("PC2 (%.2f%%)", v[2]),
           colour = "Temperature") +
      theme_classic(base_size = 10, base_family = font) +
      theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
            axis.ticks = element_line(linewidth = 0.35, colour = "black"),
            axis.text = element_text(size = 9, colour = "black"),
            axis.title = element_text(size = 10),
            legend.text = element_text(size = 9), legend.title = element_text(size = 10),
            legend.key.height = grid::unit(5.5, "mm"),
            plot.margin = margin(8, 8, 8, 8),
            plot.background = element_rect(fill = "white", colour = NA))
    if (labeled) p <- p + ggrepel::geom_text_repel(aes(label = sample_id),
      family = font, colour = "#303030", size = 2.8, seed = 20261010,
      nudge_x = ifelse(d$sample_id == "T8_1", -0.08 * diff(range(d$PC1)), 0),
      max.overlaps = Inf, max.iter = 100000, max.time = Inf,
      box.padding = 0.65, point.padding = 0.6, min.segment.length = 0,
      segment.colour = "#888888", segment.size = 0.25)
    p
  }
  plots <- list(RNA_PCA = make_plot("top"), RNA_PCA_labeled = make_plot("top", TRUE),
                 RNA_PCA_all_genes = make_plot("all"))
  check(identical(sample_source_md5, tools::md5sum(sample_source_files)),
        "Input quantification files or sample logs changed during analysis; rerun before export")
  out <- opt$output_dir
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  check(dir.exists(out) && file.access(out, 2) == 0L, "Output directory is not writable")
  write_tsv <- function(x, name) write.table(x, file.path(out, name), sep = "\t",
                                            quote = FALSE, row.names = FALSE, na = "NA")
  for (name in names(plots)) {
    dims <- if (name == "RNA_PCA_labeled") c(180, 130) else c(170, 120)
    grDevices::cairo_pdf(file.path(out, paste0(name, ".pdf")), width = dims[1] / 25.4,
                         height = dims[2] / 25.4, family = font, bg = "white")
    print(plots[[name]])
    grDevices::dev.off()
    grDevices::png(file.path(out, paste0(name, ".png")), width = dims[1], height = dims[2],
                   units = "mm", res = 600, type = "cairo", bg = "white")
    print(plots[[name]])
    grDevices::dev.off()
  }
  write_tsv(samples, "RNA_PCA_samples.tsv")
  write_tsv(gene_table, "RNA_PCA_gene_selection.tsv")
  for (key in names(pcs)) {
    suffix <- if (key == "top") "" else "_all_genes"
    write_tsv(scores[[key]], paste0("RNA_PCA", suffix, "_scores.tsv"))
    write_tsv(variances[[key]], paste0("RNA_PCA", suffix, "_variance.tsv"))
    for (metric in c("distances", "correlations"))
      write_tsv(data.frame(sample_id = expected, qc[[key]][[metric]], check.names = FALSE),
                 paste0("RNA_PCA", suffix, "_", metric, ".tsv"))
  }
  write_tsv(do.call(rbind, lapply(qc, function(x) x$summary)), "RNA_PCA_sample_QC.tsv")
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  script <- if (length(script_arg)) norm(sub("^--file=", "", script_arg[1])) else NA_character_
  inputs <- unique(c(unname(files), cmd_files, meta_files, index_file, gff_path, script))
  inputs <- inputs[!is.na(inputs)]
  write_tsv(data.frame(path = inputs, md5 = unname(tools::md5sum(inputs))),
              "RNA_PCA_input_checksums.tsv")
  if (!is.na(script) && norm(script) != normalizePath(file.path(out, "5.8.pca.R"),
      mustWork = FALSE)) check(file.copy(script, file.path(out, "5.8.pca.R"), overwrite = TRUE),
                               "Failed to copy standalone R script")
  saveRDS(list(samples = samples, gene_counts = counts, effective_lengths = txi$length,
                transformed_expression = mat, selected_genes = selected,
                eligible_genes = eligible, pca = pcs, parameters = opt),
            file.path(out, "RNA_PCA_analysis.rds"), compress = "gzip")
  writeLines(capture.output(sessionInfo()), file.path(out, "RNA_PCA_sessionInfo.txt"))
  writeLines(vst_messages, file.path(out, "RNA_PCA_VST_messages.txt"))
  qc_notes <- unlist(lapply(names(qc), function(key) {
    q <- qc[[key]]$summary
    ord <- order(-q$mean_within_group_distance)
    c(paste0(key, ": largest mean within-temperature distances (descriptive ranking): ",
              paste(sprintf("%s=%.3f", q$sample_id[head(ord, 3)],
                              q$mean_within_group_distance[head(ord, 3)]), collapse = "; ")),
      paste0(key, ": nearest sample from another temperature: ",
              if (all(q$nearest_sample_same_temperature)) "none" else
                paste(q$sample_id[!q$nearest_sample_same_temperature], collapse = ", ")))
  }))
  caption_main <- c(
    "图 X 不同温度处理下拟南芥 Col-0 转录组的主成分分析",
    sprintf(paste0("拟南芥 Col-0 经 4、8、12、16、20、24、28 和 32 ℃处理 3 d 后进行 RNA-seq 分析，",
      "其中 24 ℃为对照。每个点代表一个独立的生物学重复，每个温度处理包含 3 个重复，共 24 个样本；",
      "不同颜色表示不同温度。采用 tximport 直接导入 Salmon 基因水平定量文件 quant.genes.sf 中的 NumReads 估计计数，",
      "保留至少 3 个样本中取整后的估计计数≥10的基因，共 %s 个，并使用 DESeq2 进行方差稳定化转换",
      "（VST，blind = TRUE）。选取转换后方差最高的 %s 个基因，使用 R 的 prcomp() 进行 PCA，",
      "数据按基因中心化，不进行单位方差标准化。PC1 和 PC2 分别解释总变异的 %.2f%% 和 %.2f%%，",
      "累计解释 %.2f%%。图中保留全部样本，未添加置信椭圆；样本间相对位置用于描述表达模式的相似性，",
      "不作为组间差异显著性的检验。"), format(sum(keep), big.mark = ",", trim = TRUE),
      format(length(selected), big.mark = ",", trim = TRUE), variances$top$percent[1],
      variances$top$percent[2], sum(variances$top$percent[1:2])))
  caption_all <- c(
    "图 X 使用全部合格基因的拟南芥 Col-0 转录组主成分分析",
    sprintf(paste0("样本设计、计数过滤和 VST 方法同主图。本图使用全部 %s 个通过过滤且转换后方差非零的基因",
      "进行中心化但不作单位方差标准化的 PCA，用于检查结果对高变基因筛选的敏感性。",
      "每点代表一个生物学重复，每个温度 n = 3，共 24 个样本；颜色表示处理温度，24 ℃为对照。",
      "PC1 和 PC2 分别解释总变异的 %.2f%% 和 %.2f%%，累计解释 %.2f%%。"),
      format(length(eligible), big.mark = ",", trim = TRUE), variances$all$percent[1],
      variances$all$percent[2], sum(variances$all$percent[1:2])))
  captions <- c(caption_main, "", caption_all, "",
    "样本编号辅助图说明：在主图相同的 PCA 坐标上标注当前样本编号，用于检查生物学重复的一致性。",
    "编号格式为 T温度_重复编号；为避免遮挡，仅移动文字标签，样本坐标未作平移或抖动。",
    "", "排版说明：将图 X 替换为论文中的实际图号；图注作为正文文字置于图片下方。")
  writeLines(enc2utf8(captions), file.path(out, "RNA_PCA_caption.txt"), useBytes = TRUE)
  methods <- c(
    "Arabidopsis Col-0 RNA-seq PCA: methods and reproducibility record",
    "Design: 3 d; 4, 8, 12, 16, 20, 24, 28, 32 degrees Celsius; n=3 biological samples each.",
    "Biological replication and treatment duration follow the supplied experimental design.",
    "Sample identity and temperature follow current directory names.",
    "No cross-temperature pairing, batch correction, imputation or sample exclusion.",
    paste("Input directory:", norm(opt$input_dir)), paste("Output directory:", norm(out)),
    paste("Gene map recorded by Salmon (provenance only; no PCA aggregation):", recorded_gene_maps), paste("Salmon index:", index_paths),
    paste("Salmon version:", paste(unique(vapply(cmds, function(x) x$salmon_version,
                                                               character(1))), collapse = ", ")),
    paste("Genes quantified per sample:", length(ids)),
    "Input: quant.genes.sf; Name=gene ID; NumReads=estimated gene counts; EffectiveLength=gene effective length.",
    "Gene IDs are unique, matched across all samples, and verified against the reference GFF.",
    "tximport(type='salmon', txIn=FALSE, txOut=FALSE, geneIdCol='Name', countsFromAbundance='no').",
    "No transcript quant.sf input or tx2gene aggregation is used in this PCA analysis.",
    "DESeqDataSetFromTximport; design=~1. DESeq2 rounds the imported gene estimated counts.",
    "EffectiveLength values from quant.genes.sf are retained for normalization with library-size correction.",
    "Low-expression filter: rounded gene counts >=10 in at least 3 of 24 samples.",
    paste("Genes before filtering:", nrow(counts)), paste("Genes after count filtering:", sum(keep)),
    paste("Nonzero-variance genes in all-gene PCA:", length(eligible)),
    "varianceStabilizingTransformation(blind=TRUE, fitType='parametric'); unsupervised sample QC.",
    paste("Actual dispersion fit type:", fit_type),
    "Blind VST ignores temperature labels when estimating the dispersion trend; widespread treatment",
    "effects can be treated as dispersion. PCA is exploratory and is not a group-significance test.",
    "VST models count mean-variance dependence with library/length normalization. log2(TPM+1)",
    "is an abundance pseudocount transform without this count-dispersion model; it was not used.",
    paste("Top-variable genes requested/used:", opt$top_n, "/", length(selected)),
    "Rank: decreasing sample variance after VST, ties resolved by gene ID; zero variance excluded.",
    "prcomp(t(expression), center=TRUE, scale.=FALSE); samples are observations, genes are features.",
    "PC signs fixed so the largest absolute loading is positive; this does not change distances.",
    "All returned PCs are exported, including the numerically zero last PC expected after centering.",
    sprintf("Top-variable PCA: PC1=%.8f%%; PC2=%.8f%%.", variances$top$percent[1],
              variances$top$percent[2]),
    sprintf("All-eligible PCA: PC1=%.8f%%; PC2=%.8f%%.", variances$all$percent[1],
              variances$all$percent[2]),
    "QC: Pearson sample correlation and Euclidean distance in the full selected VST feature space.",
    "Within-group distances and nearest samples are descriptive, not formal outlier tests.", qc_notes,
    "No ellipse, inferential P value, or automatic outlier removal; review both PCA feature sets.",
    paste("Font:", font, "(Arial preferred; compatible sans-serif fallback when absent)."),
    "Main/all-gene plots: 170 x 120 mm; labeled: 180 x 130 mm; Cairo vector PDF and 600 dpi PNG.",
    "Points: size=3 mm, alpha=0.85; coordinates are not jittered; Cartesian axes fit the page layout.",
    "Legend entries contain temperature only; the control condition is documented in the caption.",
    "Text labels use ggrepel seed=20261010, max.iter=100000, max.time=Inf; points remain fixed.",
    "T8_1 text starts 8% of the PC1 range to its left to prevent a label/leader crossing.",
    paste("Palette in ascending temperature order:", paste(palette, collapse = ", ")),
    "Single-panel figures: multi-panel alignment audit is not applicable.",
    "",
    "论文方法：对24个拟南芥Col-0样本，使用tximport直接导入Salmon基因水平定量文件quant.genes.sf，",
    "设置txIn=FALSE、geneIdCol='Name'，使用NumReads估计计数，不再汇总转录本。保留至少3个样本中计数经取整后",
    "达到10的基因，使用DESeq2进行blind=TRUE的方差稳定化转换（VST）。按转换后的基因",
    sprintf("方差选择前%d个基因，以样本为观察对象，使用prcomp进行中心化且不标准化的PCA；", length(selected)),
    "同时使用全部通过过滤且方差非零的基因进行对照分析。未删除任何生物学重复。",
    "",
    caption_main,
    "主图、全部合格基因对照图和样本编号辅助图的图注见 RNA_PCA_caption.txt。",
    "",
    "Re-run: activate the existing conda environment named R, then Rscript 5.8.pca.R.",
    "Use --input-dir, --output-dir, --top-n to override defaults; --validate-only writes nothing.",
    "Versions: see RNA_PCA_sessionInfo.txt; input MD5: RNA_PCA_input_checksums.tsv."
  )
  writeLines(enc2utf8(methods), file.path(out, "RNA_PCA_methods.txt"), useBytes = TRUE)
  number_zh <- function(n) format(n, big.mark = " ", scientific = FALSE, trim = TRUE)
  method_zh <- c(
    "# RNA-seq 基因表达数据预处理与主成分分析", "", "## 论文方法正文", "",
    sprintf(paste0("以经 4、8、12、16、20、24、28 和 32 ℃处理 3 d 的拟南芥 Col-0 为研究材料，",
      "其中 24 ℃为对照，每个温度处理设置 3 个生物学重复，共 24 个 RNA-seq 样本。采用 Salmon（v%s）",
      "进行表达定量，使用其输出的基因水平定量文件 quant.genes.sf 开展后续分析。在 R（v%s）中，",
      "利用 tximport 包（v%s）直接导入各样本的基因水平估计计数（estimated counts，NumReads），",
      "设置 type = \"salmon\"、txIn = FALSE、txOut = FALSE、geneIdCol = \"Name\" 和 countsFromAbundance = \"no\"。",
      "以 Name 列中的基因 ID 对齐全部样本，直接读取每个基因的计数和 EffectiveLength；",
      "本次 PCA 不使用转录本水平的 quant.sf，也不在 R 中再次进行转录本至基因的汇总。",
      "使用 DESeq2 包（v%s）构建表达数据对象，并在导入时对基因估计计数取整。",
      "保留至少 3 个样本中取整后计数≥10 的基因；经过滤，%s 个基因中共有 %s 个用于后续分析。"),
      cmds[[1]]$salmon_version, as.character(getRversion()), as.character(packageVersion("tximport")),
      as.character(packageVersion("DESeq2")), number_zh(nrow(counts)), number_zh(sum(keep))), "",
    sprintf(paste0("利用 DESeq2 对过滤后的计数进行标准化，在校正测序文库规模差异的同时，",
      "保留 quant.genes.sf 中的基因有效长度信息以校正样本间长度差异。随后调用 ",
      "varianceStabilizingTransformation() 进行方差稳定化转换（variance stabilizing transformation，VST），",
      "以减弱表达均值与方差之间的依赖。为开展探索性样本质量评估，采用不包含温度分组信息的设计",
      "（design = ~1），设置 blind = TRUE，实际离散度趋势拟合方式为 %s。根据每个基因在全部 24 个样本中的",
      "VST 转换后表达方差进行降序排列，选择方差最高的 %s 个基因，以样本为观察对象、基因为特征，",
      "调用 R 的 prcomp() 函数进行主成分分析（principal component analysis，PCA），",
      "设置 center = TRUE、scale. = FALSE，即按基因中心化而不进行单位方差缩放。",
      "同时使用全部 %s 个通过过滤且转换后方差非零的基因进行对照 PCA，以检查结果对高变基因筛选的敏感性。"),
      if (fit_type == "parametric") "参数化拟合（fitType = \"parametric\"）" else
        paste0("局部回归拟合（fitType = \"", fit_type, "\"）"),
      number_zh(length(selected)), number_zh(length(eligible))), "",
    sprintf(paste0("采用 ggplot2 包（v%s）绘制 PC1 与 PC2 散点图，在坐标轴标注相应主成分的解释方差百分比。",
      "不同颜色表示不同温度处理，每个点代表一个生物学重复，以评估不同温度处理下样本整体表达模式及生物学重复的一致性。",
      "基于 VST 转换后的表达矩阵计算样本间 Pearson 相关系数和欧氏距离，辅助检查样本一致性。",
      "本分析保留全部 24 个样本，不进行组内均值替代，不添加置信椭圆，也不依据 PCA 分布自动删除样本。",
      "图片导出为矢量 PDF 和 600 dpi PNG。"), as.character(packageVersion("ggplot2"))), "",
    "## 实际参数核对表", "",
    "| 项目 | 本次实际采用的设置 |", "| --- | --- |",
    "| 输入文件 | 24 个 Salmon quant.genes.sf；每行是一个基因 |",
    "| 输入数值 | NumReads 基因 estimated counts；不以 TPM 或 log2(TPM + 1) 作为 PCA 输入 |",
    "| 基因导入 | tximport(type = \"salmon\", txIn = FALSE, txOut = FALSE, geneIdCol = \"Name\", countsFromAbundance = \"no\")；不再次汇总转录本 |",
    "| 计数处理 | DESeqDataSetFromTximport()；DESeq2 对导入的估计计数取整 |",
    "| 低表达过滤 | 至少 3/24 个样本的取整后基因计数≥10；在标准化和 VST 前执行 |",
    sprintf("| 基因数量 | 过滤前 %s 个；过滤后 %s 个；其中非零方差基因 %s 个 |",
            number_zh(nrow(counts)), number_zh(sum(keep)), number_zh(length(eligible))),
    "| 标准化 | DESeq2 校正文库规模差异，并保留 quant.genes.sf 的 EffectiveLength 信息构建基因及样本特异的标准化因子 |",
    sprintf("| 方差稳定化转换（VST） | varianceStabilizingTransformation()；design = ~1，blind = TRUE，实际 fitType = \"%s\" |", fit_type),
    sprintf("| 主分析（PCA） | VST 后方差最高的 %s 个基因；方差并列时按基因 ID 排序 |", number_zh(length(selected))),
    sprintf("| 对照分析（PCA） | 全部 %s 个通过过滤且转换后方差非零的基因 |", number_zh(length(eligible))),
    "| PCA 设置 | prcomp(t(expression_matrix), center = TRUE, scale. = FALSE) |",
    "| 样本处理 | 保留全部 24 个样本，采用当前目录编号；不进行组内取均值、缺失值填补、批次校正或自动剔除 |",
    "", "## 写作与解释说明", "",
    "quant.genes.sf 是本次 PCA 的直接输入；quant.sf 为转录本水平定量文件。本次流程从基因水平文件直接读取计数，不需要再以 tx2gene 汇总。样本对应关系、基因 ID 和参考索引已核验。",
    "DESeq2 标准化与 PCA 阶段的单位方差缩放是不同的处理步骤。低表达过滤在计数数据上执行，高变基因筛选在 VST 转换后的数据上执行；本次采用 VST，不采用 log2(TPM + 1)。",
    "PCA、样本相关性和距离用于探索性描述和质量评估；图中分离或接近不能单独作为组间差异显著性的证据。潜在异常样本需结合实验记录和其他质量控制信息复核。",
    "", "本文档由分析脚本按实际输入、保留基因数和软件版本生成；软件环境详见 RNA_PCA_sessionInfo.txt。")
  writeLines(enc2utf8(method_zh), file.path(out, "RNA_PCA_methods_zh.md"), useBytes = TRUE)
  report_lines <- startsWith(methods, "Genes ") | grepl("PCA: PC1", methods, fixed = TRUE) |
    grepl("fit type", methods, fixed = TRUE)
  message(paste(methods[report_lines], collapse = "\n"))
  message("Output complete: ", norm(out))
}
main()
