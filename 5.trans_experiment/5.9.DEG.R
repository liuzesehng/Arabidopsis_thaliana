#!/usr/bin/env Rscript
# Arabidopsis RNA-seq: seven temperature contrasts against 24 degrees Celsius.
# Run: Rscript 5.9.DEG.R [--input-dir DIR] [--output-dir DIR] [--validate-only]
# Inputs are read-only. Three independent biological replicates per temperature.
# Figure contract: describe expression changes relative to 24 C and the selected
# genes' sample expression patterns; no assumed temperature separation or mechanism.
# Independent single-panel figures; R is the exclusive rendering backend.
# DejaVu Sans is the installed sans-serif fallback.

options(stringsAsFactors = FALSE, warn = 1)
check <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
read_tsv <- function(path) read.delim(path, check.names = FALSE, quote = "",
                                     comment.char = "", stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE,
                                          row.names = FALSE, na = "NA")
write_matrix <- function(x, path) write_tsv(data.frame(gene_id = rownames(x), x,
                                                      check.names = FALSE), path)
write_json <- function(x, path) jsonlite::write_json(x, path, auto_unbox = TRUE,
                                                   pretty = TRUE, na = "null", digits = 16)
sha256 <- function(paths) vapply(paths, function(p) digest::digest(file = p,
                                algo = "sha256", serialize = FALSE), character(1))
classify_de <- function(pvalue, padj, lfc, alpha = 0.05, threshold = 1) {
  state <- rep("Not_DEG", length(lfc))
  state[is.na(padj)] <- "Adjusted_p_unavailable"
  state[is.na(pvalue)] <- "P_value_unavailable"
  state[!is.finite(lfc)] <- "Effect_unavailable"
  eligible <- is.finite(padj) & is.finite(lfc) & padj < alpha
  state[eligible & lfc >= threshold] <- "Up"
  state[eligible & lfc <= -threshold] <- "Down"
  state
}
rank_de <- function(x) x[order(x$padj, -abs(x$log2FoldChange), x$gene_id,
                                na.last = TRUE, method = "radix"), , drop = FALSE]

parse_options <- function(args) {
  root <- "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana"
  opt <- list(input_dir = file.path(root, "list/rnadata"),
              output_dir = file.path(root, "list/RCA/DEG"), validate_only = FALSE)
  i <- 1L
  while (i <= length(args)) {
    k <- args[i]
    if (k == "--help") {
      cat("Usage: Rscript 5.9.DEG.R [--input-dir DIR] [--output-dir DIR] [--validate-only]\n",
          "Defaults: 24 C control; padj < 0.05 and abs(log2FC) >= 1; top 10 DEGs per contrast for heatmap.\n")
      return(NULL)
    }
    if (k == "--validate-only") opt$validate_only <- TRUE else {
      check(k %in% c("--input-dir", "--output-dir") && i < length(args),
            paste("Unknown option or missing value:", k))
      i <- i + 1L
      opt[[gsub("-", "_", substring(k, 3L))]] <- args[i]
    }
    i <- i + 1L
  }
  opt
}

validate_inputs <- function(input_dir) {
  expected <- c("T4_2", "T4_4", "T4_5", "T8_1", "T8_2", "T8_4",
                "T12_2", "T12_3", "T12_6", "T16_3", "T16_5", "T16_6",
                "T20_1", "T20_2", "T20_4", "T24_1", "T24_2", "T24_4",
                "T28_1", "T28_4", "T28_6", "T32_1", "T32_3", "T32_4")
  files <- list.files(input_dir, pattern = "^quant\\.genes\\.sf$", recursive = TRUE,
                      full.names = TRUE)
  check(length(files) == length(expected),
        paste("Expected 24 quant.genes.sf files; found", length(files)))
  check(all(basename(dirname(files)) == "salmon"), "Expected sample/salmon/quant.genes.sf")
  sample_ids <- basename(dirname(dirname(files)))
  check(!anyDuplicated(sample_ids) && setequal(sample_ids, expected),
        "Sample IDs differ from the verified 24-sample design")
  files <- files[match(expected, sample_ids)]
  names(files) <- expected
  samples <- data.frame(sample_id = expected,
    temperature_C = as.integer(sub("^T([0-9]+)_.*", "\\1", expected)),
    replicate_id = as.integer(sub("^T[0-9]+_", "", expected)),
    quant_file = normalizePath(files, mustWork = TRUE))
  samples$is_control <- samples$temperature_C == 24L
  check(all(table(samples$temperature_C) == 3L), "Each temperature must have three samples")
  cmd_files <- file.path(dirname(files), "cmd_info.json")
  meta_files <- file.path(dirname(files), "aux_info/meta_info.json")
  source_files <- c(unname(files), cmd_files, meta_files)
  check(all(file.exists(source_files)), "Missing quantification provenance files")
  input_hashes <- sha256(source_files)
  cmds <- lapply(cmd_files, jsonlite::fromJSON)
  metas <- lapply(meta_files, jsonlite::fromJSON)
  for (j in seq_along(files)) {
    for (key in c("mates1", "mates2", "output")) {
      value <- cmds[[j]][[key]]
      check(length(value) >= 1L && all(grepl(paste0("(^|/)", expected[j], "(/|$)"), value)),
            paste("Salmon sample identity mismatch:", expected[j], key))
    }
    check(length(metas[[j]]$quant_errors) == 0L, paste("Salmon quantification error:", expected[j]))
  }
  indices <- unique(vapply(cmds, function(x) normalizePath(x$index, mustWork = TRUE), character(1)))
  maps <- unique(vapply(cmds, function(x) x$geneMap, character(1)))
  check(length(indices) == 1L && length(maps) == 1L, "Inconsistent reference indexes or gene maps")
  index_file <- file.path(indices, "info.json")
  gff_path <- file.path(dirname(indices), "refgen/GCF_000001735.4_TAIR10.1_genomic.gff")
  check(all(file.exists(c(index_file, gff_path))), "Missing reference index metadata or GFF")
  source_files <- c(source_files, index_file, gff_path)
  input_hashes <- c(input_hashes, sha256(c(index_file, gff_path)))
  idx <- jsonlite::fromJSON(index_file)
  for (m in metas) check(identical(m$index_seq_hash, idx$SeqHash) &&
    identical(m$index_name_hash, idx$NameHash), "Salmon reference hash mismatch")
  gff <- read.delim(gff_path, header = FALSE, comment.char = "#", quote = "",
                    colClasses = "character")
  attrs <- gff[[9]][grepl("(^|;)ID=gene-", gff[[9]])]
  annotated <- unique(sub(".*(?:^|;)ID=gene-([^;]+).*", "\\1", attrs, perl = TRUE))
  quant <- lapply(files, read_tsv)
  ids <- sort(quant[[1]]$Name, method = "radix")
  check(length(ids) > 0L && !anyDuplicated(ids) && all(grepl("^AT[1-5CM]G[0-9]+$", ids)) &&
          all(ids %in% annotated), "Invalid or unannotated Arabidopsis gene identifiers")
  for (j in seq_along(quant)) {
    q <- quant[[j]]
    check(identical(names(q), c("Name", "Length", "EffectiveLength", "TPM", "NumReads")),
          paste("Invalid quantification schema:", expected[j]))
    check(!anyDuplicated(q$Name) && setequal(q$Name, ids), paste("Gene set mismatch:", expected[j]))
    check(all(vapply(q[-1], is.numeric, logical(1))) && all(is.finite(as.matrix(q[-1]))) &&
      all(as.matrix(q[-1]) >= 0) && all(q$Length > 0) && all(q$EffectiveLength > 0) &&
      sum(q$NumReads) > 0, paste("Invalid quantification values:", expected[j]))
    quant[[j]] <- q[match(ids, q$Name), , drop = FALSE]
  }
  samples$percent_mapped <- vapply(metas, function(m) m$percent_mapped, numeric(1))
  samples$estimated_fragments <- vapply(quant, function(q) sum(q$NumReads), numeric(1))
  message("Input validation passed: ", length(files), " samples; ", length(ids), " genes per sample.")
  list(files = files, samples = samples, ids = ids, quant = quant,
       source_files = source_files, hashes_before = input_hashes,
       salmon_version = unique(vapply(metas, function(m) m$salmon_version, character(1))))
}

# Every finite gene is drawn; the PDFs retain vector points, axes and text.
export_figures <- function(results, summary, z, samples, out, font) {
  suppressPackageStartupMessages(library(ggplot2))
  figures <- file.path(out, "figures")
  sources <- file.path(out, "source_data")
  colors <- c(Down = "#0072B2", Not_DEG = "#B7B7B7", Up = "#D55E00", Unavailable = "#666666")
  state_levels <- c("Not_DEG", "Unavailable", "Down", "Up")
  base_theme <- theme_classic(base_size = 9, base_family = "DejaVu Sans") +
    theme(axis.text = element_text(size = 8, colour = "black"),
      axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.3),
      legend.position = "bottom", legend.title = element_blank(),
      legend.text = element_text(size = 7.5), legend.key.width = grid::unit(4, "mm"),
      legend.key.height = grid::unit(3.5, "mm"), legend.margin = margin(0, 0, 0, 0),
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8), plot.caption = element_text(size = 6.5, hjust = 0),
      plot.margin = margin(8, 10, 8, 8), plot.background = element_rect(fill = "white", colour = NA))
  save_plot <- function(draw, name, width_mm, height_mm) {
    prefix <- file.path(figures, name)
    grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = width_mm / 25.4,
                         height = height_mm / 25.4, family = font, bg = "white")
    tryCatch(draw(), finally = grDevices::dev.off())
    grDevices::png(paste0(prefix, ".png"), width = width_mm, height = height_mm,
                   units = "mm", res = 600, type = "cairo", bg = "white")
    tryCatch(draw(), finally = grDevices::dev.off())
    message("Exported ", name)
  }
  for (comparison in names(results)) {
    d <- results[[comparison]]
    d$plot_group <- ifelse(d$direction %in% c("Up", "Down", "Not_DEG"), d$direction, "Unavailable")
    d$plot_group <- factor(d$plot_group, levels = state_levels)
    d <- d[order(d$plot_group, d$gene_id), ]
    d$volcano_included <- is.finite(d$padj) & is.finite(d$log2FoldChange)
    d$neglog10_padj <- NA_real_
    d$neglog10_padj[d$volcano_included] <- -log10(pmax(d$padj[d$volcano_included], .Machine$double.xmin))
    source_name <- paste0(comparison, ".plot_data.tsv")
    write_tsv(d, file.path(sources, source_name))
    s <- summary[summary$comparison == comparison, ]
    title <- paste0(s$treatment_C, " °C vs 24 °C")
    subtitle <- sprintf("Up: %s   Down: %s", format(s$up, big.mark = ","), format(s$down, big.mark = ","))
    caption <- "DEG: BH-adjusted P < 0.05 and abs(log2FC) ≥ 1; n = 3 per temperature."
    v <- d[d$volcano_included, ]
    x_extent <- max(1.5, abs(v$log2FoldChange)) * 1.06
    p <- ggplot(v, aes(log2FoldChange, neglog10_padj, colour = plot_group)) +
      geom_point(size = 0.55, alpha = 0.7, stroke = 0) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3, colour = "#555555") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, colour = "#555555") +
      scale_colour_manual(values = colors, breaks = c("Down", "Not_DEG", "Up"),
                          labels = c("Down", "Not DEG", "Up"), drop = FALSE) +
      scale_x_continuous(limits = c(-x_extent, x_extent), expand = expansion(mult = 0)) +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
      labs(title = title, subtitle = subtitle, x = "log2 fold change", y = "−log10 adjusted P",
           caption = paste0(caption, "\nUnavailable adjusted P/effect: ", sum(!d$volcano_included), "; omitted from coordinates.")) + base_theme
    save_plot(function() print(p), paste0(comparison, ".volcano"), 125, 110)
  }
  bars <- rbind(data.frame(temperature_C = summary$treatment_C, direction = "Up", count = summary$up),
                data.frame(temperature_C = summary$treatment_C, direction = "Down", count = summary$down))
  bars$temperature <- factor(paste0(bars$temperature_C, " °C"), levels = paste0(summary$treatment_C, " °C"))
  bars$direction <- factor(bars$direction, levels = c("Up", "Down"))
  write_tsv(bars, file.path(sources, "DEG_counts.plot_data.tsv"))
  p <- ggplot(bars, aes(temperature, count, fill = direction)) +
    geom_col(position = position_dodge(width = 0.76), width = 0.68) +
    geom_text(aes(label = format(count, big.mark = ",", trim = TRUE)),
              position = position_dodge(width = 0.76), vjust = -0.5, size = 2.7, family = font) +
    scale_fill_manual(values = colors[c("Up", "Down")], drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(x = "Treatment temperature", y = "Number of genes",
      caption = "BH-adjusted P < 0.05 and abs(log2FC) ≥ 1; n = 3 independent biological samples per temperature.\nBars count classified genes; they are not replicate means and have no error bars.") + base_theme
  save_plot(function() print(p), "DEG_counts", 183, 115)
  if (nrow(z)) {
    temp_colors <- setNames(c("#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377", "#BBBBBB", "#000000"), seq(4, 32, 4))
    limit <- max(abs(z))
    row_cluster <- if (nrow(z) > 1L) hclust(dist(z, method = "euclidean"), method = "complete") else FALSE
    row_order <- if (isFALSE(row_cluster)) 1L else row_cluster$order
    write_matrix(z[row_order, , drop = FALSE], file.path(sources, "heatmap_zscore_display_order.tsv"))
    ha <- ComplexHeatmap::HeatmapAnnotation(Temperature = as.character(samples$temperature_C),
      col = list(Temperature = temp_colors), show_legend = FALSE, show_annotation_name = FALSE,
      simple_anno_size = grid::unit(3, "mm"))
    make_heatmap <- function() ComplexHeatmap::Heatmap(z, name = "Row Z-score",
      col = circlize::colorRamp2(c(-limit, 0, limit), c("#2166AC", "#F7F7F7", "#B2182B")),
      cluster_rows = row_cluster, cluster_columns = FALSE, row_dend_reorder = FALSE,
      show_row_names = TRUE, show_column_names = TRUE, top_annotation = ha,
      column_labels = samples$sample_id, column_names_rot = 45,
      column_names_gp = grid::gpar(fontsize = 7.5, fontfamily = font),
      row_names_gp = grid::gpar(fontsize = 7, fontfamily = font),
      row_names_max_width = grid::unit(28, "mm"), row_dend_width = grid::unit(10, "mm"),
      column_title = NULL,
      rect_gp = grid::gpar(col = NA), use_raster = FALSE,
      heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 8, fontfamily = font),
        labels_gp = grid::gpar(fontsize = 7, fontfamily = font),
        direction = "vertical", title_position = "topleft", legend_height = grid::unit(35, "mm")))
    draw_ht <- function() {
      previous_legend_padding <- ComplexHeatmap::ht_opt("HEATMAP_LEGEND_PADDING")
      on.exit(ComplexHeatmap::ht_opt(HEATMAP_LEGEND_PADDING = previous_legend_padding), add = TRUE)
      ComplexHeatmap::ht_opt(HEATMAP_LEGEND_PADDING = grid::unit(8, "mm"))
      grid::pushViewport(grid::viewport(gp = grid::gpar(fontfamily = font)))
      ht <- make_heatmap()
      ComplexHeatmap::draw(ht, newpage = FALSE, heatmap_legend_side = "right",
                          align_heatmap_legend = "heatmap_top",
                          padding = grid::unit(c(15, 4, 5, 4), "mm"))
      grid::grid.text(sprintf("Top DEGs across temperatures (%d genes)\nTop 10 DEGs per contrast, union; all 24 biological samples.\nVST (blind = FALSE), row Z-score; rows: Euclidean / complete linkage.", nrow(z)),
        x = grid::unit(4, "mm"), y = grid::unit(4, "mm"), just = c("left", "bottom"),
        gp = grid::gpar(fontsize = 7, fontfamily = font))
      grid::popViewport()
    }
    save_plot(draw_ht, "DEG_heatmap", 183, max(145, 70 + 2.7 * nrow(z)))
  } else {
    writeLines("No eligible significant genes for heatmap; no expression values were fabricated.",
               file.path(figures, "DEG_heatmap.NOT_GENERATED.txt"))
  }
  invisible(NULL)
}

main <- function() {
  opt <- parse_options(commandArgs(trailingOnly = TRUE))
  if (is.null(opt)) return(invisible(NULL))
  required <- c("tximport", "DESeq2", "jsonlite", "digest")
  if (!opt$validate_only) required <- c(required, "ggplot2", "ComplexHeatmap", "circlize")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  check(!length(missing), paste("Missing R packages:", paste(missing, collapse = ", ")))
  input <- validate_inputs(opt$input_dir)
  if (opt$validate_only) return(invisible(input$samples))
  check(capabilities("cairo"), "R Cairo support is required for figure exports")
  out <- opt$output_dir
  check(!normalizePath(out, mustWork = FALSE) %in% normalizePath(c(opt$input_dir, dirname(input$files)), mustWork = TRUE),
        "Output directory must differ from input directories")
  for (d in c(out, file.path(out, c("tables", "figures", "source_data"))))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  check(dir.exists(out) && file.access(out, 2) == 0L, "Output directory is not writable")
  started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  write_json(list(status = "RUNNING", started = started), file.path(out, "run_status.json"))
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  script_path <- if (length(script_arg)) normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE) else ""
  importer <- function(f) {
    j <- match(normalizePath(f, mustWork = TRUE), normalizePath(input$files, mustWork = TRUE))
    check(!is.na(j), "Importer received an unknown file")
    input$quant[[j]]
  }
  txi <- tximport::tximport(input$files, type = "salmon", txIn = FALSE, txOut = FALSE,
    geneIdCol = "Name", countsFromAbundance = "no", importer = importer)
  check(identical(rownames(txi$counts), input$ids) &&
        identical(colnames(txi$counts), input$samples$sample_id), "Import changed gene or sample order")
  check(all(txi$counts == do.call(cbind, lapply(input$quant, function(q) q$NumReads))) &&
        all(txi$length == do.call(cbind, lapply(input$quant, function(q) q$EffectiveLength))),
        "Imported counts or lengths differ from the source")
  coldata <- data.frame(row.names = input$samples$sample_id,
    temperature = relevel(factor(as.character(input$samples$temperature_C),
                               levels = as.character(seq(4, 32, 4))), ref = "24"))
  dds_all <- DESeq2::DESeqDataSetFromTximport(txi, colData = coldata, design = ~temperature)
  rounded <- DESeq2::counts(dds_all)
  keep <- rowSums(rounded >= 10L) >= 3L
  check(sum(keep) >= 2L, "Too few genes pass the low-count filter")
  filter_table <- data.frame(gene_id = rownames(rounded), samples_count_ge10 = rowSums(rounded >= 10L),
    passes_count_filter = keep, reason = ifelse(keep, "Retained", "Count_below_10_in_at_least_22_samples"))
  write_tsv(filter_table, file.path(out, "tables/gene_filter.tsv"))
  write_tsv(input$samples, file.path(out, "tables/samples.tsv"))
  write_matrix(txi$counts, file.path(out, "source_data/estimated_counts_all_genes.tsv"))
  write_matrix(txi$length, file.path(out, "source_data/effective_lengths_all_genes.tsv"))
  message("Fitting DESeq2: ", sum(keep), " retained genes; all 24 samples; control = 24 C.")
  dds <- DESeq2::DESeq(dds_all[keep, ], test = "Wald", fitType = "parametric",
    sfType = "ratio", betaPrior = FALSE, minReplicatesForReplace = Inf, parallel = FALSE)
  check(identical(DESeq2::counts(dds), rounded[keep, , drop = FALSE]), "Fitting modified observed counts")
  normalization <- DESeq2::normalizationFactors(dds)
  check(all(is.finite(normalization)) && all(normalization > 0), "Invalid normalization factors")
  normalized <- DESeq2::counts(dds, normalized = TRUE)
  write_matrix(DESeq2::counts(dds), file.path(out, "tables/rounded_counts_tested_genes.tsv"))
  write_matrix(normalized, file.path(out, "tables/normalized_counts.tsv"))
  write_matrix(normalization, file.path(out, "source_data/normalization_factors.tsv"))
  model_info <- as.data.frame(S4Vectors::mcols(dds))
  write_tsv(data.frame(gene_id = rownames(dds), model_info, check.names = FALSE),
            file.path(out, "tables/gene_model_diagnostics.tsv"))
  write_matrix(SummarizedExperiment::assay(dds, "cooks"), file.path(out, "source_data/cooks_distances.tsv"))
  p <- ncol(model.matrix(~temperature, coldata))
  cooks_cutoff <- qf(0.99, p, ncol(dds) - p)
  contrasts <- c(4L, 8L, 12L, 16L, 20L, 28L, 32L)
  results <- list(); summaries <- list(); selections <- list(); comparison_metadata <- list()
  for (temp in contrasts) {
    name <- paste0("T", temp, "_vs_T24")
    res <- DESeq2::results(dds, contrast = c("temperature", as.character(temp), "24"),
      alpha = 0.05, pAdjustMethod = "BH", independentFiltering = TRUE,
      cooksCutoff = cooks_cutoff, lfcThreshold = 0, altHypothesis = "greaterAbs", parallel = FALSE)
    d <- data.frame(gene_id = rownames(res), as.data.frame(res), check.names = FALSE)
    d$direction <- classify_de(d$pvalue, d$padj, d$log2FoldChange)
    d$beta_converged <- model_info$betaConv
    d$cooks_exceeds_cutoff <- !is.na(model_info$maxCooks) & model_info$maxCooks > cooks_cutoff
    d$treatment_mean_normalized <- rowMeans(normalized[, input$samples$temperature_C == temp, drop = FALSE])
    d$control_mean_normalized <- rowMeans(normalized[, input$samples$temperature_C == 24L, drop = FALSE])
    d <- rank_de(d)
    significant <- rank_de(d[d$direction %in% c("Up", "Down"), , drop = FALSE])
    for (which in c("all", "DEG", "up", "down")) {
      table <- switch(which, all = d, DEG = significant,
        up = significant[significant$direction == "Up", , drop = FALSE],
        down = significant[significant$direction == "Down", , drop = FALSE])
      write_tsv(table, file.path(out, "tables", paste0(name, ".", which, ".tsv")))
    }
    results[[name]] <- d
    summaries[[name]] <- data.frame(comparison = name, treatment_C = temp, control_C = 24L,
      n_treatment = 3L, n_control = 3L, input_genes = length(input$ids), retained_genes = nrow(d),
      genes_with_pvalue = sum(is.finite(d$pvalue)), genes_with_padj = sum(is.finite(d$padj)),
      up = sum(d$direction == "Up"), down = sum(d$direction == "Down"), total_DEG = nrow(significant),
      tested_not_DEG = sum(d$direction == "Not_DEG"), unavailable = sum(!d$direction %in% c("Up", "Down", "Not_DEG")))
    selected <- head(significant, 10L)
    selections[[name]] <- data.frame(comparison = rep(name, nrow(selected)), rank = seq_len(nrow(selected)),
      gene_id = selected$gene_id, padj = selected$padj, log2FoldChange = selected$log2FoldChange,
      direction = selected$direction)
    meta <- S4Vectors::metadata(res)
    comparison_metadata[[name]] <- list(contrast = c("temperature", as.character(temp), "24"),
      alpha = 0.05, lfc_null = 0, BH_family = "Genes passing independent filtering within this contrast",
      filter_threshold = unname(meta$filterThreshold), filter_quantile = unname(meta$filterTheta),
      cooks_cutoff = cooks_cutoff)
    message(name, ": Up=", summaries[[name]]$up, "; Down=", summaries[[name]]$down)
  }
  summary <- do.call(rbind, summaries); rownames(summary) <- NULL
  selection <- do.call(rbind, selections); rownames(selection) <- NULL
  write_tsv(summary, file.path(out, "DEG_summary.tsv"))
  write_tsv(selection, file.path(out, "source_data/heatmap_selection.tsv"))
  vst <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  mat <- SummarizedExperiment::assay(vst)
  check(all(is.finite(mat)), "Non-finite VST expression values")
  write_matrix(mat, file.path(out, "tables/vst_expression.tsv"))
  ids <- sort(unique(selection$gene_id), method = "radix")
  chosen <- mat[ids, , drop = FALSE]
  sds <- if (length(ids)) apply(chosen, 1, sd) else numeric()
  valid <- is.finite(sds) & sds > 0
  write_tsv(data.frame(gene_id = ids, row_sd = sds, included = valid),
            file.path(out, "source_data/heatmap_row_checks.tsv"))
  chosen <- chosen[valid, , drop = FALSE]
  z <- if (nrow(chosen)) sweep(sweep(chosen, 1, rowMeans(chosen), "-"), 1, sds[valid], "/") else chosen
  write_matrix(chosen, file.path(out, "source_data/heatmap_vst.tsv"))
  write_matrix(z, file.path(out, "source_data/heatmap_zscore.tsv"))
  font_cache <- file.path(tempdir(), "font-cache")
  dir.create(font_cache, showWarnings = FALSE)
  Sys.setenv(XDG_CACHE_HOME = font_cache)
  font <- "DejaVu Sans"
  export_figures(results, summary, z, input$samples, out, font)
  actual_fit <- attr(DESeq2::dispersionFunction(dds), "fitType")
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  versions <- vapply(required, function(p) as.character(packageVersion(p)), character(1))
  config <- list(started = started, input_dir = opt$input_dir, output_dir = opt$output_dir,
    script_path = script_path, script_sha256 = if (nzchar(script_path)) unname(sha256(script_path)) else NA_character_,
    design = "~ temperature", reference = "24", biological_replicates_per_temperature = 3L,
    temperatures = seq(4, 32, 4), contrasts = contrasts, count_filter = "rounded count >= 10 in >= 3 samples",
    import = list(type = "salmon", txIn = FALSE, txOut = FALSE, geneIdCol = "Name", countsFromAbundance = "no"),
    test = "two-sided Wald; H0: log2FC = 0", alpha = 0.05, lfc_threshold_post_filter = 1,
    lfc_shrinkage = FALSE, automatic_count_replacement = FALSE, dispersion_fit = actual_fit,
    comparison_metadata = comparison_metadata, R_version = R.version.string, packages = as.list(versions),
    heatmap = list(top_per_comparison = 10L, unique_selected = length(ids), plotted = nrow(z),
      transform = "VST blind=FALSE; row Z-score using sample SD", clustering = "Euclidean distance; complete linkage",
      column_order = input$samples$sample_id), plotting_backend = "R", font = font)
  write_json(config, file.path(out, "analysis_parameters.json"))
  methods <- c("# 温度处理相对 24°C 的差异表达分析", "",
    sprintf("使用 %s、tximport %s 和 DESeq2 %s 分析拟南芥 RNA-seq 数据。", R.version.string, versions['tximport'], versions['DESeq2']),
    "依据已有实验设计记录，4、8、12、16、20、24、28、32°C 各含 3 个独立生物学重复，共 24 个样本。24°C 为对照。样本分组按当前目录名确定，不进行跨温度配对、批次校正或样本剔除；重复编号不作为配对或批次变量。", "",
    "直接读取各样本 Salmon quant.genes.sf 的基因水平 NumReads 和 EffectiveLength，按 Name 对齐基因。tximport 设置 type='salmon'、txIn=FALSE、txOut=FALSE、geneIdCol='Name'、countsFromAbundance='no'，不重复进行转录本至基因汇总，也不对 TPM 进行差异检验。DESeqDataSetFromTximport 对估计计数取整，并保留有效长度信息用于基因和样本特异的标准化因子。",
    sprintf("预过滤保留至少 3 个样本中取整计数≥10 的基因：共 %d 个输入基因，%d 个通过过滤。", length(input$ids), sum(keep)), "",
    sprintf("使用全部 24 个样本统一拟合负二项广义线性模型 design=~temperature；温度作为分类变量，24°C 为参考水平。离散度拟合实际使用 %s。分别提取 4、8、12、16、20、28、32°C 相对 24°C 的对比，进行双侧 Wald 检验，原假设为 log2FC=0。", actual_fit),
    "每个比较分别启用独立过滤（alpha=0.05）和 Benjamini–Hochberg 校正；校正家族为该比较中进入多重检验校正的基因，不是七个比较联合的全局 FDR。显著 DEG 定义为 padj<0.05 且未经收缩的 abs(log2FC)≥1。倍数门槛是检验之后的筛选规则，不代表检验原假设为 abs(log2FC)≤1。正 log2FC 为处理组相对对照上调，负值为下调。",
    "保留默认含义的 Cook 距离检查（F 分布 99% 分位数阈值），关闭自动异常计数替换；未删除样本。缺失 P 值或 padj 保持 NA，并单独标记，不能解释为无差异。模型诊断和各比较的独立过滤阈值另行保存。", "",
    "火山图展示未经收缩的 log2FC 与 −log10(padj)；零 padj 仅在作图坐标计算中采用双精度最小正规正数下限，源表不改写。无法计算坐标的基因数记录于图注和绘图源表。每个散点代表一个基因，全部有限坐标均绘制，未抽样；散点、文字和坐标轴保留为 PDF 矢量。上下调数量图统计基因数，不是生物学重复均值，因此不绘制误差条。", "",
    sprintf("各比较在显著 DEG 内依次按 padj 升序、abs(log2FC)降序及基因 ID 升序排序，取前 10 个后合并去重。选中 %d 个基因，最终热图展示 %d 个有限且行标准差非零的基因。", length(ids), nrow(z)),
    "热图使用 varianceStabilizingTransformation(blind=FALSE)，再对每个基因跨 24 个样本按样本标准差进行 Z-score 标准化。列按温度及重复编号排序，不对列聚类；行使用欧氏距离和 complete linkage 聚类。热图为筛选后的描述性表达模式，不是独立验证或组间分离检验。", "",
    "图件使用 R 绘制，导出矢量 PDF 与 600 dpi PNG；热图色标覆盖完整 Z-score 范围。仅输出上下调数量图、热图及各比较的火山图，绘图源数据保存在 source_data 目录。未知的随机化、盲法信息不作补充推断。")
  writeLines(methods, file.path(out, "DEG_methods_zh.md"), useBytes = TRUE)
  writeLines(c("Temperature DEG analysis", "", "DEG_summary.tsv: seven contrast summaries.",
    "tables/: full retained-gene results, DEG/up/down tables, sample/filter/model records and expression matrices.",
    "figures/: seven volcano plots, one DEG count plot and one heatmap, each as PDF and 600 dpi PNG.",
    "source_data/: values used by every figure.",
    "analysis_parameters.json and sessionInfo.txt record parameters and software versions.",
    "NA statistics remain unavailable, not non-significant. The seven BH correction families are separate.",
    "Run the script with --validate-only to validate inputs without creating results.",
    "No raw files, sample labels or biological replicates are modified."), file.path(out, "README.txt"))
  hashes_after <- sha256(input$source_files)
  check(identical(input$hashes_before, hashes_after), "Input files changed during analysis")
  check(all(summary$total_DEG == summary$up + summary$down), "DEG count mismatch")
  check(all(summary$retained_genes == summary$total_DEG + summary$tested_not_DEG + summary$unavailable),
        "Gene classification is not exhaustive")
  check(!nrow(z) || (max(abs(rowMeans(z))) < 1e-10 && max(abs(apply(z, 1, sd) - 1)) < 1e-10),
        "Heatmap Z-score calculation failed")
  write_json(list(status = "ANALYSIS_COMPLETE", started = started,
    finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), file.path(out, "run_status.json"))
  print(summary[, c("comparison", "up", "down", "total_DEG")], row.names = FALSE)
  message("Analysis complete: ", out)
  invisible(summary)
}

if (sys.nframe() == 0L) main()
