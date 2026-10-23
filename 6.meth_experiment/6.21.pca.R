#!/usr/bin/env Rscript
# Arabidopsis WGBS PCA. Raw cytosine reports are read only.
# Rscript 6.21.pca.R --output-dir DIR [--workers 3] [--validate-only]
# Streaming awk performs count-preserving aggregation; all PCA and drawing use R.
options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})

main <- function() {
  base <- "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng"
  opt <- list(input_dir=file.path(base,"Arabidopsis_thaliana/list/methdata/bismark"),
    output_dir=file.path(base,"Arabidopsis_thaliana/list/meth_experiment/PCA"),
    reference_sizes=file.path(base,"ref/Arabidopsis_thaliana/faidx/Arabidopsis_thaliana.chrSize.txt"),
    cg_windows="1000,2000,5000", chg_windows="1000,2000,5000", chh_windows="5000,10000,20000",
    min_depth=5, high_quantile=0.999, min_sites=5, top_n=2000, workers=3,
    validate_only=FALSE, aggregate_only=FALSE,
    qa_scripts="/datapool/home/2023102768/.agents/skills/nature-figure/scripts")
  args <- commandArgs(TRUE); i <- 1L
  while(i <= length(args)) {
    key <- gsub("-","_",sub("^--","",args[i]))
    if(key=="help") {cat("WGBS PCA options:\n",paste(names(opt),collapse="\n"),"\n");return(invisible(NULL))}
    if(!key %in% names(opt)) stop("Unknown option: ",args[i])
    if(key %in% c("validate_only","aggregate_only")) opt[[key]] <- TRUE else {
      i <- i+1L; if(i>length(args)) stop("Missing option value"); opt[[key]] <- args[i]
    }; i <- i+1L
  }
  for(k in c("min_depth","high_quantile","min_sites","top_n","workers")) opt[[k]] <- as.numeric(opt[[k]])
  stopifnot(opt$high_quantile>0,opt$high_quantile<=1,
    all(unlist(opt[c("min_depth","min_sites","top_n","workers")])>=1),opt$top_n>=2)
  widths <- lapply(opt[c("cg_windows","chg_windows","chh_windows")],function(x) sort(unique(as.integer(strsplit(x,",")[[1]]))))
  names(widths) <- c("CG","CHG","CHH")
  stopifnot(all(vapply(widths,function(w) all(is.finite(w)&w>0&w%%min(w)==0),logical(1))))
  for(p in c("ggrepel","jsonlite","systemfonts")) if(!requireNamespace(p,quietly=TRUE)) stop("Missing package ",p)
  stopifnot(capabilities("cairo"),nzchar(Sys.which("awk")),nzchar(Sys.which("pigz")))
  setDTthreads(2)
  check <- function(ok,msg) if(!isTRUE(ok)) stop(msg,call.=FALSE)
  temps <- seq(4,32,4)
  ids <- c("T4-1","T4-3","T4-6","T8-3","T8-5","T8-6","T12-4","T12-7",
    "T16-1","T16-2","T16-4","T20-3","T20-4","T20-5","T24-3","T24-5","T24-6",
    "T28-2","T28-3","T28-5","T32-5","T32-6")
  files <- list.files(opt$input_dir,pattern="^[^.]+[.]nameSorted[.]deduplicated[.]CX_report[.]txt[.]gz$",recursive=TRUE,full.names=TRUE)
  found <- basename(dirname(files)); check(length(files)==22&&!anyDuplicated(found)&&setequal(found,ids),"Input samples differ from verified 22-sample manifest")
  files <- files[match(ids,found)]
  samples <- data.table(sample_id=ids,temperature_C=as.integer(sub("^T([0-9]+)-.*","\\1",ids)),
    replicate_id=as.integer(sub(".*-","",ids)),input_file=normalizePath(files))
  samples[,is_control:=temperature_C==24]
  check(identical(as.integer(table(factor(samples$temperature_C,levels=temps))),c(3L,3L,2L,3L,3L,3L,3L,2L)),"Group sizes differ")
  ref <- fread(opt$reference_sizes,col.names=c("chr","length"))
  expected_chr <- c("NC_003070.9","NC_003071.7","NC_003074.8","NC_003075.7","NC_003076.8","NC_037304.1","NC_000932.1")
  check(identical(ref$chr,expected_chr),"Unexpected reference dictionary")
  ref[,scope:=c(rep("nuclear",5),"mitochondrion","chloroplast")]
  ref[,chr_order:=.I]
  bam_dictionary <- function(path) {
    h <- gzfile(path,"rb");on.exit(close(h));stopifnot(rawToChar(readBin(h,"raw",4))=="BAM\001")
    n <- readBin(h,"integer",1,size=4,endian="little");header<-readBin(h,"raw",n)
    nr <- readBin(h,"integer",1,size=4,endian="little")
    rbindlist(lapply(seq_len(nr),function(i) {
      n<-readBin(h,"integer",1,size=4,endian="little");z<-readBin(h,"raw",n)
      data.table(chr=rawToChar(z[-n]),length=readBin(h,"integer",1,size=4,endian="little"))
    }))
  }
  for(j in seq_along(ids)) {
    d<-dirname(files[j]);sid<-ids[j]
    b<-bam_dictionary(file.path(d,paste0(sid,".nameSorted.deduplicated.bam")))
    check(identical(b,ref[,.(chr,length)]),paste("Reference mismatch",sid))
    report<-readLines(file.path(d,paste0(sid,".R1.fit_bismark_bt2_PE_report.txt")))
    check(grepl(sid,report[1],fixed=TRUE)&&grepl("ref/Arabidopsis_thaliana/bismark/",report[2],fixed=TRUE),paste("Mapping log mismatch",sid))
    samples[j,mapping_percent:=as.numeric(sub(".*\\t([0-9.]+).*","\\1",grep("^Mapping efficiency:",report,value=TRUE)))]
    dedup<-readLines(file.path(d,paste0(sid,".nameSorted.deduplicated_splitting_report.txt")))
    check(any(grepl("No overlapping methylation calls specified",dedup,fixed=TRUE)),"Overlapping mates handling differs")
  }
  message("Input validation: 22 samples, expected temperature groups, identical TAIR10.1 BAM dictionaries.")
  if(opt$validate_only) return(invisible(samples))
  out<-opt$output_dir;dir.create(out,recursive=TRUE,showWarnings=FALSE)
  check(file.access(out,2)==0,"Output directory not writable")
  for(d in c("cache","matrices","sensitivity","QC","QA")) dir.create(file.path(out,d),showWarnings=FALSE)
  wt <- function(x,name) fwrite(x,file.path(out,name),sep="\t",na="NA",quote=FALSE)
  wt(samples,"WGBS_PCA_samples.tsv");wt(ref,"WGBS_PCA_reference.tsv")
  jsonlite::write_json(opt,file.path(out,"WGBS_PCA_parameters.json"),pretty=TRUE,auto_unbox=TRUE)
  writeLines(capture.output(sessionInfo()),file.path(out,"WGBS_PCA_sessionInfo.txt"))
  writeLines(c("Bismark extraction/alignment: v0.24.2 (sample reports)",
    system2("awk","--version",stdout=TRUE,stderr=TRUE)[1],
    system2("pigz","--version",stdout=TRUE,stderr=TRUE)[1]),file.path(out,"WGBS_PCA_command_line_versions.txt"))
  self <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
  file.copy(self,file.path(out,"6.21.pca.R"),overwrite=TRUE)
  # Inspect every source once for provenance; later verify the same digests again.
  message("Computing input checksums...")
  checksums<-data.table(path=c(files,opt$reference_sizes),md5=unname(tools::md5sum(c(files,opt$reference_sizes))))
  wt(checksums,"WGBS_PCA_input_checksums.tsv")
  sizes_path<-file.path(out,"cache/reference_sizes.tsv");fwrite(ref[,.(chr,length)],sizes_path,sep="\t",col.names=FALSE)
  # Complete input validation and exact coverage histograms, with bounded memory.
  histogram_awk <- r"---(
BEGIN { FS=OFS="\t" }
FNR==NR { sizes[$1]=$2; next }
{
 if(NF!=7 || !($1 in sizes) || $2 !~ /^[0-9]+$/ || $2<1 || $2>sizes[$1] ||
    ($3!="+" && $3!="-") || $4 !~ /^[0-9]+$/ || $5 !~ /^[0-9]+$/ ||
    ($6!="CG" && $6!="CHG" && $6!="CHH" && $6!="CN")) {
   print "Invalid cytosine record at line " FNR > "/dev/stderr"; exit 2
 }
 if($1!=prevchr) { if(seen[$1]++) {print "Noncontiguous chromosome" > "/dev/stderr";exit 2}; prevpos=0 }
 if($2<=prevpos) { print "Duplicate or unsorted cytosine coordinate" > "/dev/stderr";exit 2 }
 prevchr=$1;prevpos=$2
 hist[$1 OFS $6 OFS ($4+$5)]++
}
END { for(k in hist) print k,hist[k] }
)---"
  aggregate_awk <- r"---(
BEGIN { FS=OFS="\t"; w["CG"]=cg;w["CHG"]=chg;w["CHH"]=chh }
FNR==NR { cap[$1 OFS $2]=$3;next }
($6 in w) && ($4+$5>=mindepth) {
 k=$1 OFS $6 OFS int(($2-1)/w[$6]); dep=$4+$5
 n[k]++;m[k]+=$4;t[k]+=dep
 if(dep<=cap[$1 OFS $6]) {nc[k]++;mc[k]+=$4;tc[k]+=dep}
}
END { for(k in n) { print k,"uncapped",n[k],m[k]+0,t[k];if(nc[k]>0)print k,"capped",nc[k],mc[k]+0,tc[k] } }
)---"
  hprog<-file.path(out,"cache/coverage.awk");aprog<-file.path(out,"cache/aggregate.awk")
  writeLines(histogram_awk,hprog);writeLines(aggregate_awk,aprog)
  run_stream <- function(input,program,aux,target,vars=character()) {
    # bash pipefail detects both decompression and awk failures. Every path is shell quoted.
    command<-paste("set -o pipefail;",shQuote(Sys.which("pigz")),"-dc -p 1",shQuote(input),"|",
      shQuote(Sys.which("awk")),paste(vars,collapse=" "),"-f",shQuote(program),shQuote(aux),"- >",shQuote(paste0(target,".partial")))
    status<-system2("bash",c("-c",shQuote(command)))
    check(status==0,paste("Streaming input failed",input))
    check(file.rename(paste0(target,".partial"),target),"Cannot finalize stream output")
  }
  quantile_hist <- function(h,q) {setorder(h,depth); h$depth[which(cumsum(h$n)>=q*sum(h$n))[1]]}
  process_sample <- function(j) {
    sid<-ids[j]; message(format(Sys.time())," ",sid," coverage pass")
    dir<-file.path(out,"cache",sid);dir.create(dir,showWarnings=FALSE)
    signature<-list(md5=checksums$md5[j],ref=checksums$md5[nrow(checksums)],
      widths=widths,min_depth=opt$min_depth,high_quantile=opt$high_quantile,version=3L)
    rds<-file.path(dir,"aggregated.rds")
    if(file.exists(rds)) {z<-readRDS(rds);if(identical(z$signature,signature)){message(sid," using validated cache");return(z)}}
    histfile<-file.path(dir,"coverage_histogram.tsv")
    run_stream(files[j],hprog,sizes_path,histfile)
    h<-fread(histfile,col.names=c("chr","context","depth","n"));h<-merge(h,ref[,.(chr,scope)],by="chr")
    check(setequal(h$chr,ref$chr),paste("Missing reference chromosome",sid))
    thresholds<-h[context %in% names(widths)&depth>=opt$min_depth,
      .(upper_depth=quantile_hist(.SD,opt$high_quantile)),by=.(scope,context)]
    check(nrow(thresholds)==9,paste("Missing context/scope",sid))
    caps<-merge(ref[,.(chr,scope)],thresholds,by="scope",allow.cartesian=TRUE)
    capfile<-file.path(dir,"upper_depth.tsv");fwrite(caps[,.(chr,context,upper_depth)],capfile,sep="\t",col.names=FALSE)
    h<-merge(h,thresholds,by=c("scope","context"),all.x=TRUE)
    qc<-h[,.(report_sites=sum(n),covered_sites=sum(n[depth>0]),
      depth_ge_min_sites=sum(n[depth>=opt$min_depth]),
      retained_sites=sum(n[depth>=opt$min_depth & depth<=upper_depth]),
      high_depth_removed=sum(n[depth>upper_depth]),
      mean_covered_depth=sum(depth*n)/sum(n[depth>0]),
      median_covered_depth=quantile_hist(.SD[depth>0],0.5),
      p99_covered_depth=quantile_hist(.SD[depth>0],0.99),
      max_depth=max(depth),upper_depth=unique(upper_depth)),by=.(scope,context)]
    qc[,sample_id:=sid];thresholds[,sample_id:=sid]
    message(format(Sys.time())," ",sid," window pass")
    aggfile<-file.path(dir,"window_counts.tsv")
    vars<-c(paste0("-v mindepth=",opt$min_depth),paste0("-v cg=",min(widths$CG)),
      paste0("-v chg=",min(widths$CHG)),paste0("-v chh=",min(widths$CHH)))
    run_stream(files[j],aprog,capfile,aggfile,vars)
    agg<-fread(aggfile,col.names=c("chr","context","bin","filter","sites","meth","total"))
    agg[,sample_id:=sid];check(all(agg$meth>=0&agg$meth<=agg$total),"Invalid window counts")
    z<-list(signature=signature,aggregate=agg,qc=qc,thresholds=thresholds)
    saveRDS(z,rds,compress=FALSE);message(format(Sys.time())," ",sid," aggregation complete");z
  }
  results<-parallel::mclapply(seq_along(ids),process_sample,mc.cores=opt$workers,mc.preschedule=FALSE)
  check(!any(vapply(results,inherits,logical(1),"try-error")),"At least one sample aggregation failed")
  siteqc<-rbindlist(lapply(results,`[[`,"qc"))
  readqc<-rbindlist(lapply(seq_along(ids),function(j) {
    h<-fread(file.path(out,"cache",ids[j],"coverage_histogram.tsv"),col.names=c("chr","context","depth","n"))
    h<-merge(h,ref[,.(chr,scope)],by="chr")
    h<-merge(h,results[[j]]$thresholds,by=c("scope","context"))
    h[,.(total_reads_before=sum(as.double(depth)*n),
      reads_after_min=sum(as.double(depth[depth>=opt$min_depth])*n[depth>=opt$min_depth]),
      retained_reads=sum(as.double(depth[depth>=opt$min_depth & depth<=upper_depth])*n[depth>=opt$min_depth & depth<=upper_depth]),
      high_depth_removed_reads=sum(as.double(depth[depth>upper_depth])*n[depth>upper_depth])),by=.(sample_id,scope,context)]
  }))
  siteqc<-merge(siteqc,readqc,by=c("sample_id","scope","context"),all.x=TRUE)
  siteqc[,mean_retained_depth:=retained_reads/retained_sites]
  wt(siteqc,"WGBS_PCA_site_coverage_QC.tsv")
  check(all(siteqc[,.(n_samples=.N,n_distinct_report_counts=uniqueN(report_sites)),by=.(scope,context)][,
    n_samples==22 & n_distinct_report_counts==1]),
    "Full cytosine reports contain inconsistent reference-site totals across samples")
  wt(rbindlist(lapply(results,`[[`,"thresholds")),"WGBS_PCA_depth_thresholds.tsv")
  aggregated<-rbindlist(lapply(results,`[[`,"aggregate"));rm(results);gc()
  if(opt$aggregate_only) {message("Aggregation complete.");return(invisible(NULL))}
  
  # Main and sensitivity matrices are constructed from exactly the same count sums.
  window_qc<-list();summary_rows<-list();analysis<-list();matrix_meta<-list()
  build_matrix <- function(ctx,width,mode="capped",scope="nuclear",record=TRUE) {
    chruse<-if(scope=="nuclear") ref[scope=="nuclear",chr] else ref$chr
    grid<-rbindlist(lapply(chruse,function(ch) {
      len<-ref[chr==ch,length];start<-seq(0,len-1,by=width)
      data.table(chr=ch,start0=start,end0=pmin(start+width,len),bin=seq_along(start)-1L)
    }))
    grid[,window_id:=paste(chr,start0,end0,sep=":")]
    d<-aggregated[context==ctx & filter==mode & chr %in% chruse,
      .(sites=sum(sites),meth=sum(meth),total=sum(total)),
      by=.(sample_id,chr,bin=floor(bin*min(widths[[ctx]])/width))]
    d[,p:=fifelse(sites>=opt$min_sites,meth/total,NA_real_)]
    d<-merge(d,grid,by=c("chr","bin"))
    wide<-dcast(d,window_id~sample_id,value.var="p",fill=NA_real_)
    wide<-merge(grid,wide,by="window_id",all.x=TRUE,sort=FALSE)
    setorder(wide,chr,start0)
    m<-as.matrix(wide[,..ids]);rownames(m)<-wide$window_id
    complete<-rowSums(is.na(m))==0;variances<-rep(NA_real_,nrow(m))
    transformed<-asin(sqrt(m[complete,,drop=FALSE]))
    variances[complete]<-apply(transformed,1,var)
    good<-complete & is.finite(variances)&variances>0
    idx<-which(good);ranked<-idx[order(-variances[idx],match(wide$chr[idx],ref$chr),wide$start0[idx])]
    selected<-head(ranked,opt$top_n)
    meta<-wide[,.(window_id,chr,start0,end0)]
    meta[,`:=`(n_reliable_samples=rowSums(!is.na(m)),complete=complete,transformed_variance=variances,
      eligible=good,variance_rank=match(seq_len(.N),ranked),selected_top=seq_len(.N)%in%selected)]
    key<-paste(ctx,width,mode,scope,sep="_")
    if(record) {
      q<-d[,.(covered_windows=.N,reliable_windows=sum(!is.na(p)),
        mean_retained_sites=mean(sites),median_window_total=as.numeric(median(total))),by=sample_id]
      q[,`:=`(context=ctx,width_bp=width,filter=mode,scope=scope,
        possible_windows=nrow(grid),reliable_fraction=reliable_windows/nrow(grid))]
      window_qc[[key]]<<-q
      summary_rows[[key]]<<-data.table(context=ctx,width_bp=width,filter=mode,scope=scope,
        possible_windows=nrow(grid),any_covered_windows=nrow(unique(d[,.(chr,bin)])),
        complete_windows=sum(complete),eligible_windows=sum(good),selected_windows=length(selected))
    }
    list(key=key,meta=meta,raw=m,transformed=asin(sqrt(m[good,,drop=FALSE])),
      selected=meta$window_id[selected],counts=d)
  }
  run_pc <- function(mat) {check(nrow(mat)>=2&&all(is.finite(mat)),"PCA matrix incomplete or too small");prcomp(t(mat),center=TRUE,scale.=FALSE)}
  pc_variance <- function(p) {v<-p$sdev^2;data.table(PC=colnames(p$x),eigenvalue=v,percent=100*v/sum(v),cumulative_percent=100*cumsum(v)/sum(v))}
  scores <- function(p) cbind(samples[,.(sample_id,temperature_C,replicate_id,is_control)],as.data.table(p$x))
  save_analysis <- function(b,ctx,prefix,top=TRUE) {
    m<-if(top)b$transformed[b$selected,,drop=FALSE] else b$transformed
    p<-run_pc(m);wt(scores(p),paste0(prefix,"_scores.tsv"));wt(pc_variance(p),paste0(prefix,"_variance.tsv"))
    list(pca=p,matrix=m,context=ctx,key=b$key,prefix=prefix)
  }
  main_width<-setNames(integer(3),names(widths))
  for(ctx in names(widths)) {
    for(width in widths[[ctx]]) {
      b<-build_matrix(ctx,width)
      message(ctx," ",width," bp: ",nrow(b$transformed)," complete nonconstant windows")
      if(length(b$selected)>=opt$top_n && main_width[ctx]==0L) {
        main_width[ctx]<-width
        wt(b$meta,paste0(ctx,"_window_selection.tsv"))
        fwrite(data.table(window_id=rownames(b$raw),b$raw),file.path(out,"matrices",paste0(ctx,"_window_methylation.tsv.gz")),sep="\t",na="NA")
        fwrite(b$counts,file.path(out,"matrices",paste0(ctx,"_window_counts.tsv.gz")),sep="\t",na="NA")
        saveRDS(list(metadata=b$meta,raw=b$raw,transformed=b$transformed,selected=b$selected),file.path(out,"matrices",paste0(ctx,"_matrices.rds")))
        analysis[[ctx]]<-save_analysis(b,ctx,paste0(ctx,"_PCA"))
        analysis[[paste0(ctx,"_all")]]<-save_analysis(b,ctx,paste0("sensitivity/",ctx,"_PCA_all_windows"),FALSE)
      }
      if(nrow(b$transformed)>=2) analysis[[paste0(ctx,"_w",width)]]<-save_analysis(b,ctx,paste0("sensitivity/",ctx,"_PCA_",width,"bp"))
    }
    if(main_width[ctx]==0L) {
      wt(rbindlist(summary_rows),"WGBS_PCA_window_filtering.tsv")
      stop(ctx,": fewer than ",opt$top_n," complete nonconstant windows at every configured width; no implicit relaxation")
    }
    for(kind in c("uncapped","organelles")) {
      b<-build_matrix(ctx,main_width[ctx],mode=if(kind=="uncapped")"uncapped" else "capped",
        scope=if(kind=="organelles")"all" else "nuclear")
      analysis[[paste0(ctx,"_",kind)]]<-save_analysis(b,ctx,paste0("sensitivity/",ctx,"_PCA_",kind))
    }
  }
  wt(rbindlist(summary_rows),"WGBS_PCA_window_filtering.tsv")
  winqc<-rbindlist(window_qc);wt(winqc,"WGBS_PCA_sample_window_QC.tsv")
  opt$actual_main_windows<-as.list(main_width)
  jsonlite::write_json(opt,file.path(out,"WGBS_PCA_parameters.json"),pretty=TRUE,auto_unbox=TRUE)
  saveRDS(analysis,file.path(out,"WGBS_PCA_analysis.rds"))

  # Descriptive diagnostics do not classify or delete samples automatically.
  qc_rows<-list();associations<-list();sensitivity<-list()
  for(name in names(analysis)) {
    a<-analysis[[name]];m<-a$matrix;dd<-as.matrix(dist(t(m)));cc<-cor(m)
    wt(data.table(sample_id=ids,dd),paste0(a$prefix,"_distances.tsv"))
    wt(data.table(sample_id=ids,cc),paste0(a$prefix,"_correlations.tsv"))
    qc_rows[[name]]<-rbindlist(lapply(seq_along(ids),function(j) {
      other<-setdiff(seq_along(ids),j);within<-other[samples$temperature_C[other]==samples$temperature_C[j]]
      nearest<-other[which.min(dd[j,other])]
      data.table(analysis=name,context=a$context,sample_id=ids[j],temperature_C=samples$temperature_C[j],
        mean_within_distance=mean(dd[j,within]),mean_between_distance=mean(dd[j,setdiff(other,within)]),
        min_within_correlation=min(cc[j,within]),nearest_sample=ids[nearest],
        nearest_same_temperature=samples$temperature_C[j]==samples$temperature_C[nearest])
    }))
    main<-analysis[[a$context]]
    sensitivity[[name]]<-data.table(analysis=name,context=a$context,windows=nrow(m),
      PC1_percent=pc_variance(a$pca)$percent[1],PC2_percent=pc_variance(a$pca)$percent[2],
      distance_spearman_vs_main=cor(as.vector(dist(t(main$matrix))),as.vector(dist(t(m))),method="spearman"))
  }
  qc<-rbindlist(qc_rows);wt(qc,"QC/WGBS_PCA_sample_diagnostics.tsv")
  wt(rbindlist(sensitivity),"sensitivity/WGBS_PCA_sensitivity_comparison.tsv")
  for(ctx in names(widths)) {
    d<-merge(samples,siteqc[context==ctx & scope=="nuclear"],by="sample_id")
    d<-merge(d,winqc[context==ctx & scope=="nuclear" & filter=="capped" & width_bp==main_width[ctx],.(sample_id,reliable_fraction)],by="sample_id")
    d<-d[match(ids,d$sample_id)]
    for(pc in 1:2) for(metric in c("mean_covered_depth","retained_sites","reliable_fraction","temperature_C")) {
      associations[[paste(ctx,pc,metric)]]<-data.table(context=ctx,PC=paste0("PC",pc),metric=metric,
        spearman_rho=cor(analysis[[ctx]]$pca$x[,pc],d[[metric]],method="spearman"))
    }
  }
  assoc<-rbindlist(associations);wt(assoc,"QC/WGBS_PCA_PC_QC_associations.tsv")
  integrated<-merge(qc[analysis %in% names(widths)],siteqc[scope=="nuclear"],
    by=c("context","sample_id"))
  chosen_qc<-winqc[scope=="nuclear" & filter=="capped"]
  chosen_qc<-chosen_qc[width_bp==unname(main_width[context])]
  integrated<-merge(integrated,chosen_qc[,.(context,sample_id,reliable_fraction,reliable_windows)],
    by=c("context","sample_id"))
  wt(integrated,"QC/WGBS_PCA_integrated_sample_QC.tsv")

  palette<-setNames(c("#4477AA","#66CCEE","#228833","#CCBB44","#EE6677","#AA3377","#BBBBBB","#000000"),temps)
  fonts<-systemfonts::system_fonts()$family
  font<-if("Arial"%in%fonts)"Arial" else if("DejaVu Sans"%in%fonts)"DejaVu Sans" else "sans"
  make_plot <- function(a,labeled=FALSE,combined=FALSE) {
    d<-scores(a$pca);d[,temperature:=factor(temperature_C,levels=temps)]
    v<-pc_variance(a$pca)$percent
    p<-ggplot(d,aes(PC1,PC2,colour=temperature))+geom_point(size=2.5,alpha=0.9)+
      scale_colour_manual(values=palette,breaks=as.character(temps),labels=paste0(temps," °C"),drop=FALSE)+
      scale_x_continuous(expand=expansion(mult=if(labeled)0.23 else 0.10),breaks=scales::breaks_pretty(n=3))+
      scale_y_continuous(expand=expansion(mult=if(labeled)0.23 else 0.10),breaks=scales::breaks_pretty(n=3),labels=function(x)sprintf("% .1f",x))+
      labs(title=a$context,x=sprintf("PC1 (%.2f%%)",v[1]),y=sprintf("PC2 (%.2f%%)",v[2]),colour="Temperature")+
      theme_classic(base_size=9,base_family=font)+theme(
        axis.line=element_line(linewidth=0.3),axis.ticks=element_line(linewidth=0.3),
        axis.text=element_text(size=8,colour="black"),axis.title=element_text(size=9),
        plot.title=element_text(size=9,face="bold",hjust=0.5,margin=margin(b=8)),
        plot.title.position="panel",
        aspect.ratio=0.85,legend.position="right",
        legend.text=element_text(size=8),legend.title=element_text(size=9),
        plot.tag=element_text(size=10,face="bold",hjust=0,vjust=1),
        plot.tag.position=c(0,1),plot.tag.location="plot",
        legend.key.height=grid::unit(5,"mm"),plot.margin=margin(6,7,6,6),
        plot.background=element_rect(fill="white",colour=NA))
    if(labeled) p<-p+ggrepel::geom_text_repel(aes(label=sample_id),seed=20261010,size=2.8,
      family=font,colour="#303030",max.overlaps=Inf,max.iter=100000,max.time=Inf,
      box.padding=0.8,point.padding=0.65,min.segment.length=0,segment.size=0.2,segment.colour="#888888",
      nudge_x=ifelse((a$context=="CG" & d$sample_id=="T16-4") |
        (a$context=="CHG" & d$sample_id=="T4-3"),-0.1*diff(range(d$PC1)),0))
    p
  }
  # Use the same Cairo font metrics for the alignment probe and final PDF.
  # The generic helper's base pdf() probe cannot resolve system DejaVu fonts.
  alignment_source<-readLines(file.path(opt$qa_scripts,"panel_alignment.R"))
  panel_row_line<-grep("  rows <- grob$layout[grepl(",alignment_source,fixed=TRUE)
  check(length(panel_row_line)==1L,"Unexpected panel-row helper interface")
  alignment_source[panel_row_line]<-'  rows <- grob$layout[grepl("^panel(?:-[0-9]+|;.*-[0-9]+)?$", grob$layout$name, perl = TRUE), , drop = FALSE]'
  probe_line<-grep("grDevices::pdf(probe_path",alignment_source,fixed=TRUE)
  check(length(probe_line)==1L,"Unexpected panel alignment helper interface")
  alignment_source[probe_line]<-"  grDevices::cairo_pdf(probe_path, width = width_in, height = height_in, family = font)"
  # null-unit gtable columns convert to zero outside their rendered viewport.
  # Measure the drawn panel viewports in device coordinates instead.
  measure_start<-grep("  widths_pt <- grid::convertWidth",alignment_source,fixed=TRUE)
  measure_end<-grep("  manifest <- list(",alignment_source,fixed=TRUE)
  check(length(measure_start)==1L && length(measure_end)==1L,"Unexpected alignment measurement interface")
  measure_code<-r"---(
  panels <- lapply(seq_len(nrow(panel_rows)), function(index) {
    row <- panel_rows[index, ]
    vp <- sprintf("%s.%d-%d-%d-%d", row$name, row$t, row$l, row$b, row$r)
    grid::seekViewport(vp)
    if (startsWith(row$name, "panel;")) grid::downViewport("panel.1-1-1-1")
    lo <- grid::deviceLoc(grid::unit(0, "npc"), grid::unit(0, "npc"), valueOnly = TRUE)
    hi <- grid::deviceLoc(grid::unit(1, "npc"), grid::unit(1, "npc"), valueOnly = TRUE)
    grid::upViewport(0)
    list(id = panel_ids[index], bbox_pt = c(lo$x, lo$y, hi$x, hi$y) * 72,
      grid_id = "patchwork-grid-1", row_start = as.integer(row$t - 1),
      row_stop = as.integer(row$b), col_start = as.integer(row$l - 1),
      col_stop = as.integer(row$r))
  })
)---"
  alignment_source<-c(alignment_source[seq_len(measure_start-1L)],
    strsplit(measure_code,"\n",fixed=TRUE)[[1]],alignment_source[measure_end:length(alignment_source)])
  # Open Cairo before patchwork builds text grobs, not only before grid.draw().
  probe_start<-grep("  probe_path <- tempfile",alignment_source,fixed=TRUE)
  probe_end<-grep("  grid::grid.newpage()",alignment_source,fixed=TRUE)-1L
  probe_setup<-alignment_source[probe_start:probe_end]
  alignment_source<-alignment_source[-(probe_start:probe_end)]
  grob_start<-grep("  grob <- patchwork::patchworkGrob",alignment_source,fixed=TRUE)
  alignment_source<-append(alignment_source,probe_setup,after=grob_start-1L)
  alignment_helper<-file.path(out,"QA/panel_alignment_cairo.R")
  writeLines(alignment_source,alignment_helper)
  source(alignment_helper,local=TRUE)
  export_plot <- function(p,prefix,width=160,height=120,multi=FALSE) {
    if(multi) require_patchwork_panel_alignment(p,
      manifest_path=file.path(out,"QA",paste0(basename(prefix),".alignment-layout.json")),
      report_path=file.path(out,"QA",paste0(basename(prefix),".alignment.json")),
      width_in=width/25.4,height_in=height/25.4,
      audit_script=file.path(opt$qa_scripts,"audit_panel_alignment.py"),strict=TRUE)
    grDevices::cairo_pdf(file.path(out,paste0(prefix,".pdf")),width=width/25.4,height=height/25.4,family=font,bg="white")
    print(p);invisible(dev.off())
    grDevices::png(file.path(out,paste0(prefix,".png")),width=width,height=height,units="mm",res=600,type="cairo",bg="white")
    print(p);invisible(dev.off())
  }
  for(name in names(analysis)) export_plot(make_plot(analysis[[name]]),analysis[[name]]$prefix)
  for(ctx in names(widths)) {
    export_plot(make_plot(analysis[[ctx]],TRUE),paste0(ctx,"_PCA_labeled"),185,140)
    for(ext in c("pdf","png")) file.copy(file.path(out,paste0(ctx,"_PCA_labeled.",ext)),
      file.path(out,"QC",paste0(ctx,"_PCA_labeled.",ext)),overwrite=TRUE)
  }
  make_combined <- function(labeled=FALSE) {
    panels<-lapply(seq_along(widths),function(j) {
      ctx<-names(widths)[j]
      make_plot(analysis[[ctx]],labeled)+labs(tag=paste0("(",letters[j],")"))
    })
    wrap_plots(panels,nrow=1,guides="keep")
  }
  combined<-make_combined()
  combined_labeled<-make_combined(TRUE)
  export_plot(combined,"WGBS_PCA_combined",480,120,TRUE)
  export_plot(combined_labeled,"WGBS_PCA_combined_labeled",555,140,TRUE)
  # R-only preview, distinct from the 600 dpi manuscript exports.
  grDevices::png(file.path(out,"QA/WGBS_PCA_preview.png"),width=2400,height=600,type="cairo",res=127)
  print(combined);invisible(dev.off())

  stats<-rbindlist(summary_rows)
  method<-c("# 全基因组 DNA 甲基化主成分分析方法", "",
    "## 可用于硕士学位论文的方法正文",
    "拟南芥（Arabidopsis thaliana）Col-0 植株经 3 d 温度处理，设置 4、8、12、16、20、24、28 和 32 ℃共 8 个温度梯度，其中 24 ℃为对照。各组有效 WGBS 生物学重复数依次为 3、3、2、3、3、3、3 和 2，共 22 个样本。分析保留全部有效重复，未按温度求平均，未依据 PCA 结果删除样本。",
    "采用 Bismark v0.24.2 输出的完整胞嘧啶报告文件，分别读取 CG、CHG 和 CHH 背景的染色体、坐标、链方向、甲基化及未甲基化 reads 数。上游已去除重复比对，甲基化提取指定不重复计算配对 reads 重叠区域。所有样本使用 TAIR10.1 参考基因组；分析前核对 BAM 参考序列名称与长度，并检查胞嘧啶坐标与计数。主分析覆盖全部 5 条核染色体，叶绿体和线粒体单独进行覆盖度 QC，并用于包含细胞器的敏感性检查。",
    sprintf("按染色体建立固定、不重叠的基因组窗口，CG、CHG、CHH 主分析的实际窗口大小分别为 %d、%d、%d bp。胞嘧啶报告使用 1-based 坐标，窗口表使用 0-based 半开区间。末端窗口截断至染色体末端，正负链胞嘧啶均参与窗口汇总。",main_width[1],main_width[2],main_width[3]),
    sprintf("保留单个位点总覆盖度不低于 %d× 的记录。为减少极端高覆盖位点影响，按样本和序列背景，将核基因组满足最低覆盖度记录的覆盖度第 %.3f 百分位作为上限，剔除严格超过该阈值的记录；细胞器单独计算阈值。离散分位数采用累计频数首次达到目标比例的深度值，不插值、不对 reads 计数缩尾。另以不设置覆盖度上限进行敏感性检查。",opt$min_depth,100*opt$high_quantile),
    sprintf("每个样本每个窗口至少包含 %d 个通过覆盖度筛选的胞嘧啶记录。窗口甲基化水平计算为有效位点甲基化 reads 数之和除以甲基化与未甲基化 reads 总和，不使用位点百分比的算术平均值。可靠覆盖下甲基化 reads 为零的窗口保留为真实零值；未达到可靠性标准的窗口记为缺失值 NA，不进行零填充或其他插补。只保留在全部 22 个样本均具有可靠数据的窗口。",opt$min_sites),
    sprintf("将窗口甲基化比例 p 转换为 asin(sqrt(p))，根据转换后在 22 个样本中的方差排序，去除零方差窗口，分别选取方差最高的 %d 个窗口用于主分析。采用 R 的 prcomp() 函数，以样本为行、窗口为列，设置 center = TRUE、scale. = FALSE。主成分解释方差由各背景的特征值独立计算。",opt$top_n),
    paste0("敏感性检查包括全部合格非恒定窗口、不同窗口尺度（CG：",paste(widths$CG,collapse="/")," bp；CHG：",paste(widths$CHG,collapse="/")," bp；CHH：",paste(widths$CHH,collapse="/")," bp）、不设覆盖度上限及纳入细胞器。不同尺度不以组间分离程度进行选择。"),
    "采用 ggplot2 绘制 PC1–PC2 散点图，每个点表示一个生物学重复，颜色与 RNA-seq PCA 的温度映射一致。图中不绘制置信椭圆，不抖动或人为改变 PCA 坐标。输出矢量 PDF 和 600 dpi PNG，并提供带样本编号的 QC 图。组内距离、相关性、PC 与覆盖度及温度的 Spearman 相关仅作为描述性诊断，不作为差异显著性检验或删除样本的依据。", "",
    "## 本次实际保留窗口与解释方差")
  caption<-c("图：不同温度处理下拟南芥全核基因组 DNA 甲基化 PCA。(a)，CG；(b)，CHG；(c)，CHH。Col-0 经 3 d 温度处理，24 ℃为对照。22 个点对应 22 个有效生物学重复；12 ℃和 32 ℃各 2 个重复，其余温度各 3 个重复。")
  qc_text<-c("# WGBS PCA 结果与质量检查", "", "所有结果为探索性描述，不进行样本剔除或温度组显著性推断。")
  for(ctx in names(widths)) {
    s<-stats[context==ctx & width_bp==main_width[ctx] & filter=="capped" & scope=="nuclear"]
    v<-pc_variance(analysis[[ctx]]$pca)$percent
    line<-sprintf("%s：%d bp；理论窗口 %d；22 样本完整覆盖窗口 %d；非恒定合格窗口 %d；主分析窗口 %d；PC1 = %.2f%%，PC2 = %.2f%%。",ctx,main_width[ctx],s$possible_windows,s$complete_windows,s$eligible_windows,s$selected_windows,v[1],v[2])
    depth_summary<-siteqc[context==ctx & scope=="nuclear"]
    depth_line<-sprintf("%s 各样本采用的位点覆盖度上限为 %d–%d×；覆盖位点的平均深度范围为 %.2f–%.2f×。",ctx,
      min(depth_summary$upper_depth),max(depth_summary$upper_depth),
      min(depth_summary$mean_covered_depth),max(depth_summary$mean_covered_depth))
    method<-c(method,line,depth_line);caption<-c(caption,line);qc_text<-c(qc_text,"",paste0("## ",ctx),line,depth_line)
    q<-qc[analysis==ctx];setorder(q,-mean_within_distance)
    qc_text<-c(qc_text,sprintf("最近邻来自同温度组的样本：%d/22。组内平均距离最大的样本为 %s（%.4f），此排序不等于异常检验。",sum(q$nearest_same_temperature),q$sample_id[1],q$mean_within_distance[1]))
    iq<-integrated[context==ctx];setorder(iq,-mean_within_distance)
    qc_text<-c(qc_text,"优先人工复核组内距离最大的三个样本，以下指标用于区分覆盖度问题与生物学差异；不据此删除样本：",
      vapply(seq_len(min(3,nrow(iq))),function(j)sprintf("%s：组内平均距离 %.4f；组内最低相关性 %.4f；平均覆盖度 %.2f×；可靠窗口比例 %.2f%%。",
        iq$sample_id[j],iq$mean_within_distance[j],iq$min_within_correlation[j],iq$mean_covered_depth[j],100*iq$reliable_fraction[j]),character(1)))
    groups<-q[,.(mean_within_distance=mean(mean_within_distance),
      mean_between_distance=mean(mean_between_distance),min_within_correlation=min(min_within_correlation)),by=temperature_C]
    setorder(groups,temperature_C)
    qc_text<-c(qc_text,vapply(seq_len(nrow(groups)),function(j)sprintf("%d ℃：平均组内距离 %.4f；平均组间距离 %.4f；最低组内相关性 %.4f。",
      groups$temperature_C[j],groups$mean_within_distance[j],groups$mean_between_distance[j],groups$min_within_correlation[j]),character(1)))
    aa<-assoc[context==ctx & metric %in% c("mean_covered_depth","reliable_fraction","temperature_C")]
    qc_text<-c(qc_text,apply(as.data.frame(aa),1,function(z)sprintf("%s 与 %s 的 Spearman ρ = %.4f。",z[2],z[3],as.numeric(z[4]))))
  }
  qc_text<-c(qc_text,"","## 跨背景比较与敏感性分析",
    "PC1–PC2 图只展示总变异的一部分，二维图上相邻或分离不能替代全维样本距离和组内相关性检查。主成分符号本身没有生物学意义。")
  sens<-rbindlist(sensitivity)
  for(ctx in names(widths)) {
    v<-pc_variance(analysis[[ctx]]$pca)$percent
    ss<-sens[context==ctx]
    qc_text<-c(qc_text,sprintf("%s 的前两个主成分累计解释 %.2f%%；全部合格窗口、不设深度上限、纳入细胞器的样本距离与主分析的 Spearman ρ 依次为 %.4f、%.4f、%.4f。",ctx,sum(v[1:2]),
      ss[analysis==paste0(ctx,"_all"),distance_spearman_vs_main],
      ss[analysis==paste0(ctx,"_uncapped"),distance_spearman_vs_main],
      ss[analysis==paste0(ctx,"_organelles"),distance_spearman_vs_main]),
      sprintf("%s 不同窗口尺度的距离一致性范围为 %.4f–%.4f。",ctx,
        min(ss[grepl("_w[0-9]+$",analysis),distance_spearman_vs_main]),
        max(ss[grepl("_w[0-9]+$",analysis),distance_spearman_vs_main])))
  }
  qc_text<-c(qc_text,
    "观察到主成分与覆盖度相关时，应同时考虑测序精度和生物学差异；仅凭本数据中的相关性不能确定异常原因。建议结合测序批次、DNA 质量和独立实验信息核实偏离重复。全部样本均保留。")
  method<-c(method,"","## 软件与复现",R.version.string,
    paste(vapply(c("data.table","ggplot2","patchwork","ggrepel"),function(p)paste(p,packageVersion(p)),character(1)),collapse="；"),
    "完整参数、样本清单、输入校验值及 sessionInfo 随结果提供。全部原始数据保持只读。")
  writeLines(method,file.path(out,"WGBS_PCA_methods_zh.md"),useBytes=TRUE)
  pandoc<-Sys.which("pandoc")
  if(!nzchar(pandoc)) pandoc<-file.path(R.home(),"../../bin/pandoc")
  check(file.exists(pandoc),"Pandoc is required to export the Chinese methods as Word")
  status<-system2(pandoc,c(shQuote(file.path(out,"WGBS_PCA_methods_zh.md")),
    "--from=markdown","--to=docx","--standalone","-o",shQuote(file.path(out,"WGBS_PCA_methods_zh.docx"))))
  check(status==0,"Chinese Word export failed")
  writeLines(caption,file.path(out,"WGBS_PCA_caption_zh.txt"),useBytes=TRUE)
  writeLines(qc_text,file.path(out,"QC/WGBS_PCA_results_zh.md"),useBytes=TRUE)
  check(identical(checksums$md5,unname(tools::md5sum(checksums$path))),"Raw input checksum changed during analysis")
  message("Completed: ",out)
}
main()
