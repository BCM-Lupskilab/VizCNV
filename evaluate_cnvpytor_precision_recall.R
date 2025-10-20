#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(dplyr)
})

# ---------- user input ----------
cnvpytor_dir <- "1kg38_cnvpytor_out/tsv/"      # folder containing per-sample CNVpytor TSVs
truth_dir   <- "1kg38/" # folder containing true call
sample_list  <- "sample_list.txt"   # one sample name per line
min_size     <- 2000
out_prefix   <- "cnvpytor_eval"

# ---------- load truth set ----------
truth <- fread(truth_file, col.names = c("chr","start","end","id"))
truth.gr <- makeGRangesFromDataFrame(truth, keep.extra.columns = TRUE)

# ---------- overlap helper ----------
reciprocal_hits <- function(query, subject, frac = 0.5){
  hits <- findOverlaps(query, subject)
  qw <- width(query)[queryHits(hits)]
  sw <- width(subject)[subjectHits(hits)]
  iw <- width(pintersect(query[queryHits(hits)], subject[subjectHits(hits)]))
  keep <- iw/qw >= frac & iw/sw >= frac
  hits[keep]
}

# ---------- mean ± 95 % CI helper ----------
mean_ci <- function(x){
  m <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))
  ci <- 1.96*se
  c(mean = m, lower = m - ci, upper = m + ci)
}

# ---------- per-sample evaluation ----------
#samples <- fread(sample_list, header = FALSE)[[1]]
res <- list()
autosomes = paste0("chr", 1:22)
for (s in control.samples) {
  f <- file.path(cnvpytor_dir, paste0(s, ".100.raw.tsv"))
  f2 <- file.path(truth_dir, paste0(s, "_truth_autosomal_het_DEL_DUP_over2kb_nosegdup.bed"))
 
  if (!file.exists(f)) { warning("Missing ", f); next }
  message("Processing ", s)
  if (!file.exists(f2)) { warning("Missing ", f2); next }
  message("Processing ", s)
  truth <- fread(f2,header = T)
  truth.gr <- makeGRangesFromDataFrame(truth,keep.extra.columns = T)
  call <- fread(f)%>%
    tidyr::separate(V2,into = c("chr","start","end"),sep=":|-")%>%
    filter(V11>1000)%>% # distance from closest large (>1000 bp) gap in reference genome.
    filter(V10<0.5)%>% # fraction of reference genome gaps (Ns) in call region less than 0.5
    filter(V9<0.5)%>% # fraction of reads mapped with q0 quality in call region over than 0.5
    filter(chr %in% autosomes)%>%
    filter(
      (V4 >= 1*0.9/2 & V4 <= 1*1.1/2) |
        (V4 >= 3*0.85/2 & V4 <= 3*1.15/2)
    )%>%
    mutate(svtype=ifelse(V1=="deletion","DEL","DUP"),sample=s)
  call <- call %>% filter(V3>min_size)
  call.gr <- makeGRangesFromDataFrame(call,keep.extra.columns = T)
  
  #hits <- reciprocal_hits(call.gr, truth.gr, frac = 0.5) 
  
  
  recall_raw <- rbindlist(lapply(c("DEL","DUP"),function(type){
    type.gr <- call.gr[call.gr$svtype==type]
    truth_type.gr <- truth.gr[truth.gr$svtype==type]
    hits <- findOverlaps(type.gr, truth_type.gr ) # for any overlap
    type.dt <- as.data.table(type.gr)
    truth_type.dt <- as.data.table(truth_type.gr)
    type.dt[,cnvpytorcnvID:=paste0(seqnames,"_",start,"_",end)]
    truth_type.dt$predictor <- 0
    truth_type.dt$predictor_region <- "NA"
    truth_type.dt$predictor_segmean <- "NA"
    truth_type.dt$predictor[subjectHits(hits)] <- 1
    truth_type.dt$predictor_region[subjectHits(hits)] <- type.dt[queryHits(hits),cnvpytorcnvID]
    truth_type.dt$predictor_segmean[subjectHits(hits)] <- type.dt[queryHits(hits),V4]
    return(truth_type.dt)
  }))
  precision_raw <- rbindlist(lapply(c("DEL","DUP"),function(type){
    type.gr <- call.gr[call.gr$svtype==type]
    truth_type.gr <- truth.gr[truth.gr$svtype==type]
    hits <- findOverlaps(type.gr, truth_type.gr ) # for any overlap
    type.dt <- as.data.table(type.gr)
    truth_type.dt <- as.data.table(truth_type.gr)
    truth_type.dt[,truethID:=paste0(seqnames,"_",start,"_",end)]
    type.dt[,cnvpytorcnvID:=paste0(seqnames,"_",start,"_",end)]
    type.dt$predictor <- 0
    type.dt$predictor_region <- "NA"
    type.dt$predictor_segmean <- "NA"
    type.dt$predictor[queryHits(hits)] <- 1
    type.dt$predictor_region[queryHits(hits)] <- truth_type.dt[subjectHits(hits),id]
    type.dt$predictor_segmean[queryHits(hits)] <- truth_type.dt[subjectHits(hits),truethID]
    return(type.dt)
  }))
  bysize <- rbindlist(lapply(c("DEL","DUP"),function(type){
    type.gr <- call.gr[call.gr$svtype==type]
    truth_type.gr <- truth.gr[truth.gr$svtype==type]
    hits <- findOverlaps(type.gr, truth_type.gr ) # for any overlap
    TP <- length(unique(subjectHits(hits)))
    FP <- length(type.gr) - length(unique(queryHits(hits)))
    FN <- length(truth_type.gr) - TP
    
    precision <- ifelse((TP + FP) > 0, TP/(TP + FP), NA)
    recall    <- ifelse((TP + FN) > 0, TP/(TP + FN), NA)
    F1        <- ifelse(is.na(precision) | is.na(recall) | (precision+recall)==0,
                        NA, 2*(precision*recall)/(precision+recall))
    df <- data.frame(sample = s, TP, FP, FN,
                     precision, recall, F1, svtype=type)
    return(df)
  }))
  res[[s]] <- list(recall_raw,precision_raw,bysize)
}
cnvpytor_eval <- rbindlist(lapply(res,function(x){
  return(x[[3]])
}))
cnvpytor_eval_recall_mtx <- rbindlist(lapply(res,function(x){
  return(x[[1]])
}))
cnvpytor_eval_precision_mtx <- rbindlist(lapply(res,function(x){
  return(x[[2]])
}))

df_cnvpytor_recall <- cnvpytor_eval_recall_mtx%>%
  mutate(Size=cut(width,breaks = c(2000,10000,Inf)))%>%
  filter(!is.na(Size))%>%
  group_by(Size,svtype,sample)%>%
  summarise(Value=sum(predictor)/length(predictor),Size=unique(Size))%>%
  group_by(Size,svtype)%>%
  summarise(Size=unique(Size),mValue=mean(Value),ssd=sd(Value),count=n())%>%
  mutate(se=ssd/sqrt(count),
         lower_ci=mValue-qt(0.975,df=count-1)*se,
         upper_ci=mValue+qt(0.975,df=count-1)*se,
         Metrics="recall",
         caller="cnvpytor")%>%
  dplyr::select(-ssd,-count,-se)

df_cnvpytor_precision <- cnvpytor_eval_precision_mtx%>%
  mutate(Size=cut(width,breaks = c(2000,10000,Inf)))%>%
  filter(!is.na(Size))%>%
  group_by(Size,svtype,sample)%>%
  summarise(Value=sum(predictor)/length(predictor),Size=unique(Size))%>%
  group_by(Size,svtype)%>%
  summarise(Size=unique(Size),mValue=mean(Value),ssd=sd(Value),count=n())%>%
  mutate(se=ssd/sqrt(count),
         lower_ci=mValue-qt(0.975,df=count-1)*se,
         upper_ci=mValue+qt(0.975,df=count-1)*se,
         Metrics="precision",
         caller="cnvpytor")%>%
  dplyr::select(-ssd,-count,-se)

df_two_methods <- rbindlist(list(df_cnvpytor_recall,df_cnvpytor_precision,df_size_precision,df_size_precision_dup,df_size_recall,df_size_recall_dup),fill=TRUE)%>%
  filter(is.na(mapbility)|mapbility=="q30")%>%
  mutate(
    caller   = factor(caller),
    type     = factor(svtype, levels = c("DEL","DUP")),
    size_bin = factor(Size, levels = c("(2e+03,1e+04]", "(1e+04,Inf]"))
  )
library(ggplot2)
library(scales)
p_two_method <- df_two_methods%>%ggplot(., aes(x=size_bin,y=mValue,fill = caller)) +
  geom_col(position = "dodge2", width = 0.6, show.legend = T) +
  facet_grid(Metrics ~ type , scales = "free_y") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    x = NULL, y = NULL,
    title = "Caller comparison across size and CNV types",
    subtitle = "Metrics: Precision / Recall; Size bins: 2kb–10kb and >10kb"
  ) +
  geom_text(aes(label = percent(mValue, accuracy = 1)),
            position = position_dodge(width = 0.6),
            vjust = -0.3, size = 3) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )
ggsave(filename = "cnvpytor_vizcnv_comparison.pdf",device = "pdf",width = 7,height = 7,units = "in")
# ---------- summary across samples ----------
prec_ci <- mean_ci(perf$precision)
rec_ci  <- mean_ci(perf$recall)
f1_ci   <- mean_ci(perf$F1)

summary.df <- data.frame(
  metric = c("precision","recall","F1"),
  mean   = c(prec_ci["mean"], rec_ci["mean"], f1_ci["mean"]),
  lower  = c(prec_ci["lower"], rec_ci["lower"], f1_ci["lower"]),
  upper  = c(prec_ci["upper"], rec_ci["upper"], f1_ci["upper"])
)

fwrite(perf, paste0(out_prefix, "_per_sample.tsv"), sep = "\t")
fwrite(summary.df, paste0(out_prefix, "_summary.tsv"), sep = "\t")

message("✅ Results written to: ", out_prefix, "_*.tsv")
