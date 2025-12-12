source("./mod/mod_plot_output.R")
source("./mod/mod_dnCNV.R")
source("./mod/mod_hmzcnv.R")
source("./mod/mod_findCNV.R")
source("./mod/mod_upload.R")
source("./mod/mod_UCSC.R")
source("./helper/wg_plot.R")
source("./helper/misc.R")
source("./mod/mod_allele_imbalance.R")
server <- function(input, output,session) {
  # Reavtive Values --------------------------
  values <- reactiveValues()
  values$pr_sv <- data.frame(stringsAsFactors = F)
  values$m_sv <- data.frame(stringsAsFactors = F)
  values$f_sv <- data.frame(stringsAsFactors = F)
  values$pr_sv_fil <- data.frame(stringsAsFactors = F)
  values$m_sv_fil <- data.frame(stringsAsFactors = F)
  values$f_sv_fil <- data.frame(stringsAsFactors = F)
  values$pr_rd <- data.frame(stringsAsFactors = F)
  values$m_rd <- data.frame(stringsAsFactors = F)
  values$f_rd <- data.frame(stringsAsFactors = F)
  values$snp_gvcf_file_ref <- vector()
  values$snp_gvcf_file_path <- vector()
  values$snp_gvcf_file <- data.frame(stringsAsFactors = F)
  values$local_file_paths <- data.frame(file=c("CNV call file",
                                               "Proband read depth file",
                                               "Mom's read depth file",
                                               "Dad's read depth file",
                                               "Joint SNP vcf file"),datapath=rep("None",5),stringsAsFactors = F)
  values$selected_record <- data.frame(stringsAsFactors = F)
  
  values$ref_info <- data.frame(stringsAsFactors = F)
  values$anno_rect <- data.frame(stringsAsFactors = F)
  
  plots <- reactiveValues()
  plots$pr_rd <- data.frame(stringsAsFactors = F)
  plots$baf <- data.frame(stringsAsFactors = F)
  plots$m_rd <- data.frame(stringsAsFactors = F)
  plots$f_rd <- data.frame(stringsAsFactors = F)
  plots$pr_seg <- data.frame(stringsAsFactors = F)
  plots$m_seg <- data.frame(stringsAsFactors = F)
  plots$f_seg <- data.frame(stringsAsFactors = F)
  plots$snp_chr <- data.frame(stringsAsFactors = F)
  plots$xlabel <- character()
  plots$SNPcols <- vector(length = 3) ## placeholder for color in SNP plot
  
  names <- reactiveValues()
  names$pr_rd <- NULL
  names$m_rd <- NULL
  names$f_rd <- NULL
  
  
  volumes <- c(Home="~/Downloads/","R installation" = R.home(),shinyFiles::getVolumes()())
  mod_rd_upload_Server("pr_rd",volumes=volumes,values, names) 
  mod_rd_upload_Server("m_rd",volumes=volumes,values, names) 
  mod_rd_upload_Server("f_rd",volumes=volumes,values, names) 
  mod_sv_upload_Server("pr_sv",volumes=volumes,values) 
  mod_sv_upload_Server("m_sv",volumes=volumes,values) 
  mod_sv_upload_Server("f_sv",volumes=volumes,values) 
  mod_snp_upload_Server("snp_file",volumes=volumes,values)
  
  
  
  
  ## dynamic UI -------------
  
  output$blt_baf_seg <- renderUI({
    if (input$baf_seg)
      numericInput("target_spacing", "Minimum distance between adjacent SNV(<10kbp):", 150, min = 0, max = 10000)
  })
  
  
  output$blt_dnSNV_ui <- shiny::renderUI({
    if(length(values$snp_gvcf_file_path)!=0){
      shiny::tagList(
        p(HTML("<b>Show de novo SNV ?</b>"),
          span(shiny::icon("info-circle"), id = "info_nor"),
          checkboxGroupInput(inputId="include_dnSNV",label = NULL,c("Show"="TRUE"))))
    }
    else{
      return(NULL)
    }
  })
  output$ui_chkbox_RefSeq <- renderUI({
    if(!is.null(values$p1_file)){
      mod_checkbox_UI("RefSeq")
    }else{NULL}

  })
  output$ui_chkbox_IDR <- renderUI({
    if(!is.null(values$p2_file)){
      mod_checkbox_UI("IDR",value = F)
    }else{NULL}
  })
  output$ui_chkbox_SegDup <- renderUI({
    if(!is.null(values$p3_file)){
      mod_checkbox_UI("SegDup")
    }else{NULL}
  })
  output$ui_chkbox_OMIM <- renderUI({
    if(!is.null(values$p4_file)){
      mod_checkbox_UI("OMIM")
    }else{NULL}
  })
  output$ui_chkbox_gnomAD <- renderUI({
    if(!is.null(values$p5_file)){
      mod_checkbox_UI("gnomAD",value = F)
    }else{NULL}
  })
  output$ui_chkbox_RMSK <- renderUI({
    if(!is.null(values$p6_file)){
      mod_checkbox_UI("RMSK",value = F)
    }else{NULL}
  })

  ## ref file depending on genome build  ------
  observeEvent(input$ref,{
    if(input$ref=="hg38"){
      # reserve the blacklist for now Mar1, 2024
      #blacklist <- data.table::fread("GRCh38_unified_blacklist.bed.gz")%>%
      #  regioneR::toGRanges()
      values$ref_info <- data.table::fread("data/hg38.info.txt")
      values$gaps <- data.table::fread("data/hg38_gaps.bed")%>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)
      values$p1_file <-  "data/hg38_MANE.v1.0.refseq.parquet"
      values$p2_file <-  NULL
      values$p3_file <-  "data/hg38_ucsc_sugdups.parquet"
      values$p4_file <-  "data/OMIM_gene2_hg38_MANE_all.bed"
      values$p5_file <-  NULL
      values$p6_file <-  "data/hg38_rmsk.parquet"
      values$SegDup_merge <- "data/SegDup_hg38_UCSC_sorted_merged_1k.bed"
    }else if(input$ref=="hg19"){
     # blacklist <- data.table::fread("ENCFF001TDO.bed.gz")%>%
     #   regioneR::toGRanges()
      values$ref_info <- data.table::fread("data/hg19.info.txt")
      values$gaps <- data.table::fread("data/hg19_gaps.bed")%>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)
      values$p1_file <-  "data/NCBI_RefSeq_hg19_clean.bed.parquet"
      values$p2_file <-  "data/Claudia_hg19_MergedInvDirRpts_sorted.bed"
      values$p3_file <-  "data/hg19_ucsc_sugdups.parquet"
      values$p4_file <-  "data/OMIM_gene2_hg19_UCSC_all.bed"
      values$p5_file <-  "data/gnomAD_allSV_hg19_UCSC.bed"
      values$p6_file <-  "data/hg19_rmsk.parquet"
      values$SegDup_merge <- "data/SegDup_hg19_UCSC_no_PAR_merged_1k.bed"
    }else if(input$ref=="T2T"){
      values$ref_info <- data.table::fread("data/CHM13v2.0.info.txt")
      values$p1_file <-  NULL
      values$p2_file <-  NULL
      values$p3_file <-  NULL
      values$p4_file <-  NULL
      values$p5_file <-  NULL
      values$p6_file <-  NULL
    }
  })
  
  # button to filter range---------
 
  w <- waiter::Waiter$new(html = spin_3(), 
                          color = transparent(.5))
  observeEvent(input$btn_filter,{
    seg_option <- input$seg_option
    norm_option <- input$norm_options
    chr <- input$chr
    if(nrow(values$pr_rd)==0){return(NULL)
    }else{
      w$show()
      showNotification("Normalize and segment proband read depth", duration = 3, type = "message")
      pr_rd.gr <- normalization_method(values$pr_rd, chr, norm_option)%>%
        setDT()%>%
        setnames(.,c("V1","V2","V3"),c("chrom","start","end"))%>%
        makeGRangesFromDataFrame(keep.extra.columns = T)
      # remove regions that mapped to gaps
      ov <- findOverlaps(pr_rd.gr,values$gaps)
      if (isTRUE(input$mask_option)){
        pr_rd.gr <- pr_rd.gr[-queryHits(ov)]
      }
      plots$pr_rd <- pr_rd.gr%>%as.data.frame()%>%
        dplyr::select(-c("width","strand"))%>%
        setDT()%>%
        setnames(.,c("seqnames","start","end"),c("V1","V2","V3"))
      plots$pr_seg <- SegNormRD(plots$pr_rd,id="Proband",seg.method = seg_option)
      w$hide()
    }
    if(nrow(values$m_rd)==0){return(NULL)
    }else{
      w$show()
      showNotification("Normalize and segment mother read depth", duration = 3, type = "message")
      m_rd.gr <- normalization_method(values$m_rd, chr, norm_option)%>%
        setDT()%>%
        setnames(.,c("V1","V2","V3"),c("chrom","start","end"))%>%
        makeGRangesFromDataFrame(keep.extra.columns = T)
      # remove regions that mapped to gaps
      ov <- findOverlaps(m_rd.gr,values$gaps)
      if (isTRUE(input$mask_option)){
        m_rd.gr <- m_rd.gr[-queryHits(ov)]
      }
      plots$m_rd <- m_rd.gr%>%as.data.frame()%>%
        dplyr::select(-c("width","strand"))%>%
        setDT()%>%
        setnames(.,c("seqnames","start","end"),c("V1","V2","V3"))
      plots$m_seg <- SegNormRD(plots$m_rd,id="Mother",seg.method = seg_option)
      w$hide()
    }
    if(nrow(values$f_rd)==0){return(NULL)
    }else{
      w$show()
      showNotification("Normalize and segment father read depth", duration = 3, type = "message")
      f_rd.gr <- normalization_method(values$f_rd, chr, norm_option)%>%
        setDT()%>%
        setnames(.,c("V1","V2","V3"),c("chrom","start","end"))%>%
        makeGRangesFromDataFrame(keep.extra.columns = T)
      # remove regions that mapped to gaps
      ov <- findOverlaps(f_rd.gr,values$gaps)
      if (isTRUE(input$mask_option)){
        f_rd.gr <- f_rd.gr[-queryHits(ov)]
      }
      plots$f_rd <- f_rd.gr%>%as.data.frame()%>%
        dplyr::select(-c("width","strand"))%>%
        setDT()%>%
        setnames(.,c("seqnames","start","end"),c("V1","V2","V3"))
      plots$f_seg <- SegNormRD(plots$f_rd,id="Father",seg.method = seg_option)
      w$hide()
    }
  })
  observeEvent(input$btn_filter,{
    snp_gvcf_file_path=values$snp_gvcf_file_path
    ref_genome <- input$ref
    chr <- input$chr
    if(!chr%in%values$snp_gvcf_file_ref){
      if(chr%in%chrom_id){
        chr <- names(chrom_id)[which(chrom_id==chr)]
      }else{
        chr <- chrom_id[which(names(chrom_id)==chr)][[1]]
      }
    }
    if(length(snp_gvcf_file_path)==0){
      return(NULL)
    }else{
      w$show()
      showNotification("Filtering GATK VCF file", duration = 5, type = "message")
      loc.start <- 0
      loc.end <- values$ref_info%>%
        filter(chrom==chr | chrom == paste0("chr", chr))%>%
        dplyr::select(seqlengths)%>%unlist
      range.gr <- GenomicRanges::GRanges(chr,ranges = IRanges(loc.start,loc.end))
      if (isTRUE(input$mask_option)){
        if (chr %in% c(1:22, "X", "Y")){
          range.gr <- GenomeInfoDb::renameSeqlevels(range.gr, paste0("chr", GenomeInfoDb::seqlevels(range.gr)))
          range.gr <- unlist(GenomicRanges::subtract(range.gr, values$gaps))
          range.gr <- GenomeInfoDb::renameSeqlevels(range.gr, gsub("chr", "", GenomeInfoDb::seqlevels(range.gr)))
        } else {
          range.gr <- unlist(GenomicRanges::subtract(range.gr, values$gaps))
        }
      }
      plots$snp_chr <- ReadGVCF(snp_gvcf_file_path,ref_genome=input$ref,param = range.gr,target_spacing=input$target_spacing)%>%
        as.data.frame()%>%
        filter(pr_count>7,p1_count>7,p2_count>7,width<2) ## set the minimum for the read coverage
      InhFrom <- unique(plots$snp_chr$B_InhFrom)
      if(length(InhFrom)==3){
        names(plots$SNPcols) <- InhFrom
        plots$SNPcols[names(plots$SNPcols)!="Notphased"] <- SNPCOLOR2
        plots$SNPcols["Notphased"] <- c("#999999")
      }
      w$hide()
      
      anno_table_Server("gatk_table", plots$snp_chr, ranges, chr)
    }
  })
  
  
  # Basic Plots -----
  ### "Global" reactive values
  wg_ranges <- reactiveValues(x = NULL, pr = NULL, m = NULL, f = NULL, pr_ploidy=NULL)
  wg_dnCNV_table <- reactiveValues(t = data.frame(start_cum = c(0), end_cum = c(0), stringsAsFactors = F))
  wg_hmzCNV_table <- reactiveValues(t = data.frame(start_cum = c(0), end_cum = c(0), stringsAsFactors = F))
  
  ## reset plots upon changing chr
  observeEvent(input$btn_wg_rd, {
    wg_ranges$x <-  NULL
    wg_dnCNV_table$t <-  data.frame(start_cum = c(0), end_cum = c(0), stringsAsFactors = F)
  })
  
  ## WG Plot section
  observeEvent(input$btn_wg_rd, {
    req(nrow(values$pr_rd) != 0)
    chr_list <- paste0("chr", c(1:22,"X","Y"))
    showNotification("Normalize and segment proband read depth", duration = 3, type = "message")
    withProgress(message = 'Making plot', value = 0, {
      n=length(chr_list)
      tmplist <- lapply(chr_list, function(chr){
        pr_rd.gr <- normalization_method(values$pr_rd, chr, input$wg_norm_options)%>%
          setDT()%>%
          setnames(.,c("V1","V2","V3"),c("chrom","start","end"))%>%
          makeGRangesFromDataFrame(keep.extra.columns = T)
        # remove regions that mapped to gaps
        ov <- findOverlaps(pr_rd.gr,values$gaps)
        pr_rd <- pr_rd.gr[-queryHits(ov)]%>%as.data.frame()%>%
          dplyr::select(-c("width","strand"))%>%
          setDT()%>%
          setnames(.,c("seqnames","start","end"),c("V1","V2","V3"))
        pr_seg <- SegNormRD(pr_rd,id="Proband",seg.method = "slm")
        incProgress(1/n, detail = paste("Segment", chr))
        return(pr_seg)
      })
      wg_ranges$pr <- rbindlist(tmplist)
    })
    wg_pr <- wg_seg2plot(wg_ranges$pr,input$ref)
    wg_ranges$pr_ploidy <- wg_ranges$pr%>%
      group_by(chrom)%>%
      summarise(chromosome_counts=2*median(seg.mean))
    wg_ranges <- mod_plot_wg_Server("wg_pr_rd", wg_pr, wg_ranges, wg_dnCNV_table, wg_hmzCNV_table)
  })
  observeEvent(input$btn_wg_rd, {
    req(nrow(values$m_rd) != 0)
    chr_list <- paste0("chr", c(1:22,"X","Y"))
    showNotification("Normalize and segment mother's read depth", duration = 3, type = "message")
    withProgress(message = 'Making plot', value = 0, {
      n=length(chr_list)
      tmplist <- lapply(chr_list, function(chr){
        m_rd.gr <- normalization_method(values$m_rd, chr, input$wg_norm_options)%>%
          setDT()%>%
          setnames(.,c("V1","V2","V3"),c("chrom","start","end"))%>%
          makeGRangesFromDataFrame(keep.extra.columns = T)
        # remove regions that mapped to gaps
        ov <- findOverlaps(m_rd.gr,values$gaps)
        m_rd <- m_rd.gr[-queryHits(ov)]%>%as.data.frame()%>%
          dplyr::select(-c("width","strand"))%>%
          setDT()%>%
          setnames(.,c("seqnames","start","end"),c("V1","V2","V3"))
        m_seg <- SegNormRD(m_rd,id="Mother",seg.method = "slm")
        incProgress(1/n, detail = paste("Segment", chr))
        return(m_seg)
      })
      wg_ranges$m <- rbindlist(tmplist)
    })
    wg_m <- wg_seg2plot(wg_ranges$m,input$ref)
    wg_ranges <- mod_plot_wg_Server("wg_m_rd", wg_m, wg_ranges, wg_dnCNV_table,wg_hmzCNV_table)
  })
  observeEvent(input$btn_wg_rd, {
    req(nrow(values$f_rd) != 0)
    chr_list <- paste0("chr", c(1:22,"X","Y"))
    showNotification("Normalize and segment father's read depth", duration = 3, type = "message")
    withProgress(message = 'Making plot', value = 0, {
      n=length(chr_list)
      tmplist <- lapply(chr_list, function(chr){
        f_rd.gr <- normalization_method(values$f_rd, chr, input$wg_norm_options)%>%
          setDT()%>%
          setnames(.,c("V1","V2","V3"),c("chrom","start","end"))%>%
          makeGRangesFromDataFrame(keep.extra.columns = T)
        # remove regions that mapped to gaps
        ov <- findOverlaps(f_rd.gr,values$gaps)
        f_rd <- f_rd.gr[-queryHits(ov)]%>%as.data.frame()%>%
          dplyr::select(-c("width","strand"))%>%
          setDT()%>%
          setnames(.,c("seqnames","start","end"),c("V1","V2","V3"))
        f_seg <- SegNormRD(f_rd,id="Father",seg.method = "slm")
        incProgress(1/n, detail = paste("Segment", chr))
        return(f_seg)
      })
      wg_ranges$f <- rbindlist(tmplist)
    })
    wg_f <- wg_seg2plot(wg_ranges$f,input$ref)
    wg_ranges <- mod_plot_wg_Server("wg_f_rd", wg_f, wg_ranges, wg_dnCNV_table,wg_hmzCNV_table)
  })

  ## dnCNV table
  observeEvent(input$btn_wg_dnCNV, {
    req(!is.null(wg_ranges$pr))
    req(!is.null(wg_ranges$m))
    req(!is.null(wg_ranges$f))
    seg_data <- mod_dnCNV_Server("wg_dnCNV",wg_ranges$pr, wg_ranges$m, wg_ranges$f)
    names(seg_data)[1] = "chrom"
    print(seg_data)
    temp <- wg_ranges$pr %>% 
      group_by(chrom) %>% 
      summarise(max_end = max(loc.end)) %>% 
      mutate(across("chrom", str_replace, "chr", "")) %>% 
      arrange(as.numeric(chrom)) %>% 
      mutate(loc_add = lag(cumsum(as.numeric(max_end)), default = 0)) %>% 
      mutate(chrom = paste0("chr", chrom))
    seg_data <- seg_data %>% 
      inner_join(temp, by = "chrom") %>% 
      mutate(end_cum = loc_add + end) 
    seg_data <- seg_data %>% 
      mutate(start_cum = end_cum- width)
    wg_dnCNV_table$t <- seg_data
  })
  

  
  
  ## Plots section
  ranges <- reactiveValues(x = NULL, y = NULL,cur = NULL, click = NULL, pb = NULL)
  dnCNV_table <- reactiveValues(t = data.frame(start = c(0), end = c(0), stringsAsFactors = F),
                                hl = data.frame(start = c(0), end = c(0)), 
                                hl_col = c("white"))
  hmzCNV_table <- reactiveValues(t = data.frame(start = c(0), end = c(0), stringsAsFactors = F),
                                 hl = data.frame(start = c(0), end = c(0)), 
                                 hl_col = c("white"))
  findCNV_table <- reactiveValues(t = data.frame(start = c(0), end = c(0), stringsAsFactors = F),
                                 hl = data.frame(start = c(0), end = c(0)), 
                                 hl_col = c("white"))
  
  ## reset plots upon changing chr
  observeEvent(input$btn_plot, {
    ranges$x <-  NULL
    ranges$cur <-  NULL
    ranges$click <-  NULL
    ranges$pb <-  NULL
    dnCNV_table$t <-  data.frame(start = c(0), end = c(0), stringsAsFactors = F)
    hmzCNV_table$t <- data.frame(start = c(0), end = c(0), stringsAsFactors = F)
  })

  ## dynamic highlight
  mod_col_pick_Server("highlight", dnCNV_table, ranges)
  mod_col_pick_Server("highlight", hmzCNV_table, ranges)
  
  ## RD plots
  observeEvent(input$btn_plot,{
    req(nrow(plots$pr_rd) != 0)
    
    showNotification("Plotting Read Depth Plot", duration = 12, type = "message")
    include_seg <- input$include_seg
    df <- rbindlist(list(plots$pr_seg,plots$m_seg,plots$f_seg))%>%
      filter(ID%in%include_seg)%>%
      mutate(ID=as.factor(ID))%>%
      mutate(seg.mean=ifelse(seg.mean < -2.5,-2,seg.mean))
    plots$xlabel=input$chr
    rd <- ggplot(plots$pr_rd, aes(x=V2, y=log2(ratio+0.00001))) +
      geom_point(shape=".")+
      #scattermore::geom_scattermore(shape=".",pixels=c(1024,1024))+
      geom_point(data = subset(plots$pr_rd, ratio < 0.7),aes(V2,log2(ratio+0.00001)),shape=".",color="green")+
      geom_point(data = subset(plots$pr_rd, ratio > 1.3),aes(V2,log2(ratio+0.00001)),shape=".",color="red")+
      geom_segment(data = df,aes(x=loc.start,y=seg.mean,xend=loc.end,yend=seg.mean,color=factor(ID)),size=1)+
      scale_color_manual(name="Segment",values = c("Proband"="blue","Mother"="#E69F00","Father"="#39918C"),)+
      xlab(plots$xlabel)+
      scale_rd+
      style_rd
    
 
    btnValrds <- mod_checkbox_Server("RD-static")
    ranges <- mod_plot_switch_Server("RD-static", btnValrds$box_state, rd, ranges, dnCNV_table,hmzCNV_table, zoom= F)
    btnValrdd <- mod_checkbox_Server("RD-dynamic")
    ranges <- mod_plot_switch_Server("RD-dynamic", btnValrdd$box_state, rd, ranges, dnCNV_table,hmzCNV_table)

  })
  
  #Baf-B plot
  observeEvent(input$btn_plot,{
    
    req(nrow(plots$snp_chr) != 0)
    
    noti_id <- showNotification("Plotting B-allele frequency plots", type = "message", duration = NULL)
    df <- plots$snp_chr%>%filter(likelyDN%in%c(input$include_dnSNV,"FALSE"))
    matsnp.df <-  plots$snp_chr%>%filter(!is.na(P1_phased_BAF))
    patsnp.df <- plots$snp_chr%>%filter(!is.na(P2_phased_BAF))
    bafseg.obj <- DNAcopy::CNA(plots$snp_chr$pr_absBAF, plots$snp_chr$seqnames,plots$snp_chr$start, data.type = "binary",sampleid = "Index")%>%
      DNAcopy::segment(verbose = 1)
    matseg.obj = DNAcopy::CNA(matsnp.df$P1_phased_BAF, matsnp.df$seqnames,matsnp.df$start, data.type = "binary",sampleid = "Index")%>%
      DNAcopy::segment(verbose = 1)
    patseg.obj = DNAcopy::CNA(patsnp.df$P2_phased_BAF, patsnp.df$seqnames,patsnp.df$start, data.type = "binary",sampleid = "Index")%>%
      DNAcopy::segment(verbose = 1)
    snp_out <- as.data.table(bafseg.obj$output)%>%
      mutate(len=round(as.numeric(loc.end)-as.numeric(loc.start),0),
             loc.end=as.numeric(loc.end),
             loc.start=as.numeric(loc.start),
             seg.mean=seg.mean+0.55)%>%
      filter(len>minseg,seg.mean>minsegmean)
    
    mat_out <- as.data.table(matseg.obj$output)%>%
      mutate(len=round(as.numeric(loc.end)-as.numeric(loc.start),0),
             loc.end=as.numeric(loc.end),
             loc.start=as.numeric(loc.start))%>%
      pivot_longer(cols = c("loc.start","loc.end"),names_to = "source",values_to = "pos")
    pat_out <- as.data.table(patseg.obj$output)%>%
      mutate(len=round(as.numeric(loc.end)-as.numeric(loc.start),0),
             loc.end=as.numeric(loc.end),
             loc.start=as.numeric(loc.start))%>%
      pivot_longer(cols = c("loc.start","loc.end"),names_to = "source",values_to = "pos")
    
    cols <- plots$SNPcols
    xlabel=unique(df$chrom)[1]
    
    infromb <- df$B_InhFrom %>% 
      unique()
    snpb_cols <- vector(length = length(infromb))
    for (i in 1:length(snpb_cols)){
      snpb_cols[i] <- df$B_col[which(df$B_InhFrom == infromb[i])[1]]
      names(snpb_cols)[i] <- infromb[i]
    }
    
    #infroma <- df$A_InhFrom %>% 
    #  unique()
    #snpa_cols <- vector(length = length(infroma))
    #for (i in 1:length(snpa_cols)){
    #  snpa_cols[i] <- df$A_col[which(df$A_InhFrom == infroma[i])[1]]
    #  names(snpa_cols)[i] <- infroma[i]
    #}
    
        
    # snp_a <- ggplot(df, aes(x=start,y=pr_ALT_Freq,col=A_InhFrom))+
    #   geom_point(shape=20, size = 1.5)+
    #   #scattermore::geom_scattermore(shape=".",pixels=c(1024,1024))+
    #   geom_point(data = subset(df, likelyDN %in%c("TRUE")),size = 2,shape=8,color="red")+
    #   scale_fill_manual("LikelyDN",limits=c("dnSNV"),values = "red")+
    #   xlab(xlabel)+
    #   scale_snp+
    #   style_snp+
    #   scale_colour_manual(values = snpa_cols)+
    #   guides(color = guide_legend(override.aes = list(size = 4)))+
    #   scale_x_continuous(labels = scales::label_number())
    # 
    # btnVala <- mod_checkbox_Server("Baf-A_allele")
    # ranges <- mod_plot_switch_Server("Baf-A_allele", btnVala$box_state, snp_a, ranges, dnCNV_table,hmzCNV_table)
    # 
    snp_b <- ggplot(df, aes(x=start,y=pr_ALT_Freq,col=B_InhFrom))+
      geom_point(shape=20, size = 1.5)+
     # shape = . will significant improve the plot speed
      #geom_point(shape=".",alpha=1)+
      geom_point(data = subset(df, likelyDN %in%c("TRUE")),size = 2,shape=8,color="red")+
      # geom_line(data=pat_out,aes(x=pos,y=seg.mean),col="#39918C",size=1.5)+
      # geom_line(data=mat_out,aes(x=pos,y=seg.mean),col="#E69F00",size=1.5)+
      scale_fill_manual("LikelyDN",limits=c("dnSNV"),values = "red")+
      guides(color = guide_legend(override.aes = list(size = 4)))+
      xlab(xlabel)+
      scale_snp+
      style_snp+
      scale_colour_manual(values = snpb_cols)+
      scale_x_continuous(labels = scales::label_number())
    
    if (!input$is_trio){
      snp_b <- snp_b+
        scale_color_manual(values = c("#999999", "#999999", "#999999","#999999", "#999999"))+
        theme(legend.position = "none")
    }
    
    if (input$baf_seg && input$is_trio){
      snp_b <- snp_b +
        geom_segment(data=snp_out,aes(x=loc.start,xend=loc.end,y=seg.mean,yend=seg.mean),col="purple",size=1.5)+
        geom_line(data=pat_out,aes(x=pos,y=seg.mean),col="#39918C",size=1.5)+
        geom_line(data=mat_out,aes(x=pos,y=seg.mean),col="#E69F00",size=1.5)
    }
    
    removeNotification(noti_id)
    btnValb <- mod_checkbox_Server("Baf-B_allele")
    ranges <- mod_plot_switch_Server("Baf-B_allele", btnValb$box_state, snp_b, ranges, dnCNV_table,hmzCNV_table)
    
  })
  
  ##Anno tracks-----
  
  ## switch chr of sv table 
  observeEvent(input$btn_filter,{
    chr <- input$chr
    if(nrow(values$pr_sv) == 0){
      return(NULL)
    }else{
      values$pr_sv_fil <- values$pr_sv%>%filter(CHROM==chr)
      return(values$pr_sv_fil)
    }
  })
  observeEvent(input$btn_filter, {
    chr <- input$chr
    if(nrow(values$m_sv) == 0){
      return(NULL)
    }else{
      values$m_sv_fil <- values$m_sv%>%filter(CHROM==chr)
      return(values$m_sv_fil)
    }
  })
  observeEvent(input$btn_filter, {
    chr <- input$chr
    if(nrow(values$f_sv) == 0){
      return(NULL)
    }else{
      values$f_sv_fil <- values$f_sv%>%filter(CHROM==chr)
      return(values$f_sv_fil)
    }
  })
  
 
  ## create plots
  observeEvent(input$btn_anno,{
    id <- showNotification("Pulling Data", type = "message", duration = NULL)
    chrn = input$chr
    path = "./data/"
    
    loc.end <- values$ref_info%>%
      filter(chrom==chrn | chrom == paste0("chr", chrn))%>%
      dplyr::select(seqlengths)%>%unlist
    if (nrow(values$pr_sv_fil) != 0){
      pr_sv <- process_sv(values$pr_sv_fil)
      pr_sv <- pr_sv %>% 
        add_row(CHROM = "dummy", start = 0, end = loc.end, color = "white")
      pr_sv_plot <- ggplot(pr_sv, aes(x = start, y = idx)) +
        annotate("rect", xmin = pr_sv$start, xmax = pr_sv$end, ymin = pr_sv$idx, ymax = pr_sv$idx+0.0001, color = pr_sv$color)+
        style_anno+
        scale_anno+
        ylab("pr_sv")
      btnVal_pr_sv <- mod_checkbox_Server("pr_sv")
      ranges <- mod_plot_switch_Server("pr_sv", btnVal_pr_sv$box_state, pr_sv_plot, ranges, dnCNV_table, hmzCNV_table)
      anno_table_Server("pr_sv", pr_sv, ranges, chrn)
    }
    if (nrow(values$m_sv_fil) != 0){
      m_sv <- process_sv(values$m_sv_fil)
      m_sv <- m_sv %>% 
        add_row(CHROM = "dummy", start = 0, end = loc.end, color = "white")
      m_sv_plot <- ggplot(m_sv, aes(x = POS, y = idx)) +
        annotate("rect", xmin = m_sv$start, xmax = m_sv$end, ymin = m_sv$idx, ymax = m_sv$idx+0.0001, color = m_sv$color)+
        style_anno+
        scale_anno+
        ylab("m_sv")
      btnVal_m_sv <- mod_checkbox_Server("m_sv")
      ranges <- mod_plot_switch_Server("m_sv", btnVal_m_sv$box_state, m_sv_plot, ranges, dnCNV_table, hmzCNV_table)
      anno_table_Server("m_sv", m_sv, ranges, chrn)
    }
    if (nrow(values$f_sv_fil) != 0){
      f_sv <- process_sv(values$f_sv_fil)
      f_sv <- f_sv %>% 
        add_row(CHROM = "dummy", start = 0, end = loc.end, color = "white")
      f_sv_plot <- ggplot(f_sv, aes(x = POS, y = idx)) +
        annotate("rect", xmin = f_sv$start, xmax = f_sv$end, ymin = f_sv$idx, ymax = f_sv$idx+0.0001, color = f_sv$color)+
        style_anno+
        scale_anno+
        ylab("f_sv")
      btnVal_f_sv <- mod_checkbox_Server("f_sv")
      ranges <- mod_plot_switch_Server("f_sv", btnVal_f_sv$box_state, f_sv_plot, ranges, dnCNV_table, hmzCNV_table)
      anno_table_Server("f_sv", f_sv, ranges, chrn)
    }
    
    if(!is.null(values$p1_file)){
      RefSeq_data <- read_parquet(values$p1_file,as_data_frame = F)
      RefSeq <- RefSeq_data %>% 
        filter(seqname ==  chrn,type=="exon") %>% 
        collect() %>%
        mutate(gene_num=round(as.numeric(as.factor(gene_id)),3)/100,
               strand=as.factor(strand))
      gene_x <- RefSeq%>%
        group_by(gene_num)%>%
        summarise(start=(min(start)+max(end))/2,
                  end=(min(start)+max(end))/2,
                  gene_id=unique(gene_id))
      
      p1 <- RefSeq %>% 
        ggplot(aes(xstart = start, xend = end, y = gene_num))+
        #geom_range(aes(fill = strand),height = 0.02) +
        geom_intron(data = to_intron(RefSeq, "gene_num"), arrow.min.intron.length = 1000)+
        geom_text(data=gene_x,aes(x=start,label=gene_id),vjust = -1.2,check_overlap = T,fontface="italic")+
        ylab("RefSeq")+
        style_anno+
        scale_genes+
        scale_fill_manual(values = c("+"="#E69F00","-"="#39918C"))+
        scale_x_continuous(labels = scales::label_number())
      anno_table_Server("RefSeq", RefSeq, ranges, chrn)
    }
    if(!is.null(values$p2_file)){
      IDR<-data.table::fread(values$p2_file) %>% 
        dplyr::select(V1,V2,V3,V4,V5,V6) %>% 
        dplyr::rename("chrom" = V1, "start" = V2, "end" = V3, "name" = V4, "score" = V5, "strand" = V6) %>%
        dplyr::filter(chrom == chrn)
      IDR <- IDR %>%
        mutate(idx = sample(1:1, size = dim(IDR)[1], replace = T)/1000)
      IDR_pos <- IDR %>%
        filter(strand == "+")
      IDR_neg <- IDR %>%
        filter(strand == "-")
      p2 <- ggplot(IDR, aes(x = start, y = idx)) +
        geom_segment(data = IDR_pos, aes(x = start, y = idx, xend = end, yend = idx), arrow = arrow(length = unit(0.05, "inches")), color = "dodgerblue")+
        geom_segment(data = IDR_neg, aes(x = end, y = idx, xend = start, yend = idx), arrow = arrow(length = unit(0.05, "inches")), color = "dodgerblue3")+
        style_anno+
        ylab("IDR")
      anno_table_Server("IDR", IDR, ranges, chrn)
    }
    if(!is.null(values$p3_file)){
      SegDup_data <- arrow::read_parquet( values$p3_file,as_data_frame = F)
      SegDup <- SegDup_data %>%
        dplyr::select(all_of(segdups.keep.col))%>%
        filter(chrom ==chrn)%>%
        dplyr::rename("start" = "chromStart", "end" = "chromEnd") %>%
        collect()
        
      SegDup <- SegDup %>%
        mutate(idx = sample(1:100, size = dim(SegDup)[1], replace = T)/1000) %>%
        mutate(color = case_when(level == "SD_low" ~ "gray40",
                                 level == "SD_mid" ~ "yellow2",
                                 level == "SD_high" ~ "darkorange"))
      SegDup_pos <- SegDup %>%
        filter(strand == "+")
      SegDup_neg <- SegDup %>%
        filter(strand == "-")
      p3 <- ggplot(SegDup, aes(x = start, y = idx)) +
        geom_segment(data = SegDup_pos, aes(x = start, y = idx, xend = end, yend = idx), arrow = arrow(length = unit(0.05, "inches")), color = SegDup_pos$color)+
        geom_segment(data = SegDup_neg, aes(x = end, y = idx, xend = start, yend = idx), arrow = arrow(length = unit(0.05, "inches")), color = SegDup_neg$color)+
        style_anno+
        scale_anno+
        ylab("SegDup")
      anno_table_Server("SegDup", SegDup, ranges, chrn)
    }
    if(!is.null(values$p4_file)){
      OMIM <- data.table::fread(values$p4_file)
      OMIM <- OMIM %>%
        filter(chrom == chrn) %>%
        mutate_at(vars(pheno_key), as.factor)
      OMIM <- OMIM %>%
        mutate(idx = sample(1:100, size = dim(OMIM)[1], replace = T)/1000) %>%
        mutate(color = case_when(pheno_key == "0" ~ "grey",
                                 pheno_key == "1" ~ "lightgreen",
                                 pheno_key == "2" ~ "green3",
                                 pheno_key == "3" ~ "green4",
                                 pheno_key == "4" ~ "purple", ))
      
      OMIM_label <- OMIM %>%
        filter(pheno_key  %in% c("3","4"))
      p4 <- ggplot(OMIM, aes(x = start, y = idx)) +
        annotate("rect", xmin = OMIM$start, xmax = OMIM$end, ymin = OMIM$idx, ymax = OMIM$idx+0.0001, color = OMIM$color)+
        annotate("text", x = OMIM_label$start, y = OMIM_label$idx, label = OMIM_label$gene_symbol, size = 3.5, hjust = 1.1)+
        style_anno+
        scale_anno+
        ylab("OMIM")
      anno_table_Server("OMIM", OMIM, ranges, chrn)
    }
    if(!is.null(values$p5_file)){
      gnomAD <- data.table::fread(values$p5_file)
      gnomAD <- gnomAD %>%
        dplyr::rename("chrom" = "#chrom", "start" = "chromStart", "end" = "chromEnd") %>%
        filter(chrom == chrn)
      gnomAD <- gnomAD %>%
        mutate(idx = sample(1:100, size = dim(gnomAD)[1], replace = T)/1000) %>%
        mutate(color = case_when(svtype == "BND" ~ "grey",
                                 svtype == "DEL" ~ "red",
                                 svtype == "DUP" ~ "blue",
                                 svtype == "INS" ~ "darkorange",
                                 TRUE ~ "cyan4"))
      p5 <- ggplot(gnomAD, aes(x = start, y = idx)) +
        annotate("rect", xmin = gnomAD$start, xmax = gnomAD$end, ymin = gnomAD$idx, ymax = gnomAD$idx+0.0001, color = gnomAD$color)+
        style_anno+
        scale_anno+
        ylab("gnomAD")
      anno_table_Server("gnomAD", gnomAD, ranges, chrn)
      
    }
    if(!is.null(values$p6_file)){
      rmsk_data <- read_parquet(values$p6_file,as_data_frame = F)
      rmsk <- rmsk_data %>% 
        filter(chrom ==  chrn) %>% 
        collect()%>%
        dplyr::rename("start" = "chromStart", "end" = "chromEnd")
      
      rmsk <- rmsk %>% 
        mutate(idx = case_when(repClass == "SINE" ~ 0.014*7,
                               repClass == "LINE" ~ 0.014*6,
                               repClass == "LTR" ~ 0.014*5,
                               repClass == "DNA" ~ 0.014*4,
                               repClass == "Simple_repeat" ~ 0.014*3,
                               repClass == "Low_complexity" ~ 0.014*2,
                               TRUE ~ 0.014))
      rmsk_pos <- rmsk %>%
        filter(strand == "+")
      rmsk_neg <- rmsk %>%
        filter(strand == "-")
      print(rmsk_pos)
      print(rmsk_neg)
      p6 <- ggplot(rmsk, aes(x = start, y = idx)) + 
        geom_segment(data = rmsk_pos, aes(x = start, y = idx, xend = end, yend = idx), arrow = arrow(length = unit(0.05, "inches")), color = "black")+
        geom_segment(data = rmsk_neg, aes(x = end, y = idx, xend = start, yend = idx), arrow = arrow(length = unit(0.05, "inches")), color = "black")+
        style_anno+
        scale_anno+
        ylab("RMSK")
      anno_table_Server("rmsk", rmsk, ranges, chrn)
    }
    
    
    showNotification("Annotating", type = "message", duration = 8)
    btnVal1 <- mod_checkbox_Server("RefSeq")
    btnVal2 <- mod_checkbox_Server("IDR")
    btnVal3 <- mod_checkbox_Server("SegDup")
    btnVal4 <- mod_checkbox_Server("OMIM")
    btnVal5 <- mod_checkbox_Server("gnomAD")
    btnVal6 <- mod_checkbox_Server("RMSK")
    ranges <- mod_plot_switch_Server("RefSeq", btnVal1$box_state, p1, ranges, dnCNV_table, hmzCNV_table)
    ranges <- mod_plot_switch_Server("IDR", btnVal2$box_state, p2, ranges, dnCNV_table, hmzCNV_table)
    ranges <- mod_plot_switch_Server("Segdup", btnVal3$box_state, p3, ranges, dnCNV_table, hmzCNV_table)
    ranges <- mod_plot_switch_Server("OMIM", btnVal4$box_state, p4, ranges, dnCNV_table, hmzCNV_table)
    ranges <- mod_plot_switch_Server("gnomAD", btnVal5$box_state, p5, ranges, dnCNV_table, hmzCNV_table)
    ranges <- mod_plot_switch_Server("RMSK", btnVal6$box_state, p6, ranges, dnCNV_table, hmzCNV_table)
  
    removeNotification(id)
  })
  

  ## Misc functionalities----
  ## find CNV table
  observeEvent(input$btn_findCNV, {
    req(nrow(values$pr_rd)!=0)
    path <- "./data/"
    if(!is.null(values$SegDup_merge)){
      SegDup_merge <- data.table::fread(values$SegDup_merge)
      names(SegDup_merge) <- c("chrom", "start", "end")
    }
    # if(!is.null(values$p1_file)){
    #   getmode <- function(v) {
    #     uniqv <- unique(v)
    #     uniqv[which.max(tabulate(match(v, uniqv)))]
    #   }
    #   RefSeq_data <- read_parquet(paste0(path,values$p1_file))
    #   RefSeq_gr <- RefSeq_data %>%
    #     group_by(gene_id) %>%
    #     dplyr::summarise(chr = getmode(seqname),
    #               start = min(start),
    #               end = max(end)) %>%
    #     dplyr::select(2,3,4,1) %>%
    #     makeGRangesFromDataFrame(keep.extra.columns = T)
      
    # }
    if(!is.null(values$p4_file)){
      OMIM <- data.table::fread(values$p4_file)
      OMIM <- OMIM %>%
        mutate_at(vars(pheno_key), as.factor)
      OMIM <- OMIM %>%
        mutate(idx = sample(1:100, size = dim(OMIM)[1], replace = T)/1000) %>%
        mutate(color = case_when(pheno_key == "0" ~ "grey",
                                 pheno_key == "1" ~ "lightgreen",
                                 pheno_key == "2" ~ "green3",
                                 pheno_key == "3" ~ "green4",
                                 pheno_key == "4" ~ "purple", ))

      OMIM_label <- OMIM %>%
        filter(pheno_key  %in% c("3","4"))
    }
    w$show()
    if((nrow(values$m_rd)==0) | (nrow(values$f_rd) == 0)){
    	findCNV_table$t <- mod_findCNV_pr_only_Server("findCNV",values$pr_rd, SegDup_merge, RefSeq_gr, OMIM)
    } else {
    findCNV_table$t <- mod_findCNV_Server("findCNV",values$pr_rd, values$m_rd, values$f_rd, SegDup_merge, RefSeq_gr, OMIM)
    }
    updateTabItems(session, "tabs", "table")
    w$hide()
  })
  
  ## Show current ranges
  observe({
    output$cur_range <- renderText({
      req(!is.null(ranges$pb))
      paste0("range: ", round(ranges$pb[1]), "-", round(ranges$pb[2]), "    width: ", round(ranges$pb[2])-round(ranges$pb[1])+1)
    })
  })
  
  ## Show cur location
  observe({
    output$cur_loc <- renderText({
      req(!is.null(ranges$cur))
      paste0("location: ", round(ranges$cur))
    })
  })
  
  ## copy clicked location
  # observe({
  #   if(is.null(ranges$click)){ranges$click <- 0}
  #   clipr::write_clip(str_remove_all(as.character(round(ranges$click)), "[\r\n]"), allow_non_interactive = T)
  # })

  ## btn_goto
  observeEvent(input$btn_go,{
    req(!is.null(input$goto_reg))
    str <- stringr::str_trim(input$goto_reg)
    str <- strsplit(str,"-|_")
    if (length(str[[1]]) == 1) {
      if (!is.na(as.numeric(str[[1]]))){
        showNotification("Jumping to coordinates", type = "message")
        from <- as.numeric(str[[1]])-500000
        if (from < 0){from <- 0}
        to <- as.numeric(str[[1]])+500000
        ranges$x <- c(from, to)
      }else {
        showNotification("Looking up gene name", type = "message")
        path <- "./data/"
        if (input$ref=="hg38"){
          p1_file <- "hg38_MANE.v1.0.refseq.parquet"
        } else if (input$ref == "hg19"){
          p1_file <- "NCBI_RefSeq_hg19_clean.bed.parquet"
        }
        
        RefSeq <- read_parquet(paste0(path,p1_file))
        search <- as.character(str[[1]])
        found <- RefSeq %>% 
          filter(seqname == input$chr) %>% 
          filter(grepl(search, gene_id, ignore.case = T))
        if (nrow(found) != 0){
          print(ranges$x)
          from <- min(found$start)-500000
          if (from < 0){
            from <- 0
          }
          to <- max(found$end)+500000
          ranges$x <- c(as.numeric(from),
                        as.numeric(to))
        }else{
          showNotification("Gene not found", type = "error")
        }

      }
    } else if (length(str[[1]]) == 2){
      showNotification("Jumping to coordinates", type = "message")
      from <- as.numeric(str[[1]][1])
      to <- as.numeric(str[[1]][2])
      ranges$x <- c(from, to)
    }
  }, ignoreInit = T)
  
  ##UCSC
  observe({
    mod_UCSC_Server("UCSC", input$ref, input$chr, ranges)
  })



  ## Show cur input
  observe({
    output$pr_name <- renderText({
      req(!is.null(names$pr_rd))
      fname <- unlist(strsplit(names$pr_rd, "\\."))[1]
      paste0("Proband: ", fname)
    })
    output$pr_name2 <- renderText({
      req(!is.null(names$pr_rd))
      fname <- unlist(strsplit(names$pr_rd, "\\."))[1]
      paste0("Proband: ", fname)
    })
  })
  observe({
    output$m_name <- renderText({
      req(!is.null(names$m_rd))
      fname <- unlist(strsplit(names$m_rd, "\\."))[1]
      paste0("Mother: ", fname)
    })
    output$m_name2 <- renderText({
      req(!is.null(names$m_rd))
      fname <- unlist(strsplit(names$m_rd, "\\."))[1]
      paste0("Mother: ", fname)
    })
  })
  observe({
    output$f_name <- renderText({
      req(!is.null(names$f_rd))
      fname <- unlist(strsplit(names$f_rd, "\\."))[1]
      paste0("Father: ", fname)
    })
    output$f_name2 <- renderText({
      req(!is.null(names$f_rd))
      fname <- unlist(strsplit(names$f_rd, "\\."))[1]
      paste0("Father: ", fname)
    })
  })
  onStop(function() {
    rm(list=ls())
    cat("Session stopped\n")
    })
}
