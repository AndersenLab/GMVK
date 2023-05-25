library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(showtext)
library(ape)
library(ggrepel)

plotNormGeneFeatures <- function(sub_gff) {
  
  showtext::showtext_auto()
  
  GFFbyL2 <- sub_gff %>% 
    dplyr::arrange(source,desc(mAlias)) %>%
    dplyr::group_by(desc(source),desc(mAlias),desc(L2)) %>%
    dplyr::mutate(ngroup=cur_group_id()) %>%
    dplyr::ungroup()
  
  neg_str <- GFFbyL2 %>% dplyr::filter(strand=="-")
  neg_L3 <- neg_str %>% 
    dplyr::filter(type=="exon" | type =="intron" | type =="CDS") %>%
    dplyr::group_by(L2) %>%
    dplyr::mutate(start=abs(start-(max(end)))) %>%
    dplyr::mutate(end=abs(end-max(end))) %>%
    dplyr::mutate(temp=start) %>%
    dplyr::mutate(start=end) %>%
    dplyr::mutate(end=temp) %>%
    dplyr::select(-temp) %>%
    dplyr::ungroup()
  neg_L1L2 <- neg_str %>% dplyr::filter(!(type=="exon") & !(type =="intron") & !(type =="CDS"))
  neg_corr <- rbind(neg_L1L2,neg_L3) %>% dplyr::arrange(start)
  
  pos_str <- GFFbyL2 %>% dplyr::filter(strand=="+")
  
  grouped_gff <- rbind(neg_corr,pos_str)
  
  CDS_start <- grouped_gff %>% 
    dplyr::group_by(L2) %>%
    dplyr::filter(type=="CDS") %>% 
    dplyr::filter(start==(min(start))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(start==max(start))
  
  maxCDS <- CDS_start %>% dplyr::select(start,ngroup)
  
  reference_tran <- grouped_gff %>% dplyr::filter(ngroup==maxCDS$ngroup)
  sub_tran <- grouped_gff %>% 
    dplyr::filter(!(ngroup==maxCDS$ngroup)) %>%
    dplyr::group_by(L2) %>%
    dplyr::mutate(CDSstr=ifelse(type=="CDS",start,1e9)) %>%
    dplyr::mutate(hjust_cds=min(CDSstr)) %>%
    dplyr::mutate(shift=maxCDS$start-hjust_cds) %>%
    dplyr::mutate(start=start+(maxCDS$start-hjust_cds)) %>%
    dplyr::mutate(end=end+(maxCDS$start-hjust_cds)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-CDSstr,-hjust_cds,-shift)
  
  plottable_df <- rbind(reference_tran,sub_tran)
  
  geneL <- plottable_df %>%
    dplyr::filter(type=="gene") %>%
    dplyr::mutate(len=end-start)
  
  span <- max(geneL$len)
  
  restE <- plottable_df %>% 
    dplyr::filter(type=="exon")
  
  lastE <- restE %>% dplyr::group_by(L2) %>%
    dplyr::filter(end==max(end))  %>%
    dplyr::ungroup() %>%
    dplyr::select(L2,end) 
    
  
  labels <- grouped_gff %>%
    dplyr::filter(type=="gene") %>%
    dplyr::select(L2,mAlias,ngroup,strand) %>%
    dplyr::group_by(mAlias) %>%
    dplyr::mutate(iso=1:n()) %>%
    #dplyr::rowwise() %>%
    dplyr::mutate(verif=ifelse(any(iso > 1),"Y","N")) %>% #mAlias," (",as.character(iso),")"), mAlias)) 
    dplyr::ungroup() %>%
    dplyr::mutate(label=ifelse(verif=="Y",paste0(mAlias," (",L2,")"),mAlias)) %>%
    dplyr::left_join(lastE, by="L2")
  
  cdsE <- plottable_df %>% 
    dplyr::filter(type=="CDS")
  
  ggplot() + geom_rect(data = restE, aes(xmin = start/1e3 ,xmax = end/1e3 ,ymin=ngroup-0.25 , ymax=ngroup+0.25),fill="grey") +
    geom_rect(data = cdsE, aes(xmin = start/1e3 ,xmax = end/1e3 ,ymin=ngroup-0.25 , ymax=ngroup+0.25,fill=source)) +
    geom_segment(data = plottable_df %>% dplyr::filter(type=="intron"), aes(x=start/1e3,xend=(start+((end-start)/2))/1e3,y=ngroup,yend=ngroup+0.25),color="darkgrey") +
    geom_segment(data = plottable_df %>% dplyr::filter(type=="intron"), aes(x=(start+((end-start)/2))/1e3,xend=end/1e3,y=ngroup+0.25,yend=ngroup),color="darkgrey") +
    annotate('text',x=0,y=labels$ngroup+0.45,label=paste0((labels$label)),parse=F,fontface='italic',hjust=0) +
    annotate('text',x=labels$end/1e3+0.015*(span/1e3),y=labels$ngroup,label=paste0("(",labels$strand,")"),parse=F,hjust=0) +
    theme(axis.title.y  = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_line(),
          panel.background = element_blank(),
          legend.position=c(0.95,0.10),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0,size = 12)) + xlab("Transcript length (kb)") + scale_fill_manual(values=c("#e69f00","#cc799d","#009e73"))
}

tair <-read_tsv("/projects/b1059/projects/Nicolas/shiny_gmvt/data/master_DB.tsv")
orthos <-read_tsv("/projects/b1059/projects/Nicolas/shiny_gmvt/data/masterOrthoDB_wAlias.tsv")
spp <- c("N2","QX1410","NIC58") 
geneid <- "ben-1"

test <- orthos %>% 
  dplyr::select(c(spp,c("seqname", "WB_id", "WB_alias"))) %>%
  dplyr::filter(if_any(all_of(c("seqname", "WB_id", "WB_alias")), ~grepl(geneid, .)))

filter <- c(unlist(strsplit(test$WB_id,", ")),unlist(strsplit(test$QX1410,", ")),unlist(strsplit(test$NIC58,", ")))

plot_df <- tair %>% dplyr::filter(L1 %in% filter)
sub_gff <- tair %>% dplyr::filter(L1 %in% filter)

plotNormGeneFeatures(plot_df)
