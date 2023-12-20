library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(showtext)
library(ape)
library(ggrepel)

plotNormGeneFeatures <- function(sub_gff,mode,sep,labeller,cexpt,lgy,lgx) {
  
  showtext::showtext_auto()
  
  GFFbyL2 <- sub_gff %>% 
    dplyr::arrange(source,dplyr::desc(mAlias)) %>%
    dplyr::group_by(dplyr::desc(source),dplyr::desc(mAlias),dplyr::desc(L2)) %>%
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
  
  neg_L1L2 <- neg_str %>% 
    dplyr::filter(!(type=="exon") & !(type =="intron") & !(type =="CDS"))
  
  neg_corr <- rbind(neg_L1L2,neg_L3) %>% 
    dplyr::arrange(start)
  
  pos_str <- GFFbyL2 %>% dplyr::filter(strand=="+")
  
  grouped_gff <- rbind(neg_corr,pos_str)
  
  if (mode == "lon") {
    longestCDS <- grouped_gff %>%
      dplyr::filter(type=="CDS") %>%
      dplyr::mutate(len=end-start) %>%
      dplyr::group_by(L2) %>%
      dplyr::mutate(maxLen=sum(len)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(L2,.keep_all = T) %>%
      dplyr::group_by(L1) %>%
      dplyr::mutate(longest=ifelse(maxLen==max(maxLen),"Y","N")) %>%
      dplyr::ungroup() %>%
      dplyr::filter(longest=="Y") %>%
      dplyr::select(L2,maxLen) %>%
      dplyr::rename(CDSlen=maxLen)
    
    longestTran <- grouped_gff %>%
      dplyr::filter(type=="mRNA") %>%
      dplyr::left_join(longestCDS,by="L2") %>%
      dplyr::filter(!is.na(CDSlen)) %>%
      dplyr::mutate(len=end-start) %>%
      dplyr::group_by(L1) %>%
      dplyr::mutate(longest=ifelse(len==max(len),"Y","N")) %>%
      dplyr::ungroup() %>%
      dplyr::filter(longest=="Y") #%>%
    # dplyr::group_by(L1) %>%
    # dplyr::mutate(size=n())
    
    temp <- grouped_gff %>% 
      dplyr::filter(L2 %in% longestTran$L2) %>%
      dplyr::arrange(source,dplyr::desc(mAlias)) %>%
      dplyr::group_by(dplyr::desc(source),dplyr::desc(mAlias),dplyr::desc(L2)) %>%
      dplyr::mutate(ngroup=cur_group_id()) %>%
      dplyr::ungroup()
    
    grouped_gff <- temp
  }
  
  if(sep=="gsep" | sep == "gspp"){
    temp <- grouped_gff %>%
      dplyr::arrange(Orthogroup,source,dplyr::desc(mAlias)) %>%
      dplyr::group_by(Orthogroup,dplyr::desc(source),dplyr::desc(mAlias),dplyr::desc(L2)) %>%
      dplyr::mutate(ngroup=cur_group_id()) %>%
      dplyr::ungroup()
    
    grouped_gff <- temp
  }
  
  CDS_start <- grouped_gff %>% 
    dplyr::group_by(L2) %>%
    dplyr::filter(type=="CDS") %>% 
    dplyr::filter(start==(min(start))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(start==max(start))
  
  maxCDS <- CDS_start %>% dplyr::select(start,ngroup)
  
  reference_tran <- grouped_gff %>% dplyr::filter(ngroup==maxCDS[1,]$ngroup)
  sub_tran <- grouped_gff %>% 
    dplyr::filter(!(ngroup==maxCDS[1,]$ngroup)) %>%
    dplyr::group_by(L2) %>%
    dplyr::mutate(CDSstr=ifelse(type=="CDS",start,1e9)) %>%
    dplyr::mutate(hjust_cds=min(CDSstr)) %>%
    dplyr::mutate(shift=maxCDS[1,]$start-hjust_cds) %>%
    dplyr::mutate(start=start+(maxCDS[1,]$start-hjust_cds)) %>%
    dplyr::mutate(end=end+(maxCDS[1,]$start-hjust_cds)) %>%
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
    dplyr::select(Orthogroup,L2,mAlias,ngroup,strand,source) %>%
    dplyr::group_by(mAlias) %>%
    dplyr::mutate(iso=1:n()) %>%
    #dplyr::rowwise() %>%
    dplyr::mutate(verif=ifelse(any(iso > 1),"Y","N")) %>% #mAlias," (",as.character(iso),")"), mAlias)) 
    dplyr::ungroup() %>%
    dplyr::mutate(label=paste0("italic(",mAlias,")","~(","italic(",L2,"))")) %>%
    dplyr::left_join(lastE, by="L2") 
  
  cdsE <- plottable_df %>% 
    dplyr::filter(type=="CDS")
  
  
  plot <- ggplot() + geom_rect(data = restE, aes(xmin = start/1e3 ,xmax = end/1e3 ,ymin=ngroup-0.25 , ymax=ngroup+0.25),fill="grey") +
    geom_rect(data = cdsE, aes(xmin = start/1e3 ,xmax = end/1e3 ,ymin=ngroup-0.25 , ymax=ngroup+0.25,fill=source)) +
    geom_segment(data = plottable_df %>% dplyr::filter(type=="intron"), aes(x=start/1e3,xend=(start+((end-start)/2))/1e3,y=ngroup,yend=ngroup+0.25),color="darkgrey") +
    geom_segment(data = plottable_df %>% dplyr::filter(type=="intron"), aes(x=(start+((end-start)/2))/1e3,xend=end/1e3,y=ngroup+0.25,yend=ngroup),color="darkgrey") +
    geom_text(data = labels, aes(label=label,x=0,y=ngroup+0.45,hjust=0),parse = T,size=cexpt)+
    geom_text(data = labels, aes(label=paste0("(",strand,")"),x=end/1e3+0.015*(span/1e3),y=ngroup,hjust=0),size=cexpt)+
    #annotate('text',x=0,y=labels$ngroup+0.45,label=paste0((labels$label)),parse=F,fontface='italic',hjust=0) +
    #annotate('text',x=labels$end/1e3+0.015*(span/1e3),y=labels$ngroup,label=paste0("(",labels$strand,")"),parse=F,hjust=0) +
    theme(axis.title.y  = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_line(),
          panel.background = element_blank(),
          legend.position=c(lgx,lgy),
          legend.title = element_blank(),
          legend.key.size = unit(cexpt/10,'cm'),
          legend.text = element_text(size=(cexpt*2.5)),
          axis.text.x = element_text(size=cexpt*2.5),
          axis.title.x = element_text(size=cexpt*3),
          strip.text = element_text(size=cexpt*4),
          plot.title = element_text(hjust = 0,size = 12)) + xlab("Transcript length (kb)") + scale_fill_manual(values=c("#e69f00","#cc799d","#009e73")) 
  
  if (sep=="spp"){
    finalPlot <- plot + facet_wrap(~source,ncol = 1,scales=labeller)
  } else if(sep=="gsep"){
    finalPlot <- plot + facet_wrap(~Orthogroup,ncol = 1,scales=labeller)
  } else if(sep=="gspp"){
    finalPlot <- plot + facet_wrap(~Orthogroup+source,ncol = 1,scales=labeller)
  } else {
    finalPlot <- plot 
  }
  
  
  return(finalPlot)
}