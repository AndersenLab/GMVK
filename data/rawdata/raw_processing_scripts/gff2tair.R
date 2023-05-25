library(dplyr)
library(tidyr)
library(readr)
library(ape)


#showtext::showtext_auto()


processGFF <- function(path,intron) {
  
  gff <- ape::read.gff(file = path)
  
  geneF <- gff %>% 
    dplyr::filter(type =="gene") %>%
    tidyr::separate(attributes, into = c("ID", "post"),sep = ";Name=") %>%
    dplyr::mutate(ID=gsub("ID=","",ID)) %>%
    dplyr::mutate(ID=gsub("gene:","",ID)) %>%
    dplyr::mutate(ID=gsub("Gene:","",ID)) %>%
    tidyr::separate(post,into=c("pre","aliases"), sep=";Alias=") %>%
    tidyr::separate(aliases,into=c("mAlias","rest"), sep=",",extra = "merge") %>%
    dplyr::mutate(L1=ID)
  
  L1 <- geneF %>% 
    dplyr::select(seqid,start,end,type,strand,ID,L1)
  
  geneAlias <- geneF %>% 
    dplyr::select(ID,mAlias) 
  
  L2 <- gff %>% 
    dplyr::filter(type=="mRNA") %>%
    tidyr::separate(attributes,into = c("ID","post"),sep=";Parent=") %>%
    dplyr::mutate(ID=gsub("ID=","",ID)) %>%
    dplyr::mutate(ID=gsub("transcript:","",ID)) %>%
    dplyr::mutate(ID=gsub("Transcript:","",ID)) %>%
    tidyr::separate(post,into = c("L1","rest"),sep=";Name=") %>%
    dplyr::mutate(L1=gsub("gene:","",L1)) %>%
    dplyr::mutate(L1=gsub("Gene:","",L1)) %>%
    dplyr::select(seqid,start,end,type,strand,ID,L1)
  
  if(any(grepl("biotype",L2$L1))) {
    temp <- L2 %>% 
      tidyr::separate(L1,into=c("L1_fix","rest"),sep=";biotype", remove = T) %>%
      dplyr::rename(L1=L1_fix) %>%
      dplyr::select(-rest)
    L2 <- temp
  }
  
  parentL2 <- L2 %>% 
    dplyr::select(ID,L1)
  
  children <- gff %>% 
    dplyr::filter(!(type=="mRNA") & !(type=="gene")) %>%
    dplyr::filter(!(type=="three_prime_UTR") & !(type=="five_prime_UTR")) %>%
    tidyr::separate(attributes,into=c("L3","L2"),sep=";") %>%
    dplyr::mutate(L2=ifelse(is.na(L2) | grepl("Note=",L2),L3,L2)) %>%
    dplyr::mutate(L3=gsub("ID=","",L3)) %>%
    dplyr::mutate(L3=gsub("CDS:","",L3)) %>%
    dplyr::mutate(L2=gsub("Parent=","",L2)) %>%
    dplyr::mutate(L2=gsub("transcript:","",L2)) %>%
    dplyr::mutate(L2=gsub("Transcript:","",L2)) %>%
    dplyr::mutate(L3=gsub("Parent=","",L3)) %>%
    dplyr::mutate(L3=gsub("transcript:","",L3)) %>%
    dplyr::mutate(L3=gsub("Transcript:","",L3)) 
  
  if(any(grepl("biotype",children$L1))) {
    temp <- children %>% 
      tidyr::separate(L1,into=c("L1_fix","rest"),sep=";biotype", remove = T) %>%
      dplyr::rename(L1=L1_fix) %>%
      dplyr::select(-rest)
    children <- temp
  }
  
  if(any(grepl(",",children$L2))) {
    children <- children %>%
      tidyr::separate_rows(L2,sep = ",") %>%
      dplyr::mutate(L3=L2) %>%
      dplyr::left_join(parentL2,by = c("L2"="ID")) %>%
      dplyr::mutate(L3=paste0(L3,"_",type)) %>%
      dplyr::group_by(L3) %>%
      dplyr::mutate(L3=paste0(L3,row_number()))%>%
      dplyr::ungroup()
  } else {
    children <- children %>%
      dplyr::left_join(parentL2,by = c("L2"="ID")) 
  }
  
  parentL3 <- children %>% 
    dplyr::select(L2,L3)
  
  L3 <- children %>% 
    dplyr::rename(ID=L3) %>%
    dplyr::select(seqid,start,end,type,strand,ID,L1)
  
  tabGFF <- BiocGenerics::rbind(L1,L2,L3) %>%
    dplyr::arrange(seqid,start) %>%
    dplyr::left_join(geneAlias,by = c("L1"="ID")) #%>%
    #dplyr::group_by(L1) %>%
    #dplyr::mutate(end=end-min(start)) %>%
    #dplyr::mutate(start=start-min(start))
  
  geneMain <- tabGFF %>% 
    dplyr::filter(type=="gene")
  
  subHead <- tabGFF %>% 
    dplyr::filter(type=="mRNA")
  
  geneHead <- subHead %>%
    dplyr::select(ID,L1) %>%
    dplyr::left_join(geneMain,by ="L1") %>%
    dplyr::select(-ID.y) %>%
    dplyr::rename(L2="ID.x") %>%
    dplyr::mutate(ID=L2) %>%
    dplyr::select(seqid,start,end,type,strand,ID,L1,mAlias,L2)
  
  remFeatures <- tabGFF %>% 
    dplyr::filter(!(type=="gene")) %>%
    dplyr::left_join(parentL3,by=c("ID"="L3")) %>%
    dplyr::mutate(L2=ifelse(type=="mRNA",ID,L2)) 
  
  if (intron) {
    finalFeatures <- rbind(geneHead,remFeatures) %>%
      dplyr::group_by(L2) %>%
      dplyr::arrange(start) %>%
      #dplyr::left_join(geneAlias,by = c("L1"="ID")) %>%
      dplyr::mutate(startPos=min(start)) %>%
      dplyr::mutate(end=end-min(start)) %>%
      dplyr::mutate(start=start-min(start))
    
  } else {
    introns <- remFeatures %>% 
      dplyr::filter(type=="exon") %>% 
      dplyr::group_by(L2) %>%
      dplyr::mutate(istart=end+1,iend=lead(start)-1) %>%
      dplyr::select(-start,-end) %>%
      dplyr::rename(start=istart,end=iend) %>%
      dplyr::filter(!is.na(end)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(type="intron") %>%
      dplyr::select(seqid,start,end,type,strand,ID,L1,mAlias,L2) %>%
      dplyr::mutate(ID=gsub("exon","intron",ID))
    
    finalFeatures <- rbind(geneHead,remFeatures,introns) %>%
      dplyr::group_by(L2) %>%
      dplyr::arrange(start) %>%
      #dplyr::left_join(geneAlias,by = c("L1"="ID"))%>%
      dplyr::mutate(startPos=min(start)) %>%
      dplyr::mutate(end=end-min(start)) %>%
      dplyr::mutate(start=start-min(start)) 
    
  }

  #features <- dplyr::group_split(finalFeatures)
  return(finalFeatures)
}

path <- "data/rawdata/gff3/Curation-VF-230214.PC.clean.renamed.csq.gff3"
#gff_cbrig <- ape::read.gff(file = path)
intron <- FALSE
tair_cbrig <- processGFF(path,intron)
tair_cbrig$source <- "QX1410"

path <- "data/rawdata/gff3/NIC58.final_annotation.fixed.CSQ.gff3"
#gff_ctrop <- ape::read.gff(file = path)
intron <- FALSE
tair_ctrop <- processGFF(path,intron)
tair_ctrop$source <- "NIC58"

path <- "data/rawdata/gff3/c_elegans.PRJNA13758.WS279.protein_coding.gff3"
#gff_cele <- ape::read.gff(file = path)
intron <- TRUE
tair_celeg <- processGFF(path,intron)
tair_celeg$source <- "N2"

masterDB <- write.table(rbind(tair_cbrig,tair_celeg,tair_ctrop),file = "data/master_DB.tsv",quote = F,row.names = F,sep = '\t')

