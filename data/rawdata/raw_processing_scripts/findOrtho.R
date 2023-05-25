library(dplyr)
library(tidyr)
library(readr)
library(ape)
library(stringr)

path <- "data/rawdata/c_elegans.PRJNA13758.WS279.protein_coding.gff"
gff <- ape::read.gff(file = path)

geneF <- gff %>% 
  dplyr::filter(type =="gene") %>%
  tidyr::separate(attributes, into = c("ID", "post"),sep = ";Name=") %>%
  dplyr::mutate(ID=gsub("ID=","",ID)) %>%
  dplyr::mutate(ID=gsub("gene:","",ID)) %>%
  dplyr::mutate(ID=gsub("Gene:","",ID)) %>%
  tidyr::separate(post,into=c("pre","aliases"), sep=";Alias=") %>%
  tidyr::separate(aliases,into=c("mAlias","rest"), sep=",",extra = "merge") %>%
  dplyr::mutate(L1=ID) %>%
  tidyr::separate(pre,into=c("prep","seqname"),sep = ";sequence_name=") %>%
  tidyr::separate(seqname,into=c("seqname_clean","rest2"), sep= ";biotype=")

geneAlias <- geneF %>% 
  dplyr::select(ID,mAlias,seqname_clean) 

tranF <- gff %>% 
  dplyr::filter(type =="mRNA") %>%
  tidyr::separate(attributes, into = c("L2", "post"),sep = ";Parent=") %>%
  tidyr::separate(post, into=c("L1","post2"), sep= ";Name") %>%
  dplyr::select(L1,L2) %>%
  dplyr::mutate(L1=gsub("Gene:","",L1)) %>%
  dplyr::mutate(L2=gsub("ID=Transcript:","",L2)) %>%
  dplyr::left_join(geneAlias,by=c("L1"="ID"))

orthogroups <- readr::read_tsv("data/masterOrthoDB.tsv")
colnames(orthogroups) <- c("Orthogroup","QX1410","NIC58","N2")

orthos_wAlias <- orthogroups %>% 
  dplyr::select(Orthogroup,N2,QX1410,NIC58) %>%
  dplyr::mutate(N2=gsub("Transcript_","",N2)) %>%
  dplyr::mutate(QX1410=gsub("transcript_","",QX1410)) %>%
  tidyr::separate_rows(N2,sep=", ") %>%
  dplyr::left_join(tranF, by=c("N2"="L2")) %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::mutate(L1_group=paste(L1,collapse = ", ")) %>%
  dplyr::mutate(mAlias_group=paste(mAlias,collapse = ", ")) %>%
  dplyr::mutate(seqname_group=paste(seqname_clean,collapse = ", ")) %>%
  dplyr::mutate(N2_group=paste(N2,collapse=", ")) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(Orthogroup, .keep_all = T) %>%
  dplyr::select(-seqname_clean,-L1,-mAlias,-N2) %>%
  tidyr::separate_rows(QX1410,sep=", ") %>%
  tidyr::separate(QX1410,into=c("prefix","Ngene","Ntran"),sep="\\.") %>%
  dplyr::select(-Ntran) %>%
  tidyr::unite("QX1410",prefix,Ngene,sep = ".",remove = T) %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::mutate(QX1410_group=paste(QX1410,collapse = ", ")) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(Orthogroup, .keep_all = T) %>%
  dplyr::select(-QX1410) %>%
  tidyr::separate_rows(NIC58,sep=", ") %>%
  tidyr::separate(NIC58,into=c("prefix","Ngene","Ntran"),sep="\\.") %>%
  dplyr::select(-Ntran) %>%
  tidyr::unite("NIC58",prefix,Ngene,sep = ".",remove = T) %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::mutate(NIC58_group=paste(NIC58,collapse = ", ")) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(Orthogroup, .keep_all = T) %>%
  dplyr::select(-NIC58) %>%
  dplyr::mutate(QX1410_group=ifelse(grepl("NA.NA",QX1410_group),NA,QX1410_group)) %>%
  dplyr::mutate(NIC58_group=ifelse(grepl("NA.NA",NIC58_group),NA,NIC58_group))
  
colnames(orthos_wAlias) <- c("Orthogroup","WB_id","WB_alias","seqname","N2","QX1410","NIC58")

write.table(orthos_wAlias,file="data/masterOrthoDB_wAlias.tsv",quote = F,row.names = F,sep = "\t")  
