################### Metascape Analysis to extract biological proccesses
################### Packages ################### 
install.packages("readxl")
install.packages('DescTools')
install.packages('tidyverse')
install.packages("dplyr")
library(readxl)
library(data.table)
library(tidyverse)
library(DescTools)
library(matrixStats)

#-----------------------------------
################### Pathways ################### 
#These are the most biologically significant pathways chosen by looking at Metascape Go heatmaps of expression 
################### Actin ################### 

cluster_group_dataframe <- read_xlsx(file.choose(),sheet = 2) #import excel file output from metascape containing all pulled cluster pathways manually
dirname(rstudioapi::getActiveDocumentContext()$path) #setting working directory to source file location for data import
dir.create("GenePathwayLists") #create an output folder to export merged datasets to

#indexing search for pathways containing actin

cluster_group_dataframe_actin <- cluster_group_dataframe[grep("\\bactin\\b", cluster_group_dataframe$Description), ]

#calling for column containing gene symbol entries and splitting them into a list 
actin_list_symbols <- strsplit(as.character(cluster_group_dataframe_actin$Symbols[1]),",",fixed = T) 
actin_list_symbols <- do.call(rbind.data.frame, actin_list_symbols)
actin_list_symbols <- t(actin_list_symbols)
rownames(actin_list_symbols) <- c(1:length(row_number(actin_list_symbols)))
colnames(actin_list_symbols) = 'Gene.Symbols'

IKDNprBINDALL <- read.table(file.choose(), header=TRUE, sep="\t",fill = TRUE) #import master sheet of expression 

colnames(IKDNprBINDALL)[8] <- 'Gene.Symbols' #rename column because merge() function can only merge by a shared column

actin_merged_df <- merge(IKDNprBINDALL,actin_list_symbols,by='Gene.Symbols') #merged dataframe created 

actin_merged_df_cluster <- actin_merged_df[c(1,9:24)]  #subsetting columns necessary for clustering 
actin_merged_df_fc <- actin_merged_df[c(1,9:24,25,28,31,34,37,40,43)]  #subsetting columns necessary 

colnames(actin_merged_df_cluster) <- c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2")
colnames(actin_merged_df_fc) <- c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2","D1.FC","D2.FC","D3.FC","D4.FC","D6.FC","D12.FC","D18.FC")

write.table(actin_merged_df_cluster, "GenePathwayLists/actin_merged_df.txt",sep="\t",row.names = F)


################### Vasculature ################### 

#indexing search for pathways containing vasculature
cluster_group_dataframe_vasc <- cluster_group_dataframe[grep("\\vasculature\\b", cluster_group_dataframe$Description), ]

#calling for column containing gene symbol entries and splitting them into a list 
vasc_list_symbols <- strsplit(as.character(cluster_group_dataframe_vasc$Symbols[1]),",",fixed = T)
vasc_list_symbols <- do.call(rbind.data.frame, vasc_list_symbols)
vasc_list_symbols <- t(vasc_list_symbols)
rownames(vasc_list_symbols) <- c(1:length(row_number(vasc_list_symbols)))
colnames(vasc_list_symbols) = 'Gene.Symbols'
colnames(IKDNprBINDALL)[8] <- 'Gene.Symbols'
vasc_merged_df <- merge(IKDNprBINDALL,vasc_list_symbols,by='Gene.Symbols')
colnames(vasc_merged_df) 
vasc_merged_cluster_df <- vasc_merged_df[c(1,9:24)] 
vasc_merged_df_fc <- actin_merged_df[c(1,9:24,25,28,31,34,37,40,43)] 
colnames(vasc_merged_cluster_df) <- c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2")
colnames(vasc_merged_df_fc) <-  c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2","D1.FC","D2.FC","D3.FC","D4.FC","D6.FC","D12.FC","D18.FC")

#exporting 
write.table(vasc_merged_cluster_df, "GenePathwayLists/vasc_merged_cluster_df.txt",sep="\t",row.names = F)

################### Response to Wounding ################### 

#indexing search for pathways containing wounding
cluster_group_dataframe_wound <- cluster_group_dataframe[grep("\\wounding\\b", cluster_group_dataframe$Description), ]
summary(cluster_group_dataframe_wound)

#calling for column containing gene symbol entries and splitting them into a list 
wound_list_symbols <- strsplit(as.character(cluster_group_dataframe_wound$Symbols[1]),",",fixed = T)
wound_list_symbols <- do.call(rbind.data.frame, wound_list_symbols)
wound_list_symbols <- t(wound_list_symbols)
rownames(wound_list_symbols) <- c(1:length(row_number(wound_list_symbols)))
colnames(wound_list_symbols) = 'Gene.Symbols'
colnames(IKDNprBINDALL)[8] <- 'Gene.Symbols'
wound_merged_df <- merge(IKDNprBINDALL,wound_list_symbols,by='Gene.Symbols')
wound_merged_cluster_df <- wound_merged_df[c(1,9:24)] 
wound_merged_df_fc <- wound_merged_df[c(1,9:24,25,28,31,34,37,40,43)] 
colnames(wound_merged_cluster_df) <- c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2")
colnames(wound_merged_df_fc) <-  c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2","D1.FC","D2.FC","D3.FC","D4.FC","D6.FC","D12.FC","D18.FC")
write.table(wound_merged_cluster_df, "GenePathwayLists/wound_merged_df.txt",sep="\t",row.names = F)

################### Regulation of GTPase activity ################### 

cluster_group_dataframe_gtpase <- cluster_group_dataframe[grep("\\GTPase\\b", cluster_group_dataframe$Description), ]

cluster_group_dataframe_gtpase$Description
gtpase_list_symbols <- strsplit(as.character(cluster_group_dataframe_gtpase$Symbols[1]),",",fixed = T)
gtpase_list_symbols <- do.call(rbind.data.frame, gtpase_list_symbols)
gtpase_list_symbols <- t(gtpase_list_symbols)
rownames(gtpase_list_symbols) <- c(1:length(row_number(gtpase_list_symbols)))
colnames(gtpase_list_symbols) = 'Gene.Symbols'
colnames(IKDNprBINDALL)[8] <- 'Gene.Symbols'
gtpase_merged_df <- merge(IKDNprBINDALL,gtpase_list_symbols,by='Gene.Symbols')
gtpase_merged_cluster_df <- gtpase_merged_df[c(1,9:24)] 
gtpase_merged_df_fc <- gtpase_merged_df[c(1,9:24,25,28,31,34,37,40,43)] 
colnames(gtpase_merged_cluster_df) <- c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2")
colnames(gtpase_merged_df_fc) <-  c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2","D1.FC","D2.FC","D3.FC","D4.FC","D6.FC","D12.FC","D18.FC")
write.table(gtpase_merged_cluster_df, "GenePathwayLists/gtpase_merged_df.txt",sep="\t",row.names = F)

################### Axon development ################### 

cluster_group_dataframe_axon <- cluster_group_dataframe[grep("axon", cluster_group_dataframe$Description), ] #find when contained within strings too 

(cluster_group_dataframe_axon$Description)
axon_list_symbols <- strsplit(as.character(cluster_group_dataframe_axon$Symbols[1]),",",fixed = T)
axon_list_symbols <- do.call(rbind.data.frame, axon_list_symbols)
axon_list_symbols <- t(axon_list_symbols)
rownames(axon_list_symbols) <- c(1:length(row_number(axon_list_symbols)))
colnames(axon_list_symbols) = 'Gene.Symbols'
colnames(IKDNprBINDALL)[8] <- 'Gene.Symbols'
axon_merged_df <- merge(IKDNprBINDALL,axon_list_symbols,by='Gene.Symbols')
axon_merged_cluster_df <- axon_merged_df[c(1,9:24)] 
axon_merged_df_fc <- axon_merged_df[c(1,9:24,25,28,31,34,37,40,43)] 
colnames(axon_merged_cluster_df) <- c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2")
colnames(axon_merged_df_fc) <-  c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2","D1.FC","D2.FC","D3.FC","D4.FC","D6.FC","D12.FC","D18.FC")
write.table(axon_merged_cluster_df, "GenePathwayLists/axon_merged_df.txt",sep="\t",row.names = F)

################### Cell junction organisation ################### 

cluster_group_dataframe_cj <- cluster_group_dataframe[grep("cell junction", cluster_group_dataframe$Description), ]

(cluster_group_dataframe_cj$Description)
cj_list_symbols <- strsplit(as.character(cluster_group_dataframe_cj$Symbols[1]),",",fixed = T)
cj_list_symbols <- do.call(rbind.data.frame, cj_list_symbols)
cj_list_symbols <- t(cj_list_symbols)
rownames(cj_list_symbols) <- c(1:length(row_number(cj_list_symbols)))
colnames(cj_list_symbols) = 'Gene.Symbols'
colnames(IKDNprBINDALL)[8] <- 'Gene.Symbols'
cj_merged_df <- merge(IKDNprBINDALL,cj_list_symbols,by='Gene.Symbols')
cj_merged_cluster_df <- cj_merged_df[c(1,9:24)] 
cj_merged_df_fc <- cj_merged_df[c(1,9:24,25,28,31,34,37,40,43)] 
colnames(cj_merged_cluster_df) <- c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2")
colnames(cj_merged_df_fc) <-  c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2","D1.FC","D2.FC","D3.FC","D4.FC","D6.FC","D12.FC","D18.FC")
write.table(cj_merged_cluster_df, "GenePathwayLists/cj_merged_df.txt",sep="\t",row.names = F)

################### Cell-cell adhesion ################### 

cluster_group_dataframe_cell_ad <- cluster_group_dataframe[grep("cell-cell adhesion", cluster_group_dataframe$Description), ]

(cluster_group_dataframe_cell_ad$Description)
cell_ad_list_symbols <- strsplit(as.character(cluster_group_dataframe_cell_ad$Symbols[1]),",",fixed = T)
cell_ad_list_symbols <- do.call(rbind.data.frame, cell_ad_list_symbols)
cell_ad_list_symbols <- t(cell_ad_list_symbols)
rownames(cell_ad_list_symbols) <- c(1:length(row_number(cell_ad_list_symbols)))
colnames(cell_ad_list_symbols) = 'Gene.Symbols'
colnames(IKDNprBINDALL)[8] <- 'Gene.Symbols'
cell_ad_merged_df <- merge(IKDNprBINDALL,cell_ad_list_symbols,by='Gene.Symbols')
cell_ad_merged_cluster_df <- cell_ad_merged_df[c(1,9:24)] 
cell_ad_merged_df_fc <- cell_ad_merged_df[c(1,9:24,25,28,31,34,37,40,43)] 
colnames(cell_ad_merged_cluster_df) <- c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2")
colnames(cell_ad_merged_df_fc) <-  c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2","D1.FC","D2.FC","D3.FC","D4.FC","D6.FC","D12.FC","D18.FC")

write.table(cell_ad_merged_cluster_df, "GenePathwayLists/cell_ad_merged_df.txt",sep="\t",row.names = F)


################### Transemmbrane receptor protein tyrosine kinase ################### 

cluster_group_dataframe_tyrk <- cluster_group_dataframe[grep("tyrosine kinase", cluster_group_dataframe$Description), ]

(cluster_group_dataframe_tyrk$Description)
tyrk_list_symbols <- strsplit(as.character(cluster_group_dataframe_tyrk$Symbols[1]),",",fixed = T)
tyrk_list_symbols <- do.call(rbind.data.frame, tyrk_list_symbols)
tyrk_list_symbols <- t(tyrk_list_symbols)
rownames(tyrk_list_symbols) <- c(1:length(row_number(tyrk_list_symbols)))
colnames(tyrk_list_symbols) = 'Gene.Symbols'
colnames(IKDNprBINDALL)[8] <- 'Gene.Symbols'
tyrk_merged_df <- merge(IKDNprBINDALL,tyrk_list_symbols,by='Gene.Symbols')
tyrk_merged_cluster_df <- tyrk_merged_df[c(1,9:24)] 
tyrk_merged_df_fc <- tyrk_merged_df[c(1,9:24,25,28,31,34,37,40,43)] 
colnames(tyrk_merged_cluster_df) <- c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2")
colnames(tyrk_merged_df_fc) <-  c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2","D1.FC","D2.FC","D3.FC","D4.FC","D6.FC","D12.FC","D18.FC")
write.table(tyrk_merged_cluster_df, "GenePathwayLists/tyrk_merged_df.txt",sep="\t",row.names = F)

################### Head development ################### 

cluster_group_dataframe_head <- cluster_group_dataframe[grep("head development", cluster_group_dataframe$Description), ]
(cluster_group_dataframe_head$Description)
head_list_symbols <- strsplit(as.character(cluster_group_dataframe_head$Symbols[1]),",",fixed = T)
head_list_symbols <- do.call(rbind.data.frame, head_list_symbols)
head_list_symbols <- t(head_list_symbols)
rownames(head_list_symbols) <- c(1:length(row_number(head_list_symbols)))
colnames(head_list_symbols) = 'Gene.Symbols'
colnames(IKDNprBINDALL)[8] <- 'Gene.Symbols'
head_merged_df <- merge(IKDNprBINDALL,head_list_symbols,by='Gene.Symbols')
head_merged_cluster_df <- head_merged_df[c(1,9:24)] 
head_merged_df_fc <- head_merged_df[c(1,9:24,25,28,31,34,37,40,43)] 
colnames(head_merged_cluster_df) <- c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2")
colnames(head_merged_df_fc) <-  c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2","D1.FC","D2.FC","D3.FC","D4.FC","D6.FC","D12.FC","D18.FC")
write.table(head_merged_cluster_df, "GenePathwayLists/head_merged_df.txt",sep="\t",row.names = F)

################### MAP cascade ################### 

cluster_group_dataframe_mapk <- cluster_group_dataframe[grep("MAPK cascade", cluster_group_dataframe$Description), ]


(cluster_group_dataframe_mapk$Description)
mapk_list_symbols <- strsplit(as.character(cluster_group_dataframe_mapk$Symbols[1]),",",fixed = T)
mapk_list_symbols <- do.call(rbind.data.frame, mapk_list_symbols)
mapk_list_symbols <- t(mapk_list_symbols)
rownames(mapk_list_symbols) <- c(1:length(row_number(mapk_list_symbols)))
colnames(mapk_list_symbols) = 'Gene.Symbols'
colnames(IKDNprBINDALL)[8] <- 'Gene.Symbols'
mapk_merged_df <- merge(IKDNprBINDALL,mapk_list_symbols,by='Gene.Symbols')
mapk_merged_cluster_df <- mapk_merged_df[c(1,9:24)] 
mapk_merged_df_fc <- mapk_merged_df[c(1,9:24,25,28,31,34,37,40,43)] 
colnames(mapk_merged_cluster_df) <- c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2")
colnames(mapk_merged_df_fc) <-  c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2","D1.FC","D2.FC","D3.FC","D4.FC","D6.FC","D12.FC","D18.FC")
write.table(mapk_merged_cluster_df, "GenePathwayLists/mapk_merged_df.txt",sep="\t",row.names = F)


################### Transcription Regulators David List ################### 

transcription_regulators_list <- read.table(file.choose(),header=TRUE, sep="\t",fill = TRUE)#import file downloaded from david containing gene lists associated with transcription regulators
colnames(transcription_regulators_list)[1] <- 'Gene.Symbols'
transcription_reg_merged <- merge(transcription_regulators_list,IKDNprBINDALL,by='Gene.Symbols') #merged dataframe created 
transcription_reg_merged <- transcription_reg_merged[c(1,11:26)] 
colnames(transcription_reg_merged) <- c("Gene.Symbols","DMSO1","DMSO2","4OHT.D1.1","4OHT.D1.2","4OHT.D2.1","4OHT.D2.2","4OHT.D3.1","4OHT.D3.2","4OHT.D4.1","4OHT.D4.2","4OHT.D6.1","4OHT.D6.2","4OHT.D12.1","4OHT.D12.2","4OHT.D18.1","4OHT.D18.2")


sort(IKDNprBINDALL$Gene.Symbols)

dir.create("DavidLists") #create an output folder to export merged datasets to

dim(transcription_reg_merged)

write.table(transcription_reg_merged, "DavidLists/transcription_reg_merged.txt",sep="\t",row.names = F) #file will be used to cluster and produce heatmaps 



#-----------------------------------
################### Overlap ################### 

# For Unique Values
symbol_lst <- list(
  mapk_list_symbols,
  
  head_list_symbols,
  
  tyrk_list_symbols,
  
  cell_ad_list_symbols,
  
  cj_list_symbols,
  
  axon_list_symbols,
  
  gtpase_list_symbols,
  
  wound_list_symbols,
  
  vasc_list_symbols)

unique(unlist(lapply(symbol_lst, function(x) unique(x[,1]))))


################### Range Bar Plots ################### 

# Actin

actin_matrix <- as.matrix(actin_merged_df_fc[,2:17]) # create a matrix to calculate range 
actin_range <- data.frame(rowRanges(actin_matrix)) # row ranges function calculates the minimum and maximum values of each row, in this case the min and max expression level per gene 
colnames(actin_range) <- c('Min','Max') #label for clarity 
actin_range_df <- cbind(actin_merged_df_fc,actin_range,actin_range[2]-actin_range[1]) #find the range by subtracting min from max 
colnames(actin_range_df)[ncol(actin_range_df)] <- "Range" #label for clarity 

actin_range

#Not average but differential 

# Plotingt high time points 

#log 2 for expression and fold change 

actin.d12 <- data.frame(cbind(log2(actin_range_df$Range),actin_range_df$D12.FC)) #ploting range against fc 


colnames(actin.d12) <- c("Range","D12.FC")

actin.d12 <- mutate(actin.d12,Colour = ifelse(D12.FC > 1 & Range < 0.05,"red", 
                                              ifelse(D12.FC < -1 & Range < 0.05, "blue","grey"))) #color coding values by threshold for plotting

plot(actin.d12$D12.FC,actin.d12$Range)

ggplot(actin.d12, aes(x=D12.FC, y=Range)) +
  geom_point(aes(color=color))+
  scale_color_manual(values=c("blue", "light grey", "red"))+
  theme(panel.background = element_rect(fill="white"),
        axis.line=element_line(colour="black"),
        legend.position="none",
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 

actin.d18 <-data.frame(cbind(log2(actin_range_df$Range),actin_range_df$D18.FC))
colnames(actin.d18) <- c("Range","D18.FC")

actin.d18 <- mutate(actin.d18,Colour = ifelse(D18.FC > 1 & Range < 0.05,"red", 
                                              ifelse(D18.FC < -1 & Range < 0.05, "blue","grey")))

plot(actin.d18$D18.FC,actin.d18$Range)


ggplot(actin.d18, aes(x=D18.FC, y=Range)) +
  geom_point(aes(color=color))+
  scale_color_manual(values=c("blue", "light grey", "red"))+
  theme(panel.background = element_rect(fill="white"),
        axis.line=element_line(colour="black"),
        legend.position="none",
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 


# Vascular

vasc_matrix <- as.matrix(vasc_merged_df_fc[,2:17])
vasc_range <- data.frame(rowRanges(vasc_matrix))
colnames(vasc_range) <- c('Min','Max')
vasc_df_range <- cbind(vasc_merged_df_fc,vasc_range,vasc_range[2]-vasc_range[1])
colnames(vasc_df_range)[ncol(vasc_df_range)] <- "Range"

vasc.d12 <-data.frame(cbind(log10(vasc_df_range$Range),vasc_df_range$D12.FC))
colnames(vasc_df_range) <- c("Range","D12.FC")

# Wound 


wound_matrix <- as.matrix(wound_merged_df_fc[,2:17])
wound_range <- data.frame(rowRanges(wound_matrix))
colnames(wound_range) <- c('Min','Max')
wound_df_range <- cbind(wound_merged_df_fc,wound_range,wound_range[2]-wound_range[1])
colnames(wound_df_range)[ncol(wound_df_range)] <- "Range"

wound.d12 <-data.frame(cbind(log10(wound_df_range$Range),wound_df_range$D12.FC))
colnames(wound.d12) <- c("Range","D12.FC")

# Axon

axon_matrix <- as.matrix(axon_merged_df_fc[,2:17])
axon_range <- data.frame(rowRanges(axon_matrix))
colnames(axon_range) <- c('Min','Max')
axon_df_range <- cbind(axon_merged_df_fc,axon_range,axon_range[2]-axon_range[1])
colnames(axon_df_range)[ncol(axon_df_range)] <- "Range"

axon.d12 <-data.frame(cbind(log10(axon_df_range$Range),axon_df_range$D12.FC))
colnames(axon.d12) <- c("Range","D12.FC")

# Cell junction

cj_matrix <- as.matrix(cj_merged_df_fc[,2:17])
cj_range <- data.frame(rowRanges(cj_matrix))
colnames(cj_range) <- c('Min','Max')
cj_df_range <- cbind(cj_merged_df_fc,cj_range,cj_range[2]-cj_range[1])
colnames(cj_df_range)[ncol(cj_df_range)] <- "Range"

cj.d12 <-data.frame(cbind(log10(cj_df_range$Range),cj_df_range$D12.FC))
colnames(cj.d12) <- c("Range","D12.FC")

# Cell Adhesion

cell_ad_matrix <- as.matrix(cell_ad_merged_df_fc[,2:17])
cell_ad_range <- data.frame(rowRanges(cell_ad_matrix))
colnames(cell_ad_range) <- c('Min','Max')
cell_ad_df_range <- cbind(cell_ad_merged_df_fc,cell_ad_range,cell_ad_range[2]-cell_ad_range[1])
colnames(cell_ad_df_range)[ncol(cell_ad_df_range)] <- "Range"

cell_ad.d12 <-data.frame(cbind(log10(cell_ad_df_range$Range),cell_ad_df_range$D12.FC))
colnames(cell_ad.d12) <- c("Range","D12.FC")

# Mapk

mapk_matrix <- as.matrix(mapk_merged_df_fc[,2:17])
mapk_range <- data.frame(rowRanges(mapk_matrix))
colnames(mapk_range) <- c('Min','Max')
mapk_df_range <- cbind(mapk_merged_df_fc,mapk_range,mapk_range[2]-mapk_range[1])
colnames(mapk_df_range)[ncol(mapk_df_range)] <- "Range"

mapk.d12 <-data.frame(cbind(log10(mapk_df_range$Range),mapk_df_range$D12.FC))
colnames(mapk.d12) <- c("Range","D12.FC")


# Tyrk

tyrk_matrix <- as.matrix(tyrk_merged_df_fc[,2:17])
tyrk_range <- data.frame(rowRanges(tyrk_matrix))
colnames(tyrk_range) <- c('Min','Max')
tyrk_df_range <- cbind(tyrk_merged_df_fc,tyrk_range,tyrk_range[2]-tyrk_range[1])
colnames(tyrk_df_range)[ncol(tyrk_df_range)] <- "Range"

tyrk.d12 <-data.frame(cbind(log10(tyrk_df_range$Range),tyrk_df_range$D12.FC))
colnames(tyrk.d12) <- c("Range","D12.FC")

# Head 

head_matrix <- as.matrix(head_merged_df_fc[,2:17])
head_range <- data.frame(rowRanges(head_matrix))
colnames(head_range) <- c('Min','Max')
head_df_range <- cbind(head_merged_df_fc,head_range,head_range[2]-head_range[1])
colnames(head_df_range)[ncol(head_df_range)] <- "Range"


# Bar plots 


pathway_ranges <- cbind(mean(actin_range_df$Range),    
                        mean(vasc_df_range$Range),
                        mean(wound_df_range$Range),
                        mean(axon_df_range$Range),
                        mean(cj_df_range$Range),
                        mean(cell_ad_df_range$Range),
                        mean(mapk_df_range$Range),
                        mean(tyrk_df_range$Range),
                        mean(head_df_range$Range))


pathway_minimum <- cbind(mean(actin_range_df$Min),    
                         mean(vasc_df_range$Min),
                         mean(wound_df_range$Min),
                         mean(axon_df_range$Min),
                         mean(cj_df_range$Min),
                         mean(cell_ad_df_range$Min),
                         mean(mapk_df_range$Min),
                         mean(tyrk_df_range$Min),
                         mean(head_df_range$Min))

pathway_maximum <- cbind(mean(actin_range_df$Max),    
                         mean(vasc_df_range$Max),
                         mean(wound_df_range$Max),
                         mean(axon_df_range$Max),
                         mean(cj_df_range$Max),
                         mean(cell_ad_df_range$Max),
                         mean(mapk_df_range$Max),
                         mean(tyrk_df_range$Max),
                         mean(head_df_range$Max))

pathway_min_max <- rbind(log10(pathway_minimum),log10(pathway_maximum))

colnames(pathway_min_max) <- c("Actin","Vascular","Wound","Axon","C Junction","C Adhesion","Mapk","Tyr. Kinase","Head Dev") 


barplot(pathway_min_max,
        main = "Expression levels of Metascape enriched Pathways",
        ylab = 'Log 10 Expression Range',
        xlab = "Pathway",
        col = c("yellow","navy"), yaxp=c(0, 6, 6))


barplot(pathway_ranges,
        main = "Expression levels of Metascape enriched Pathways",
        ylab = 'Log 10 Expression Range',
        xlab = "Pathway", yaxp=c(0, 1300, 4))

#Transcription Reg 

#Range

transcription_reg_df <- as.matrix(transcription_reg_merged[2:17])
transcription_reg_range <- data.frame(rowRanges(transcription_reg_df))
colnames(transcription_reg_range) <- c('Min','Max')
transcription_reg_range_df <- cbind(transcription_reg_range,transcription_reg_range[2]-transcription_reg_range[1])
colnames(transcription_reg_range_df)[ncol(transcription_reg_range_df)] <- "Range"

transcription_reg_range_df <- cbind(transcription_reg_merged,transcription_reg_range_df)


table1 <- t(transcription_reg_range_df$Min)
table2 <- t(transcription_reg_range_df$Max)
table.gene.symbols <- t(transcription_reg_range_df$Gene.Symbols)

transcription_table <- rbind(table1,table2)
colnames(transcription_table) <- table.gene.symbols

#Mean


transcription_reg_df <- as.matrix(transcription_reg_merged[2:17])
transcription_reg_means <- data.frame(rowMeans(transcription_reg_df))

colnames(transcription_reg_means)[ncol(transcription_reg_means)] <- "Mean"

transcription_reg_range_df2 <- cbind(transcription_reg_merged,transcription_reg_means)


table.means <- t(transcription_reg_means)
table.gene.symbols <- t(transcription_reg_range_df$Gene.Symbols)

colnames(table.means) <- table.gene.symbols


table.means.sorted <- t(apply(table.means, 1, sort,decreasing = TRUE))


top20 <- t(data.frame(table.means.sorted[,1:20]))

barplot(top20,
        main = "Mean Expression levels of David Transcription Regulator Gene List",
        ylab = 'Mean Expression',
        xlab = "Pathway",
        col = c("pink","navy"),las=2,yaxp=c(0, 15000, 3))



