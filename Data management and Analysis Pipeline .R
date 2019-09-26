################# Pacakge installation #################
install.packages("readxl")
install.packages('VennDiagram')
install.packages('tidyr')
install.packages('dplyr')
install.packages('stringr')
install.packages('ggplot2')
install.packages('plotly')
install.packages('heatmaply')
install.packages('gplots')
install.packages('cluster')
install.packages('pheatmap')
pkgs <- c("factoextra",  "NbClust")
install.packages(pkgs)
library(factoextra)
library(NbClust)
library(readxl)
library(ggplot2)
library(reshape2)
library(heatmaply)
library(gplots)
library(dplyr)
library(plotly)
require(gridExtra)
library(dplyr)
library(plyr)
library(tidyverse)
library(pheatmap)
library(stringr)
library(tidyr)
library(VennDiagram)
library(cluster)
library(cluster)
library(reshape2)
library(plotly)
library(devtools)
library(ggbiplot)

################# Data Import ################# 
current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) #setting working directory to source file location for data import
filenames <- list.files(current_dir, pattern="*.txt", full.names=TRUE) #read all "txt" files in chosen directory path and add to combined list

file_list_df <- plyr::llply(filenames, read.table,header=FALSE, fill = TRUE) #fill=true for differing sizes, fills gaps with blank spaces

res <- plyr::llply(file_list_df, summary)
names(file_list_df) <- basename(filenames)
names_df <- basename(filenames)
for (file in filenames) {
  
  merged_dataset <- read.table(file, header=TRUE, sep="\t",fill = TRUE)
} # loop to merge into a single dataframe

names(file_list_df) <- gsub("-","_",names(file_list_df)) #replace all "-" instances

names_df <- names(file_list_df)
for (i in 1:length(file_list_df)) {
  assign(paste0(names_df[i]), as.data.frame(file_list_df[[i]]))
} # loop to convert to separate dataframes each named using the original filename 

complete_read_tables <- as.data.frame(names_df[58:59])
DESEQ_tables <- as.data.frame(names_df[1:57])
up_tables <- as.data.frame(names_df[30:57])


#-----------------------------------


################# Up Regulated Gene Lists ################# 

#DOWN Symbols 

IKDNprBIND.Down_4OHT_d1_vs_DMSO_V7 <- strsplit(as.character(IKDNprBIND.Down_4OHT_d1_vs_DMSO.txt$V7),"|",fixed = T)
IKDNprBIND.Down_4OHT_d1_vs_DMSO_V7 <- do.call(rbind.data.frame, IKDNprBIND.Down_4OHT_d1_vs_DMSO_V7)
IKDNprBIND.Down_4OHT_d1_vs_DMSO_GeneSymbols <- IKDNprBIND.Down_4OHT_d1_vs_DMSO_V7[1]
IKDNprBIND.Down_4OHT_d1_vs_DMSO_GeneSymbols <- as.data.frame(IKDNprBIND.Down_4OHT_d1_vs_DMSO_GeneSymbols[-1,])
colnames(IKDNprBIND.Down_4OHT_d1_vs_DMSO_GeneSymbols) = 'Gene.Symbols'


IKDNprBIND.Down_4OHT_d2_vs_DMSO_V7 <- strsplit(as.character(IKDNprBIND.Down_4OHT_d2_vs_DMSO.txt$V7),"|",fixed = T)
IKDNprBIND.Down_4OHT_d2_vs_DMSO_V7 <- do.call(rbind.data.frame, IKDNprBIND.Down_4OHT_d2_vs_DMSO_V7)
IKDNprBIND.Down_4OHT_d2_vs_DMSO_GeneSymbols <- IKDNprBIND.Down_4OHT_d2_vs_DMSO_V7[1]
IKDNprBIND.Down_4OHT_d2_vs_DMSO_GeneSymbols <- as.data.frame(IKDNprBIND.Down_4OHT_d2_vs_DMSO_GeneSymbols[-1,])
colnames(IKDNprBIND.Down_4OHT_d2_vs_DMSO_GeneSymbols) = 'Gene.Symbols'


IKDNprBIND.Down_4OHT_d3_vs_DMSO_V7 <- strsplit(as.character(IKDNprBIND.Down_4OHT_d3_vs_DMSO.txt$V7),"|",fixed = T)
IKDNprBIND.Down_4OHT_d3_vs_DMSO_V7 <- do.call(rbind.data.frame, IKDNprBIND.Down_4OHT_d3_vs_DMSO_V7)
IKDNprBIND.Down_4OHT_d3_vs_DMSO_GeneSymbols <- IKDNprBIND.Down_4OHT_d3_vs_DMSO_V7[1]
IKDNprBIND.Down_4OHT_d3_vs_DMSO_GeneSymbols <- as.data.frame(IKDNprBIND.Down_4OHT_d3_vs_DMSO_GeneSymbols[-1,])
colnames(IKDNprBIND.Down_4OHT_d3_vs_DMSO_GeneSymbols) = 'Gene.Symbols'

IKDNprBIND.Down_4OHT_d4_vs_DMSO_V7 <- strsplit(as.character(IKDNprBIND.Down_4OHT_d4_vs_DMSO.txt$V7),"|",fixed = T)
IKDNprBIND.Down_4OHT_d4_vs_DMSO_V7 <- do.call(rbind.data.frame, IKDNprBIND.Down_4OHT_d4_vs_DMSO_V7)
IKDNprBIND.Down_4OHT_d4_vs_DMSO_GeneSymbols <- IKDNprBIND.Down_4OHT_d4_vs_DMSO_V7[1]
IKDNprBIND.Down_4OHT_d4_vs_DMSO_GeneSymbols <- as.data.frame(IKDNprBIND.Down_4OHT_d4_vs_DMSO_GeneSymbols[-1,])
colnames(IKDNprBIND.Down_4OHT_d4_vs_DMSO_GeneSymbols) = 'Gene.Symbols'

#

IKDNprBIND.Down_4OHT_d6_vs_DMSO_V7 <- strsplit(as.character(IKDNprBIND.Down_4OHT_d6_vs_DMSO.txt$V7),"|",fixed = T)
IKDNprBIND.Down_4OHT_d6_vs_DMSO_V7 <- do.call(rbind.data.frame, IKDNprBIND.Down_4OHT_d6_vs_DMSO_V7)
IKDNprBIND.Down_4OHT_d6_vs_DMSO_GeneSymbols <- IKDNprBIND.Down_4OHT_d6_vs_DMSO_V7[1]
IKDNprBIND.Down_4OHT_d6_vs_DMSO_GeneSymbols <- as.data.frame(IKDNprBIND.Down_4OHT_d6_vs_DMSO_GeneSymbols[-1,])
colnames(IKDNprBIND.Down_4OHT_d6_vs_DMSO_GeneSymbols) = 'Gene.Symbols'

IKDNprBIND.Down_4OHT_d12_vs_DMSO_V7 <- strsplit(as.character(IKDNprBIND.Down_4OHT_d12_vs_DMSO.txt$V7),"|",fixed = T)
IKDNprBIND.Down_4OHT_d12_vs_DMSO_V7 <- do.call(rbind.data.frame, IKDNprBIND.Down_4OHT_d12_vs_DMSO_V7)
IKDNprBIND.Down_4OHT_d12_vs_DMSO_GeneSymbols <- IKDNprBIND.Down_4OHT_d12_vs_DMSO_V7[1]
IKDNprBIND.Down_4OHT_d12_vs_DMSO_GeneSymbols <- as.data.frame(IKDNprBIND.Down_4OHT_d12_vs_DMSO_GeneSymbols[-1,])
colnames(IKDNprBIND.Down_4OHT_d12_vs_DMSO_GeneSymbols) = 'Gene.Symbols'

IKDNprBIND.Down_4OHT_d18_vs_DMSO_V7 <- strsplit(as.character(IKDNprBIND.Down_4OHT_d18_vs_DMSO.txt$V7),"|",fixed = T)
IKDNprBIND.Down_4OHT_d18_vs_DMSO_V7 <- do.call(rbind.data.frame, IKDNprBIND.Down_4OHT_d18_vs_DMSO_V7)
IKDNprBIND.Down_4OHT_d18_vs_DMSO_GeneSymbols <- IKDNprBIND.Down_4OHT_d18_vs_DMSO_V7[1]
IKDNprBIND.Down_4OHT_d18_vs_DMSO_GeneSymbols <- as.data.frame(IKDNprBIND.Down_4OHT_d18_vs_DMSO_GeneSymbols[-1,])
colnames(IKDNprBIND.Down_4OHT_d18_vs_DMSO_GeneSymbols) = 'Gene.Symbols'


#Merge Unique

MyMerge <- function(x, y){
  df <- merge(x, y, by= "Gene.Symbols", all.x= TRUE, all.y= TRUE)
  return(df)
}

down_Gene_Symbols_merged <- Reduce(MyMerge, unique(list(
  IKDNprBIND.Down_4OHT_d1_vs_DMSO_GeneSymbols,
  IKDNprBIND.Down_4OHT_d2_vs_DMSO_GeneSymbols,
  IKDNprBIND.Down_4OHT_d3_vs_DMSO_GeneSymbols,
  IKDNprBIND.Down_4OHT_d4_vs_DMSO_GeneSymbols,
  IKDNprBIND.Down_4OHT_d6_vs_DMSO_GeneSymbols,
  IKDNprBIND.Down_4OHT_d12_vs_DMSO_GeneSymbols,
  IKDNprBIND.Down_4OHT_d18_vs_DMSO_GeneSymbols)))

dim(down_Gene_Symbols_merged)

colnames(IKDNprBINDALL.txt)[7] <- 'Gene.Symbols'

down_all_merged <- merge(IKDNprBINDALL,down_Gene_Symbols_merged,by='Gene.Symbols')

dim(down_all_merged)

#UP Symbols 

IKDNprBIND.Up_4OHT_d1_vs_DMSO_V7 <- strsplit(as.character(IKDNprBIND.Up_4OHT_d1_vs_DMSO.txt$V7),"|",fixed = T)
IKDNprBIND.Up_4OHT_d1_vs_DMSO_V7 <- do.call(rbind.data.frame, IKDNprBIND.Up_4OHT_d1_vs_DMSO_V7)
IKDNprBIND.Up_4OHT_d1_vs_DMSO_GeneSymbols <- IKDNprBIND.Up_4OHT_d1_vs_DMSO_V7[1]
IKDNprBIND.Up_4OHT_d1_vs_DMSO_GeneSymbols <- as.data.frame(unique(IKDNprBIND.Up_4OHT_d1_vs_DMSO_GeneSymbols[-1,]))
colnames(IKDNprBIND.Up_4OHT_d1_vs_DMSO_GeneSymbols) = 'Gene.Symbols'

IKDNprBIND.Up_4OHT_d2_vs_DMSO_V7 <- strsplit(as.character(IKDNprBIND.Up_4OHT_d2_vs_DMSO.txt$V7),"|",fixed = T)
IKDNprBIND.Up_4OHT_d2_vs_DMSO_V7 <- do.call(rbind.data.frame, IKDNprBIND.Up_4OHT_d2_vs_DMSO_V7)
IKDNprBIND.Up_4OHT_d2_vs_DMSO_GeneSymbols <- IKDNprBIND.Up_4OHT_d2_vs_DMSO_V7[1]
IKDNprBIND.Up_4OHT_d2_vs_DMSO_GeneSymbols <- as.data.frame(unique(IKDNprBIND.Up_4OHT_d2_vs_DMSO_GeneSymbols[-1,]))
colnames(IKDNprBIND.Up_4OHT_d2_vs_DMSO_GeneSymbols) = 'Gene.Symbols'

IKDNprBIND.Up_4OHT_d3_vs_DMSO_V7 <- strsplit(as.character(IKDNprBIND.Up_4OHT_d3_vs_DMSO.txt$V7),"|",fixed = T)
IKDNprBIND.Up_4OHT_d3_vs_DMSO_V7 <- do.call(rbind.data.frame, IKDNprBIND.Up_4OHT_d3_vs_DMSO_V7)
IKDNprBIND.Up_4OHT_d3_vs_DMSO_GeneSymbols <- IKDNprBIND.Up_4OHT_d3_vs_DMSO_V7[1]
IKDNprBIND.Up_4OHT_d3_vs_DMSO_GeneSymbols <- as.data.frame(unique(IKDNprBIND.Up_4OHT_d3_vs_DMSO_GeneSymbols[-1,]))
colnames(IKDNprBIND.Up_4OHT_d3_vs_DMSO_GeneSymbols) = 'Gene.Symbols'

IKDNprBIND.Up_4OHT_d4_vs_DMSO_V7 <- strsplit(as.character(IKDNprBIND.Up_4OHT_d4_vs_DMSO.txt$V7),"|",fixed = T)
IKDNprBIND.Up_4OHT_d4_vs_DMSO_V7 <- do.call(rbind.data.frame, IKDNprBIND.Up_4OHT_d4_vs_DMSO_V7)
IKDNprBIND.Up_4OHT_d4_vs_DMSO_GeneSymbols <- IKDNprBIND.Up_4OHT_d4_vs_DMSO_V7[1]
IKDNprBIND.Up_4OHT_d4_vs_DMSO_GeneSymbols <- as.data.frame(unique(IKDNprBIND.Up_4OHT_d4_vs_DMSO_GeneSymbols[-1,]))
colnames(IKDNprBIND.Up_4OHT_d4_vs_DMSO_GeneSymbols) = 'Gene.Symbols'

#

IKDNprBIND.Up_4OHT_d6_vs_DMSO_V7 <- strsplit(as.character(IKDNprBIND.Up_4OHT_d6_vs_DMSO.txt$V7),"|",fixed = T)
IKDNprBIND.Up_4OHT_d6_vs_DMSO_V7 <- do.call(rbind.data.frame, IKDNprBIND.Up_4OHT_d6_vs_DMSO_V7)
IKDNprBIND.Up_4OHT_d6_vs_DMSO_GeneSymbols <- IKDNprBIND.Up_4OHT_d6_vs_DMSO_V7[1]
IKDNprBIND.Up_4OHT_d6_vs_DMSO_GeneSymbols <- as.data.frame(unique(IKDNprBIND.Up_4OHT_d6_vs_DMSO_GeneSymbols[-1,]))
colnames(IKDNprBIND.Up_4OHT_d6_vs_DMSO_GeneSymbols) = 'Gene.Symbols'

IKDNprBIND.Up_4OHT_d12_vs_DMSO_V7 <- strsplit(as.character(IKDNprBIND.Up_4OHT_d12_vs_DMSO.txt$V7),"|",fixed = T)
IKDNprBIND.Up_4OHT_d12_vs_DMSO_V7 <- do.call(rbind.data.frame, IKDNprBIND.Up_4OHT_d12_vs_DMSO_V7)
IKDNprBIND.Up_4OHT_d12_vs_DMSO_GeneSymbols <- IKDNprBIND.Up_4OHT_d12_vs_DMSO_V7[1]
IKDNprBIND.Up_4OHT_d12_vs_DMSO_GeneSymbols <- as.data.frame(unique(IKDNprBIND.Up_4OHT_d12_vs_DMSO_GeneSymbols[-1,]))
colnames(IKDNprBIND.Up_4OHT_d12_vs_DMSO_GeneSymbols) = 'Gene.Symbols'

IKDNprBIND.Up_4OHT_d18_vs_DMSO_V7 <- strsplit(as.character(IKDNprBIND.Up_4OHT_d18_vs_DMSO.txt$V7),"|",fixed = T)
IKDNprBIND.Up_4OHT_d18_vs_DMSO_V7 <- do.call(rbind.data.frame, IKDNprBIND.Up_4OHT_d18_vs_DMSO_V7)
IKDNprBIND.Up_4OHT_d18_vs_DMSO_GeneSymbols <- IKDNprBIND.Up_4OHT_d18_vs_DMSO_V7[1]
IKDNprBIND.Up_4OHT_d18_vs_DMSO_GeneSymbols <- as.data.frame(unique(IKDNprBIND.Up_4OHT_d18_vs_DMSO_GeneSymbols[-1,]))
colnames(IKDNprBIND.Up_4OHT_d18_vs_DMSO_GeneSymbols) = 'Gene.Symbols'


up_Gene_Symbols_merged <- Reduce(MyMerge,(list(
  IKDNprBIND.Up_4OHT_d1_vs_DMSO_GeneSymbols,
  IKDNprBIND.Up_4OHT_d2_vs_DMSO_GeneSymbols,
  IKDNprBIND.Up_4OHT_d3_vs_DMSO_GeneSymbols,
  IKDNprBIND.Up_4OHT_d4_vs_DMSO_GeneSymbols,
  IKDNprBIND.Up_4OHT_d6_vs_DMSO_GeneSymbols,
  IKDNprBIND.Up_4OHT_d12_vs_DMSO_GeneSymbols,
  IKDNprBIND.Up_4OHT_d18_vs_DMSO_GeneSymbols)))

up_Gene_Symbols_merged <- unique(up_Gene_Symbols_merged)

dim(up_Gene_Symbols_merged)

colnames(IKDNprBINDALL.txt)[7] <- 'Gene.Symbols'

UP_all_merged <- merge(IKDNprBINDALL,up_Gene_Symbols_merged,by='Gene.Symbols')

write.csv(UP_all_merged,"UP_all_merged.csv")

################# Up Regulated Venn Diagram ################# 

#D1-4
venn.diagram(main= "Up Regulated Genes From D1-D4 Samples",    main.pos= c( 0.5, 1.05),    main.just= c( 0.5, 1),    sub.pos= c( 0.5, 1.05),    sub.just= c( 0.5, 1),   lwd= 2,    lty= "solid",    col= "black",    fill= c( "red1", "orange","royalblue", "purple"),    alpha = 0.5,   rotation.degree= 0,    rotation.centre= c(0.5, 0.5),    label.col= "black",    cex= 2,    fontface= "plain",    fontfamily= "serif",    category.names= c("D1", "D2", "D3", "D4"), cat.dist= 0.1,    cat.cex= 1,    cat.col= "black",    cat.fontface= "plain",    cat.fontfamily= "serif",    cat.prompts= FALSE, ext.text= FALSE,    ext.pos= "",    ext.percent= "",    ext.line.lwd= "",    ext.line.lty= "",    ext.dist= "",    ext.length= "", euler.d= TRUE,    scaled= TRUE,    sep.dist= 0.05,    offset = 0,    height= 6,    width= 6,    resolution= 1000, description= "",    x = list(x1 = IKDNprBIND.Up_4OHT_d1_vs_DMSO_GeneSymbols, x2= IKDNprBIND.Up_4OHT_d2_vs_DMSO_GeneSymbols, x3= IKDNprBIND.Up_4OHT_d3_vs_DMSO_GeneSymbols, x4=IKDNprBIND.Up_4OHT_d4_vs_DMSO_GeneSymbols),    units= "in",    filename= "UpRegulated_D1toD4.tiff" )

#D6-18

venn.diagram(main= "Up Regulated Genes From D4-D18 Samples",    main.pos= c( 0.5, 1.05),    main.just= c( 0.5, 1),    sub.pos= c( 0.5, 1.05),    sub.just= c( 0.5, 1),   lwd= 2,    lty= "solid",    col= "black",    fill= c( "coral", "gold","blue", "lavender"),    alpha = 0.5,   rotation.degree= 0,    rotation.centre= c(0.5, 0.5),    label.col= "black",    cex= 2,    fontface= "plain",    fontfamily= "serif",    category.names= c("D4", "D6", "D12", "D18"), cat.dist= 0.1,    cat.cex= 1,    cat.col= "black",    cat.fontface= "plain",    cat.fontfamily= "serif",    cat.prompts= FALSE, ext.text= FALSE,    ext.pos= "",    ext.percent= "",    ext.line.lwd= "",    ext.line.lty= "",    ext.dist= "",    ext.length= "", euler.d= TRUE,    scaled= TRUE,    sep.dist= 0.05,    offset = 0,    height= 6,    width= 6,    resolution= 1000, description= "",    x = list(x1 = IKDNprBIND.Up_4OHT_d4_vs_DMSO_GeneSymbols, x2= IKDNprBIND.Up_4OHT_d6_vs_DMSO_GeneSymbols, x3= IKDNprBIND.Up_4OHT_d12_vs_DMSO_GeneSymbols, x4=IKDNprBIND.Up_4OHT_d18_vs_DMSO_GeneSymbols),    units= "in",    filename= "UpRegulated_D4toD18.tiff" )

################# Upregulated Overlap lists ################# 
overlap1_2 <- calculate.overlap(x = list(x1 = IKDNprBIND.Up_4OHT_d1_vs_DMSO_GeneSymbols, x2= IKDNprBIND.Up_4OHT_d2_vs_DMSO_GeneSymbols))
sum1 <- na.omit(as.data.frame(IKDNprBIND.Up_4OHT_d1_vs_DMSO_GeneSymbols))
sum2 <- na.omit(as.data.frame(IKDNprBIND.Up_4OHT_d2_vs_DMSO_GeneSymbols))

unique1_2 <- unique(overlap1_2[3])
capture.output(as.data.frame(unique1_2), file = "unique1_2.txt")
capture.output(as.data.frame(sum1), file = "sum1.txt")
capture.output(as.data.frame(sum2), file = "sum2.txt")

overlap2_3 <- calculate.overlap(x = list(x1 = IKDNprBIND.Up_4OHT_d2_vs_DMSO_GeneSymbols, x2= IKDNprBIND.Up_4OHT_d3_vs_DMSO_GeneSymbols))
sum3 <- na.omit(as.data.frame(IKDNprBIND.Up_4OHT_d3_vs_DMSO_GeneSymbols))
unique2_3 <- unique(overlap2_3[3])
capture.output(as.data.frame(unique2_3), file = "unique2_3.txt")
capture.output(as.data.frame(sum3), file = "sum3.txt")
length(sum3)
overlap3_4 <- calculate.overlap(x = list(x1 = IKDNprBIND.Up_4OHT_d3_vs_DMSO_GeneSymbols, x2= IKDNprBIND.Up_4OHT_d4_vs_DMSO_GeneSymbols))
sum4 <- na.omit(as.data.frame(IKDNprBIND.Up_4OHT_d4_vs_DMSO_GeneSymbols))
unique3_4 <- unique(overlap3_4[3])
capture.output(as.data.frame(unique3_4), file = "unique3_4.txt")
capture.output(as.data.frame(sum4), file = "sum4.txt")

overlap4_6 <- calculate.overlap(x = list(x1 = IKDNprBIND.Up_4OHT_d4_vs_DMSO_GeneSymbols, x2= IKDNprBIND.Up_4OHT_d6_vs_DMSO_GeneSymbols))
sum6 <- na.omit(as.data.frame(IKDNprBIND.Up_4OHT_d6_vs_DMSO_GeneSymbols))
overlap4_6
unique4_6  <- unique(overlap4_6 [3])
capture.output(as.data.frame(unique4_6), file = "unique4_6.txt")
capture.output(as.data.frame(sum6), file = "sum6.txt")

overlap6_12 <- calculate.overlap(x = list(x1 = IKDNprBIND.Up_4OHT_d6_vs_DMSO_GeneSymbols, x2= IKDNprBIND.Up_4OHT_d12_vs_DMSO_GeneSymbols))
sum6 <- na.omit(as.data.frame(IKDNprBIND.Up_4OHT_d6_vs_DMSO_GeneSymbols))
unique6_12 <- unique(overlap6_12[3])
capture.output(as.data.frame(unique6_12), file = "unique6_12.txt")
capture.output(as.data.frame(sum12), file = "sum12.txt")

overlap12_18 <- calculate.overlap(x = list(x1 = IKDNprBIND.Up_4OHT_d12_vs_DMSO_GeneSymbols, x2= IKDNprBIND.Up_4OHT_d18_vs_DMSO_GeneSymbols))
sum12 <- na.omit(as.data.frame(IKDNprBIND.Up_4OHT_d12_vs_DMSO_GeneSymbols))
sum18 <- na.omit(as.data.frame(IKDNprBIND.Up_4OHT_d18_vs_DMSO_GeneSymbols))
unique12_18 <- unique(overlap12_18[3])
capture.output(as.data.frame(unique12_18), file = "unique12_18.txt")
capture.output(as.data.frame(sum18), file = "sum18.txt")

################# Extracting Unique Entries ################# 

#Function to merge unique by gene symbol 



write.csv(UP_all_merged, file = "UP_all_merged.csv")


#using the already created list 

clusterall <- UP_all_merged[c(1,9:24)]



colnames(clusterall) <- c("GeneSymbols","DMSO.1","DMSO.2","D1.1","D1.2","D2.1","D2.2","D3.1","D3.2","D4.1","D4.2","D6.1","D6.2","D12.1","D12.2","D18.1","D18.2") 


#-----------------------------------
################# Clustering  ################# 

# x: A data matrix or data frame or dissimilarity matrix
# k: The desired number of clusters to be generated
# metric: Metric for calculating dissimilarities between observations
# stand: If TRUE, variables are standardized before calculating the dissimilarities

unlogged_clusterall_df <- clusterall %>% remove_rownames %>% column_to_rownames(var="GeneSymbols")
clusterall_df <- log10(unlogged_clusterall_df)
df <- clusterall_df %>%  select_if(is.numeric)
scaled_df <- as.matrix(df)
kmm <- kmeans(scaled_df,10,nstart = 50,iter.max = 15) #we keep number of iter.max=15 to ensure the algorithm converges and nstart=50 to #ensure that atleat 50 random sets are choosen  
# between_SS / total_SS =  71.2 %
scaled_df[!is.finite(scaled_df)] <- 0


fviz_nbclust(scaled_df, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")
#4 is the optimum

#K means
set.seed(123)
scaled_df[!rowSums(!is.finite(scaled_df)),]
scaled_df[!is.finite(scaled_df)] <- 0
scaled_df_kmeans <- scale(scaled_df)
k2 <- kmeans(scaled_df_kmeans, centers = 10, nstart = 25)
df %>%
  mutate(Cluster = k2$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")
dfc <- cbind(df, id=clusterall$GeneSymbols, cluster=k2$cluster)
dfc$idsort <- dfc$id[order(dfc$cluster)]
dfc$idsort <- order(dfc$idsort)
dfm <- melt(dfc, id.vars=c("cluster", "id"))
dfm$value[!is.finite(dfm$value)] <- 0

png(filename="kmeans_heatmap.png",width = 800 , height = 3000)
heatmap(as.matrix(df)[(k2$cluster),],Rowv=NA,Colv=NA,legend=c('col'))
dev.off()


any(duplicated(rownames(k2)))
k2m <- as.matrix(df)[(k2$cluster),]
k2m[!is.finite(k2m)] <- 0
rownames(k2m) <- make.names(clusterall[,1], unique = TRUE)
heatmaply(k2m,Rowv=NA,trace = "none",xlab = "Clusters", ylab = "Gene Symbols", 
          main = "K=10 K-means Cluster",
          margins = c(60,100,40,20))
k2df <- as.matrix(df)[(k2$cluster),]
rownames(k2df) <- make.names(clusterall[,1], unique = TRUE)

png(filename="K104HeatMaply.png",width = 800, height = 1000)
heatmaply(k2df$cluster,Rowv=NA,trace = "none",xlab = "Clusters", ylab = "Gene Symbols", 
          main = "K-means K=10 Cluster",
          margins = c(60,100,40,20))
dev.off()

png(filename="kmeans_centers_heatmap.png",width = 800, height = 2000)
kmeans_centers_heatmap <- heatmap.2(k2m, col=redgreen(75), trace=c("none"))
dev.off()


################### Tables for clusters 

kmeans_table <- dfm[c(2,1)]

kmeans_table <- kmeans_table %>% 
  group_by(cluster) %>% 
  mutate(grouped_id = row_number())


gene_list_1 <- filter(kmeans_table, cluster == 1)
gene_list_2 <- filter(kmeans_table, cluster == 2)
gene_list_3 <- filter(kmeans_table, cluster == 3)
gene_list_4 <- filter(kmeans_table, cluster == 4)
gene_list_5 <- filter(kmeans_table, cluster == 5)
gene_list_6 <- filter(kmeans_table, cluster == 6)
gene_list_7 <- filter(kmeans_table, cluster == 7)
gene_list_8 <- filter(kmeans_table, cluster == 8)
gene_list_9 <- filter(kmeans_table, cluster == 9)
gene_list_10 <- filter(kmeans_table, cluster == 10)


################# Pearson correlation matrix ################### 

up_and_down_df <- IKDNprBINDALL[c(1,9:24)]
colnames(up_and_down_df) <-  c("GeneSymbols","DMSO.1","DMSO.2","D1.1","D1.2","D2.1","D2.2","D3.1","D3.2","D4.1","D4.2","D6.1","D6.2","D12.1","D12.2","D18.1","D18.2") 
c <- round(cor(up_and_down_df[2:17]),2) # Pearson Correlation




# Get lower triangle of the correlation matrix
get_lower_tri<-function(c){
  c[upper.tri(c)] <- NA
  return(c)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(c){
  c[lower.tri(c)]<- NA
  return(c)
}

upper_tri <- get_upper_tri(c)

melted_cormat <- melt(upper_tri, na.rm = TRUE)

ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "navy", mid = "blue", 
                       midpoint = 0.95, limit = c(0.90,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

png(filename="Pearson.png")
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "white", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

dev.off()
################# PCA ############# 

colnames(up_and_down_df) <-  c("GeneSymbols","DMSO.1","DMSO.2","D1.1","D1.2","D2.1","D2.2","D3.1","D3.2","D4.1","D4.2","D6.1","D6.2","D12.1","D12.2","D18.1","D18.2") 

cluster_pca<-prcomp(t(up_and_down_df[2:17]))


# Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component, looking for elbow to find msot accurate PC

fviz_eig(cluster_pca) #extract and visualise eigen values of dimensions: plotting eigen values/variances against number of dimensions 
help(fviz_eig)
png(filename="PCA.png")
fviz_pca_ind(cluster_pca,
             axes = c(1, 2), #Choosing PCs 
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("blue", "navy", "black"),
             repel = TRUE     # Avoid text overlapping
)
dev.off()
scores = as.data.frame(cluster_pca$x) 

#3D plot of first three Prinicpal Components
png(filename="PCA3D.png")
plot_ly(x=scores$PC1, y=scores$PC2, z=scores$PC3, type="scatter3d", text = row.names(scores)) %>% add_markers() %>%  add_text(textposition = 'topleft') %>%  add_lines(name = ~"scatter") 
dev.off()
################# Gene count bar plots  #####################

############### TEMPORAL (Compared to DMSO)


DESEQ_tables # (list of differentially expressed tables)

change_in_genes_down_tempo <- cbind(nrow(IKDNprBIND.Down_4OHT_d1_vs_DMSO.txt),    
                                    nrow(IKDNprBIND.Down_4OHT_d2_vs_DMSO.txt),
                                    nrow(IKDNprBIND.Down_4OHT_d3_vs_DMSO.txt),
                                    nrow(IKDNprBIND.Down_4OHT_d4_vs_DMSO.txt),
                                    nrow(IKDNprBIND.Down_4OHT_d6_vs_DMSO.txt),
                                    nrow(IKDNprBIND.Down_4OHT_d12_vs_DMSO.txt),
                                    nrow(IKDNprBIND.Down_4OHT_d18_vs_DMSO.txt))



change_in_genes_up_tempo <- cbind(nrow(IKDNprBIND.Up_4OHT_d1_vs_DMSO.txt),    
                                  nrow(IKDNprBIND.Up_4OHT_d2_vs_DMSO.txt),
                                  nrow(IKDNprBIND.Up_4OHT_d3_vs_DMSO.txt),
                                  nrow(IKDNprBIND.Up_4OHT_d4_vs_DMSO.txt),
                                  nrow(IKDNprBIND.Up_4OHT_d6_vs_DMSO.txt),
                                  nrow(IKDNprBIND.Up_4OHT_d12_vs_DMSO.txt),
                                  nrow(IKDNprBIND.Up_4OHT_d18_vs_DMSO.txt))

changes_in_genes_tempo <- rbind(change_in_genes_up_tempo,change_in_genes_down_tempo)

colnames(changes_in_genes_tempo) <- c("D1vDMSO","D2vDMSO","D3vDMSO","D4vDMSO","D6vDMSO","D12vDMSO","D18vDMSO") 
colnames(changes_in_genes_tempo)

png(filename="Genecount.png")
barplot(changes_in_genes_tempo,
        main = "Number of up and down-regulated genes across DMSO/D0-D18 ",
        ylab = 'Number of Genes',
        xlab = "Time",
        col = c("red","green"), legend = colnames(change_in_genes),beside=TRUE)
legend("topleft",
       c("UP","DOWN"),
       fill = c("red","green")
)
dev.off()


################ DIFFERENTIAL plots
#create a dataframe with the change in gene expression for the 7 time points down-regulated DESEQ files
change_in_genes_down <- cbind(nrow(IKDNprBIND.Down_4OHT_d1_vs_DMSO.txt),    
                              nrow(IKDNprBIND.Down_4OHT_d2_vs_4OHT_d1.txt),
                              nrow(IKDNprBIND.Down_4OHT_d3_vs_4OHT_d2.txt),
                              nrow(IKDNprBIND.Down_4OHT_d4_vs_4OHT_d3.txt),
                              nrow(IKDNprBIND.Down_4OHT_d6_vs_4OHT_d4.txt),
                              nrow(IKDNprBIND.Down_4OHT_d12_vs_4OHT_d6.txt),
                              nrow(IKDNprBIND.Down_4OHT_d18_vs_4OHT_d12.txt))

change_in_genes_up <- cbind(nrow(IKDNprBIND.Up_4OHT_d1_vs_DMSO.txt),    
                            nrow(IKDNprBIND.Up_4OHT_d2_vs_4OHT_d1.txt),
                            nrow(IKDNprBIND.Up_4OHT_d3_vs_4OHT_d2.txt),
                            nrow(IKDNprBIND.Up_4OHT_d4_vs_4OHT_d3.txt),
                            nrow(IKDNprBIND.Up_4OHT_d6_vs_4OHT_d4.txt),
                            nrow(IKDNprBIND.Up_4OHT_d12_vs_4OHT_d6.txt),
                            nrow(IKDNprBIND.Up_4OHT_d18_vs_4OHT_d12.txt))

changes_in_genes <- rbind(change_in_genes_up,change_in_genes_down)

colnames(changes_in_genes) <- c("D1vDMSO","D2vD1","D3vD2","D4vD3","D6vD4","D12vD6","D18vD12") 
colnames(changes_in_genes)

png(filename="Rageofchange,png")
barplot(changes_in_genes,
        main = "Rate of change in number of up and down-regulated genes across DMSO/D0-D18 ",
        ylab = 'Number of Genes',
        xlab = "Time",
        col = c("red","green"), legend = colnames(change_in_genes),beside=TRUE)
legend("topleft",
       c("UP","DOWN"),
       fill = c("red","green")
)
dev.off()

#-----------------------------------


