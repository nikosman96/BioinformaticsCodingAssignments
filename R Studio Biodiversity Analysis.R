################### Packages ####################
install.packages("bold")
install.packages("rentrez")
install.packages("vegan")
install.packages("seqinr")
install.packages("stringr")
install.packages("tidyverse")
biocLite("msa")
browseVignettes("msa")
biocLite("muscle")
biocLite("DECIPHER")
install.packages("stringi")
install.packages("ape")
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)
install.packages("dplyr")
install.packages("devtools")
devtools::install_github("ropensci/rentrez")
install.packages("kmer")
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
install.packages("phangorn")
install.packages("picante")
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
install.packages("PhyloMeasures")
install.packages("seqinr")
devtools::install_github("zachcp/phylogeo")
install.packages("rich")

library(seqinr)
library(phytools)
library(mapdata)
library(picante)
library(phytools)
library(Biostrings)
library(seqinr)
library(rentrez)
library(stringr)
library(tidyverse)
library("msa")
library(Biostrings)
library(muscle)
library(msa)
library(DECIPHER)
library(vegan)
library(ape)
library(stringi)
library(stats)
library(ade4)
library("rentrez")
library("kmer")
library("dplyr")
library(phangorn)
library("bold")
library(raster)
library(sp)
library(ggplot2)
library(ape)
library(ggtree)
library("PhyloMeasures")
library(rich)
library(ggplot2)
library(gridExtra)
library(phylogeo)
library(phytools)

################### Dataset  #############################

#Genomic and regional data retrieved from Bold
#Data input 

formicidae_bold<- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=formicidae&format=tsv") #Bold api loading data from the formicidae family
write_tsv(formicidae_bold, 'formicidae.tsv') #hard writing data in directory for quick access

formicidae_bold<- read_tsv('formicidae.tsv')
################### Quality Control + Data Exploration ################################

summary(formicidae_bold) #Explore data by looking at summary
length(unique(formicidae_bold$species_name)) #How many unqieu species are present

#using select function to extract wanted columns
formicidae_bold %>% dplyr::select(processid, markercode, bin_uri, species_name, nucleotides, region, country, province_state, elev, lat, lon) -> formicidae_bold_extracted 

#sample 1 sequence per species to narrow down data using dplyr sample function
library(dplyr); formicidae_bold_extracted %>% group_by(species_name) %>% sample_n(1) -> formicidae_bold_grouped 

formicidae_bold_df <- as.data.frame(formicidae_bold_grouped) #convert to data frame

sum(is.na(formicidae_bold_df$bin_uri)) #check for missing values
formicidae_bold_df <- na.omit(formicidae_bold_df) #and remove them

#checking country data 
formicidae_bold_df %>%
  group_by(country) %>%
  summarize(count = length(processid)) %>%
  arrange(desc(count))

#Drop states with inadequate results states 
formicidae_bins <-formicidae_bold_df[which(formicidae_bold_df$country == "Madagascar" | formicidae_bold_df$country == "Canada" | formicidae_bold_df$country == "United States" ), ]

#check again to see what's left: three biggest
formicidae_bins %>%
  group_by(country) %>%
  summarize(count = length(processid)) %>%
  arrange(desc(count))

#most well represented markercode? COI-5p 188 compared to 12s - 4 and h3 - 3 

formicidae_bold_df %>%
  group_by(markercode) %>% 
  summarize(n = length(processid)) %>%
  arrange(desc(n)) %>%
  print()

################### First visualisation: Spec Accum Curves ############################## 

#Spec Accum curve for biological diversity over sampling effort (country)

formicidae_bins_comm <- formicidae_bins %>% #group by country and bin_uri
  group_by(country, bin_uri) %>%
  count(bin_uri)

formicidae_bins_comm <- spread(formicidae_bins_comm, bin_uri, n) #spread data into community object by bin uri and frequency

formicidae_bins_comm[is.na(formicidae_bins_comm)] <- 0 #set NA values to 0

formicidae_bins_comm <- formicidae_bins_comm %>% #remove reownames
  remove_rownames %>% column_to_rownames(var="country")

formicidae_bins_comm <- specaccum(formicidae_bins_comm, method = "collector") #speccaccum 

plot(formicidae_bins_comm,xlab="Number of sampling sites (countries)", ylab="Biological Diversity")

#Spec Accum curve for biological diversity over sampling effort (region)

formicidae_bins_comm2 <- formicidae_bins %>% #group by country and bin_uri
  group_by(region, bin_uri) %>%
  count(bin_uri)

formicidae_bins_comm2 <- spread(formicidae_bins_comm2, bin_uri, n) #spread data into community object by bin uri and frequency

formicidae_bins_comm2[is.na(formicidae_bins_comm2)] <- 0 #set NA values to 0

formicidae_bins_comm2 <- formicidae_bins_comm2 %>% #remove reownames
  remove_rownames %>% column_to_rownames(var="region")

formicidae_bins_comm2 <- specaccum(formicidae_bins_comm2,method = "collector") # perform speccaccum 

plot(formicidae_bins_comm2,xlab="Number of sampling sites (regions)", ylab="Biological Diversity")

#################### Second visualisation: Wolrd Climate Heat Map ############################## 

r <- getData("worldclim",var="bio",res=10) #extract data from object using getData as well as geohraphical data in raster package 
r <- r[[c(1,12)]] #2 parameters indexed
names(r) <- c("Temp","Prec") #named 

lats <- formicidae_bins$lat #assining latitude and longitude values to an object
lons <- formicidae_bins$lon 

coords <- data.frame(x=lons,y=lats)
coords <- na.omit(coords) #remove NA values

points <- SpatialPoints(coords, proj4string = r@crs) #assign points

values <- extract(r,points)

coord_df <- cbind.data.frame(coordinates(points),values) #create data frame of coordinates

plot(r[[1]]) %>%
plot(points,add=T) %>%
title("Annual Mean Temperature")

#################### Third visualisation: Biodiversity as Species Richness  ######################### 

formicidae_country_species <- formicidae_bins %>%
  group_by(region, species_name) %>%
  count(species_name)

formicidae_country_species <- spread(formicidae_country_species, species_name, n)
formicidae_country_species[is.na(formicidae_country_species)] <- 0 #setting to 0 
formicidae_country_species <- as.data.frame(formicidae_country_species)

formicidae_matrix <- data.matrix(formicidae_country_species)
formicidae_rich <- rich(formicidae_matrix, nrandom=100, verbose=TRUE)

formicidae_rich$cr # observed cumulative species richness
formicidae_rich$mr # observed mean value of species richness over the n samples

#subsetting data again for each country 
formicidae_mad<- data.matrix(formicidae_country_species[which(formicidae_bins$country == "Madagascar"), ]) formicidae_mad[is.na(formicidae_mad)] <- 0 
formicidae_can <- data.matrix(formicidae_country_species[which(formicidae_bins$country ==  "Canada"), ])
formicidae_can[is.na(formicidae_can)] <- 0 
formicidae_us <- data.matrix(formicidae_country_species[which(formicidae_bins$country == "United States" ), ])
formicidae_us[is.na(formicidae_us)] <- 0 

#species richness using rarc
mad_rarc <- (rarc(formicidae_mad))
can_rarc<- (rarc(formicidae_can))
us_rarc <- (rarc(formicidae_us))

#mean species rirchness
mean(mad_rarc[["out"]][["mean.richness"]])
mean(can_rarc[["out"]][["mean.richness"]])
mean(us_rarc[["out"]][["mean.richness"]])

#plotting rarc function against number of regions in each country
plot(mad_rarc$out[,6],mad_rarc$out[,1], type="b", ylim=range(c(mad_rarc$out[,2],mad_rarc$out[,3])), xlab="Number of sampling units/regions in Madagascar", ylab=" Species richness")
plot(can_rarc$out[,6],can_rarc$out[,1], type="b", ylim=range(c(can_rarc$out[,2],can_rarc$out[,3])), xlab="Number of sampling units (regions) in Canada", ylab="Species richness")
plot(us_rarc$out[,6],us_rarc$out[,1], type="b", ylim=range(c(us_rarc$out[,2],us_rarc$out[,3])), xlab="Number of sampling units (regions) in US", ylab="Species richness")



#################### Fourth visualisation: Species Richness against Climate  ######################### 

#subsetting country data and assigning them to new variables 
formicidae_mad_lat<- formicidae_bins[which(formicidae_bins$country == "Madagascar"), ]
formicidae_can_lat <-formicidae_bins[which(formicidae_bins$country ==  "Canada"), ] 
formicidae_us_lat <- formicidae_bins[which(formicidae_bins$country == "United States" ), ]

#mean annual temperatures 
mean(us_coord_df$Temp)
mean(can_coord_df$Temp)
mean(mad_coord_df$Temp)

##### US #####
uslats <- formicidae_us_lat$lat #assining latitude and longitude values to new object
uslons <- formicidae_us_lat$lon 

uscoords <- data.frame(x=uslons,y=uslats )
uscoords <- na.omit(uscoords) #remove NA values

uspoints <- SpatialPoints(uscoords, proj4string = r@crs) #assign points

usvalues <- extract(r,uspoints)

us_coord_df <- cbind.data.frame(coordinates(uspoints),usvalues) 

us_coord_df$Temp #subset temperature data retrieved through worldclim 

#plot with mean richness on y axis against temperature on x-axis against 
us_plot <- plot(x=us_coord_df$Temp, y=us_rarc[["out"]][["mean.richness"]],xlab="Mean Annual Temperature in US", ylab="Species richness")


##### CAN #####
canlats <- formicidae_can_lat$lat #assining latitude and longitude values to an object
canlons <- formicidae_can_lat$lon 

cancoords <- data.frame(x=canlons,y=canlats )
cancoords <- na.omit(cancoords) #remove NA values

canpoints <- SpatialPoints(cancoords, proj4string = r@crs) #assign points

canvalues <- extract(r,canpoints)

can_coord_df <- cbind.data.frame(coordinates(canpoints),canvalues) 

can_coord_df$Temp 

can_plot <- plot(x=can_coord_df$Temp, y=can_rarc[["out"]][["mean.richness"]],xlab="Mean Annual Temperature in Canada", ylab="Species richness")




##### MAD #####
madlats <- formicidae_mad_lat$lat #assining latitude and longitude values to an object
madlons <- formicidae_mad_lat$lon 

madcoords <- data.frame(x=madlons,y=madlats )
madcoords <- na.omit(madcoords) #remove NA values

madpoints <- SpatialPoints(madcoords, proj4string = r@crs) #assign points

madvalues <- extract(r,madpoints)

mad_coord_df <- cbind.data.frame(coordinates(madpoints),madvalues) 

mad_coord_df$Temp 

mad_plot <- plot(x=mad_coord_df$Temp, y=mad_rarc[["out"]][["mean.richness"]],xlab="Mean Annual Temperature in Madagascar", ylab="Species richness")


###################### Phylogenetic Analysis ######################

#Sample species for clustering
formicidae_sample <- sample_n(formicidae_bins,15)

formicidae_sample$nucleotides <- DNAStringSet(formicidae_sample$nucleotides)

MSA <- DNAStringSet(muscle::muscle(formicidae_sample$nucleotides,  maxiters = 2, diags = True))

dnaBin_formicidae <- as.DNAbin(MSA)

distance_matrix <- dist.dna(dnaBin_formicidae, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE) #creating a a pairwise distance matrix using the dist.dna function

rownames(distance_matrix) <- c(formicidae_sample$species_name) #To add tip labels rename columns as species_names from sample 

# Phylo <- ladderize(nj(distance_matrix))
# plot(Phylo, cex = 1,tip.color=rainbow(5)) 
# title("Muscle Alignment Dendrogram with coloured sets")

Phylo2 <-untangle(upgma(distance_matrix))
plot(Phylo2, cex = 1,tip.color=rainbow(5)) 
title("Muscle Alignment Dendrogram with coloured sets")

###################### Bonus ######################

#coordinates
phylo_lats <- formicidae_sample$lat #assigning latitude and longitude values to an object
phylo_lons <- formicidae_sample$lon 

phylo_coords <- data.frame(x=phylo_lons,y=phylo_lats)
phylo_coords <- na.omit(phylo_coords) #remove NA values


plot(phylo.to.map(Phylo2, phylo_coords, rotate=FALSE))

obj<-phylo.to.map(Phylo2, phylo_coords ,plot=FALSE)



