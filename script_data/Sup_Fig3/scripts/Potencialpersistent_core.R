# By Maryam Mahmoudi , maryam.mahmoudi@uni-tuebingen.de

# Taxonomy based phylogenetic tree of persistent core microbe over 3 years 

library(cluster)
library(ape)
require(tidyverse)

packages <- c("readxl" , "vegan" , "tidyr" , "ggplot2" , "rowr" , "dplyr" , "reshape" , "BiodiversityR" , "phyloseq" , "cowplot" , "otuSummary" , "RColorBrewer" , "forcats") # list of packages that are needed 
lapply(packages, require, character.only = TRUE)


table = read.table("data/OccurenceData.txt" , sep = "\t" , header = T , row.names = 1)

#make tree for otus that have more than 95 percent occurrence for Bac and PV9 and more than 98 percent for Fungi at least in one year
table2 = table
table2$otu = row.names(table2)
otu_taxa <- table2 %>% separate("otu" , c("Genera", "otu1" ,"NumOtu") , "_" , remove = FALSE)
otu_taxa = otu_taxa %>% unite(otunum, c("otu1" , "NumOtu") , sep = "" , remove = FALSE)
otu_taxa = otu_taxa %>% unite(otu_genera, c("otunum" , "Genus") , sep = "_" , remove = FALSE)
otuBac = otu_taxa[otu_taxa$Kingdom == "k__Bacteria" ,]
#Identified potential core microbe
otu_Bac1 = otuBac[otuBac$Ocurrence.y1 >=0.98 | otuBac$Ocurrence.y2 >=0.98 |otuBac$Ocurrence.y3 >=0.98, ]
otuFungi = otu_taxa[otu_taxa$Kingdom == "k__Fungi",]
otuFungi1 = otuFungi[otuFungi$Ocurrence.y1 >=0.95 | otuFungi$Ocurrence.y2 >=0.95 |otuFungi$Ocurrence.y3 >=0.95, ]
otuoomy = otu_taxa[otu_taxa$Kingdom == "k__Stramenopila",]
otuoomy1 = otuoomy[otuoomy$Ocurrence.y1 >=0.95 | otuoomy$Ocurrence.y2 >=0.95 | otuoomy$Ocurrence.y3 >=0.95, ]
otu2 = rbind(otuFungi1,otu_Bac1,otuoomy1)
rownames(otu2) <- otu2$otu_genera
taxafinal = otu2[,c(13:17, 19:20)]
colnames(taxafinal) = c("K" , "P" , "C" , "O" , "F" , "G" , "S")
cols <- c("K" , "P" , "C" , "O" , "F" , "G" , "S")
taxafinal[cols] <- lapply(taxafinal[cols], factor)  ## as.factor() could also be used

dist <- daisy(taxafinal, metric = "gower")
dendo <- hclust(dist )
pdf("sup_fig3_dendo_treeOccurenceOver3years.pdf" , width = 14)
plot(dendo )
dev.off()
dendo_tree <- as.phylo(dendo) 
write.tree(phy=dendo_tree, file="dendo_treeOccurenceOver3years.newick") 
