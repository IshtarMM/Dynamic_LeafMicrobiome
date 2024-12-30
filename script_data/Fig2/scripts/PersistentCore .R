# By Maryam Mahmoudi , maryam.mahmoudi@uni-tuebingen.de

# Taxonomy based phylogenetic tree of persistent core microbe 
library(cluster)
library(ape)
require(tidyverse)

packages <- c("readxl" , "vegan" , "tidyr" , "ggplot2" , "rowr" , "dplyr" , "reshape" , "BiodiversityR" , "phyloseq" , "cowplot" , "otuSummary" , "RColorBrewer" , "forcats" , "funrar") # list of packages that are needed 
lapply(packages, require, character.only = TRUE)


table = read.table("data/OccurenceData.txt" , sep = "\t" , header = T , row.names = 1)
table$otu = row.names(table)
table1 <- table %>% separate("otu" , c("Genera", "otu1" ,"NumOtu") , "_" , remove = FALSE)
table1 = table1 %>% unite(otunum, c("otu1" , "NumOtu") , sep = "" , remove = FALSE)
table1 = table1 %>% unite(otu_genera, c("otunum" , "Genus") , sep = "_" , remove = FALSE)
rownames(table1) <- table1$otu_genera
core = subset(table1,table1$Core=="yes")
core = core[,c(13:17,19:20)]
colnames(core) = c("K" , "P" , "C" , "O" , "F" , "G" , "S")
cols <- c("K" , "P" , "C" , "O" , "F" , "G" , "S")
core[cols] <- lapply(core[cols], factor)  ## as.factor() could also be used

dist <- daisy(core, metric = "gower")
dendo <- hclust(dist )
pdf("fig2a_dendo_treeCores.pdf" , width = 14)
plot(dendo )
dev.off()
dendo_tree <- as.phylo(dendo) 
write.tree(phy=dendo_tree, file="dendo_treeCores.newick") 


# Changes in the relative abundance of core taxa over time

coreRA = subset(table , Core == "yes")
coreRA = coreRA[c(1,20)]
data = c("data/BV5.txt" , "data/FITS2.txt" , "data/PV9.txt")
listP = list()
tab = data.frame()
for (i in 1:3){
  if (i < 3){
    otu = read.table(data[i] , header = TRUE)
    row.names(otu) = otu$Row.names
    otu1 = otu[c(17:ncol(otu))]
    otu2 = merge(y = otu1, x = otu[,c(11,12)] , by = "row.names")
    row.names(otu2) = otu2$Row.names
    otu2 = otu2[-c(1,2)]
    agg =  aggregate(. ~  Month, data = otu2, sum) 
    row.names(agg) = agg$Month
    aggRA = make_relative(as.matrix(agg[-c(1)]))          #abundance matrix, with sites in rows and species in columns
    aggmelt = melt(data = aggRA)
    aggmelt = aggmelt[(aggmelt$X2 %in% coreRA$otu),]
    colnames(aggmelt) = c("Month" , "OTU" , "TotalRelativeAbundance")
    aggmelt$Month = factor(aggmelt$Month, levels = c("Nov","Dec" ,"Jan","Feb","Mar"))
    tab = rbind(tab,aggmelt)
    p = ggplot(aggmelt, aes(Month,y=TotalRelativeAbundance ,group=OTU, fill = OTU , position="stack")) + 
      geom_area( colour="black", size=.2, alpha=.8) + 
      theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 8), 
            axis.text.y = element_text(angle=90, hjust=1 , size = 8 , colour="black" ) , 
            axis.title.y =  element_text(angle=90, size = 8) ,
            axis.title.x = element_blank() ,
            panel.background = element_rect(colour = "black" , fill = "NA"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylim(0,0.75)
    listP[[i]] = p
  }
  else{
    
    otu = read.table(data[i] , header = TRUE)
    row.names(otu) = otu$Row.names
    otu1 = otu[c(17:ncol(otu))]
    otu2 = merge(y = otu1, x = otu[,c(11,12)] , by = "row.names")
    row.names(otu2) = otu2$Row.names
    otu2 = otu2[-c(1,2)]
    agg =  aggregate(. ~  Month, data = otu2, sum) 
    row.names(agg) = agg$Month
    aggRA = make_relative(as.matrix(agg[-c(1)]))          #abundance matrix, with sites in rows and species in columns
    aggmelt = melt(data = aggRA)
    aggmelt = aggmelt[(aggmelt$X2 %in% coreRA$otu),]
    colnames(aggmelt) = c("Month" , "OTU" , "TotalRelativeAbundance")
    aggmelt$Month = factor(aggmelt$Month, levels = c("Nov","Dec" ,"Jan","Feb","Mar"))
    tab = rbind(tab , aggmelt)
    p = ggplot(aggmelt, aes(Month,y=TotalRelativeAbundance ,group=OTU, fill = OTU , position="stack")) + 
      geom_area( colour="black", size=.2, alpha=.8) + 
      theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 8), 
            axis.text.y = element_text(angle=90, hjust=1 , size = 8 , colour="black" ) , 
            axis.title.y =  element_text(angle=90, size = 8) ,
            axis.title.x = element_blank() ,
            panel.background = element_rect(colour = "black" , fill = "NA"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylim(0,0.75)
    listP[[i]] = p
  }
}
f = do.call("grid.arrange", c(listP, ncol=3 , nrow = 1))
ggsave("Fig2b_RelativeAbundanceCores.pdf" , f, width = 12 , height = 3 , units = "in")

