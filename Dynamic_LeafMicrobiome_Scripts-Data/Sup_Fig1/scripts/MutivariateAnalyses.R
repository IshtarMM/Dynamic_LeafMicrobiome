# By Maryam Mahmoudi , maryam.mahmoudi@uni-tuebingen.de

# NMDS and permanova


# Load libreries
packages <- c("readxl" , "vegan" , "tidyr" , "ggplot2" , "rowr" , "dplyr" ,
              "reshape" , "BiodiversityR" , "phyloseq" , "cowplot" , 
              "otuSummary" , "RColorBrewer" , "forcats" ,"ggforce" , "plyr"  , "ggpubr" , "gridExtra") # list of packages that are needed 
lapply(packages, require, character.only = TRUE)


library(ggplot2)             
library(agricolae)
library(gridExtra) # for grid.arrange
library(grid) 

# NMDS plot
pltListnmds = list()
listAfilesRA = c("data/BV5RAlog10.txt","data/FITS2RAlog10.txt" , "data/PV9RAlog10.txt")
factor = c("Month","Experiment","Ecotype")
color = c("#bb2200","#00bb3b","#009cbb"  , "#5e00bb" , "#bbaf00" , "#0072B2", "#D55E00", "#CC79A7")
scale = function(x) sprintf("%.1f",x)
x = 0
for (i in 1:3){
  for (j in factor){ 
        x = x + 1
        file = read.table(listAfilesRA[i], header = T)
        file$Month = factor(file$Month, levels = c("Nov","Dec" ,"Jan","Feb","Mar"))
        row.names(file) = file$Replicate
        sample_otu <- file[-c(1:17)]  
        sampleinfo <- file[1:17]  #samples information
        rownames(sampleinfo) <- sampleinfo$Replicate
        sample_otu <- file[-c(1:17)]  
        otu_mat <- t(sample_otu) # otu matrix as a input for phyloseq (rownames otu and sample in the columns)
        sampleinfo <- file[1:17]  #samples information
        rownames(sampleinfo) <- sampleinfo$Replicate
        otu <- otu_table(otu_mat , taxa_are_rows = TRUE,errorIfNULL = TRUE)
        samples <- sample_data(sampleinfo)
        filePh <- phyloseq(otu,samples)
        file.ord <- ordinate(filePh, "NMDS", "bray")
        ##############
        nmds = plot_ordination(filePh, file.ord, type="Samples" , color = j) +
          theme_bw() + 
          theme(legend.text = element_text(colour="black", size = 12),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                text = element_text(size = 12) ,  
                panel.border = element_rect(colour = "black", fill=NA, size=0.7) ,
                axis.title.x = element_text(colour = "black" , size = 12) , 
                axis.text.x = element_text(colour = "black" , size = 12) , 
                axis.text.y = element_text(colour = "black" , size = 12) , 
                axis.title.y = element_text(colour = "black" ,size = 12),
                strip.text = element_text(colour = "black" ,size = 12)) + geom_point( size = 0.3) +
          stat_ellipse(aes( size = 3),geom = "polygon", type="t", alpha=0 , size = 0.6) +
          scale_y_continuous(labels =scale) + scale_color_manual(values = color) 
        pltListnmds[[x]] = nmds
    
  }
  
}


x = grid.arrange(grobs = lapply(pltListnmds[1:3], "+", theme(plot.margin=margin(10,10,10,10))))
y = grid.arrange(grobs = lapply(pltListnmds[4:6], "+", theme(plot.margin=margin(10,10,10,10))))
z = grid.arrange(grobs = lapply(pltListnmds[7:9], "+", theme(plot.margin=margin(10,10,10,10))))
f = grid.arrange(x,y,z , ncol = 3)
ggsave("sup_fig2a.pdf" , f, width = 14, height = 10, units = "in")


# PERMANOVA on relative abundance log10 transformation of data

# Bacteria
permanovaresults <- data.frame()
Bac <- read.table("data/BV5RAlog10.txt" , header = T)
rownames(Bac) <- Bac$Replicate
sampleinfo <- Bac[1:17]
Bac <- Bac[-c(1:17)]

ado <- adonis(Bac ~ sampleinfo$Experiment*sampleinfo$Month*sampleinfo$Ecotype, permutations = 10000, method = 'bray') 
permanovaresults <- rbind(permanovaresults , data.frame(ado$aov.tab))
permanovaresults$Taxa <- "Bacteria"

# Fungi
Fungi <- read.table("data/FITS2RAlog10.txt" , header = T)
rownames(Fungi) <- Fungi$Replicate
sampleinfo <- Fungi[1:17]
Fungi <- Fungi[-c(1:17)]

ado <- adonis(Fungi ~ sampleinfo$Experiment*sampleinfo$Month*sampleinfo$Ecotype, permutations = 10000, method = 'bray') 
adores <- data.frame(ado$aov.tab)
adores$Taxa <- "Fungi"
permanovaresults <- rbind(permanovaresults , adores)

# oomycete
Oomycete <- read.table("data/PV9RAlog10.txt" , header = T)
rownames(Oomycete) <- Oomycete$Replicate
sampleinfo <- Oomycete[1:17]
Oomycete <- Oomycete[-c(1:17)]

ado <- adonis(Oomycete ~ sampleinfo$Experiment*sampleinfo$Month*sampleinfo$Ecotype, permutations = 10000, method = 'bray') 
adores <- data.frame(ado$aov.tab)
adores$Taxa <- "Oomycete"
permanovaresults <- rbind(permanovaresults , adores)

write.table(permanovaresults , "sup_Fig2b_Permanova.txt" ,row.names = T , col.names = NA , sep = "\t")
