packages <- c("readxl" , "vegan" , "tidyr" , "ggplot2" , "rowr" , "dplyr" ,
              "reshape" , "BiodiversityR" , "phyloseq" , "cowplot" , 
              "otuSummary" , "RColorBrewer" , "forcats" ,"ggforce" , "plyr" , "microbiome") # list of packages that are needed 
library(agricolae)
library(gridExtra) # for grid.arrange
library(grid) 

lapply(packages, require, character.only = TRUE)
setwd("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/Dataset2019/")

listAfiles = c("BV5.txt","FITS2.txt" , "PV9.txt")
pltList <- list()
color = c("#009E73","#F70C0C","#0072B2")
for (i in 1:3){
  file = read.table(listAfiles[i], header = T)
  file$Month = factor(file$Month, levels = c("Nov","Dec" ,"Jan","Feb","Mar"))
  row.names(file) = file$Replicate
  sample_otu <- file[-c(1:16)]  
  otu_mat <- t(sample_otu) # otu matrix as a input for phyloseq (rownames otu and sample in the columns)
  sampleinfo <- file[1:16]  #samples information
  sampleinfo$Experiment <- sub("^", "Year", sampleinfo$Experiment)
  rownames(sampleinfo) <- sampleinfo$Replicate
  
  otu <- otu_table(otu_mat , taxa_are_rows = TRUE,errorIfNULL = TRUE)
  samples <- sample_data(sampleinfo)
  filePh <- phyloseq(otu,samples)
  ############################################################################################################################################################
  rich = estimate_richness(filePh,measures = c("Observed", "Chao1", "Shannon"))
  rich = select(rich , -c(3))
  rich1 = merge(rich,data.frame(sample_data(filePh)),by=0)
  rich1 = select(rich1 , -c(1))
  print(shapiro.test(rich$Shannon))
  #################################################
  p1 <- ggplot(rich1 , aes(x = Month , y = Shannon))  + 
    geom_boxplot( fill=color[i], color="black"  , outlier.shape=NA ,outlier.size=-1,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                  position = position_dodge(width=6))  + 
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill=color[i] ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
    geom_jitter(alpha = 1 , color = "black" , position = position_jitter(width = 0.3) , size = 0.2) + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 10), 
          axis.text.y = element_text(angle=90, hjust=1 , size = 10 , colour="black" ) , 
          axis.title.y =  element_text(angle=90, size = 10) ,
          axis.title.x = element_blank() ,
          panel.background = element_rect(colour = "black" , fill = "NA"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylim(0,6)
  pltList[[i]] = p1
}

#grid.arrange(pltList[[1]], pltList[[2]], pltList[[3]], ncol = 1) # say you have 4 plots
x = 3
listAfilesRA = c("BV5RAlog.txt","FITS2RAlog.txt" , "PV9RAlog.txt")
color = c("#009E73","#F70C0C","#0072B2")
for (j in 1:3){
  x = x + 1
  file = read.table(listAfilesRA[j], header = T)
  file$Month = factor(file$Month, levels = c("Nov","Dec" ,"Jan","Feb","Mar"))
  row.names(file) = file$Replicate
  sample_otu <- file[-c(1:17)]  
  sampleinfo <- file[1:17]  #samples information
  sampleinfo$Experiment <- sub("^", "Year", sampleinfo$Experiment)
  rownames(sampleinfo) <- sampleinfo$Replicate
  info_distance = betadisper(vegdist(sample_otu, method="bray"), sampleinfo$Month, type = c("centroid"), bias.adjust = TRUE)
  dis_to_center = data.frame(info_distance$distances)
  colnames(dis_to_center) = "distance"
  dis_to_center1  = merge(dis_to_center , sampleinfo  , by = "row.names")
  print(shapiro.test(dis_to_center1$distance))
  
  #################################################
  scale = function(x) sprintf("%.1f",x)
  p2 <- ggplot(dis_to_center1 , aes(x = Month , y = distance))  + 
    geom_boxplot( fill=color[j], color="black"  , outlier.shape=NA ,outlier.size=-1,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                  position = position_dodge(width=6))  + 
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill=color[j] ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
    geom_jitter(alpha = 1 , color = "black" , position = position_jitter(width = 0.3) , size = 0.2) + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 10), 
          axis.text.y = element_text(angle=90, hjust=1 , size = 10 , colour="black" ) , 
          axis.title.y =  element_text(angle=90, size = 10) ,
          axis.title.x = element_blank() ,
          panel.background = element_rect(colour = "black" , fill = "NA"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_y_continuous(labels =scale,limits = c(0,1) ) 
  pltList[[x]] = p2
}

#grid.arrange(pltList[[1]], pltList[[2]], pltList[[3]], ncol = 1)
x = grid.arrange(grobs = lapply(pltList[1:3], "+", theme(plot.margin=margin(10,10,10,10))))
y = grid.arrange(grobs = lapply(pltList[4:6], "+", theme(plot.margin=margin(10,10,10,10))))

f = grid.arrange(x,y , ncol = 2)
f = grid.arrange(x,y,z , ncol = 3) #z come form Bc_time_compa_110221.R
ggsave("results/PlotsPaper/div_dist4.pdf" , f, width = 18 , height = 18 , units = "cm")



####NMDS
pltListnmds = list()
listAfilesRA = c("BV5RAlog.txt","FITS2RAlog.txt" , "PV9RAlog.txt")
factor = c("Month","Experiment","Ecotype")
color = c("#bb2200","#00bb3b","#009cbb"  , "#5e00bb" , "#bbaf00" , "#0072B2", "#D55E00", "#CC79A7")
scale = function(x) sprintf("%.1f",x)
x = 0
for (i in 1:3){
  for (j in factor){
        x = x + 1
        file = read.table(listAfilesRA[3], header = T)
        file$Month = factor(file$Month, levels = c("Nov","Dec" ,"Jan","Feb","Mar"))
        row.names(file) = file$Replicate
        sample_otu <- file[-c(1:17)]  
        sampleinfo <- file[1:17]  #samples information
        #sampleinfo$Experiment <- sub("^", "Year", sampleinfo$Experiment)
        rownames(sampleinfo) <- sampleinfo$Replicate
        sample_otu <- file[-c(1:17)]  
        otu_mat <- t(sample_otu) # otu matrix as a input for phyloseq (rownames otu and sample in the columns)
        sampleinfo <- file[1:17]  #samples information
        #sampleinfo$Experiment <- sub("^", "Year", sampleinfo$Experiment)
        rownames(sampleinfo) <- sampleinfo$Replicate
        otu <- otu_table(otu_mat , taxa_are_rows = TRUE,errorIfNULL = TRUE)
        samples <- sample_data(sampleinfo)
        filePh <- phyloseq(otu,samples)
        file.ord <- ordinate(filePh, "NMDS", "bray")
        ##############
        nmds = plot_ordination(filePh, file.ord, type="Samples", color=j) +
          theme_bw() + 
          theme(legend.position = "none" ,legend.text = element_text(colour="black", size = 12),
                legend.title = element_text(face = NULL),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                text = element_text(size = 12) ,  
                panel.border = element_rect(colour = "black", fill=NA, size=0.7) ,
                axis.title.x = element_text(colour = "black" , size = 12) , 
                axis.text.x = element_text(colour = "black" , size = 12) , 
                axis.text.y = element_text(colour = "black" , size = 12) , 
                axis.title.y = element_text(colour = "black" ,size = 12),
                strip.text = element_text(colour = "black" ,size = 12)) + geom_point( size = 0.3) +
                 stat_ellipse(aes(fill=j , size = 3),geom = "polygon", type="t", alpha=0 , size = 0.6) +
          scale_y_continuous(labels =scale) + scale_color_manual(values = color) 
          pltListnmds[[x]] = nmds
        
        }
          
  }
  

x = grid.arrange(grobs = lapply(pltListnmds[1:3], "+", theme(plot.margin=margin(10,10,10,10))))
y = grid.arrange(grobs = lapply(pltListnmds[4:6], "+", theme(plot.margin=margin(10,10,10,10))))
z = grid.arrange(grobs = lapply(pltListnmds[7:9], "+", theme(plot.margin=margin(10,10,10,10))))
f = grid.arrange(x,y,z , ncol = 3)
ggsave("results/PlotsPaper/nmdslegendno.pdf" , f, width = 25 , height = 20 , units = "cm")

#############################################################################################################

####NMDS
BV5 = read.table("BV5RAlog.txt", header = T)
BV51 = subset(BV5,BV5$Experiment=="Year1")
BV52 = subset(BV5,BV5$Experiment=="Year2")
BV53 = subset(BV5,BV5$Experiment=="Year3")

FI = read.table("FITS2RAlog.txt", header = T)
FI1 = subset(FI,FI$Experiment=="Year1")
FI2 = subset(FI,FI$Experiment=="Year2")
FI3 = subset(FI,FI$Experiment=="Year3")

PV = read.table("PV9RAlog.txt", header = T)
PV1 = subset(PV,PV$Experiment=="Year1")
PV2 = subset(PV,PV$Experiment=="Year2")
PV3 = subset(PV,PV$Experiment=="Year3")







pltListnmds = list()
listAfilesRA = c("BV51","BV52" , "BV53", "FI1","FI2","FI3" ,"PV1","PV2","PV3")
factor = c("Month")
color = c("#bb2200","#00bb3b","#009cbb"  , "#5e00bb" , "#bbaf00" , "#0072B2", "#D55E00", "#CC79A7")
scale = function(x) sprintf("%.1f",x)
x = 0
for (i in 1:9){
  for (j in factor){
    x = x + 1
    file = get(listAfilesRA[i])
    file$Month = factor(file$Month, levels = c("Nov","Dec" ,"Jan","Feb","Mar"))
    row.names(file) = file$Replicate
    sample_otu <- file[-c(1:17)]  
    sampleinfo <- file[1:17]  #samples information
    #sampleinfo$Experiment <- sub("^", "Year", sampleinfo$Experiment)
    rownames(sampleinfo) <- sampleinfo$Replicate
    sample_otu <- file[-c(1:17)]  
    otu_mat <- t(sample_otu) # otu matrix as a input for phyloseq (rownames otu and sample in the columns)
    sampleinfo <- file[1:17]  #samples information
    #sampleinfo$Experiment <- sub("^", "Year", sampleinfo$Experiment)
    rownames(sampleinfo) <- sampleinfo$Replicate
    otu <- otu_table(otu_mat , taxa_are_rows = TRUE,errorIfNULL = TRUE)
    samples <- sample_data(sampleinfo)
    filePh <- phyloseq(otu,samples)
    file.ord <- ordinate(filePh, "CCA", "bray")
    ##############
    nmds = plot_ordination(filePh, file.ord, type="Samples", color=j) +
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
      stat_ellipse(aes(fill=j , size = 3),geom = "polygon", type="t", alpha=0 , size = 0.6) +
      scale_y_continuous(labels =scale) + scale_color_manual(values = color) 
    pltListnmds[[x]] = nmds
    
  }
  
}


x = grid.arrange(grobs = lapply(pltListnmds[1:3], "+", theme(plot.margin=margin(10,10,10,10))))
y = grid.arrange(grobs = lapply(pltListnmds[4:6], "+", theme(plot.margin=margin(10,10,10,10))))
z = grid.arrange(grobs = lapply(pltListnmds[7:9], "+", theme(plot.margin=margin(10,10,10,10))))
f = grid.arrange(x,y,z , ncol = 3)
ggsave("results/PlotsPaper/1CCAMonthYear.pdf" , f, width = 35 , height = 20 , units = "cm")




