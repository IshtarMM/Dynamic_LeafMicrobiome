######################################Create the plot SparCC
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(plyr)
library(data.table)
library(rowr)
library(grid)
library(gridExtra)
###########for one file sparcc
setwd("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/NetworkResult/ResultSparcc/Otu/Month/Pval/NetFeature/")

###################plot all together sparcc
#listfiles = c("2.csv")
listfiles = c("11.csv" , "12.csv" , "1.csv" , "2.csv" , "3.csv")
#listfiles = c("1.csv" , "2.csv" , "3.csv" , "All.csv")
plotsCCBC<- list()
for (i in 1:length(listfiles)){
  data<-read.csv(listfiles[i], sep=",", dec = ".", stringsAsFactors = F, header=T, row.names = NULL )
  data["Genus_Otu"] = paste(data$Genus,"_",data$name)
  data$kingdom = gsub("\\s*\\([^\\)]+\\)","",as.character(data$Kingdom))
  fname = sapply(basename(listfiles[i]), function(x) gsub(".csv", "", x))
  print(fname)
  n <- 5
  top5_BC<-data[data$BetweennessCentrality > quantile(data$BetweennessCentrality,prob=1-n/100),]$Genus_Otu
  
  top5_BC_min<-min(data[data$BetweennessCentrality > quantile(data$BetweennessCentrality,prob=1-n/100),]$BetweennessCentrality)
  
  top5_CC<-data[data$ClosenessCentrality > quantile(data$ClosenessCentrality,prob=1-n/100),]$Genus_Otu
  
  top5_CC_min<-min(data[data$ClosenessCentrality > quantile(data$ClosenessCentrality,prob=1-n/100),]$ClosenessCentrality)
  
  
  p5_Intersection <- top5_BC[which((top5_BC %in% top5_CC) == TRUE)]
  
  p5_Intersection1 = as.data.frame(p5_Intersection)
  
  expression <- ggplot(data, aes(x = ClosenessCentrality, y= BetweennessCentrality, size= TotalRelativeabundance, shape=Core, label=Genus_Otu, ymax=max(BetweennessCentrality)*1.1)) +
    theme(text = element_text(size=15))
  
  z <- expression + geom_point(aes(colour=factor(kingdom), shape=factor(Core)), stat= "identity",position="identity",alpha= 0.7) +
    labs(title = fname) +
    scale_size_area(max_size = 15)+
    geom_vline(xintercept = top5_CC_min, linetype="dashed") +
    geom_hline(yintercept = top5_BC_min, linetype= "dashed") +
    scale_x_continuous("Closeness Centrality") +
    scale_y_continuous("Betweenness Centrality") +
    scale_shape_manual(values=c(21, 19)) +
    geom_text_repel(data=subset(data, Core =="yes"),size = 2,
                    box.padding   = 1.5,
                    point.padding = 0.5,
                    force         = 100,
                    segment.size  = 0.2,
                    segment.color = "grey50") +
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  plotsCCBC[[i]] = z
}
resultplot = marrangeGrob(grobs=plotsCCBC, nrow=5, ncol=1)

ggsave("MonthOtu.jpg" , resultplot , dpi = 1000 , width = 12, height = 30, units = "in")
ggsave("MonthOtu.pdf" , resultplot , dpi = 1000 , width = 12, height = 30, units = "in")
dev.off()



##########################################################################################################################################

#CoNet all file plot Otu

setwd("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/NetworkResult/CoNetResults/Otu/Edges/NetFeature/Month/")

listfiles = c("11.csv" , "12.csv" , "1.csv" , "2.csv" , "3.csv")
#listfiles = c("11p.csv" , "12p.csv" , "1p.csv" , "2p.csv" , "3p.csv" , "allp.csv")
#listfiles = c("1.csv" , "2.csv" , "3.csv" , "All.csv")
#listfiles = c("3p.csv")
##########Draw box plot for all files in one page
plotsCCBC<- list()
for (i in 1:length(listfiles)){
  data<-read.csv(listfiles[i], sep=",", dec = ".", stringsAsFactors = F, header=T, row.names = NULL )
  data["Genus_Otu"] = paste(data$Genus,"_",data$name)
  data$kingdom = gsub("\\s*\\([^\\)]+\\)","",as.character(data$Kingdom))
  fname = sapply(basename(listfiles[i]), function(x) gsub(".csv", "", x))
  print(fname)
  n <- 5
  top5_BC<-data[data$BetweennessCentrality > quantile(data$BetweennessCentrality,prob=1-n/100),]$Genus_Otu
  
  top5_BC_min<-min(data[data$BetweennessCentrality > quantile(data$BetweennessCentrality,prob=1-n/100),]$BetweennessCentrality)
  
  top5_CC<-data[data$ClosenessCentrality > quantile(data$ClosenessCentrality,prob=1-n/100),]$Genus_Otu
  
  top5_CC_min<-min(data[data$ClosenessCentrality > quantile(data$ClosenessCentrality,prob=1-n/100),]$ClosenessCentrality)
  
  
  p5_Intersection <- top5_BC[which((top5_BC %in% top5_CC) == TRUE)]
  
  p5_Intersection1 = as.data.frame(p5_Intersection)
  
  expression <- ggplot(data, aes(x = ClosenessCentrality, y= BetweennessCentrality, size= TotalRelativeabundance, shape=Core, label=Genus_Otu, ymax=max(BetweennessCentrality)*1.1)) +
    theme(text = element_text(size=15))
  
  z <- expression + geom_point(aes(colour=factor(kingdom), shape=factor(Core)), stat= "identity",position="identity",alpha= 0.7) +
    labs(title = fname) +
    scale_size_area(max_size = 15)+
    geom_vline(xintercept = top5_CC_min, linetype="dashed") +
    geom_hline(yintercept = top5_BC_min, linetype= "dashed") +
    scale_x_continuous("Closeness Centrality") +
    scale_y_continuous("Betweenness Centrality") +
    scale_shape_manual(values=c(21, 19)) +
    geom_text_repel(data=subset(data, Core =="yes"),size = 2,
                    box.padding   = 1.5,
                    point.padding = 0.5,
                    force         = 100,
                    segment.size  = 0.2,
                    segment.color = "grey50") +
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  plotsCCBC[[i]] = z
}
resultplot = marrangeGrob(grobs=plotsCCBC, nrow=5, ncol=1)

ggsave("MonthOtu.jpg" , resultplot , dpi = 1000 , width = 12, height = 30, units = "in")
ggsave("MonthOtu.pdf" , resultplot , dpi = 1000 , width = 12, height = 30, units = "in")
dev.off()





##############################Sparcc Summerized node



library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(plyr)
library(data.table)
library(rowr)
library(grid)
library(gridExtra)
###########for one file sparcc Overall_avg centrality
setwd("/Users/mahmoudi/Documents/Doc/project/SparCC/SparccOutputs/1000/")

data<-read.csv("summerizednodetable.csv", sep=",", dec = ".", stringsAsFactors = F, header=T, row.names = NULL )
data$Overall_avgRank = rank(-data[,27] , ties.method = "first")
n <- 5

top5_OA<-data[data$Overall_avg < quantile(data$Overall_avg,prob=1-95/100),]$OTU_taxo
top5_OA_min<-max(data[data$Overall_avg < quantile(data$Overall_avg,prob=1-95/100),]$Overall_avg)

top5_OCC<-data[data$Ocurrence.total > quantile(data$Ocurrence.total,prob=1-n/100),]$OTU_taxo

top5_OCC_min<-min(data[data$Ocurrence.total > quantile(data$Ocurrence.total,prob=1-n/100),]$Ocurrence.total)


p5_Intersection <- top5_OA[which((top5_OA %in% top5_OCC) == TRUE)]

p5_Intersection1 = as.data.frame(p5_Intersection)

expression <- ggplot(data, aes(x = Ocurrence.total, y= Overall_avg, size= Relative.abundance, shape=Core., label=OTU_Genera, ymax=max(Overall_avg)*-0.1)) +
  theme(text = element_text(size=15))

z <- expression + geom_point(aes(colour=factor(Kingdom), shape=factor(Core.)), stat= "identity",position="identity",alpha= 0.7) +
  scale_size_area(max_size = 15)+
  geom_vline(xintercept = top5_OCC_min, linetype="dashed") +
  geom_hline(yintercept = top5_OA_min, linetype= "dashed") +
  scale_x_continuous("Ocurrence.total") +
  scale_shape_manual(values=c(21, 19)) +
  geom_text_repel(data=subset(data, Core. =="yes"),size = 2,
                  box.padding   = 1.5,
                  point.padding = 0.5,
                  force         = 100,
                  segment.size  = 0.2,
                  segment.color = "grey50") +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +scale_y_reverse()
z
ggsave("CentralitAvg_OCC.pdf" , z , dpi = 1000 , width = 15, height = 10, units = "in")
dev.off()

###################################################################################################################################

######################################Create the plot SparCC finall
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(plyr)
library(data.table)
library(rowr)
library(grid)
library(gridExtra)
###########for one file sparcc
setwd("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/NetworkResult/ResultSparcc/Otu/Month/Pval/NetFeature/")

###################plot all together sparcc and each month
#listfiles = c("All.csv")
listfiles = c("11.csv" , "12.csv" , "1.csv" , "2.csv" , "3.csv")
#listfiles = c("1.csv" , "2.csv" , "3.csv" , "All.csv")
#listfiles = c("All.csv")
plotsCCBC<- list()
for (i in 1:length(listfiles)){
  data<-read.csv(listfiles[i], sep=",", dec = ".", stringsAsFactors = F, header=T, row.names = NULL )
  rownames(data) = data$name
  #data["Otrad_Otu_00001", "Core"] <- 'yes'
  data["Genus_Otu"] = paste(data$Genus,"_",data$name)
  data$kingdom = gsub("\\s*\\([^\\)]+\\)","",as.character(data$Kingdom))
  fname = sapply(basename(listfiles[i]), function(x) gsub(".csv", "", x))
  print(fname)
  n <- 5
  top5_BC<-data[data$BetweennessCentrality > quantile(data$BetweennessCentrality,prob=1-n/100),]$Genus_Otu
  
  top5_BC_min<- min(data[data$BetweennessCentrality > quantile(data$BetweennessCentrality,prob=1-n/100),]$BetweennessCentrality)
  
  top5_CC<-data[data$ClosenessCentrality > quantile(data$ClosenessCentrality,prob=1-n/100),]$Genus_Otu
  
  top5_CC_min<-min(data[data$ClosenessCentrality > quantile(data$ClosenessCentrality,prob=1-n/100),]$ClosenessCentrality)
  
  
  p5_Intersection <- top5_BC[which((top5_BC %in% top5_CC) == TRUE)]
  
  p5_Intersection1 = as.data.frame(p5_Intersection)
  
  expression <- ggplot(data, aes(x = ClosenessCentrality, y= BetweennessCentrality,shape=Core, label=shared.name, ymax=max(BetweennessCentrality))) +
    theme(text = element_text(size=2))
  
  z <- expression + geom_point(aes(colour = factor(kingdom ),size=Degree),stat= "identity",position="identity",alpha= 0.7) +
    labs(title = fname) +
    scale_size_area(max_size = 5) +
    geom_vline(xintercept = top5_CC_min, linetype="dashed") +
    geom_hline(yintercept = top5_BC_min, linetype= "dashed") +
    scale_x_continuous("Closeness Centrality") +
    scale_y_continuous("Betweenness Centrality") +
    scale_shape_manual(values=c(21, 19))  +
    theme_bw()+ theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black",size = 0.5), 
    axis.text.x = element_text( colour="black", size = 10), 
    axis.text.y = element_text(  size = 10 , colour="black")) + 
    scale_color_manual(values=c("#009E73", "#D55E00","#56A5CC")) +
    geom_text_repel(data=subset(data,shared.name=="BacV5_Otu_000012" | shared.name=="BacV5_Otu_000011" | shared.name=="BacV5_Otu_000004"  ),size = 2,
                    box.padding   = 5,
                    point.padding = 0.1,
                    segment.size  = 0.4,
                    segment.color = "black") + ylim(0,0.05) + xlim(0,0.7)
  
  ggsave(paste("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/Dataset2019/results/PlotsPaper/Network/Sparcc/",listfiles[i],".pdf" , sep="") , z  ,  width = 15 , height =10 , units = "cm")
  }

z
#resultplot = marrangeGrob(grobs=plotsCCBC, nrow=5, ncol=1) | shared.name=="Otrad_Otu_00001"

 #ggsave("all.pdf" , z , dpi = 1000 , width = 16, height = 16, units = "cm")
ggsave("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/Dataset2019/results/PlotsPaper/AllASparCC.pdf" , z  ,  width = 15 , height =10 , units = "cm")
 
dev.off()
#ggsave("11.pdf" , z , dpi = 1000 , width = 8, height = 7, units = "in")
#dev.off()

###################################################################################################################################


######################################Create the plot CoNet fainal
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(plyr)
library(data.table)
library(rowr)
library(grid)
library(gridExtra)
###########for one file sparcc
setwd("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/NetworkResult/CoNetResults/Otu/Edges/NetFeature/Month/")

###################plot all together sparcc
#listfiles = c("All.csv")
listfiles = c("11.csv" , "12.csv" , "1.csv" , "2.csv" , "3.csv")
#listfiles = c("1.csv" , "2.csv" , "3.csv" , "All.csv")
#listfiles = c("All.csv")
plotsCCBC<- list()

for (i in 1:length(listfiles)){
  data<-read.csv(listfiles[i], sep=",", dec = ".", stringsAsFactors = F, header=T, row.names = NULL )
  rownames(data) = data$name
  data["Otrad-Otu-00001", "Core"] <- 'yes'
  data["Genus_Otu"] = paste(data$Genus,"_",data$name)
  data$kingdom = gsub("\\s*\\([^\\)]+\\)","",as.character(data$Kingdom))
  fname = sapply(basename(listfiles[i]), function(x) gsub(".csv", "", x))
  print(fname)
  n <- 5
  top5_BC<-data[data$BetweennessCentrality > quantile(data$BetweennessCentrality,prob=1-n/100),]$Genus_Otu
  
  top5_BC_min<-min(data[data$BetweennessCentrality > quantile(data$BetweennessCentrality,prob=1-n/100),]$BetweennessCentrality)
  
  top5_CC<-data[data$ClosenessCentrality > quantile(data$ClosenessCentrality,prob=1-n/100),]$Genus_Otu
  
  top5_CC_min<-min(data[data$ClosenessCentrality > quantile(data$ClosenessCentrality,prob=1-n/100),]$ClosenessCentrality)
  
  
  p5_Intersection <- top5_BC[which((top5_BC %in% top5_CC) == TRUE)]
  
  p5_Intersection1 = as.data.frame(p5_Intersection)
  
  expression <- ggplot(data, aes(x = ClosenessCentrality, y= BetweennessCentrality,shape=Core, label=shared.name, ymax=max(BetweennessCentrality))) +
    theme(text = element_text(size=2))
  
  z <- expression + geom_point(aes(colour = factor(kingdom),size=Degree),stat= "identity",position="identity",alpha= 0.7) +
    labs(title = fname) +
    scale_size_area(max_size = 5)+
    geom_vline(xintercept = top5_CC_min, linetype="dashed") +
    geom_hline(yintercept = top5_BC_min, linetype= "dashed") +
    scale_x_continuous("Closeness Centrality") +
    scale_y_continuous("Betweenness Centrality") +
    scale_shape_manual(values=c(21, 19))  +
    theme_bw()+ theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      axis.line = element_line(colour = "black",size = 0.5), 
                      axis.text.x = element_text( colour="black", size = 10), 
                      axis.text.y = element_text(  size = 10 , colour="black")) + 
    scale_color_manual(values=c("#009E73", "#D55E00","#56A5CC")) +
    geom_text_repel(data=subset(data,shared.name=="BacV5-Otu-000012" | shared.name=="BacV5-Otu-000011" | shared.name=="BacV5-Otu-000004" ),size = 3,
                    box.padding   = 2,
                    point.padding = 0.1,
                    segment.size  = 0.4,
                    segment.color = "black") + ylim(0,0.06) + xlim(0,0.7)
  
  
  ggsave(paste("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/Dataset2019/results/PlotsPaper/Network/CoNet/",listfiles[i],".pdf" , sep="") , z  ,  width = 15 , height =10 , units = "cm")
}

z
#resultplot = marrangeGrob(grobs=plotsCCBC, nrow=5, ncol=1)

#ggsave("All.pdf" , z , dpi = 1000 , width = 16, height = 10, units = "cm")
ggsave("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/Dataset2019/results/PlotsPaper/AllACoNet.pdf" , z  ,  width = 15 , height =10 , units = "cm")

dev.off()
#ggsave("11.pdf" , z , dpi = 1000 , width = 8, height = 7, units = "in")
#dev.off()

#############################################################################################
######################################Create the plot Sparcc fainal
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(plyr)
library(data.table)
library(rowr)
library(grid)
library(gridExtra)
###########for one file sparcc year
setwd("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/NetworkResult/ResultSparcc/Otu/Year/Pval/NetFeature/")

###################plot all together sparcc
#listfiles = c("All.csv")
#listfiles = c("11.csv" , "12.csv" , "1.csv" , "2.csv" , "3.csv")
listfiles = c("1.csv" , "2.csv" , "3.csv")
#listfiles = c("All.csv")
plotsCCBC<- list()
for (i in 1:length(listfiles)){
  data<-read.csv(listfiles[i], sep=",", dec = ".", stringsAsFactors = F, header=T, row.names = NULL )
  rownames(data) = data$name
  data["Otrad_Otu_00001", "Core"] <- 'yes'
  data["Genus_Otu"] = paste(data$Genus,"_",data$name)
  data$kingdom = gsub("\\s*\\([^\\)]+\\)","",as.character(data$Kingdom))
  fname = sapply(basename(listfiles[i]), function(x) gsub(".csv", "", x))
  print(fname)
  n <- 5
  top5_BC<-data[data$BetweennessCentrality > quantile(data$BetweennessCentrality,prob=1-n/100),]$Genus_Otu
  
  top5_BC_min<- min(data[data$BetweennessCentrality > quantile(data$BetweennessCentrality,prob=1-n/100),]$BetweennessCentrality)
  
  top5_CC<-data[data$ClosenessCentrality > quantile(data$ClosenessCentrality,prob=1-n/100),]$Genus_Otu
  
  top5_CC_min<-min(data[data$ClosenessCentrality > quantile(data$ClosenessCentrality,prob=1-n/100),]$ClosenessCentrality)
  
  
  p5_Intersection <- top5_BC[which((top5_BC %in% top5_CC) == TRUE)]
  
  p5_Intersection1 = as.data.frame(p5_Intersection)
  
  expression <- ggplot(data, aes(x = ClosenessCentrality, y= BetweennessCentrality,shape=Core, label=shared.name, ymax=max(BetweennessCentrality))) +
    theme(text = element_text(size=2))
  
  z <- expression + geom_point(aes(colour = factor(kingdom ),size=Degree),stat= "identity",position="identity",alpha= 0.7) +
    labs(title = fname) +
    scale_size_area(max_size = 5) +
    geom_vline(xintercept = top5_CC_min, linetype="dashed") +
    geom_hline(yintercept = top5_BC_min, linetype= "dashed") +
    scale_x_continuous("Closeness Centrality") +
    scale_y_continuous("Betweenness Centrality") +
    scale_shape_manual(values=c(21, 19))  +
    theme_bw()+ theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      axis.line = element_line(colour = "black",size = 0.5), 
                      axis.text.x = element_text( colour="black", size = 10), 
                      axis.text.y = element_text(  size = 10 , colour="black")) + 
    scale_color_manual(values=c("#009E73", "#D55E00","#56A5CC"))
  
  ggsave(paste("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/Dataset2019/results/PlotsPaper/",listfiles[i],".png" , sep="") , z  ,  width = 15 , height =10 , units = "cm")
}






