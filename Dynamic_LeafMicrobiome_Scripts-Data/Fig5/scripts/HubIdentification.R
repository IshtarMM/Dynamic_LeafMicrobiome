# By Maryam Mahmoudi , maryam.mahmoudi@uni-tuebingen.de

library(ggplot2)
library(ggrepel)
library(gridExtra)
library(plyr)
library(data.table)
library(rowr)
library(grid)
library(ggpubr)
library(grid)

setwd("data/")

#Node tables of each month 11=Nov,12=Dec,1=Jan,2=Feb,3=Mar
listfiles = c("11.csv" , "12.csv" , "1.csv" , "2.csv" , "3.csv") 

# Hub identification 
plotsCCBC<- list()
for (i in 1:length(listfiles)){
  data<-read.csv(listfiles[i], sep=",", dec = ".", stringsAsFactors = F, header=T, row.names = NULL )
  rownames(data) = data$name
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
    geom_text_repel(data=subset(data,shared.name=="BacV5_Otu_000012" |  shared.name=="BacV5_Otu_000004"  ),size = 2,
                    box.padding   = 5,
                    point.padding = 0.1,
                    segment.size  = 0.4,
                    segment.color = "black") 
  plotsCCBC[[i]] = z
  
}

# Changes in the connectivity of core taxa

file = read.table("Hubs_Sparcc.txt" , sep = "\t" , header = T)
file$HubCC_BCC = file$ClosenessCentrality * file$BetweennessCentrality ###
file$Genus_Otu = gsub("\\s*\\([^\\)]+\\)","",as.character(file$Genus_Otu))
file$Genus_Otu = gsub("-","_",as.character(file$Genus_Otu))
file$HubCore = ifelse(file$Hub=="yes" , "Yes", "No")


all = file[c(file$Month !="Allmonths"),]
all$Month  = factor(all$Month ,levels = c("Nov" , "Dec" , "Jan" , "Feb" , "Mar") )
all$name1 = all$Genus_Otu
all$name1 = gsub("BacV5_","",as.character(all$name1))
all$name1 = gsub("Ftrad_","",as.character(all$name1))
all$name1 = gsub("Otrad_","",as.character(all$name1))
all$name1 = gsub("Otu_","Otu",as.character(all$name1))
all$name1 = gsub("_unclassified","_un.",as.character(all$name1))
all$name1 = gsub(" _ "," ",as.character(all$name1))

library(RColorBrewer)
n <- 19
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
col1=sample(col_vector, n)
col1=sample(col, n)
f = col1
pie(rep(1,n), col=col1)
all2 = ggplot(all, aes(x = Month,y=HubCC_BCC ,group=name1 , position="stack" , colour=name1 , label = name1 )) + 
  geom_line(size=1, alpha=.8) +  geom_point(aes(shape=factor(HubCore) , size = 4)) +
 scale_shape_manual(values=c(21, 19)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 12), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 12 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 12) ,
        axis.title.x = element_blank() ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        strip.background = element_rect(colour="black",fill="white")) + 
  geom_text_repel(data = subset( all , HubCore == "Yes")  , size = 3  , force = 100 , max.overlaps = Inf,show.legend = FALSE , vjust = -0.3 , hjust = -0.3) +
  scale_color_manual(values = col1) + ylab("Hub(BC.CC)")  + ylim(0,0.03) 

plotsCCBC[[6]] = all2

resultplot = ggarrange(plotlist = plotsCCBC[c(1,2,3,4,5)],nrow=1 , ncol=5 )
ggsave("../Fig5a_HubIdentification.pdf" , resultplot  , width = 25, height = 4, units = "in")
resultplotb = ggarrange(plotlist = plotsCCBC[6],nrow=1 )
ggsave("../Fig5b_HubIdentification.pdf" , resultplotb  , width = 10, height = 7, units = "in")

