# By Maryam Mahmoudi , maryam.mahmoudi@uni-tuebingen.de
#Changes in phyllosphere microbial interaction networks

library(multcompView)
library(ggplot2)
library(tidyr)
library(wesanderson)
library(ggpubr)
library(EnvStats)
library(tidyr)
library(wesanderson)
library(ggpubr)
library(EnvStats)
library(rstatix)
library(agricolae)
library(plyr)
library(FSA)
library(rcompanion)
library(igraph)
library(poweRlaw)


# extrct data for figure 4B and 4D
filelist = c("11-12node.csv","12-1node.csv","1-2node.csv","2-3node.csv")
#Monthname = c("Novemnber", "December" ,"January" , "February" , "March")
Monthname1 = c("Novemnber-December", "December-January" ,"January-February" , "February-March")
data <- data.frame()
for (i in 1:4){
  Node <- read.csv(filelist[i])
  colnames(Node) <- c("SUID", "Firstnodefile" ,"Secondnodefile" ,"DnScoredegreeCorrected","RewiringDnScore","EdgeCount","Name","Selected","SharedName")
  Firstmonthnode <- Node[Node$Firstnodefile=="true",]
  Secondmonthnode <- Node[Node$Secondnodefile=="true",]
  #first <- Node[Node$Secondnodefile=="true",]
  overlap <- Node[(Node$Secondnodefile=="true" & Node$Firstnodefile == "true"),]
  #second <- Node[Node$Firstnodefile == "true",]
  #print(paste(basename(Node) , nrow(first)-nrow(overlap) , nrow(overlap) , nrow(second)-nrow(overlap)))
  
  Node$Month <- Monthname1[i]
  print(Monthname1[i])
  print(basename(filelist[i]))
  #print(paste(nrow(Firstmonthnode) - nrow(overlap) ,nrow(overlap) ,nrow(Secondmonthnode) - nrow(overlap)))
  print(paste(nrow(Firstmonthnode)  ,nrow(overlap) ,nrow(Secondmonthnode)))
  print(paste("overlap/Secondmonthnode= " ,  nrow(overlap)/nrow(Secondmonthnode)))
}

#put printed number in following table , and use the table to plot Fig4Band FD

data <- data.frame(inheritededges=c(NA,0.06,0.16,0.34,0.084),
                   edges = c(8814,12748,9879,2342,5919),
                   nodes=c(355,495,314,236,327),
                   inheritednodes=c(NA,0.51,0.89,0.83,0.58),
                   Month=c("November","December","January","February","March"))

data$Month <- factor(data$Month, levels = c("November", "December" ,"January" , "February","March"))

############
#56A5CC
#D55E00
#0.89/0.34 = 2.61

# Fig 4D
#inheritededges and inheritednodes inone plot
scale = function(x) sprintf("%.1f",x)
plot1 <- ggplot(data, aes(x = Month , group=1) )  +  
  geom_point(aes(y = inheritededges ),size=3) +
  geom_line(aes(y = inheritededges ) , linetype = "twodash" , size=1.5)  +
  geom_point(aes(y = inheritednodes/2.61 ),size=3) +
  geom_line(aes(y =inheritednodes/2.61 ),size=1.5) +
  scale_y_continuous(sec.axis = sec_axis(~.*2.61, name = "Inherited nodes %" , labels = scale)) + 
  scale_color_manual(values=c("#0072B2","#D55E00")) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black",size = 12), 
        axis.text.y = element_text(angle=90, hjust=1 ,size = 12 , colour="black" ) , 
        axis.title.y =  element_text(angle=90,size = 12) ,
        axis.title.x = element_blank() ,
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(y=" Inherited edges %" )

  ggsave("Fig4_D_InheritedNodeEdges.pdf" , plot1, width = 5 , height =5 )


#-------------------------------------------------------
#Numberofnodesedges
#56A5CC
#D55E00
#12748/495 = 25.75

#  Figure "4b" Number of nodes and edges per month
plot2 <- ggplot(data, aes(x = Month , group=1) )  +  
  geom_line(aes(y = edges ),size=1.5 ,linetype = "twodash") + 
  geom_point(aes(y = edges ) ,size=3) +
  geom_line(aes(y =nodes*25.75),size=1.5) +
  geom_point(aes(y = nodes*25.75) ,size=3) + 
  scale_y_continuous(sec.axis = sec_axis(~./25.8, name = "Number of nodes")) + 
  scale_color_manual(values=c("#0072B2","#D55E00")) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black",size = 12), 
        axis.text.y = element_text(angle=90, hjust=1 ,size = 12 , colour="black" ) , 
        axis.title.y =  element_text(angle=90,size = 12) ,
        axis.title.x = element_blank() ,
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(y="Number of edges",colour = "Parameter" ) 

  ggsave("Fig4_B_InheritedNodeEdges.pdf" , plot2, width = 5 , height =5 )





