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

pltList = list()

# Figure "e" rewiring score of nodes box plots 
#
setwd("data/")
#file names shows month numbers 11=Nov,12=Dec,1=Jan,2=Feb,3=Mar
filelist = c("11-12node.csv","12-1node.csv","1-2node.csv","2-3node.csv")
Monthname1 = c("November_December", "December_January" ,"January_February" , "February_March")
data <- data.frame()
for (i in 1:4){
  Node <- read.csv(filelist[i])
  colnames(Node) <- c("SUID", "Firstnodefile" ,"Secondnodefile" ,"DnScoredegreeCorrected","RewiringDnScore","EdgeCount","Name","Selected","SharedName")
  Node$Month <- Monthname1[i]
  print(Monthname1[i])
  print(basename(filelist[i]))
  data <- rbind(data,Node)
}



########################################################################################################################################
# Rewiring  dunn test withe letter 
data$MonthLev <- factor(data$Month, levels = c("November_December","December_January","January_February","February_March"))


DT = dunnTest(RewiringDnScore ~ MonthLev,data=data,method="bh") 
PT = DT$res
PT
rownames(PT) = PT$Comparison
PvalAdj = PT[,4]
names(PvalAdj) = c("December_January-February_March",
                   "December_January-January_February",
                   "February_March-January_February",
                   "December_January-November_December",
                   "February_March-November_December",
                   "January_February-November_December")

multcompLetters2(RewiringDnScore ~ Month, PvalAdj, data , threshold = 0.05)

group_data <- cldList(P.adj ~ Comparison,
                      data = PT,
                      threshold = 0.05 , reverse=TRUE)
#find outliers


maxdata = ddply(data,~Month,summarise,max=min(max(RewiringDnScore),quantile(RewiringDnScore, probs = 0.75) + (IQR(RewiringDnScore) * 1.5)))

group_data1 <-merge(x = maxdata, y = group_data[,1:2] , by.x= "Month" , by.y = "Group") 

plotdunnletterRewairing <- ggplot(data = data, aes(x = MonthLev , y =RewiringDnScore) )  + 
  geom_boxplot( fill='#C0C0C0', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#C0C0C0' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) +
  geom_jitter(alpha = 1 , color = "black" , position = position_jitter(width = 0.3) , size = 0.2) + 
  theme(axis.text.x = element_text(angle = 63, hjust = 1 , colour="black",size = 12), 
        axis.text.y = element_text(angle=90, hjust=1 ,size = 12 , colour="black" ) , 
        axis.title.y =  element_text(angle=90,size = 12) ,
        axis.title.x = element_blank() ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(data=group_data1,aes(x=Month,y=max,label=Letter) 
            ,position = position_dodge(width = 0.9), size = 7 , vjust=-0.5 ) + ylim(c(0,125))
pltList[[1]] = plotdunnletterRewairing


# Figure "d" Inherited edges and inherited nodes over time

# this info comes from top codes
data <- data.frame(inheritededges=c(NA,0.06,0.16,0.34,0.084),
                   edges = c(8814,12748,9879,2342,5919),
                   nodes=c(355,495,314,236,327),
                   inheritednodes=c(NA,0.51,0.89,0.83,0.58),
                   Month=c("November","December","January","February","March"))

data$Month <- factor(data$Month, levels = c("November", "December" ,"January" , "February","March"))

###############
#56A5CC
#D55E00
#0.89/0.34 = 2.61

#inheritededges and inheritednodes inone plot
scale = function(x) sprintf("%.1f",x)
num <- ggplot(data, aes(x = Month , group=1) )  +  
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

pltList[[2]] = num
########################################################################################################################################################################
#Numberofnodesedges
#56A5CC
#D55E00
#12748/495 = 25.75

#  Figure "b" Nombure of nodes and edges per month
num1 <- ggplot(data, aes(x = Month , group=1) )  +  
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


pltList[[3]] = num1


####Figure "c" degree of nodes per month

#plots for dynet features rewiring and ...
filelist = c("11.csv","12.csv","1.csv","2.csv","3.csv")
Monthname = c("November", "December" ,"January" , "February" , "March")
data <- data.frame()
for (i in 1:5){
  Node <- read.csv(filelist[i])
  Node$Month <- Monthname[i]
  print(Monthname[i])
  print(nrow(Node))
  print(basename(filelist[i])) 
  data <- rbind(data,Node)
}

########################################################################################################################################
## finall tukey test degree letter
data$MonthLev <- factor(data$Month, levels = c("November","December","January","February","March"))
DT = dunnTest(Degree ~ MonthLev,data=data,method="bh") 
PT = DT$res
PT
PT = DT$res
PT
rownames(PT) = PT$Comparison
PvalAdj = PT[,4]
names(PvalAdj) = c("December-February","December-January","February-January"
                   ,"December-March","February-March","January-March"
                   ,"December-November","February-November","January-November","March-November")

multcompLetters2(Degree ~ Month, PvalAdj, data , threshold = 0.05)


group_data <- cldList(P.adj ~ Comparison,
                      data = PT,
                      threshold = 0.05 , reverse=TRUE)
#find outliers


maxdata = ddply(data,~Month,summarise,max=min(max(Degree),quantile(Degree, probs = 0.75) + (IQR(Degree) * 1.5)))

group_data1 <-merge(x = maxdata, y = group_data[,1:2] , by.x= "Month" , by.y = "Group") 


plotdunnletterdegree <- ggplot(data = data, aes(x = MonthLev , y =Degree) )  + 
  geom_boxplot( fill='#C0C0C0', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#C0C0C0' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) +
  geom_jitter(alpha = 1 , color = "black" , position = position_jitter(width = 0.3) , size = 0.2) + 
  theme(axis.text.x = element_text(angle = 63, hjust = 1 , colour="black",size = 12), 
        axis.text.y = element_text(angle=90, hjust=1 ,size = 12 , colour="black" ) , 
        axis.title.y =  element_text(angle=90,size = 12) ,
        axis.title.x = element_blank() ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(data=group_data1,aes(x=Month,y=max,label=Letter) 
            ,position = position_dodge(width = 0.9), size = 8 , vjust=-0.5 ) + ylim(0,200)

pltList[[4]] = plotdunnletterdegree

f = ggarrange(plotlist = pltList[c(3,4,2,1)],nrow=1 , ncol=4 )

ggsave("../Figb-e_NetworkPropertis.pdf" , f, width = 20 , height =5 )





