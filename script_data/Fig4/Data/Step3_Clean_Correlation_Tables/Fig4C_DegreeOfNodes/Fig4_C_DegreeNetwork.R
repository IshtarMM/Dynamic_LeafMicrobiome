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

####Figure 4 "c" degree of nodes per month

filelist = c("11.csv","12.csv","1.csv","2.csv","3.csv") # Degree of nodes are calculated in Cytoscape
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
## finall tukey test 
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
ggsave("Fig4_C_Degree.pdf" , f, width = 5 , height =5 )