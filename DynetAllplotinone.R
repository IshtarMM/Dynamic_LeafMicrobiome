#f
library(multcompView)
library(ggplot2)
library(tidyr)
# Install
#install.packages("wesanderson")
# Load
library(wesanderson)
library(ggpubr)
#install.packages("EnvStats")
library(EnvStats)
library(tidyr)
# Install
#install.packages("wesanderson")
# Load
library(wesanderson)
library(ggpubr)
#install.packages("EnvStats")
library(EnvStats)
library(rstatix)
library(agricolae)
library(plyr)
library(FSA)
library(rcompanion)
########################################################################################################################################
pltList = list()


setwd("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/NetworkResult/ResultSparcc/Otu/Month/Pval/NetFeature1/Dynet/")
#FigureF
#plots for dynet features rewiring and ...
filelist = c("11-12node.csv","12-1node.csv","1-2node.csv","2-3node.csv")
#Monthname = c("Novemnber", "December" ,"January" , "February" , "March")
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
#finall Rewiring  dunn test withe letter 
data$MonthLev <- factor(data$Month, levels = c("November_December","December_January","January_February","February_March"))
#res.kruskal <- data%>% kruskal_test(RewiringDnScore ~ MonthLev)
#dunn <- data%>% dunn_test(RewiringDnScore ~ MonthLev, p.adjust.method = "BH")
#dunn <- dunn %>% add_xy_position(x = "MonthLev")

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

#ggsave("Plots/plotdunnletterRewairing.jpg" , plotdunnletterRewairing  ,  width = 10 , height = 15 , units = "cm")
#ggsave("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/Dataset2019/results/PlotsPaper/plotdunnletterRewairing.pdf" , plotdunnletterRewairing  ,  width = 10 , height = 15 , units = "cm")

#Figure b and c
#######finall this info comes from top codes
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

##########finall inheritededges and inheritednodes inone plot
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

####Figure d
setwd("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/NetworkResult/ResultSparcc/Otu/Month/Pval/NetFeature1/Dynet/networks_features/")


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

Node <- read.csv(filelist[5])
hist(Node$Degree)

install.packages("igraph")
library(igraph)
install.packages("poweRlaw")
library(poweRlaw)
set.seed(202)
g <- static.power.law.game(500, 1000, exponent.out= 2.2, exponent.in = -1, loops = FALSE, multiple = TRUE, finite.size.correction = TRUE)
plot(g, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "")
data <- degree(g)
data <- data[data>0]
Node <- read.csv(filelist[5])
m_pl <- displ$new(data)
est_pl <- estimate_xmin(m_pl)
est_pl$distance

## load "eSetObject" containing simulated time-course data

library(splineTimeR)
f = data.frame(data(TCsimData))

## reconstruct gene association networks from time-course data
igr <- splineNetRecon(eSet = TCsimData, treatmentType = "T2", cutoff.ggm = c(0.8,0.9))

## check for scale-free properties of reconstructed networks (igraphs)
scaleFreeProp <- networkProperties(Node$Degree)
head(scaleFreeProp)

## the functional interaction pairs provided in FIs data package
library(FIs)
data(FIs)
names(FIs)
head(FIs$FIs_Reactome)
head(FIs$FIs_BioGRID)








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


####figure d

setwd("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/NetworkResult/ResultSparcc/Otu/Month/Pval/NetFeature1/Dynet/")

filelist = c("11-12edge.csv","12-1edge.csv","1-2edge.csv","2-3edge.csv")
Monthname = c("December" ,"January" , "February" , "March")
data <- data.frame()
for (i in 1:4){
  edge <- read.csv(filelist[i])
  print(Monthname[i])
  print(basename(filelist[i]))
  edge1 <- separate(edge,shared.name, c("Node1" , "Node2"), sep = "(interacts with)", remove = FALSE) #separate "intraction" part to node1 and node2 
  edge2 <- edge1[,c(4,8,16,17)] #keep edge present in network one and network2 and intraction (Node1 and node2)
  colnames(edge2) <- c("FirstMonthPresent","SecoundMonthPresent","Node1","Node2")
  edge2[] <- lapply(edge2, gsub, pattern='[( ]', replacement='')
  edge2[] <- lapply(edge2, gsub, pattern='[) ]', replacement='')
  edge2$Status <- paste(edge2$FirstMonthPresent,edge2$SecoundMonthPresent ,sep="") #concat present of edges in network1 and network2
  
  Node1 = edge2[,c(3,5)]
  Node2 = edge2[,c(4,5)]
  colnames(Node1) <- c("Node","Status") #
  colnames(Node2) <- c("Node","Status")# some nodes only are here
  Nodefile <- rbind(Node1,Node2)
  
  Nodefile1 =  Nodefile %>% group_by(Node,Status) %>% dplyr::summarise(number = n()) #number of intraction per node based on (present in only net1 truefalse, only net two falsetrue, present in both tt) 
  #convert to matrix shape
  
  Nodefile2 = with(Nodefile1, {
    out <- matrix(nrow=nlevels(factor(Node)), ncol=nlevels(factor(Status)),
                  dimnames=list(levels(factor(Node)), levels(factor(Status))))
    out[cbind(Node,Status)] <- number
    out
  })
  
  Nodefile2[is.na( Nodefile2)] = 0 #keep edges that are in both network TrueTrue
  Nodefile3 = data.frame(Nodefile2)
  Nodefile3$inherited = Nodefile3$truetrue/(Nodefile3$falsetrue + Nodefile3$truetrue) #inherited edges per node = number of edges present in both networks(TrueTrue) devided to sum(namber of TrueTrue and falsetrue )
  Nodefile3$Month <- Monthname[i]
  data <- rbind(data,Nodefile3)
}

##########################################################################
#Inherited edges per nodes dunn test withe letter plot. 
#data[is.na( data)] = 0
data = data[data$inherited!="NaN",]
data$MonthLev <- factor(data$Month, levels = c("December","January","February","March"))

DT = dunnTest(inherited ~ MonthLev,data=data,method="bh") 
PT = DT$res
PT
PvalAdj = PT[,4]
names(PvalAdj) = c("December-February" ,"December-January" , "February-January" , "December-March","February-March"  ,  "January-March"  )

rownames(PT) = PT$Comparison
multcompLetters2(inherited ~ Month, PvalAdj, data , threshold = 0.05 )


#group_data <- cldList(P.adj ~ Comparison,data = PT,threshold = 0.05,reverse=TRUE)


#find outliers

maxdata = ddply(data,~Month,summarise,max=min(max(inherited),quantile(inherited, probs = 0.75) + (IQR(inherited) * 1.5)))

group_data1 <-merge(x = maxdata, y = group_data[,1:2] , by.x= "Month" , by.y = "Group") 

data1<- data %>% 
  mutate( inheritedout = inherited > median(inherited) + 
            IQR(inherited)*1.5 | inherited < median(inherited) - IQR(inherited)*1.5)


#add sudo November 
x = data1[data1$MonthLev=="December",1:4]
x$Month = "November"
x$MonthLev = "November"
x $inheritedout = "FALSE"
colnames(x) 
data2 = rbind(data1,x)

data2$MonthLev <- factor(data2$Month, levels = c("November","December","January","February","March"))

scale = function(x) sprintf("%.1f",x)
plot <- ggplot(data = data2, aes(x = MonthLev , y =inherited) )  + 
  geom_boxplot( fill='#C0C0C0', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#C0C0C0' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) +
  geom_jitter(data = function(x) dplyr::filter(x, inheritedout=="FALSE"),alpha = 1 , color = "black" , position = position_jitter(width = 0.3) , size = 0.2) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black",size = 12), 
        axis.text.y = element_text(angle=90, hjust=1 ,size = 12 , colour="black" ) , 
        axis.title.y =  element_text(angle=90,size = 12) ,
        axis.title.x = element_blank() ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(y="Inherited edges per node") + scale_y_continuous(labels = scale , limits = c(0,1.2)) +
  geom_text(data=group_data1,aes(x=Month,y=max,label=Letter) 
            ,position = position_dodge(width = 0.9), size = 8 , vjust=-0.5 ) 
 

pltList[[5]] = plot


f = ggarrange(plotlist = pltList[c(3,4,2,5,1)],nrow=1 , ncol=5 )
f = grid.arrange(grobs = lapply(pltList[c(3,4,2,5)], "+", theme(plot.margin=margin(10,10,10,10))) , nrow=1 , ncol=4)
f1 = ggarrange(plotlist = pltList[c(1)],nrow=1 , ncol=1)

ggsave("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/Dataset2019/results/PlotsPaper/boxDznetall.pdf" , f, width = 25 , height =5 )
ggsave("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/Dataset2019/results/PlotsPaper/boxDznetall1.pdf" , f1, width = 5 , height =5 )





