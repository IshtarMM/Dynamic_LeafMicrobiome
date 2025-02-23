# By Maryam Mahmoudi , maryam.mahmoudi@uni-tuebingen.de

# Load libreries
packages <- c("agricolae" , "tidyr" , "ggplot2" , "rowr" , "dplyr" ,
              "reshape" , "BiodiversityR" , "phyloseq" , "cowplot" , 
              "otuSummary" , "RColorBrewer" , "forcats" ,"ggforce" , "plyr" , "microbiome" , "ggpubr") # list of packages that are needed 
lapply(packages, require, character.only = TRUE)
library(funrar)

# Bacteria files
otu_matorder_aggt = read.table(file = "data/BV5_RA_order.txt" , header = T)
rownames(otu_matorder_aggt) = otu_matorder_aggt$Replicate

list = c("o__Sphingomonadales" , "o__Actinomycetales" ,"o__Burkholderiales" , "o__Rhizobiales" , "o__Bacillales", "o__Pseudomonadales"
         , "o__Deinococcales" , "o__Flavobacteriales" , "o__Xanthomonadales"  , "o__Methylophilales")
listdata = data.frame(list)

subset = otu_matorder_aggt[,colnames(otu_matorder_aggt) %in% listdata$list]
for ( col in 1:ncol(subset)){
  colnames(subset)[col] <-  sub("o__*", "", colnames(subset)[col])
}
identical(rownames(subset),rownames(otu_matorder_aggt))
subset1 = cbind(otu_matorder_aggt[11] , subset)
subset1$Month = factor(subset1$Month , levels = c("Nov","Dec","Jan" , "Feb" , "Mar"))

anova_result <- aov(Sphingomonadales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Sphingomonadales),quantile(Sphingomonadales, probs = 0.75) + (IQR(Sphingomonadales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p1 <- ggplot(subset1 , aes(x = Month , y = Sphingomonadales))  + 
  geom_boxplot( fill='#3EC423', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#3EC423' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Sphingomonadales" , y = "Relative abundance") + labs(y = "Relative abundance") 


#####
anova_result <- aov(Actinomycetales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Actinomycetales),quantile(Actinomycetales, probs = 0.75) + (IQR(Actinomycetales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p2 <- ggplot(subset1 , aes(x = Month , y = Actinomycetales))  + 
  geom_boxplot( fill='#3EC423', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#3EC423' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Actinomycetales" , y = "Relative abundance") + labs(y = "Relative abundance") + ylim(0,max(subset1$Actinomycetales)+0.1)

###
anova_result <- aov(Burkholderiales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Burkholderiales),quantile(Burkholderiales, probs = 0.75) + (IQR(Burkholderiales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p3 <- ggplot(subset1 , aes(x = Month , y = Burkholderiales))  + 
  geom_boxplot( fill='#3EC423', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#3EC423' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Burkholderiales" , y = "Relative abundance") + labs(y = "Relative abundance") + ylim(0,max(subset1$Burkholderiales)+0.1)

###

anova_result <- aov(Rhizobiales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Rhizobiales),quantile(Rhizobiales, probs = 0.75) + (IQR(Rhizobiales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p4 <- ggplot(subset1 , aes(x = Month , y = Rhizobiales))  + 
  geom_boxplot( fill='#3EC423', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#3EC423' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Rhizobiales" , y = "Relative abundance") + labs(y = "Relative abundance")  + ylim(0,max(subset1$Rhizobiales)+0.1)


###
anova_result <- aov(Bacillales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Bacillales),quantile(Bacillales, probs = 0.75) + (IQR(Bacillales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p5 <- ggplot(subset1 , aes(x = Month , y = Bacillales))  + 
  geom_boxplot( fill='#3EC423', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#3EC423' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Bacillales" , y = "Relative abundance") + labs(y = "Relative abundance") + ylim(0,max(subset1$Bacillales)+0.1)




###
anova_result <- aov(Pseudomonadales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Pseudomonadales),quantile(Pseudomonadales, probs = 0.75) + (IQR(Pseudomonadales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p6 <- ggplot(subset1 , aes(x = Month , y = Pseudomonadales))  + 
  geom_boxplot( fill='#3EC423', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#3EC423' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Pseudomonadales" , y = "Relative abundance") + labs(y = "Relative abundance") 


###
anova_result <- aov(Deinococcales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Deinococcales),quantile(Deinococcales, probs = 0.75) + (IQR(Deinococcales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p7 <- ggplot(subset1 , aes(x = Month , y = Deinococcales))  + 
  geom_boxplot( fill='#3EC423', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#3EC423' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Deinococcales" , y = "Relative abundance") + labs(y = "Relative abundance") 
###

anova_result <- aov(Flavobacteriales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Flavobacteriales),quantile(Flavobacteriales, probs = 0.75) + (IQR(Flavobacteriales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p8 <- ggplot(subset1 , aes(x = Month , y = Flavobacteriales))  + 
  geom_boxplot( fill='#3EC423', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#3EC423' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Flavobacteriales" , y = "Relative abundance") + labs(y = "Relative abundance") 

###
anova_result <- aov(Xanthomonadales~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Xanthomonadales),quantile(Xanthomonadales, probs = 0.75) + (IQR(Xanthomonadales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p9 <- ggplot(subset1 , aes(x = Month , y = Xanthomonadales))  + 
  geom_boxplot( fill='#3EC423', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#3EC423' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Xanthomonadales" , y = "Relative abundance") + labs(y = "Relative abundance") 


###
anova_result <- aov(Methylophilales~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Methylophilales),quantile(Methylophilales, probs = 0.75) + (IQR(Methylophilales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p10 <- ggplot(subset1 , aes(x = Month , y = Methylophilales))  + 
  geom_boxplot( fill='#3EC423', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#3EC423' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Methylophilales" , y = "Relative abundance") + labs(y = "Relative abundance") 

##
###############################################################################################################################################
#Fungi

otu_matorder_aggt = read.table("data/FITS_RA_order.txt" , header = T)

write.table(otu_matorder_aggt , "agg_order/histogramTableFITS2.txt" , sep = "\t" ,row.names = F)

rownames(otu_matorder_aggt) = otu_matorder_aggt$Replicate

list = c("o__Cystobasidiales" ,  "o__Sporidiobolales"  ,"o__Microbotryales" , "o__Tremellales" , "o__Capnodiales" , "o__Leucosporidiales", "k__Fungi_unclassified" , "o__Filobasidiales" ,  "o__Helotiales"   , "o__Pleosporales")
listdata = data.frame(list)

subset = otu_matorder_aggt[,colnames(otu_matorder_aggt) %in% listdata$list]
for ( col in 1:ncol(subset)){
  colnames(subset)[col] <-  sub("o__*", "", colnames(subset)[col])
  colnames(subset)[col] <-  sub("k__*", "", colnames(subset)[col])
}

identical(rownames(subset),rownames(otu_matorder_aggt))
subset1 = cbind(otu_matorder_aggt[12] , subset)
subset1$Month = factor(subset1$Month , levels = c("Nov","Dec","Jan" , "Feb" , "Mar"))

anova_result <- aov(Cystobasidiales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Cystobasidiales),quantile(Cystobasidiales, probs = 0.75) + (IQR(Cystobasidiales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p11 <- ggplot(subset1 , aes(x = Month , y = Cystobasidiales))  + 
  geom_boxplot( fill='#FF7F50', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#FF7F50' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Cystobasidiales" , y = "Relative abundance") + labs(y = "Relative abundance") 


##
anova_result <- aov(Sporidiobolales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Sporidiobolales),quantile(Sporidiobolales, probs = 0.75) + (IQR(Sporidiobolales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 

p12 <- ggplot(subset1 , aes(x = Month , y = Sporidiobolales))  + 
  geom_boxplot( fill='#FF7F50', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#FF7F50' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Sporidiobolales" , y = "Relative abundance") + labs(y = "Relative abundance") + ylim(0,max(subset1$Sporidiobolales)+0.1)
##

anova_result <- aov(Microbotryales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Microbotryales),quantile(Microbotryales, probs = 0.75) + (IQR(Microbotryales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p13 <- ggplot(subset1 , aes(x = Month , y = Microbotryales))  + 
  geom_boxplot( fill='#FF7F50', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#FF7F50' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Microbotryales" , y = "Relative abundance") + labs(y = "Relative abundance") 

##

anova_result <- aov(Tremellales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Tremellales),quantile(Tremellales, probs = 0.75) + (IQR(Tremellales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p14 <- ggplot(subset1 , aes(x = Month , y = Tremellales))  + 
  geom_boxplot( fill='#FF7F50', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#FF7F50' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Tremellales" , y = "Relative abundance") + labs(y = "Relative abundance") 

##

anova_result <- aov(Capnodiales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Capnodiales),quantile(Capnodiales, probs = 0.75) + (IQR(Capnodiales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p15 <- ggplot(subset1 , aes(x = Month , y = Capnodiales))  + 
  geom_boxplot( fill='#FF7F50', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#FF7F50' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Capnodiales" , y = "Relative abundance") + labs(y = "Relative abundance") 

##

anova_result <- aov(Leucosporidiales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Leucosporidiales),quantile(Leucosporidiales, probs = 0.75) + (IQR(Leucosporidiales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p16 <- ggplot(subset1 , aes(x = Month , y = Leucosporidiales))  + 
  geom_boxplot( fill='#FF7F50', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#FF7F50' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Leucosporidiales" , y = "Relative abundance") + labs(y = "Relative abundance") 

##
anova_result <- aov(Fungi_unclassified ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Fungi_unclassified),quantile(Fungi_unclassified, probs = 0.75) + (IQR(Fungi_unclassified) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p17 <- ggplot(subset1 , aes(x = Month , y = Fungi_unclassified))  + 
  geom_boxplot( fill='#FF7F50', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#FF7F50' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Fungi_unclassified" , y = "Relative abundance") + labs(y = "Relative abundance") 

##
anova_result <- aov(Filobasidiales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Filobasidiales),quantile(Filobasidiales, probs = 0.75) + (IQR(Filobasidiales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p18 <- ggplot(subset1 , aes(x = Month , y = Filobasidiales))  + 
  geom_boxplot( fill='#FF7F50', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#FF7F50' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Filobasidiales" , y = "Relative abundance") + labs(y = "Relative abundance") 

##
anova_result <- aov(Helotiales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Helotiales),quantile(Helotiales, probs = 0.75) + (IQR(Helotiales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p19 <- ggplot(subset1 , aes(x = Month , y = Helotiales))  + 
  geom_boxplot( fill='#FF7F50', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#FF7F50' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Helotiales" , y = "Relative abundance") + labs(y = "Relative abundance") 

##
anova_result <- aov(Pleosporales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Pleosporales),quantile(Pleosporales, probs = 0.75) + (IQR(Pleosporales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p20 <- ggplot(subset1 , aes(x = Month , y = Pleosporales))  + 
  geom_boxplot( fill='#FF7F50', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#FF7F50' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Pleosporales" , y = "Relative abundance") + labs(y = "Relative abundance") + ylim(0,max(subset1$Pleosporales)+0.1)


# oomycete

otu_matorder_aggt = read.table("data/PrV9_RA_order.txt" , header = T)
rownames(otu_matorder_aggt) = otu_matorder_aggt$Replicate

list = c("o__Peronosporales" ,  "o__Pythiales"  ,"o__Albuginales")

listdata = data.frame(list)

subset = otu_matorder_aggt[,colnames(otu_matorder_aggt) %in% listdata$list]
for ( col in 1:ncol(subset)){
  colnames(subset)[col] <-  sub("o__*", "", colnames(subset)[col])
}


identical(rownames(subset),rownames(otu_matorder_aggt))
subset1 = cbind(otu_matorder_aggt[12] , subset)
subset1$Month = factor(subset1$Month , levels = c("Nov","Dec","Jan" , "Feb" , "Mar"))

##

anova_result <- aov(Peronosporales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Peronosporales),quantile(Peronosporales, probs = 0.75) + (IQR(Peronosporales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p21 <- ggplot(subset1 , aes(x = Month , y = Peronosporales))  + 
  geom_boxplot( fill='#33D4FF', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#33D4FF' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Peronosporales" , y = "Relative abundance") + labs(y = "Relative abundance") + ylim(0,max(subset1$Peronosporales)+0.1)

##

anova_result <- aov(Pythiales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Pythiales),quantile(Pythiales, probs = 0.75) + (IQR(Pythiales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p22 <- ggplot(subset1 , aes(x = Month , y = Pythiales))  + 
  geom_boxplot( fill='#33D4FF', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#33D4FF' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Pythiales" , y = "Relative abundance") + labs(y = "Relative abundance") + ylim(0,max(subset1$Pythiales)+0.1)

##
anova_result <- aov(Albuginales ~ Month, subset1)
HSD = HSD.test(anova_result, "Month", group = TRUE)
group_data <- HSD$groups[order(rownames(HSD$groups)),]
group_data$var <- rownames(group_data)
maxdata = ddply(subset1,~Month,summarise,max=min(max(Albuginales),quantile(Albuginales, probs = 0.75) + (IQR(Albuginales) * 1.5)))
group_data1 <-merge(x = maxdata, y = group_data[,2:3] , by.x= "Month" , by.y = "var") 
#
p23 <- ggplot(subset1 , aes(x = Month , y = Albuginales))  + 
  geom_boxplot( fill='#33D4FF', color="black"  , outlier.shape=NA ,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill='#33D4FF' ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 14), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 14 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 14) ,axis.title.x =  element_text(size = 14) ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,
        plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=group_data1,aes(x=Month,y=max,label=groups) ,position = position_dodge(width = 0.9), size = 4 , vjust=-0.5 )  +
  labs(title = "Albuginales" , y = "Relative abundance") + labs(y = "Relative abundance") + ylim(0,max(subset1$Albuginales)+0.1)


all = ggarrange(p1, p2, p3,p4,p5,p6,p7,p8,p9,p10,p11, p12, p13,p14,p15,p16,p17,p18,p19,p20 ,p21,p22,p23 ,  ncol = 10, nrow = 3 )

ggsave("Supp_Fig1.pdf" , all , width = 30 , height = 10 , units = "in")







