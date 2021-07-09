# By Maryam Mahmoudi , maryam.mahmoudi@uni-tuebingen.de

packages <- c("readxl" , "vegan" , "tidyr" , "ggplot2" , "rowr" , "dplyr" ,
              "reshape" , "BiodiversityR" , "phyloseq" , "cowplot" , "rcompanion",
              "otuSummary" , "RColorBrewer" , "forcats" ,"ggforce" , "plyr" , "microbiome" ,"FSA") # list of packages that are needed 
library(agricolae)
library(gridExtra) # for grid.arrange
library(grid) 
lapply(packages, require, character.only = TRUE)


# alpha diversity

listAfiles = c("data/BV5.txt","data/FITS2.txt" , "data/PV9.txt")
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

# Distance to centroid
x = 3
listAfilesRA = c("data/BV5RAlog10.txt","data/FITS2RAlog10.txt" , "data/PV9RAlog10.txt")
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


# Between month variability

# Bacteria
tab_bac<-read.table("data/BV5RAlog10.txt", h=T)
tab_bac2 <- tab_bac[-c(1:17)]
samp_bac2 <- tab_bac[1:17]
bac_dist<-vegdist(tab_bac2, method="bray")

b1 = as.matrix(bac_dist)
#### compare consecutive months

#distances between nov vs dec samples
b.nov.dec<-as.matrix(bac_dist)[samp_bac2$Month=="Nov",samp_bac2$Month=="Dec"]
b.nov.dec2<-as.vector(b.nov.dec)
b.nov.dec3<-as.data.frame(cbind(b.nov.dec2, rep("Nov_Dec", length(b.nov.dec2))))
names(b.nov.dec3) [1]<-"bc"
names(b.nov.dec3) [2]<-"cond"

#distances between dec vs jan samples
b.dec.jan<-as.matrix(bac_dist)[samp_bac2$Month=="Dec",samp_bac2$Month=="Jan"]
b.dec.jan2<-as.vector(b.dec.jan)
b.dec.jan3<-as.data.frame(cbind(b.dec.jan2, rep("Dec_Jan", length(b.dec.jan2))))
names(b.dec.jan3) [1]<-"bc"
names(b.dec.jan3) [2]<-"cond"

#distances between jan vs feb samples
b.jan.feb<-as.matrix(bac_dist)[samp_bac2$Month=="Jan",samp_bac2$Month=="Feb"]
b.jan.feb2<-as.vector(b.jan.feb)
b.jan.feb3<-as.data.frame(cbind(b.jan.feb2, rep("Jan_Feb", length(b.jan.feb2))))
names(b.jan.feb3) [1]<-"bc"
names(b.jan.feb3) [2]<-"cond"


#distances between feb vs mar samples
b.feb.mar<-as.matrix(bac_dist)[samp_bac2$Month=="Feb",samp_bac2$Month=="Mar"]
b.feb.mar2<-as.vector(b.feb.mar)
b.feb.mar3<-as.data.frame(cbind(b.feb.mar2, rep("Feb_Mar", length(b.feb.mar2))))
names(b.feb.mar3) [1]<-"bc"
names(b.feb.mar3) [2]<-"cond"

b.distances<-rbind(b.nov.dec3,b.dec.jan3,b.jan.feb3,b.feb.mar3)
b.distances$bc<-as.numeric(as.character(b.distances$bc))

DT = dunnTest(bc~cond,data=b.distances,method="bh") 
PT = DT$res
PT
rownames(PT) = PT$Comparison
group_data <- cldList(P.adj ~ Comparison,data = PT,threshold = 0.05 , reverse=TRUE)

scale = function(x) sprintf("%.1f",x)
p7 <- ggplot(b.distances , aes(x = cond , y = bc))  + 
  geom_boxplot( fill="#009E73", color="black"  , outlier.shape=NA ,outlier.size=-1,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill="#009E73" ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 10), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 10 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 10) ,
        axis.title.x = element_blank() ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_continuous(labels =scale,limits = c(0,1.05) )  +
  scale_x_discrete(labels=c("Nov","Dec","Jan", "Feb"))
pltList[[7]] = p7

#stats
library(agricolae)
print(HSD.test(aov(bc~cond, data=b.distances), "cond"))

########################################################################################################################################
###fungi

tab_fun<-read.table("data/FITS2RAlog10.txt", h=T)
tab_fun2 <- tab_fun[-c(1:17)]
samp_fun2 <- tab_fun[1:17]
fun_dist<-vegdist(tab_fun2, method="bray")

#### compare consecutive months

#distances between nov vs dec samples
f.nov.dec<-as.matrix(fun_dist)[samp_fun2$Month=="Nov",samp_fun2$Month=="Dec"]
f.nov.dec2<-as.vector(f.nov.dec)
f.nov.dec3<-as.data.frame(cbind(f.nov.dec2, rep("Nov_Dec", length(f.nov.dec2))))
names(f.nov.dec3) [1]<-"bc"
names(f.nov.dec3) [2]<-"cond"

#distances between dec vs jan samples
f.dec.jan<-as.matrix(fun_dist)[samp_fun2$Month=="Dec",samp_fun2$Month=="Jan"]
f.dec.jan2<-as.vector(f.dec.jan)
f.dec.jan3<-as.data.frame(cbind(f.dec.jan2, rep("Dec_Jan", length(f.dec.jan2))))
names(f.dec.jan3) [1]<-"bc"
names(f.dec.jan3) [2]<-"cond"

#distances between jan vs feb samples
f.jan.feb<-as.matrix(fun_dist)[samp_fun2$Month=="Jan",samp_fun2$Month=="Feb"]
f.jan.feb2<-as.vector(f.jan.feb)
f.jan.feb3<-as.data.frame(cbind(f.jan.feb2, rep("Jan_Feb", length(f.jan.feb2))))
names(f.jan.feb3) [1]<-"bc"
names(f.jan.feb3) [2]<-"cond"


#distances between feb vs mar samples
f.feb.mar<-as.matrix(fun_dist)[samp_fun2$Month=="Feb",samp_fun2$Month=="Mar"]
f.feb.mar2<-as.vector(f.feb.mar)
f.feb.mar3<-as.data.frame(cbind(f.feb.mar2, rep("Feb_Mar", length(f.feb.mar2))))
names(f.feb.mar3) [1]<-"bc"
names(f.feb.mar3) [2]<-"cond"

f.distances<-rbind(f.nov.dec3,f.dec.jan3,f.jan.feb3,f.feb.mar3)
f.distances$bc<-as.numeric(as.character(f.distances$bc))
boxplot(bc~cond, data=f.distances)

DT = dunnTest(bc~cond,data=f.distances,method="bh") 
PT = DT$res
PT
rownames(PT) = PT$Comparison

#multcompLetters2(distance ~ Month, PvalAdj, dis_to_center1 , threshold = 0.05)

group_data <- cldList(P.adj ~ Comparison,data = PT,threshold = 0.05 , reverse=TRUE)

scale = function(x) sprintf("%.1f",x)
p8 <- ggplot(f.distances , aes(x = cond , y = bc))  + 
  geom_boxplot( fill="#F70C0C", color="black"  , outlier.shape=NA ,outlier.size=-1,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill="#F70C0C" ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 10), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 10 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 10) ,
        axis.title.x = element_blank() ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_continuous(labels =scale,limits = c(0,1.05) )  + 
  scale_x_discrete(labels=c("Nov","Dec","Jan", "Feb"))
pltList[[8]] = p8

print(HSD.test(aov(bc~cond, data=f.distances), "cond"))

########################################################################################################################################
###oomyc

tab_oom <- read.table("data/PV9RAlog10.txt", h=T)
tab_oom2 <- tab_oom[-c(1:17)]
samp_oom2 <- tab_oom[1:17]
oom_dist<-vegdist(tab_oom2, method="bray")

#### compare consecutive months

#distances between nov vs dec samples
o.nov.dec<-as.matrix(oom_dist)[samp_oom2$Month=="Nov",samp_oom2$Month=="Dec"]
o.nov.dec2<-as.vector(o.nov.dec)
o.nov.dec3<-as.data.frame(cbind(o.nov.dec2, rep("Nov_Dec", length(o.nov.dec2))))
names(o.nov.dec3) [1]<-"bc"
names(o.nov.dec3) [2]<-"cond"

#distances between dec vs jan samples
o.dec.jan<-as.matrix(oom_dist)[samp_oom2$Month=="Dec",samp_oom2$Month=="Jan"]
o.dec.jan2<-as.vector(o.dec.jan)
o.dec.jan3<-as.data.frame(cbind(o.dec.jan2, rep("Dec_Jan", length(o.dec.jan2))))
names(o.dec.jan3) [1]<-"bc"
names(o.dec.jan3) [2]<-"cond"

#distances between jan vs feb samples
o.jan.feb<-as.matrix(oom_dist)[samp_oom2$Month=="Jan",samp_oom2$Month=="Feb"]
o.jan.feb2<-as.vector(o.jan.feb)
o.jan.feb3<-as.data.frame(cbind(o.jan.feb2, rep("Jan_Feb", length(o.jan.feb2))))
names(o.jan.feb3) [1]<-"bc"
names(o.jan.feb3) [2]<-"cond"


#distances between feb vs mar samples
o.feb.mar<-as.matrix(oom_dist)[samp_oom2$Month=="Feb",samp_oom2$Month=="Mar"]
o.feb.mar2<-as.vector(o.feb.mar)
o.feb.mar3<-as.data.frame(cbind(o.feb.mar2, rep("Feb_Mar", length(o.feb.mar2))))
names(o.feb.mar3) [1]<-"bc"
names(o.feb.mar3) [2]<-"cond"

o.distances<-rbind(o.nov.dec3,o.dec.jan3,o.jan.feb3,o.feb.mar3)
o.distances$bc<-as.numeric(as.character(o.distances$bc))
boxplot(bc~cond, data=o.distances)


DT = dunnTest(bc~cond,data=o.distances,method="bh") 
PT = DT$res
PT
rownames(PT) = PT$Comparison

#multcompLetters2(distance ~ Month, PvalAdj, dis_to_center1 , threshold = 0.05)

group_data <- cldList(P.adj ~ Comparison,data = PT,threshold = 0.05 , reverse=TRUE)

scale = function(x) sprintf("%.1f",x)
p9 <- ggplot(o.distances , aes(x = cond , y = bc))  + 
  geom_boxplot( fill="#0072B2", color="black"  , outlier.shape=NA ,outlier.size=-1,width=0.2 , lwd = 0.4 ,linetype = "dashed" , 
                position = position_dodge(width=6))  + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..) ,fill="#0072B2" ,outlier.shape=NA ,width=0.7 , lwd = 0.4) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..) ,width=0.3) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1 , colour="black", size = 10), 
        axis.text.y = element_text(angle=90, hjust=1 , size = 10 , colour="black" ) , 
        axis.title.y =  element_text(angle=90, size = 10) ,
        axis.title.x = element_blank() ,
        panel.background = element_rect(colour = "black" , fill = "NA"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_continuous(labels =scale,limits = c(0,1.05) )  +
  scale_x_discrete(labels=c("Nov","Dec","Jan", "Feb"))
pltList[[9]] = p9



# save plots
x = grid.arrange(grobs = lapply(pltList[1:3], "+", theme(plot.margin=margin(10,10,10,10)))) # alpha diversity
y = grid.arrange(grobs = lapply(pltList[4:6], "+", theme(plot.margin=margin(10,10,10,10)))) # distance to centroid
z = grid.arrange(grobs = lapply(pltList[7:9], "+", theme(plot.margin=margin(10,10,10,10)))) # between month variability
f = grid.arrange(x,y,z , ncol = 3) #z come form Bc_time_compa_110221.R
ggsave("Fig3_diversity.pdf" , f, width = 18 , height = 18 , units = "cm")














