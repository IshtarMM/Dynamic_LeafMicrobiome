# By Maryam Mahmoudi , maryam.mahmoudi@uni-tuebingen.de

# Load libreries
packages <- c("readxl" , "vegan" , "tidyr" , "ggplot2" , "rowr" , "dplyr" ,
              "reshape" , "BiodiversityR" , "phyloseq" , "cowplot" , 
              "otuSummary" , "RColorBrewer" , "forcats" ,"ggforce" , "plyr"  , "ggpubr" , "gridExtra") # list of packages that are needed 
lapply(packages, require, character.only = TRUE)

# open Bacteria otu and taxa files

otu <- read.table("data/BV5_OtuTable.txt" , header = T)
taxa <- read.table("data/BV5_TaxonomyTable.txt" , header = T)
otu1 <- select(otu ,-c(1:4,6:16)) #select otu table and replicate
rownames(otu1) <- otu1$Replicate  #put Replicate as a rownames of otutable
otu1 <- otu1[-1]   # remove Replicate column 

#create phyloseq file
otu_mat <- t(otu1) 
sampleinfo <- otu[1:16]  #samples information
sampleinfo$Experiment <- sub("^", "Year", sampleinfo$Experiment)
rownames(sampleinfo) <- sampleinfo$Replicate
tax_mat <- select(taxa,-c(2,3))
rownames(tax_mat) <- tax_mat$BacV5_Otu_
tax_mat <- select(tax_mat, -c(1))
colnames(tax_mat) <- c("Kingdom" , "Phylum" , "Class" , "Order" , "Family" , "Genus" ,"Species")
otu_phylo <- otu_table(otu_mat , taxa_are_rows = TRUE,errorIfNULL = TRUE)
taxa_phylo <- tax_table(as.matrix(tax_mat))
samples_phylo <- sample_data(sampleinfo)
file <- phyloseq(otu_phylo,taxa_phylo,samples_phylo)

phy <- transform_sample_counts(file, function(x) x/sum(x))  ###convert to reletive abundance

# agglomerate taxa
#glom <- tax_glom(phy, taxrank = 'Species')

# create dataframe from phyloseq object
dat <- psmelt(phy)
y = dat
dat1 = y
# convert Phylum to a character vector from a factor because R
dat1$Order <- as.character(dat1$Order)

#dat1$Order[dat1$Abundance <= 0.2] <- "Other" #rename genera with < 1% abundance
agg = aggregate(dat1[,3], list(dat1$Order), sum) ## sum of relative abundance of different orders in all samples 
colnames(agg) = c("order1" , "sumAbundance")
agg$order2 = agg$order1  # copy of this col because i want to have main col and change the order of this column to "other" category
agg$order2 = as.character(agg$order2)
agg$order2[agg$sumAbundance < 2.4] <- "Other"  ##change the order of otus that have less than a threshold to "Other"
dat2 = merge(dat1 , agg ,by.x ="Order" , by.y ="order1")  #add the information of changing orders name to main data frame
#positions = c("Nov" , "Dec" , "Jan" , "Feb" , "Mar")  # ordrt of month to seeing in the plots
dat2$Month <- factor(dat2$Month, levels=c("Nov" , "Dec" , "Jan" , "Feb" , "Mar")) # change the ordrt of months to seeing in the plots
#year1$order2 <- factor(year1$order2, levels=rev(levels(factor(year1$order2))))
dat2$order2 <- factor(dat2$order2, levels=rev(c("o__Sphingomonadales" , "o__Actinomycetales" ,"o__Burkholderiales" , "o__Rhizobiales" , "o__Bacillales", "o__Pseudomonadales"
                                                , "o__Deinococcales" , "o__Flavobacteriales" , "o__Xanthomonadales"  , "o__Methylophilales" , "Other")))  #change the sorts of order based on your sort to see in plots.

col = c("#0F85B0","#FD861E","#F9DA73" , "#16B7B2" , "#AF65B2","#905E45","#416B4C","#CE5A17","#776804","#2B2B83","#999999")
y1 = dat2
yy1 = aggregate(y1$Abundance, by=list(y1$order2  , y1$Experiment , y1$Sample,y1$Month), FUN=sum)
colnames(yy1) = c("order2" , "Year" , "Sample", "Month" , "Abundance")
y1mvals = data.frame(yy1[ yy1$order == "o__Sphingomonadales", ])
y1mxvals = y1mvals[with(y1mvals, order(Abundance)), ]$Sample
yy1$Sample <- factor(yy1$Sample, levels = y1mxvals)
year1 = yy1

scale = function(x) sprintf("%.1f",x)
ploty1 <- ggplot()  + geom_bar(data=year1, aes(x=Sample, y=Abundance, fill = order2), stat="identity", position="fill" , width =1) +
  facet_wrap(.~ Month, scales = 'free_x' , nrow = 1 , strip.position="bottom") + 
  theme_minimal_grid() +
  theme(text = element_text(size=12),
        axis.text.y = element_text(size = 12 ),
        strip.text.x = element_text(size = 12),
        panel.spacing.x=unit(0, "lines"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values=rev(col))  +
  scale_y_continuous(labels =scale,limits = c(0,1) )




# Fungi bar plot

otu <- read.table("data/FITS2_OtuTable.txt" , header = T)
taxa <- read.table("data/FITS2_TaxonomyTable.txt" , header = T)
otu1 <- select(otu ,-c(1:4,6:16)) #select otu table and replicate
rownames(otu1) <- otu1$Replicate  #put Replicate as a rownames of otutable
otu1 <- otu1[-1]   # remove Replicate column 

#create phyloseq file
otu_mat <- t(otu1) 
sampleinfo <- otu[1:16]  #samples information
sampleinfo$Experiment <- sub("^", "Year", sampleinfo$Experiment)
rownames(sampleinfo) <- sampleinfo$Replicate
tax_mat <- select(taxa,-c(2,3))
rownames(tax_mat) <- tax_mat$Ftrad_Otu_
tax_mat <- select(tax_mat, -c(1))
colnames(tax_mat) <- c("Kingdom" , "Phylum" , "Class" , "Order" , "Family" , "Genus" ,"Species")
otu_phylo <- otu_table(otu_mat , taxa_are_rows = TRUE,errorIfNULL = TRUE)
taxa_phylo <- tax_table(as.matrix(tax_mat))
samples_phylo <- sample_data(sampleinfo)
file <- phyloseq(otu_phylo,taxa_phylo,samples_phylo)

phy <- transform_sample_counts(file, function(x) x/sum(x))  ###convert to reletive abundance

# agglomerate taxa
#glom <- tax_glom(phy, taxrank = 'Species')

# create dataframe from phyloseq object
dat <- psmelt(phy)
y = dat
dat1 = y
# convert Phylum to a character vector from a factor because R
dat1$Order <- as.character(dat1$Order)

agg = aggregate(dat1[,3], list(dat1$Order), sum) ## sum of relative abundance of different orders in all samples 
colnames(agg) = c("order1" , "sumAbundance")
agg$order2 = agg$order1  # copy of this col because i want to have main col and change the order of this column to "other" category
agg$order2 = as.character(agg$order2)
agg$order2[agg$sumAbundance < 5.6] <- "Other"  ##change the order of otus that have less than a threshold to "Other"
dat2 = merge(dat1 , agg ,by.x ="Order" , by.y ="order1")  #add the information of changing orders name to main data frame
#positions = c("Nov" , "Dec" , "Jan" , "Feb" , "Mar")  # ordrt of month to seeing in the plots
dat2$Month <- factor(dat2$Month, levels=c("Nov" , "Dec" , "Jan" , "Feb" , "Mar")) # change the ordrt of months to seeing in the plots
#year1$order2 <- factor(year1$order2, levels=rev(levels(factor(year1$order2))))
dat2$order2 <- reorder(dat2$order2 , dat2$sumAbundance)
dat2$order2 <- factor(dat2$order2, levels=rev(c("o__Cystofilobasidiales" , "o__Sporidiobolales" ,"o__Microbotryales" , "o__Tremellales" , "o__Capnodiales", "o__Leucosporidiales"
                                                , "k__Fungi_unclassified" , "o__Filobasidiales" , "o__Helotiales"  , "o__Pleosporales" , "Other")))  #change the sorts of order based on your sort to see in plots.


col = c("#789B44","#FD861E","#F9DA73" , "#16B7B2" , "#AF65B2","#905E45","#416B4C","#CE5A17","#776804","#2B2B83","#999999")

y1 = dat2
yy1 = aggregate(y1$Abundance, by=list(y1$order2 ,y1$Experiment , y1$Month , y1$Sample), FUN=sum)
colnames(yy1) = c("order2" , "Year","Month" , "Sample1" , "Abundance")
y1mvals = data.frame(yy1[ yy1$order == "o__Cystofilobasidiales", ])
y1mxvals = y1mvals[with(y1mvals, order(Abundance)), ]$Sample
yy1$Sample <- factor(yy1$Sample, levels = y1mxvals)
year1 = yy1



ploty2 <- ggplot()  + geom_bar(data=year1, aes(x=Sample, y=Abundance, fill = order2), stat="identity", position="fill" , width =1) +
  facet_wrap(.~ Month, scales = 'free_x' , nrow = 1 , strip.position="bottom") + 
  theme_minimal_grid() +
  theme(text = element_text(size=12 ),
        axis.text.y = element_text(size = 12 ),
        strip.text.x = element_text(size = 12),
        panel.spacing.x=unit(0, "lines"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values=rev(col))  +
  scale_y_continuous(labels =scale,limits = c(0,1) )


# Oomycete barplot
otu <- read.table("data/PrV9_OtuTable.txt" , header = T)
taxa <- read.table("data/PrV9_TaxonomyTable.txt" , header = T)
otu1 <- select(otu ,-c(1:4,6:16)) #select otu table and replicate
rownames(otu1) <- otu1$Replicate  #put Replicate as a rownames of otutable
otu1 <- otu1[-1]   # remove Replicate column 

#create phyloseq file
otu_mat <- t(otu1) 
sampleinfo <- otu[1:16]  #samples information
sampleinfo$Experiment <- sub("^", "Year", sampleinfo$Experiment)
rownames(sampleinfo) <- sampleinfo$Replicate
tax_mat <- select(taxa,-c(2,3))
rownames(tax_mat) <- tax_mat$Otrad_Otu_
tax_mat <- select(tax_mat, -c(1))
colnames(tax_mat) <- c("Kingdom" , "Phylum" , "Class" , "Order" , "Family" , "Genus" ,"Species")
otu_phylo <- otu_table(otu_mat , taxa_are_rows = TRUE,errorIfNULL = TRUE)
taxa_phylo <- tax_table(as.matrix(tax_mat))
samples_phylo <- sample_data(sampleinfo)
file <- phyloseq(otu_phylo,taxa_phylo,samples_phylo)

phy <- transform_sample_counts(file, function(x) x/sum(x))  ###convert to reletive abundance
# create dataframe from phyloseq object
dat <- psmelt(phy)
y = dat
dat1 = y
# convert Phylum to a character vector from a factor because R
dat1$Order <- as.character(dat1$Order)
agg = aggregate(dat1[,3], list(dat1$Order), sum) ## sum of relative abundance of different orders in all samples 
colnames(agg) = c("order1" , "sumAbundance")
agg$order2 = agg$order1  # copy of this col because i want to have main col and change the order of this column to "other" category
agg$order2 = as.character(agg$order2)
agg$order2[agg$sumAbundance =< 1.160096e+00] <- "Other"  ##change the order of otus that have less than a threshold to "Other"
dat2 = merge(dat1 , agg ,by.x ="Order" , by.y ="order1")  #add the information of changing orders name to main data frame
dat2$Month <- factor(dat2$Month, levels=c("Nov" , "Dec" , "Jan" , "Feb" , "Mar")) # change the ordrt of months to seeing in the plots
dat2$order2 <- factor(dat2$order2, levels=rev(c("o__Peronosporales"  , "o__Pythiales"   ,"o__Albuginales" , "o__Chlamydomonadales" , "Other")))  #change the sorts of order based on your sort to see in plots.

col = c("#8C8CE9","#F9DA73" ,"#BD4B4B","#999999") # list of color code


y1= dat2
yy1 = aggregate(y1$Abundance, by=list(y1$order2 ,y1$Experiment , y1$Month , y1$Sample), FUN=sum)
colnames(yy1) = c("order2","Year" , "Month" , "Sample1" , "Abundance")
y1mvals = data.frame(yy1[ yy1$order == "o__Peronosporales", ])
y1mxvals = y1mvals[with(y1mvals, order(Abundance)), ]$Sample
yy1$Sample <- factor(yy1$Sample, levels = y1mxvals)
year1 = yy1

ploty3 <- ggplot()  + geom_bar(data=year1, aes(x=Sample, y=Abundance, fill = order2), stat="identity", position="fill" , width =1) +
  facet_wrap(.~ Month, scales = 'free_x' , nrow = 1 , strip.position="bottom") + 
  theme_minimal_grid() +
  theme(text = element_text(size=12 ),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        panel.spacing.x=unit(0, "lines"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values=rev(col))  +
  scale_y_continuous(labels =scale,limits = c(0,1) )

ggsave("fig1.pdf",grid.arrange(ploty1,ploty2,ploty3 , nrow = 3) , width = 7 , height = 8 , units = "in")











