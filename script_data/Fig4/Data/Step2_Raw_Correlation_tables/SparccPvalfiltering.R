#Create the Edge table sparcc otu level
library(tidyr)

corfilelist <- c("1_corsparcc.txt" ,"2_corsparcc.txt","3_corsparcc.txt","11_corsparcc.txt","12_corsparcc.txt")
pvafilelist <- c("1_pvals.two_sided.txt","2_pvals.two_sided.txt","3_pvals.two_sided.txt","11_pvals.two_sided.txt","12_pvals.two_sided.txt")

for (index in 1:length(corfilelist)){
  Pvals<-read.table(pvafilelist[index], sep="\t", header=T, row.names=1) #it's important to set the row names also (needed for coercion into as.matrix)
  Pvals1<-as.matrix(Pvals)
  Pvals1[upper.tri(Pvals, diag=TRUE)]<-NA
  Pvals2<-as.data.frame(as.table(Pvals1))
  colnames(Pvals2)<-c("Node1","Node2","Pvalue")
  
  Cors<-read.table(corfilelist[index], sep="\t", header=T, row.names=1) #it's important to set the row names also (needed for coercion into as.matrix)
  Cors1<-as.matrix(Cors)
  Cors1[upper.tri(Cors1, diag=TRUE)]<-NA
  Cors2<-as.data.frame(as.table(Cors1))
  colnames(Cors2)<-c("Node1","Node2","Cor")
  
  Edge_table<-cbind(Pvals2,Cors2$Cor, deparse.level=2)
  Edge_table_final<-Edge_table[!is.na(Edge_table$Pvalue),]
  colnames(Edge_table_final)<-c("Node1","Node2","Pvalue","Cor")#adjust headers   dimention = 496506 * 2 + 997
  #write.table(Edge_table_final, file=paste("Pval/EdgeTable_",basename(corfilelist[index])), quote=F, row.names=F , sep = "\t") #create the final edge table to use in cytoscape
  
  ####Pvalue Filter
  SparCC_Default <- Edge_table_final
  SparCC_Default_P0.001 <- SparCC_Default[which(SparCC_Default$Pvalue<=0.001),]
  write.table(SparCC_Default_P0.001, file=paste("0.001_",basename(corfilelist[index])), sep="\t", row.names=FALSE ,quote = F)
}

