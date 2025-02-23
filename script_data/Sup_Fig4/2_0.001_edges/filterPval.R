listfiles = list.files(pattern = "*.csv" , full.names = T)
for (file in listfiles){
  data1<-read.csv(file, sep=",", dec = ".", stringsAsFactors = F, header=T, row.names = NULL)
  data = data1[which(data1$pval<=0.001),]
  seprateInteraction= data%>%separate("Label" , c("N1", "N2") , "->")
  seprateInteractionNode1= seprateInteraction%>%separate("N1" , c("k1", "p1" , "c1" , "o1" , "f1" , "Node1") , ";")
  seprateInteractionNode2= seprateInteractionNode1%>%separate("N2" , c("k2", "p2" , "c2" , "o2" , "f2" , "Node2") , ";")
  finalfile = as.data.frame(cbind(seprateInteractionNode2[,c(10,16,18)]))
  finalfile= as.data.frame(sapply(finalfile, function(x) gsub("\"", "", x)))
  write.table(finalfile, paste("0.001edge/" , basename(file)), sep="\t", quote = F, row.names=FALSE)
  print(basename(file))
} 
