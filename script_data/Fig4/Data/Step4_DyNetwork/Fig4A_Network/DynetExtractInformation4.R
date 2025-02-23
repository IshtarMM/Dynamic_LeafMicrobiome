library(tidyr)
#in this Script we can create a nodeatribute table and an edge table for create a network that have information of 5 networks that are conect to each other (the input files for this script came from DyNet that is a app in the Cytoscape)
#################################################################################################################################################
###Create a node atribute table for dynanic networks contain taxonomy information, month, core 

setwd("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/NetworkResult/ResultSparcc/Otu/Month/Pval/NetFeature1/Dynet/") 
nodeatributefile <- read.table("Sparcc_OtuNodetables.txt",header = TRUE) # node atribute table for sparcc networks output that i made it before ( contains taxonomy information and cores)
nodelist = c("11-12node.csv","12-1node.csv","1-2node.csv","2-3node.csv")  # node files , (Dynet outputs)
edgelist = c("11-12edge.csv","12-1edge.csv","1-2edge.csv","2-3edge.csv") # edge files (Dynet outupt , also contains correlation and p-value for each edge)
monthname <- c("nov","dec","jan","feb","mar")   # we need it later


#################################################################################################################################################

nodeatribute<- data.frame()  
for (i in 1:4){
  Node <- read.csv(file = nodelist[i])
  colnames(Node) <- c("SUID", "Firstnodefile" ,"Secondnodefile" ,"DnScoredegreeCorrected","RewiringDnScore","EdgeCount","Name","Selected","SharedName")
  
  ###FirstMonth 
  FirstMonth <- Node[Node$Firstnodefile=="true",]  #select Firstmonth node 
  FirstMonth1 <- FirstMonth[,c("Name", "SharedName")] 
  FirstMonthB <- merge(FirstMonth1,nodeatributefile,by.x="Name", by.y = "OtuName") # merge info from node atribute table that i made it before
  FirstMonthB$Month <- monthname[i] #add info of month to dataframe (this nodes are from which month)
  FirstMonthB1 <- data.frame(NodeName = paste(monthname[i],"_",FirstMonthB$Name,sep = "") , FirstMonthB) #add an string to bigginig of each node to show the monthnem for ex:jan_Bacteriaotu1
  
  ###SecondMonth
  SecondMonth <- Node[Node$Secondnodefile=="true",]  #select Secondmonth node 
  SecondMonth1 <- SecondMonth[,c("Name", "SharedName")] 
  SecondMonthB <- merge(SecondMonth1,nodeatributefile,by.x="Name", by.y = "OtuName") # merge info from node atribute table that i made it before
  SecondMonthB$Month <- monthname[i+1] #add info of month to dataframe (this nodes are from which month)
  SecondMonthB1 <- data.frame(NodeName = paste(monthname[i+1],"_",SecondMonthB$Name,sep = "") , SecondMonthB) #add an string to bigginig of each node to show the monthnem for ex:jan_Bacteriaotu1
  
  nodeatribute <- rbind(nodeatribute,FirstMonthB1,SecondMonthB1)
}
finalnodeinfo = nodeatribute[!duplicated(nodeatribute$NodeName), ]
write.csv(finalnodeinfo,"newedgenode/nodeatributev.csv",quote = FALSE)

#################################################################################################################################################

#Create edge table contain correlation information and also psudo edge that shows connection between networks.
#we need to merge edges from different months and also define some psudo edges that are showing conection between two consecutive month 
#################################################################################################################################################

setwd("/Users/mahmoudi/Documents/Doc/2019/NewData_Octuber/NetworkResult/ResultSparcc/Otu/Month/Pval/NetFeature1/Dynet/") 
nodelist = c("11-12node.csv","12-1node.csv","1-2node.csv","2-3node.csv")  # node files , (Dynet outputs)
edgelist = c("11-12edge.csv","12-1edge.csv","1-2edge.csv","2-3edge.csv") # edge files (Dynet outupt , also contains correlation and p-value for each edge)
monthname <- c("nov","dec","jan","feb","mar")   # we need it later
counter <- 0
finaledge <- data.frame()
for (i in 1:4){
  i = 2
  counter <- counter + 1
  Node <- read.csv(file = nodelist[i])
  colnames(Node) <- c("SUID", "Firstnodefile" ,"Secondnodefile" ,"DnScoredegreeCorrected","RewiringDnScore","EdgeCount","Name","Selected","SharedName")
  overlap <- Node[(Node$Secondnodefile=="true" & Node$Firstnodefile == "true"),]
  overlap1 <- overlap[,c("Name", "SharedName")]
  colnames(overlap1) <- c("Node1","Node2")
  overlap1$Node1 <-paste(monthname[i],overlap1$Node1,sep = "_")
  overlap1$Node2 <-paste(monthname[i+1],overlap1$Node2,sep = "_")
  overlap1$binarycorrelation <- "linkeredge"
  overlap1$Correlation <- "no_info"
  ov <- overlap1
  ov$Stablenode <- "yes"
  
  edge <- read.csv(file = edgelist[i])
  #print(dim(edge))
  edge <- edge[,c(2,4,5,6,8,9,10,15)]
  colnames(edge) <- c("Corr_firstmonth","FirstedgefilePresentAbsence" ,"FirstmonthPval" ,
                      "Corr_secondmonth" , "SecondedgefilePresentAbsence" , "DyNetPairwiseComparison","SecondmonthPval" , "sharedname")
  #edgeoverlap <- edge[(edge$Firstedgefilebolian=="true" & edge$secondedgefilebolian == "true"),]
  
  ####firstmonth edge
  Firstmonthnedge <- edge[edge$FirstedgefilePresentAbsence=="true",]
  Firstmonthnedge1 <- separate(Firstmonthnedge,sharedname, c("Node1" , "Node2"), sep = "(interacts with)", remove = FALSE)
  Firstmonthnedge2 <- Firstmonthnedge1[,c(1,9,10)]
  Firstmonthnedge2[] <- lapply(Firstmonthnedge2, gsub, pattern='[( ]', replacement='') #remove "Space("
  Firstmonthnedge2[] <- lapply(Firstmonthnedge2, gsub, pattern='[) ]', replacement='')  #remove "Space)"
  Firstmonthnedge2$Node1 <-paste(monthname[i],Firstmonthnedge2$Node1,sep = "_") #add the name of month to begginig of each node
  Firstmonthnedge2$Node2 <-paste(monthname[i],Firstmonthnedge2$Node2,sep = "_")  #add the name of month to begginig of each node
  positivecor1 <- Firstmonthnedge2[Firstmonthnedge2$Corr_firstmonth>=0,]
  #print(dim(positivecor1))
  positivecor1$binarycorrelation <- "Positive"
  negativecore1 <- Firstmonthnedge2[Firstmonthnedge2$Corr_firstmonth<0,]
  #print(dim(negativecore1))
  negativecore1$binarycorrelation <- "Negative"
  FM <- rbind(positivecor1,negativecore1)
  colnames(FM) <- c("Correlation","Node1","Node2","binarycorrelation")
  FM1 <- FM[,c(2,3,4,1)]
  ####second month edges
  Secondmonthedge <- edge[edge$SecondedgefilePresentAbsence=="true",]
  Secondmonthedge1 <- separate(Secondmonthedge,sharedname, c("Node1" , "Node2"), sep = "(interacts with)", remove = FALSE)
  Secondmonthedge2 <- Secondmonthedge1[,c(4,9,10)]
  Secondmonthedge2[] <- lapply(Secondmonthedge2, gsub, pattern='[( ]', replacement='')
  Secondmonthedge2[] <- lapply(Secondmonthedge2, gsub, pattern='[) ]', replacement='')
  Secondmonthedge2$Node1 = paste(monthname[i+1],Secondmonthedge2$Node1,sep = "_")
  Secondmonthedge2$Node2 = paste(monthname[i+1],Secondmonthedge2$Node2,sep = "_")
  positivecor2 <- Secondmonthedge2[Secondmonthedge2$Corr_secondmonth>=0,]
  positivecor2$binarycorrelation <- "Positive"
  negativecore2 <- Secondmonthedge2[Secondmonthedge2$Corr_secondmonth<0,]
  negativecore2$binarycorrelation <- "Negative"
  print(dim(positivecor2))
  print(dim(negativecore2))
  SM <- rbind(positivecor2,negativecore2)
  colnames(SM) <- c("Correlation","Node1","Node2","binarycorrelation")
  SM1 <- SM[,c(2,3,4,1)]
  if (counter == 1){
    finaledge = rbind(finaledge,overlap1,FM1,SM1) 
  }else{
    finaledge = rbind(finaledge,overlap1,SM1) 
     }
  }
 
write.csv(finaledge,"newedgenode/edge.csv",quote = FALSE)

x = finaledge[finaledge$binarycorrelation=="linkeredge",]
x1 = finaledge[finaledge$binarycorrelation=="Positive",]







