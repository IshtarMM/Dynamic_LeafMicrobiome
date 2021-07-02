#A script for core microbes
#defination of core: otus that have occurance(present) in more than a "threshold" in samples   
# First open your otu tabe and sample, separeate your sample based on different year (if you have) 
# Arrange each table  to this format (OTU in row and as a row name and sample in column)
#Set a "threshold" for your occurence

threshold = 0.9 #number more than 0 
experiment = 3 #number of years
Otutable = "OTUtable.txt"
#subset samples of different years here is an example for 3 years 
year1  = subset(Otutable, year=="x")  #samples of first year
year2  = subset(Otutable, year=="y")  #samples of second year 
year3  = subset(Otutable, year=="z")  #samples of third year 
#arrange your tables to otu in rows and samples in column
#Occurrence 
exp1_Occurrence = rowSums(year1 != 0)/ncol(year1) #count non zero columns for each row and devide to number of column
exp2_Occurrence = rowSums(year2 != 0)/ncol(year2)
exp3_Occurrence = rowSums(year3 != 0)/ncol(year3)

occu_file = cbind(as.data.frame(exp1_Occurrence) , as.data.frame(exp2_Occurrence) , as.data.frame(exp3_Occurrence))
core = occu_file[rowSums(occu_file >= threshold) >= experiment , ] #occurence in more than threshold  in each year





