#Script that implements a multi-K sharing (MKS) analysis, as described in Wootton, L.M., Forest, F., Verboom, G.A., 2022. Consilience across
#multiple, independent genomic data sets reveals species in a complex with limited phenotypic variation.

#All data required to run this script can be found here:
#In section 1, data files are converted from CLUMPAK (http://clumpak.tau.ac.il/) format to a format that can be used in the MKS analysis.
#In section 2, the MKS analysis is run using the Rehmannii Clade as an example. 

#### Section 1: Data formatting ####
####################################

#Load libraries
library(pophelper)
library(tidyverse)
library(vegan)

#Converting CLUMPAK output files to R friendly format

calc.K<-list()
for (k in 2:10) { #CHANGE if using different number of Ks 
  setwd("/Path/to/Clumpak/outputs")
  calc.K[[k]]<-read.table(paste("./K=",k,"/MajorCluster/CLUMPP.files/ClumppIndFile.output",sep=""))
  rows<-paste(calc.K[[k]]$V1,calc.K[[k]]$V5,sep="")
  calc.K[[k]]<-cbind(rows,calc.K[[k]][,c(seq(from=6,len=k),4)])
}

#write files out
setwd("/Path/to/Clumpak/outputs/formatted_files")
for(k in 2:10){
  write.table(calc.K[[k]],row.names=F,col.names=F,quote=F,
              file=paste("calc.clumpp.merged.K",k,".txt",sep=""))
}

#Make a single column data frame of the population names of your samples. 
#Make sure that they are in the same order as the populations in your CLUMPAK files 

pop_names <- read.table("popmap_rehmannii_reordered", stringsAsFactors=F)

pop_names <- separate(pop_names, 
                      col = V1, 
                      into = c("letter", "Pop"),
                      sep = "A")
pop_names <- data.frame(Pop = pop_names[,2])


#Read formatted CLUMPAK files into R

calc.clumpp.list<-list.files(path="/Path/to/Clumpak/outputs/formatted_files")

#Convert to Q file 

qlist.calc <- readQ(files=calc.clumpp.list[1:9]) #CHANGE to relevant number of K

#Join population names to Q files

gen_dat <- list()
for (k in 1:9) { #CHANGE to relevant number of K
  gen_dat[[k]] <- cbind(pop_names, qlist.calc[[k]])  
}

#Now all the files should be in a usable format for further analysis

#### Section 2: MKS analysis ####
#################################

#Note: As the number of potential gene pools increases with each level of K, implementing a single loop to cycle through all Ks isn't straight forward.
#Therefore, there is a lot of repetition in the code below. When running it for your study, you can increase or decrease the number of repetitions as needed. 

## Finding average population gene pool assignments ##

#Individuals within a population are assigned to gene pools with different probabilities. Here, we find the average probability of gene pool assignment across each population for each K.  

gen_dat_reduced <- list()

#K2

gen_dat_reduced[[1]] <- gen_dat[[2]] %>% 
  group_by(Pop) %>% 
  mutate(C1 = mean(Cluster1),
         C2 = mean(Cluster2)) %>% 
  as.data.frame()

#K3

gen_dat_reduced[[2]] <- gen_dat[[3]] %>% 
  group_by(Pop) %>% 
  mutate(C1 = mean(Cluster1),
         C2 = mean(Cluster2),
         C3 = mean(Cluster3)) %>% 
  as.data.frame()

#K4

gen_dat_reduced[[3]] <- gen_dat[[4]] %>% 
  group_by(Pop) %>% 
  mutate(C1 = mean(Cluster1),
         C2 = mean(Cluster2),
         C3 = mean(Cluster3),
         C4 = mean(Cluster4)) %>% 
  as.data.frame()

#K5

gen_dat_reduced[[4]] <- gen_dat[[5]] %>% 
  group_by(Pop) %>% 
  mutate(C1 = mean(Cluster1),
         C2 = mean(Cluster2),
         C3 = mean(Cluster3),
         C4 = mean(Cluster4),
         C5 = mean(Cluster5)) %>% 
  as.data.frame()

#K6

gen_dat_reduced[[5]] <- gen_dat[[6]] %>% 
  group_by(Pop) %>% 
  mutate(C1 = mean(Cluster1),
         C2 = mean(Cluster2),
         C3 = mean(Cluster3),
         C4 = mean(Cluster4),
         C5 = mean(Cluster5),
         C6 = mean(Cluster6)) %>% 
  as.data.frame()

#K7

gen_dat_reduced[[6]] <- gen_dat[[7]] %>% 
  group_by(Pop) %>% 
  mutate(C1 = mean(Cluster1),
         C2 = mean(Cluster2),
         C3 = mean(Cluster3),
         C4 = mean(Cluster4),
         C5 = mean(Cluster5),
         C6 = mean(Cluster6),
         C7 = mean(Cluster7)) %>% 
  as.data.frame()

#K8

gen_dat_reduced[[7]] <- gen_dat[[8]] %>% 
  group_by(Pop) %>% 
  mutate(C1 = mean(Cluster1),
         C2 = mean(Cluster2),
         C3 = mean(Cluster3),
         C4 = mean(Cluster4),
         C5 = mean(Cluster5),
         C6 = mean(Cluster6),
         C7 = mean(Cluster7),
         C8 = mean(Cluster8)) %>% 
  as.data.frame()

#K9

gen_dat_reduced[[8]] <- gen_dat[[9]] %>% 
  group_by(Pop) %>% 
  mutate(C1 = mean(Cluster1),
         C2 = mean(Cluster2),
         C3 = mean(Cluster3),
         C4 = mean(Cluster4),
         C5 = mean(Cluster5),
         C6 = mean(Cluster6),
         C7 = mean(Cluster7),
         C8 = mean(Cluster8),
         C9 = mean(Cluster9)) %>% 
  as.data.frame()

#K10

gen_dat_reduced[[9]] <- gen_dat[[1]] %>% 
  group_by(Pop) %>% 
  mutate(C1 = mean(Cluster1),
         C2 = mean(Cluster2),
         C3 = mean(Cluster3),
         C4 = mean(Cluster4),
         C5 = mean(Cluster5),
         C6 = mean(Cluster6),
         C7 = mean(Cluster7),
         C8 = mean(Cluster8),
         C9 = mean(Cluster9),
         C10 = mean(Cluster10)) %>% 
  as.data.frame()

#Filter relevant columns and reorder rows 

distinct_gen_dat <- list()
for (k in 1:9) { #CHANGE to relevant number of K
  start <- k+3
  finish <- k+k+3
  distinct_gen_dat[[k]] <- gen_dat_reduced[[k]][,c(1, start:finish)] 
  distinct_gen_dat[[k]] <- distinct_gen_dat[[k]] %>% distinct()
  
  row.names(distinct_gen_dat[[k]]) <- distinct_gen_dat[[k]][,1] #rename rows
  distinct_gen_dat[[k]] <- distinct_gen_dat[[k]][,-1] #remove pop column
  
  #Reordered based on order of phylogeny
  distinct_gen_dat[[k]] <- distinct_gen_dat[[k]][c(10,11,12, 5,9,7,8,6, 1,4,3, 2),]
}

### MKS analysis ###

# Code genepools as present (1) or absent (0) in a population based on predefined cut off 

cutoff <- 0.25 

for (k in 1:9) { #CHANGE to relevant number of K
  
  distinct_gen_dat[[k]][distinct_gen_dat[[k]] >= cutoff] <- 1 
  distinct_gen_dat[[k]][distinct_gen_dat[[k]] < cutoff] <- 0 
}


#Function to evaluate whether genepools are shared between populations

compare_values <- function(x, y) {
  ifelse( x == 1 & y == 1, 1, 0)
}

#Create binary presence absence matrix for each gene pool in a given value of K

#K2 ####

k2matC1 <- matrix(nrow = 12, ncol = 12) #Nrow and ncol should equal the number of samples in your study
k2matC2 <- matrix(nrow = 12, ncol = 12)


for (i in 1:12) {
  for (j in 1:12) {
    value1 <- distinct_gen_dat[[1]][i , 1]
    value2 <- distinct_gen_dat[[1]][i , 2]
    
    k2matC1[i,j] <- compare_values(value1, distinct_gen_dat[[1]][j , 1]) 
    k2matC2[i,j] <- compare_values(value2, distinct_gen_dat[[1]][j , 2]) 
    
    
  }
}

#Add genepool matrices together
k2mat_full <- k2matC1 + 
  k2matC2

#Convert back to binary (within a K, we are not interested in how many genepools are shared across populations, just that there is at least one)
k2mat_full[k2mat_full > 1] <- 1

# K3 ####
k3matC1 <- matrix(nrow = 12, ncol = 12)
k3matC2 <- matrix(nrow = 12, ncol = 12)
k3matC3 <- matrix(nrow = 12, ncol = 12)


for (i in 1:12) {
  for (j in 1:12) {
    value1 <- distinct_gen_dat[[2]][i , 1]
    value2 <- distinct_gen_dat[[2]][i , 2]
    value3 <- distinct_gen_dat[[2]][i , 3]
    
    k3matC1[i,j] <- compare_values(value1, distinct_gen_dat[[2]][j , 1]) 
    k3matC2[i,j] <- compare_values(value2, distinct_gen_dat[[2]][j , 2]) 
    k3matC3[i,j] <- compare_values(value3, distinct_gen_dat[[2]][j , 3]) 
    
    
  }
}

k3mat_full <- k3matC1 + 
  k3matC2 +
  k3matC3 

k3mat_full[k3mat_full > 1] <- 1

# K4 ####
k4matC1 <- matrix(nrow = 12, ncol = 12)
k4matC2 <- matrix(nrow = 12, ncol = 12)
k4matC3 <- matrix(nrow = 12, ncol = 12)
k4matC4 <- matrix(nrow = 12, ncol = 12)


for (i in 1:12) {
  for (j in 1:12) {
    value1 <- distinct_gen_dat[[3]][i , 1]
    value2 <- distinct_gen_dat[[3]][i , 2]
    value3 <- distinct_gen_dat[[3]][i , 3]
    value4 <- distinct_gen_dat[[3]][i , 4]
    
    k4matC1[i,j] <- compare_values(value1, distinct_gen_dat[[3]][j , 1]) 
    k4matC2[i,j] <- compare_values(value2, distinct_gen_dat[[3]][j , 2]) 
    k4matC3[i,j] <- compare_values(value3, distinct_gen_dat[[3]][j , 3]) 
    k4matC4[i,j] <- compare_values(value4, distinct_gen_dat[[3]][j , 4])
    
    
  }
}

k4mat_full <- k4matC1 + 
  k4matC2 +
  k4matC3 +
  k4matC4 

k4mat_full[k4mat_full > 1] <- 1


# K5 ####
k5matC1 <- matrix(nrow = 12, ncol = 12)
k5matC2 <- matrix(nrow = 12, ncol = 12)
k5matC3 <- matrix(nrow = 12, ncol = 12)
k5matC4 <- matrix(nrow = 12, ncol = 12)
k5matC5 <- matrix(nrow = 12, ncol = 12)


for (i in 1:12) {
  for (j in 1:12) {
    value1 <- distinct_gen_dat[[4]][i , 1]
    value2 <- distinct_gen_dat[[4]][i , 2]
    value3 <- distinct_gen_dat[[4]][i , 3]
    value4 <- distinct_gen_dat[[4]][i , 4]
    value5 <- distinct_gen_dat[[4]][i , 5]
    
    k5matC1[i,j] <- compare_values(value1, distinct_gen_dat[[4]][j , 1]) 
    k5matC2[i,j] <- compare_values(value2, distinct_gen_dat[[4]][j , 2]) 
    k5matC3[i,j] <- compare_values(value3, distinct_gen_dat[[4]][j , 3]) 
    k5matC4[i,j] <- compare_values(value4, distinct_gen_dat[[4]][j , 4]) 
    k5matC5[i,j] <- compare_values(value5, distinct_gen_dat[[4]][j , 5])
    
    
  }
}

k5mat_full <- k5matC1 + 
  k5matC2 +
  k5matC3 +
  k5matC4 +
  k5matC5  

k5mat_full[k5mat_full > 1] <- 1

# K6 ####
k6matC1 <- matrix(nrow = 12, ncol = 12)
k6matC2 <- matrix(nrow = 12, ncol = 12)
k6matC3 <- matrix(nrow = 12, ncol = 12)
k6matC4 <- matrix(nrow = 12, ncol = 12)
k6matC5 <- matrix(nrow = 12, ncol = 12)
k6matC6 <- matrix(nrow = 12, ncol = 12)


for (i in 1:12) {
  for (j in 1:12) {
    value1 <- distinct_gen_dat[[5]][i , 1]
    value2 <- distinct_gen_dat[[5]][i , 2]
    value3 <- distinct_gen_dat[[5]][i , 3]
    value4 <- distinct_gen_dat[[5]][i , 4]
    value5 <- distinct_gen_dat[[5]][i , 5]
    value6 <- distinct_gen_dat[[5]][i , 6]
    
    k6matC1[i,j] <- compare_values(value1, distinct_gen_dat[[5]][j , 1]) 
    k6matC2[i,j] <- compare_values(value2, distinct_gen_dat[[5]][j , 2]) 
    k6matC3[i,j] <- compare_values(value3, distinct_gen_dat[[5]][j , 3]) 
    k6matC4[i,j] <- compare_values(value4, distinct_gen_dat[[5]][j , 4]) 
    k6matC5[i,j] <- compare_values(value5, distinct_gen_dat[[5]][j , 5]) 
    k6matC6[i,j] <- compare_values(value6, distinct_gen_dat[[5]][j , 6])
    
    
  }
}

k6mat_full <- k6matC1 + 
  k6matC2 +
  k6matC3 +
  k6matC4 +
  k6matC5 +
  k6matC6  

k6mat_full[k6mat_full > 1] <- 1

# K7 ####
k7matC1 <- matrix(nrow = 12, ncol = 12)
k7matC2 <- matrix(nrow = 12, ncol = 12)
k7matC3 <- matrix(nrow = 12, ncol = 12)
k7matC4 <- matrix(nrow = 12, ncol = 12)
k7matC5 <- matrix(nrow = 12, ncol = 12)
k7matC6 <- matrix(nrow = 12, ncol = 12)
k7matC7 <- matrix(nrow = 12, ncol = 12)


for (i in 1:12) {
  for (j in 1:12) {
    value1 <- distinct_gen_dat[[6]][i , 1]
    value2 <- distinct_gen_dat[[6]][i , 2]
    value3 <- distinct_gen_dat[[6]][i , 3]
    value4 <- distinct_gen_dat[[6]][i , 4]
    value5 <- distinct_gen_dat[[6]][i , 5]
    value6 <- distinct_gen_dat[[6]][i , 6]
    value7 <- distinct_gen_dat[[6]][i , 7]
    
    k7matC1[i,j] <- compare_values(value1, distinct_gen_dat[[6]][j , 1]) 
    k7matC2[i,j] <- compare_values(value2, distinct_gen_dat[[6]][j , 2]) 
    k7matC3[i,j] <- compare_values(value3, distinct_gen_dat[[6]][j , 3]) 
    k7matC4[i,j] <- compare_values(value4, distinct_gen_dat[[6]][j , 4]) 
    k7matC5[i,j] <- compare_values(value5, distinct_gen_dat[[6]][j , 5]) 
    k7matC6[i,j] <- compare_values(value6, distinct_gen_dat[[6]][j , 6]) 
    k7matC7[i,j] <- compare_values(value7, distinct_gen_dat[[6]][j , 7]) 
    
    
  }
}

k7mat_full <- k7matC1 + 
  k7matC2 +
  k7matC3 +
  k7matC4 +
  k7matC5 +
  k7matC6 +
  k7matC7  

k7mat_full[k7mat_full > 1] <- 1

# K8 ####
k8matC1 <- matrix(nrow = 12, ncol = 12)
k8matC2 <- matrix(nrow = 12, ncol = 12)
k8matC3 <- matrix(nrow = 12, ncol = 12)
k8matC4 <- matrix(nrow = 12, ncol = 12)
k8matC5 <- matrix(nrow = 12, ncol = 12)
k8matC6 <- matrix(nrow = 12, ncol = 12)
k8matC7 <- matrix(nrow = 12, ncol = 12)
k8matC8 <- matrix(nrow = 12, ncol = 12)


for (i in 1:12) {
  for (j in 1:12) {
    value1 <- distinct_gen_dat[[7]][i , 1]
    value2 <- distinct_gen_dat[[7]][i , 2]
    value3 <- distinct_gen_dat[[7]][i , 3]
    value4 <- distinct_gen_dat[[7]][i , 4]
    value5 <- distinct_gen_dat[[7]][i , 5]
    value6 <- distinct_gen_dat[[7]][i , 6]
    value7 <- distinct_gen_dat[[7]][i , 7]
    value8 <- distinct_gen_dat[[7]][i , 8]
    
    k8matC1[i,j] <- compare_values(value1, distinct_gen_dat[[7]][j , 1]) 
    k8matC2[i,j] <- compare_values(value2, distinct_gen_dat[[7]][j , 2]) 
    k8matC3[i,j] <- compare_values(value3, distinct_gen_dat[[7]][j , 3]) 
    k8matC4[i,j] <- compare_values(value4, distinct_gen_dat[[7]][j , 4]) 
    k8matC5[i,j] <- compare_values(value5, distinct_gen_dat[[7]][j , 5]) 
    k8matC6[i,j] <- compare_values(value6, distinct_gen_dat[[7]][j , 6]) 
    k8matC7[i,j] <- compare_values(value7, distinct_gen_dat[[7]][j , 7]) 
    k8matC8[i,j] <- compare_values(value8, distinct_gen_dat[[7]][j , 8]) 
    
    
  }
}

k8mat_full <- k8matC1 + 
  k8matC2 +
  k8matC3 +
  k8matC4 +
  k8matC5 +
  k8matC6 +
  k8matC7 +
  k8matC8 

k8mat_full[k8mat_full > 1] <- 1

# K9 ####
k9matC1 <- matrix(nrow = 12, ncol = 12)
k9matC2 <- matrix(nrow = 12, ncol = 12)
k9matC3 <- matrix(nrow = 12, ncol = 12)
k9matC4 <- matrix(nrow = 12, ncol = 12)
k9matC5 <- matrix(nrow = 12, ncol = 12)
k9matC6 <- matrix(nrow = 12, ncol = 12)
k9matC7 <- matrix(nrow = 12, ncol = 12)
k9matC8 <- matrix(nrow = 12, ncol = 12)
k9matC9 <- matrix(nrow = 12, ncol = 12)


for (i in 1:12) {
  for (j in 1:12) {
    value1 <- distinct_gen_dat[[8]][i , 1]
    value2 <- distinct_gen_dat[[8]][i , 2]
    value3 <- distinct_gen_dat[[8]][i , 3]
    value4 <- distinct_gen_dat[[8]][i , 4]
    value5 <- distinct_gen_dat[[8]][i , 5]
    value6 <- distinct_gen_dat[[8]][i , 6]
    value7 <- distinct_gen_dat[[8]][i , 7]
    value8 <- distinct_gen_dat[[8]][i , 8]
    value9 <- distinct_gen_dat[[8]][i , 9]
    
    k9matC1[i,j] <- compare_values(value1, distinct_gen_dat[[8]][j , 1]) 
    k9matC2[i,j] <- compare_values(value2, distinct_gen_dat[[8]][j , 2]) 
    k9matC3[i,j] <- compare_values(value3, distinct_gen_dat[[8]][j , 3]) 
    k9matC4[i,j] <- compare_values(value4, distinct_gen_dat[[8]][j , 4]) 
    k9matC5[i,j] <- compare_values(value5, distinct_gen_dat[[8]][j , 5]) 
    k9matC6[i,j] <- compare_values(value6, distinct_gen_dat[[8]][j , 6]) 
    k9matC7[i,j] <- compare_values(value7, distinct_gen_dat[[8]][j , 7]) 
    k9matC8[i,j] <- compare_values(value8, distinct_gen_dat[[8]][j , 8]) 
    k9matC9[i,j] <- compare_values(value9, distinct_gen_dat[[8]][j , 9])
    
    
  }
}

k9mat_full <- k9matC1 + 
  k9matC2 +
  k9matC3 +
  k9matC4 +
  k9matC5 +
  k9matC6 +
  k9matC7 +
  k9matC8 +
  k9matC9 

k9mat_full[k9mat_full > 1] <- 1


# K10 ####
k10matC1 <- matrix(nrow = 12, ncol = 12)
k10matC2 <- matrix(nrow = 12, ncol = 12)
k10matC3 <- matrix(nrow = 12, ncol = 12)
k10matC4 <- matrix(nrow = 12, ncol = 12)
k10matC5 <- matrix(nrow = 12, ncol = 12)
k10matC6 <- matrix(nrow = 12, ncol = 12)
k10matC7 <- matrix(nrow = 12, ncol = 12)
k10matC8 <- matrix(nrow = 12, ncol = 12)
k10matC9 <- matrix(nrow = 12, ncol = 12)
k10matC10 <- matrix(nrow = 12, ncol = 12)

for (i in 1:12) {
  for (j in 1:12) {
    value <- distinct_gen_dat[[9]][i , 1]
    value2 <- distinct_gen_dat[[9]][i , 2]
    value3 <- distinct_gen_dat[[9]][i , 3]
    value4 <- distinct_gen_dat[[9]][i , 4]
    value5 <- distinct_gen_dat[[9]][i , 5]
    value6 <- distinct_gen_dat[[9]][i , 6]
    value7 <- distinct_gen_dat[[9]][i , 7]
    value8 <- distinct_gen_dat[[9]][i , 8]
    value9 <- distinct_gen_dat[[9]][i , 9]
    value10 <- distinct_gen_dat[[9]][i ,10]
    
    k10matC1[i,j] <- compare_values(value, distinct_gen_dat[[9]][j , 1]) 
    k10matC2[i,j] <- compare_values(value2, distinct_gen_dat[[9]][j , 2]) 
    k10matC3[i,j] <- compare_values(value3, distinct_gen_dat[[9]][j , 3]) 
    k10matC4[i,j] <- compare_values(value4, distinct_gen_dat[[9]][j , 4]) 
    k10matC5[i,j] <- compare_values(value5, distinct_gen_dat[[9]][j ,5]) 
    k10matC6[i,j] <- compare_values(value6, distinct_gen_dat[[9]][j , 6]) 
    k10matC7[i,j] <- compare_values(value7, distinct_gen_dat[[9]][j , 7]) 
    k10matC8[i,j] <- compare_values(value8, distinct_gen_dat[[9]][j , 8]) 
    k10matC9[i,j] <- compare_values(value9, distinct_gen_dat[[9]][j , 9]) 
    k10matC10[i,j] <- compare_values(value10, distinct_gen_dat[[9]][j , 10])
    
    
  }
}

k10mat_full <- k10matC1 + 
  k10matC2 +
  k10matC3 +
  k10matC4 +
  k10matC5 +
  k10matC6 +
  k10matC7 +
  k10matC8 +
  k10matC9 +
  k10matC10

k10mat_full[k10mat_full > 1] <- 1


#Final step: Add all matrices together to get mulki-K sharing values

all_rehmannii <- k2mat_full + 
                 k3mat_full + 
                 k4mat_full + 
                 k5mat_full + 
                 k6mat_full + 
                 k7mat_full + 
                 k8mat_full + 
                 k9mat_full + 
                 k10mat_full

#Use this final matrix to make heatmaps, e.g. with ggtree::gheatmap()

rm(calc.K, 
   distinct_gen_dat,
   gen_dat,
   gen_dat_reduced,
   pop_names,
   qlist.calc,
   calc.clumpp.list,
   cutoff,
   finish,
   rows,
   start,
   i,j,k)
rm(list=ls(pattern="matC"))
rm(list=ls(pattern="full"))
rm(list=ls(pattern="value"))
