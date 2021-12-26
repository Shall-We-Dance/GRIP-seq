# Input: m6A site bed file
# Output: UpSetPlot
# Author: Hongjiang Liu
# Email: hongjiang.liu@ucsf.edu

#### Load necessary Packages ####
library('bedtoolsr')
library('dplyr')
library(data.table)
library(ggplot2)
library("UpSetR")
library("ComplexHeatmap")
library(gprofiler2)

setwd(".")
GRIP <- read.table("all.GRIP.sorted.merged.bed",header=F)
CIMS <- read.table("CIMS-miCLIP.bed",header=F)
CITS <- read.table("CITS-miCLIP.bed",header=F)
m6Am <- read.table("miCLIP-m6Am.bed",header=F)
all <- read.table("all.sorted.merged.bed",header=F)

all[,4:7]=0
colnames(all) <- c("chr","start","end","GRIP","CIMS","CITS","m6Am")

#Initialization
xrow=1
#Continue on breakpoint
load("UpSet.Rdata")
#This loop takes a long time, over 8 hours
for(i in xrow:nrow(all)){
  print(i)
  f_GRIP = filter(GRIP, V1 == all[i,1])
  f_CIMS = filter(CIMS, V1 == all[i,1])
  f_CITS = filter(CITS, V1 == all[i,1])
  f_m6Am = filter(m6Am, V1 == all[i,1])
  
  for(j in 1:nrow(f_GRIP)){
    if(f_GRIP[j,2]>=all[i,2] & f_GRIP[j,3]<=all[i,3]){
      all[i,4]=1
      break
    }
  }
  for(j in 1:nrow(f_CIMS)){
    if(f_CIMS[j,2]>=all[i,2] & f_CIMS[j,3]<=all[i,3]){
      all[i,5]=1
      break
    }
  }
  if(nrow(f_CITS) > 0){
  for(j in 1:nrow(f_CITS)){
    if(f_CITS[j,2]>=all[i,2] & f_CITS[j,3]<=all[i,3]){
      all[i,6]=1
      break
    }
  }}
  if(nrow(f_m6Am) > 0){
  for(j in 1:nrow(f_m6Am)){
    if(f_m6Am[j,2]>=all[i,2] & f_m6Am[j,3]<=all[i,3]){
      all[i,7]=1
      break
    }
  }}
  
  if(i%%500 == 0)
  {
    xrow = i
    save(all,xrow, file="UpSet.Rdata")
    write.table(all,file="all.csv",sep="\t",row.names = F,quote = F)
  }
  if(i%%nrow(all) == 0)
  {
    xrow = i
    save(all,xrow, file="UpSet.Rdata")
    write.csv(all,file="all.csv",sep="\t",row.names = F,quote = F)
  }
}
write.table(all,file="all.csv",sep="\t",row.names = F,quote = F)

write.csv(all,file="output_onehot.csv",sep="\t",row.names = F, col.names = T,quote = F)
df = read.csv("all.csv")
colnames(df) <- c("chr","start","end","GRIP","CIMS-miCLIP","CITS-miCLIP","m6Am-miCLIP")

upset(df,sets.bar.color = "#2d3436",
      sets = c("GRIP","CIMS-miCLIP","CITS-miCLIP","m6Am-miCLIP"),
      sets.x.label = "Interactions Called",
      queries = list(list(query = intersects, params = list("GRIP"), color = "#6c5ce7", active = T),
                     list(query = intersects, params = list("CIMS-miCLIP"), color = "#ff7675", active = T),
                     list(query = intersects, params = list("CITS-miCLIP"), color = "#fdcb6e", active = T),
                     list(query = intersects, params = list("m6Am-miCLIP"), color = "#00b894", active = T)
      ),
      
      attribute.plots = list(gridrows = 100, plots = list(list(plot = histogram,x = "chr", queries = T),
                                                         list(plot = scatter_plot, y = "chr",x = "start", queries = T)), ncols =1),
      query.legend = "bottom"
)




