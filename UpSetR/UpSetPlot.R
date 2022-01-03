# Input: m6A site bed file
# Output: UpSetPlot
# Author: Hongjiang Liu
# Email: hongjiang.liu@ucsf.edu

#### Load necessary Packages ####
library('dplyr')
library(ggplot2)
library("UpSetR")
library(progress)

#import data
print("Loading...")
args <- commandArgs(trailingOnly=TRUE)
GRIP_input <- args[1]
CIMS_input <- args[2]
CITS_input <- args[3]
m6Am_input <- args[4]
all_input <- args[5]
outputfile <- args[6]
outputpdf <- args[7]
#GRIP <- read.table("bed/all.GRIP.merged.bed",header=F)
#CIMS <- read.table("bed/CIMS-miCLIP.bed",header=F)
#CITS <- read.table("bed/CITS-miCLIP.bed",header=F)
#m6Am <- read.table("bed/miCLIP-m6Am.bed",header=F)
#all <- read.table("bed/all.merged.bed",header=F)
#outputfile <- "all.csv"
#outputpdf <- "upsetplot.pdf"

GRIP <- read.table(GRIP_input,header=F)
CIMS <- read.table(CIMS_input,header=F)
CITS <- read.table(CITS_input,header=F)
m6Am <- read.table(m6Am_input,header=F)
all <- read.table(all_input,header=F)
all[,4:7]=0
colnames(all) <- c("chr","start","end","GRIP","CIMS","CITS","m6Am")

xrow=1
p1 <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = nrow(all))
print("Making one-hot table...")
for(i in xrow:nrow(all)){
  p1$tick()
  f_GRIP = filter(GRIP, V1 == all[i,1])
  f_CIMS = filter(CIMS, V1 == all[i,1])
  f_CITS = filter(CITS, V1 == all[i,1])
  f_m6Am = filter(m6Am, V1 == all[i,1])
  
  for(j in 1:nrow(f_GRIP)){
    if(f_GRIP[j,2]>=all[i,2] & f_GRIP[j,3]<=all[i,3]){
      all[i,4]=1
      break
    }
    else{
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
    }
  }
  if(nrow(f_m6Am) > 0){
    for(j in 1:nrow(f_m6Am)){
      if(f_m6Am[j,2]>=all[i,2] & f_m6Am[j,3]<=all[i,3]){
        all[i,7]=1
        break
      }
    }
  }
}

p2 <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = nrow(all))
print("Converting...")
for(k in 1:nrow(all)){
  p2$tick()
  if(all[k,1] == "chr1"){
    all[k,1] <- 1
  }
  if(all[k,1] == "chr2"){
    all[k,1] <- 2
  }
  if(all[k,1] == "chr3"){
    all[k,1] <- 3
  }
  if(all[k,1] == "chr4"){
    all[k,1] <- 4
  }
  if(all[k,1] == "chr5"){
    all[k,1] <- 5
  }
  if(all[k,1] == "chr6"){
    all[k,1] <- 6
  }
  if(all[k,1] == "chr7"){
    all[k,1] <- 7
  }
  if(all[k,1] == "chr8"){
    all[k,1] <- 8
  }
  if(all[k,1] == "chr9"){
    all[k,1] <- 9
  }
  if(all[k,1] == "chr10"){
    all[k,1] <- 10
  }
  if(all[k,1] == "chr11"){
    all[k,1] <- 11
  }
  if(all[k,1] == "chr12"){
    all[k,1] <- 12
  }
  if(all[k,1] == "chr13"){
    all[k,1] <- 13
  }
  if(all[k,1] == "chr14"){
    all[k,1] <- 14
  }
  if(all[k,1] == "chr15"){
    all[k,1] <- 15
  }
  if(all[k,1] == "chr16"){
    all[k,1] <- 16
  }
  if(all[k,1] == "chr17"){
    all[k,1] <- 17
  }
  if(all[k,1] == "chr18"){
    all[k,1] <- 18
  }
  if(all[k,1] == "chr19"){
    all[k,1] <- 19
  }
  if(all[k,1] == "chr20"){
    all[k,1] <- 20
  }
  if(all[k,1] == "chr21"){
    all[k,1] <- 21
  }
  if(all[k,1] == "chr22"){
    all[k,1] <- 22
  }
  if(all[k,1] == "chrX"){
    all[k,1] <- 23
  }
  if(all[k,1] == "chrY"){
    all[k,1] <- 24
  }
  if(all[k,1] == "chrM"){
    all[k,1] <- 25
  }
  
}
print("Writing table...")
write.table(all,file=outputfile,sep="\t",row.names = F,quote = F)

print("Making UpSet Plot...")
colnames(all) <- c("chr","start","end","GRIP-seq","CIMS-miCLIP","CITS-miCLIP","m6Am-miCLIP")

pdf(outputpdf)
upsetplot <- upset(all,sets.bar.color = "#2d3436",
      sets = c("GRIP-seq","CIMS-miCLIP","CITS-miCLIP","m6Am-miCLIP"),
      sets.x.label = "Interactions Called",
      queries = list(list(query = intersects, params = list("GRIP-seq"), color = "#6c5ce7", active = T),
                     list(query = intersects, params = list("CIMS-miCLIP"), color = "#ff7675", active = T),
                     list(query = intersects, params = list("CITS-miCLIP"), color = "#fdcb6e", active = T),
                     list(query = intersects, params = list("m6Am-miCLIP"), color = "#00b894", active = T)
      )#,
      #attribute.plots = list(gridrows = 100, plots = list(list(plot = scatter_plot, x = "chr",y = "start", queries = T)), ncols =1),
      #query.legend = "bottom"
)
print(upsetplot)
dev.off()

