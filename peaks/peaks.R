# Input: bed file peak form clipper
# Output: real peak for GRIP-seq
# Author: Hongjiang Liu
# Email: hongjiang.liu@ucsf.edu


#load library
library(dplyr)
library(progress)
#import data
print("Loading...")
args <- commandArgs(trailingOnly=TRUE)
raw_input <- args[1]
ref_input <- args[2]
threshhold <- args[3]
output_file_name <- args[4]

raw <- read.table(raw_input,header=F)
ref <- read.table(ref_input,header=F)

colnames(raw) <- c("chr","start","end","name","pval","strain")
colnames(ref) <- c("chr","postion","depth")
prepeaks <- raw[0,]
peaks <- raw[0,]
temp <- raw
p1 <- progress_bar$new(total = nrow(temp))
print("Initializing...")
for(i in 1:nrow(temp)){
  p1$tick()
  if (temp[i,6] == '-'){
    temp[i,7] <- temp[i,3] - 3
    temp[i,8] <- temp[i,3] - 8
    temp[i,9] <- temp[i,3] + 3
    temp[i,10] <- temp[i,3] + 8
  }
  if (temp[i,6] == '+'){
    temp[i,7] <- temp[i,3] + 3
    temp[i,8] <- temp[i,3] + 8
    temp[i,9] <- temp[i,3] - 3
    temp[i,10] <- temp[i,3] - 8
  }
}
colnames(temp) <- c("chr","start","end","name","pval","strain","high3","high8","low3","low8")
print("Finding depth...")
temp <- left_join(temp,ref,by=c("chr","high3"="postion"))
temp <- left_join(temp,ref,by=c("chr","high8"="postion"))
temp <- left_join(temp,ref,by=c("chr","low3"="postion"))
temp <- left_join(temp,ref,by=c("chr","low8"="postion"))

temp[is.na(temp)] <- 5

p2 <- progress_bar$new(total = nrow(temp))
print("Finding peaks (1st trun)...")
for(i in 1:nrow(temp)){
  p2$tick()
  if (temp[i,11] / temp[i,13] > threshhold | temp[i,11] / temp[i,14] > threshhold | temp[i,12] / temp[i,13] > threshhold | temp[i,12] / temp[i,14] > threshhold){
    prepeaks <- bind_rows(prepeaks,temp[i,1:6])
  }
}

p3 <- progress_bar$new(total = nrow(prepeaks))
print("Successed, starting the 2nd turn...")
for(i in 1:nrow(prepeaks)){
  p3$tick()
  if (prepeaks[i,6] == '-'){
    prepeaks[i,7:22] <- (prepeaks[i,3] - 8):(prepeaks[i,3] + 7)
  }
  if (prepeaks[i,6] == '+'){
    prepeaks[i,7:22] <- (prepeaks[i,2] + 7):(prepeaks[i,2] - 8)
  }
}

colnames(prepeaks) <- c("chr","start","end","name","pval","strain","p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","p13","p14","p15","p16")
print("Finding depth...")
p5 <- progress_bar$new(total = 16)
prepeaks <- left_join(prepeaks,ref,by=c("chr","p1"="postion"))
p5$tick()
prepeaks <- left_join(prepeaks,ref,by=c("chr","p2"="postion"))
p5$tick()
prepeaks <- left_join(prepeaks,ref,by=c("chr","p3"="postion"))
p5$tick()
prepeaks <- left_join(prepeaks,ref,by=c("chr","p4"="postion"))
p5$tick()
prepeaks <- left_join(prepeaks,ref,by=c("chr","p5"="postion"))
p5$tick()
prepeaks <- left_join(prepeaks,ref,by=c("chr","p6"="postion"))
p5$tick()
prepeaks <- left_join(prepeaks,ref,by=c("chr","p7"="postion"))
p5$tick()
prepeaks <- left_join(prepeaks,ref,by=c("chr","p8"="postion"))
p5$tick()
prepeaks <- left_join(prepeaks,ref,by=c("chr","p9"="postion"))
p5$tick()
prepeaks <- left_join(prepeaks,ref,by=c("chr","p10"="postion"))
p5$tick()
prepeaks <- left_join(prepeaks,ref,by=c("chr","p11"="postion"))
p5$tick()
prepeaks <- left_join(prepeaks,ref,by=c("chr","p12"="postion"))
p5$tick()
prepeaks <- left_join(prepeaks,ref,by=c("chr","p13"="postion"))
p5$tick()
prepeaks <- left_join(prepeaks,ref,by=c("chr","p14"="postion"))
p5$tick()
prepeaks <- left_join(prepeaks,ref,by=c("chr","p15"="postion"))
p5$tick()
prepeaks <- left_join(prepeaks,ref,by=c("chr","p16"="postion"))
p5$tick()

prepeaks[is.na(prepeaks)] <- 5

p4 <- progress_bar$new(total = nrow(prepeaks))
print("Finding peaks (2nd trun)...")
# To find out the biggest different point
for(i in 1:nrow(prepeaks)){
  p4$tick()
  high_value <- 0
  for(j in 23:37){
    if(prepeaks[i,j] / prepeaks[i,j+1] > high_value){
      high_value <- prepeaks[i,j] / prepeaks[i,j+1]
      x_high <- j
    }
  }
  if(prepeaks[i,6] == "-"){
    temp_peak <- data.frame(prepeaks[i,1],(prepeaks[i,j-16]-1),(prepeaks[i,j+1-16]-1),prepeaks[i,4],prepeaks[i,5],prepeaks[i,6])
  }
  if(prepeaks[i,6] == "+"){
    temp_peak <- data.frame(prepeaks[i,1],prepeaks[i,j+1-16],prepeaks[i,j-16],prepeaks[i,4],prepeaks[i,5],prepeaks[i,6])
  }
  colnames(temp_peak) <- c("chr","start","end","name","pval","strain")
  peaks <- bind_rows(peaks,temp_peak)
}

print("Successed...")
print("Writing into file...")
print("Finished...")
write.table(peaks,file=output_file_name,sep="\t",row.names = F,col.names = F,quote = F)
