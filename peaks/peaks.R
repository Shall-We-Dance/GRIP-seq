# Input: bed file peak form clipper
# Output: real peak for GRIP-seq
# Author: Hongjiang Liu
# Email: hongjiang.liu@ucsf.edu

#load dplyr
library(dplyr)
#import data
args <- commandArgs(trailingOnly=TRUE)
raw_input <- args[1]
ref_input <- args[2]
threshhold <- args[3]
threshhold_peak <- args[3]
output_file_name <- args[4]


raw <- read.table(raw_input,header=F)
ref <- read.table(ref_input,header=F)


result <- raw[0,]
colnames(result) <- c("chr","start","end","strain")
colnames(ref) <- c("chr","postion","depth")

for(i in 1:nrow(raw)){
  print(i)
  ref_chr=filter(ref, chr == raw[i,1])
  #select No.i peak to identify
  if (raw[i,4] == '-'){
    temp_high_3 <- data.frame(raw[i,1],raw[i,3]-3)
    colnames(temp_high_3) <- c("chr","postion")
    temp_high_8 <- data.frame(raw[i,1],raw[i,3]-8)
    colnames(temp_high_8) <- c("chr","postion")
    temp_low_3 <- data.frame(raw[i,1],raw[i,3]+3)
    colnames(temp_low_3) <- c("chr","postion")
    temp_low_8 <- data.frame(raw[i,1],raw[i,3]+8)
    colnames(temp_low_8) <- c("chr","postion")
  }
  if (raw[i,4] == '+'){
    temp_high_3 <- data.frame(raw[i,1],raw[i,3]+3)
    colnames(temp_high_3) <- c("chr","postion")
    temp_high_8 <- data.frame(raw[i,1],raw[i,3]+8)
    colnames(temp_high_8) <- c("chr","postion")
    temp_low_3 <- data.frame(raw[i,1],raw[i,3]-3)
    colnames(temp_low_3) <- c("chr","postion")
    temp_low_8 <- data.frame(raw[i,1],raw[i,3]-8)
    colnames(temp_low_8) <- c("chr","postion")
  }
  temp <- temp_high_3
  temp <- bind_rows(temp,temp_high_8)
  temp <- bind_rows(temp,temp_low_3)
  temp <- bind_rows(temp,temp_low_8)
  temp <- left_join(temp,ref_chr,by=c("chr","postion"))
  #get depth
  #is.na
  for(k in 1:4){
    if(is.na(temp[k,3])){
    temp[k,3] <- 5
    }
  } 
  if(temp[1,3]/temp[3,3] > threshhold | temp[2,3]/temp[4,3] > threshhold | temp[1,3]/temp[4,3] > threshhold | temp[2,3]/temp[3,3] > threshhold){
    if (raw[i,4] == '-'){
      temp_peak <- data.frame(V1=raw[i,1],V2=(raw[i,3]-8):(raw[i,3]+8))
      colnames(temp_peak) <- c("chr","postion")
      temp_peak <- left_join(temp_peak,ref,by=c("chr","postion"))
      for (j in 1:(nrow(temp_peak)-1)){
        #is.na
        if(is.na(temp_peak[j,3])){
          temp_peak[j,3] <- 5
        }
        #is.na
        if(is.na(temp_peak[j+1,3])){
          temp_peak[j+1,3] <- 5
        }
        if(temp_peak[j,3]/temp_peak[j+1,3] > threshhold_peak){
          temp_result = data.frame(raw[i,1],(temp_peak[j,2]-1),(temp_peak[j+1,2]-1),raw[i,4])
          colnames(temp_result) <- c("chr","start","end","strain")
          result = bind_rows(result,temp_result)
          break
        }
      }
    }
    if (raw[i,4] == '+'){
      temp_peak <- data.frame(V1=raw[i,1],V2=(raw[i,2]-8):(raw[i,2]+8))
      colnames(temp_peak) <- c("chr","postion")
      temp_peak <- left_join(temp_peak,ref,by=c("chr","postion"))
      for (j in 1:(nrow(temp_peak)-1)){
        #is.na
        if(is.na(temp_peak[j,3])){
          temp_peak[j,3] <- 5
        }
        #is.na
        if(is.na(temp_peak[j+1,3])){
          temp_peak[j+1,3] <- 5
        }
        if(temp_peak[j+1,3]/temp_peak[j,3] > threshhold_peak){
          temp_result = data.frame(raw[i,1],temp_peak[j,2],temp_peak[j+1,2],raw[i,4])
          colnames(temp_result) <- c("chr","start","end","strain")
          result = bind_rows(result,temp_result)
          break
        }
      }
    }
  }
  if(i%%1000 ==0){
    write.table(result,file=output_file_name,sep="\t",row.names = F,col.names = F,quote = F)
  }
}
write.table(result,file=output_file_name,sep="\t",row.names = F,col.names = F,quote = F)

