# This code needs two arguments first (agrs[6]) input annotated_sample_nucleotide file and second (args[7]) output annotated_gap_region file.
#command: $ Rscript ./step2_R_code_for_gap_generation_v1.r Depth_for_SAMPLE.gap_region.txt Depth_for_SAMPLE.gap_region.txt_final.csv
args <- commandArgs()
#print(args)
Temp <- read.delim(args[6], header=FALSE, quote="")
Gap <- matrix(,nrow=length(Temp[,2]),ncol=5)
j<-1
for(i in 2:length(Temp[,2]))
{
  tryCatch({
  Gap[j,1]<-as.character(Temp[i-1,1]);
  Gap[j,2]<-Temp[i-1,2];
  Gap[j,4]<-as.character(Temp[i-1,3]);
  Gap[j,5]<-as.character(Temp[i-1,4]);
  while(as.character(Temp[i,1])==as.character(Temp[i-1,1]) && Temp[i,2]==(Temp[i-1,2]+1))
  {
    i<-i+1;
  }
  #print("this is second")
  #print(i)
  Gap[j,3]<-Temp[i-1,2];
  j<-j+1;
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
Gap<-Gap[!is.na(Gap[,1]),]
Gap[length(Gap[,3]),3]<- Temp[length(Temp[,3]),2]  # this step done coz of missing last entry
Gap_final<- matrix(,nrow=length(Gap[,2]),ncol=5)  # New this step is to pick the right informations
for (i in 1:length(Gap[,2]))
{
  if(Gap[i,2]==Gap[i,3])
  {
    Gap_final[i+1,]<-(Gap[i+1,])
  }
}
Gap_final[1,]<-Gap[1,]
Gap_final<-Gap_final[!is.na(Gap_final[,1]),]
colnames(Gap_final)<- c("Chr_no","start","end","Gene","NM_annotation")
Gap_final<-as.data.frame(Gap_final)
write.csv(x = Gap_final, file = args[7])
rm(Gap_final)
rm(Temp)
rm(Gap) 
