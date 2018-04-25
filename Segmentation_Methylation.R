library(changepoint)

options(scipen=999)

DIR="/media/behnam/Black_Seagate2/Mouse/Analysis3/K36M/Methylation/"
newname="nzi"
numberOfWindows=5      
mark="Methylation"


input_mark00<-read.csv(paste(DIR,"K36M_1kb_Methylation_CpGcount_sum_average.bedgraph",sep=""), header = FALSE,sep = "\t",col.names = c("Chr","Start","End","CpGs","Sum","Average"),colClasses = c("character","integer","integer","integer","numeric","numeric"))

#
DIR="/media/behnam/Black_Seagate2/Mouse/Analysis3/K36M/Methylation/5w/"

for (CHR in c("chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"))
{
print(CHR)
input_mark0<-subset(input_mark00,(input_mark00$Chr==CHR))

input_mark0$Average<-input_mark0$Average+0.1
input_mark0$normalized<-(input_mark0$Average/sum(input_mark0$Average))

input_mark0$scaled<-(input_mark0$normalized-min(input_mark0$normalized))/(max(input_mark0$normalized)-min(input_mark0$normalized))

input_mark<-input_mark0


cpts.cvrg=cpt.meanvar(input_mark$scaled,minseglen = numberOfWindows,penalty="SIC",method="PELT",test.stat="Normal")
#plot(x = input_mark$scaled, ylim = c(min(input_mark$scaled),max(input_mark$scaled)), type = "h", cex = .1, xaxt='n', ylab ="Scaled read count", col="black" , frame.plot=F, main ="", xlab=grange)

#segments(0,mean(input_mark$scaled[0:cpts(cpts.cvrg)[1]]),cpts(cpts.cvrg)[1],mean(input_mark$scaled[0:cpts(cpts.cvrg)[1]]),col = "red", lwd=4)

sink(file = paste(DIR,paste(newname,mark,CHR,numberOfWindows,"peaks.bedgraph",sep = "-"),sep=""), append = TRUE)



cat(paste(CHR,0,input_mark[cpts(cpts.cvrg)[1],"End"],(mean(input_mark$scaled[0:cpts(cpts.cvrg)[1]])), sep = "\t"),"\n",sep="")

#format(x$V2[cpts(cpts.cvrg)[length(cpts(cpts.cvrg))]],scientific = FALSE),"\t",format(x$V2[length(x$V7)],scientific = FALSE)
numberOfpeaks=length(cpts(cpts.cvrg))-1
for (i in 1:numberOfpeaks)
{
  #segments(cpts(cpts.cvrg)[i],mean(input_mark$scaled[cpts(cpts.cvrg)[i]:cpts(cpts.cvrg)[i+1]]),cpts(cpts.cvrg)[i+1],mean(input_mark$scaled[cpts(cpts.cvrg)[i]:cpts(cpts.cvrg)[i+1]]),col = "red", lwd=4)
  cat(paste(CHR,input_mark[cpts(cpts.cvrg)[i],"Start"],input_mark[cpts(cpts.cvrg)[i+1],"Start"],(mean(input_mark$scaled[cpts(cpts.cvrg)[i]:cpts(cpts.cvrg)[i+1]])), sep = "\t"),"\n",sep="")
}
#segments(cpts(cpts.cvrg)[length(cpts(cpts.cvrg))],mean(input_mark$scaled[cpts(cpts.cvrg)[length(cpts(cpts.cvrg))]:nrow(input_mark)]),nrow(input_mark),mean(input_mark$scaled[cpts(cpts.cvrg)[length(cpts(cpts.cvrg))]:nrow(input_mark)]),col = "red", lwd=4)

cat(paste(CHR,input_mark[cpts(cpts.cvrg)[length(cpts(cpts.cvrg))],"Start"],input_mark[cpts(cpts.cvrg)[length(cpts(cpts.cvrg))],"End"],(mean(input_mark$scaled[cpts(cpts.cvrg)[length(cpts(cpts.cvrg))]:nrow(input_mark)])), sep = "\t"),"\n",sep="")
sink()

# 
# 
# ### Round2-1 (round(scaled(zscore across chromosome)))
# 
# sink(file = paste(DIR,paste("Round2-1",mark,CHR,numberOfWindows,"peaks.bedgraph",sep = "-"),sep=""), append = FALSE)
# segmented1<-read.csv(paste(DIR,paste(newname,mark,CHR,numberOfWindows,"peaks.bedgraph",sep = "-"),sep=""), header = FALSE, col.names = c("Chr","Start","End","scaled_mark"),sep = "\t",colClasses = c("character","numeric","numeric","numeric"))
#  
#  for (g in 1:nrow(segmented1))
#  {
#    segmented1[g,"status0"]<-(segmented1[g,"scaled_mark"]-mean(segmented1$scaled_mark))/sd(segmented1$scaled_mark)
#  }
#  for (h in 1:nrow(segmented1))
#  {
#    segmented1[h,"status1"]<-round((segmented1[h,"status0"]-min(segmented1$status0))/(max(segmented1$status0)-min(segmented1$statu0)),digits = 1)
#    cat(paste(segmented1[h,"Chr"],segmented1[h,"Start"],segmented1[h,"End"],segmented1[h,"status1"], sep = "\t"),"\n",sep="")
#  }
# sink()

### Round2-2 (round(scaled(zscore in three segments)))

segmented2<-read.csv(paste(DIR,paste(newname,mark,CHR,numberOfWindows,"peaks.bedgraph",sep = "-"),sep=""), header = FALSE, col.names = c("Chr","Start","End","scaled_mark"),sep = "\t",colClasses = c("character","numeric","numeric","numeric"))
sink(file = paste(DIR,paste("Round2-2",mark,CHR,numberOfWindows,"peaks.bedgraph",sep = "-"),sep=""), append = FALSE)
for (j in 1:nrow(segmented2))
{
  if(j>1 && j < nrow(segmented2)){
    segmented2[j,"status02"]<-(segmented2[j,"scaled_mark"]-mean(segmented2[j-1:j+1,"scaled_mark"]))/(sd(segmented2[j-1:j+1,"scaled_mark"]))  
  }
  else if (j == 1)
  {
    segmented2[j,"status02"]<-(segmented2[j,"scaled_mark"]-mean(segmented2[1:3,"scaled_mark"]))/(sd(segmented2[1:3,"scaled_mark"]))  
  }
  else if (j == nrow(segmented2)){
    segmented2[j,"status02"]<-(segmented2[j,"scaled_mark"]-mean(segmented2[j-2:j,"scaled_mark"]))/(sd(segmented2[j-2:j,"scaled_mark"]))
  }
}
for (l in 1:nrow(segmented2))
{
  segmented2[l,"status2"]<-(segmented2[l,"status02"]-min(segmented2$status02))/(max(segmented2$status02)-min(segmented2$status02))
  cat(paste(segmented2[l,"Chr"],segmented2[l,"Start"],segmented2[l,"End"],segmented2[l,"status2"], sep = "\t"),"\n",sep="")
}
sink()

}
