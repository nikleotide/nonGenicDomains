## Before running this script, you need to calculate the average of methylation based on the CpG count and the ratio of methylated Cs in each window of CpGs
# intersectBed -a /media/behnam/Black_Seagate2/Mouse/Data/mm10_1kb_intervals.bed -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -c > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1kb_CpG_counts.dat
# intersectBed -a /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1kb_CpG_counts.dat -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -wa -wb | groupBy -g 1,2,3,4 -c 8 -o sum | awk '{print $0"\t"$5/$4}' > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_Methylation_CpGcount_sum_average.bedgraph


library(changepoint)

options(scipen=999)

DIR="/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/"
newname="nzi"
numberOfWindows=5      
mark="Methylation"

input_mark00<-read.csv(paste(DIR,"Parental_C3H10T_BS_Methylation_CpGcount_sum_average.bedgraph",sep=""), header = FALSE,sep = "\t",col.names = c("Chr","Start","End","CpGs","Sum","Average"),colClasses = c("character","integer","integer","integer","numeric","numeric"))

#
DIR=paste(DIR,"Segmentation_",numberOfWindows,"/",sep="")
system(paste("mkdir -p ",DIR,sep=""))

for (CHR in c("chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"))
{
  print(CHR)
  input_mark0<-subset(input_mark00,(input_mark00$Chr==CHR))
  
  input_mark0$Average<-input_mark0$Average+0.1
  input_mark0$normalized<-(input_mark0$Average/sum(input_mark0$Average))
  
  input_mark0$scaled<-(input_mark0$normalized-min(input_mark0$normalized))/(max(input_mark0$normalized)-min(input_mark0$normalized))
  
  input_mark<-input_mark0
  
  
  cpts.cvrg=cpt.meanvar(input_mark$scaled,minseglen = numberOfWindows,penalty="SIC",method="PELT",test.stat="Normal")

  
  sink(file = paste(DIR,paste(newname,mark,CHR,numberOfWindows,"peaks.bedgraph",sep = "-"),sep=""), append = TRUE)
  
  
  
  cat(paste(CHR,0,input_mark[cpts(cpts.cvrg)[1],"End"],(mean(input_mark$scaled[0:cpts(cpts.cvrg)[1]])), sep = "\t"),"\n",sep="")
  
  numberOfpeaks=length(cpts(cpts.cvrg))-1
  for (i in 1:numberOfpeaks)
  {
    cat(paste(CHR,input_mark[cpts(cpts.cvrg)[i],"Start"],input_mark[cpts(cpts.cvrg)[i+1],"Start"],(mean(input_mark$scaled[cpts(cpts.cvrg)[i]:cpts(cpts.cvrg)[i+1]])), sep = "\t"),"\n",sep="")
  }

  cat(paste(CHR,input_mark[cpts(cpts.cvrg)[length(cpts(cpts.cvrg))],"Start"],input_mark[cpts(cpts.cvrg)[length(cpts(cpts.cvrg))],"End"],(mean(input_mark$scaled[cpts(cpts.cvrg)[length(cpts(cpts.cvrg))]:nrow(input_mark)])), sep = "\t"),"\n",sep="")
  sink()


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


samplename=gsub(".*/","",gsub("/$","",DIR))

system(paste("cat ",DIR,newname,"-",mark,"-chr*",numberOfWindows,"-peaks.bedgraph  | sed 's/^chr//' | sed 's/X/23/' | sort -nk1,1 -nk2,2 | sed 's/^/chr/' | sed 's/chr23\\s/chrX\t/' | grep -v -e '^chrY' | awk '{print $1\"\t\"$2\"\t\"$3-1\"\t\"$4}' > ",DIR,samplename,"-",numberOfWindows,"kb-",mark,"-1stComplete_Genome-peaks.bedgraph", sep=""))


system(paste("grep -v Inf ",DIR,samplename,"-",numberOfWindows,"kb-",mark,"-1stComplete_Genome-peaks.bedgraph > ",DIR,samplename,"-",numberOfWindows,"kb-",mark,"-1stComplete_Genome-peaks.corrected.bedgraph",sep=""))


system(paste("rm ",DIR,"/",newname,"*",sep=""))
system(paste("rm ",DIR,"/Round","*",sep=""))
system(paste("rm ",DIR,"/","*Genome-peaks.bedgraph",sep=""))

