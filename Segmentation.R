library(changepoint)
options(scipen=999)


DIR="WHERE YOU KEEP YOUR 1kb WINDOW COVERAGE OF INPUT AND MARK READ COUNTS"
samplename=gsub(".*/","",gsub("/$","",DIR))
intermediatename="nzi"
finalname="Round2-2"
numberOfWindows=50     
mark="H3K36me2"

input_mark00<-read.csv(paste(DIR,"4252_4295_C3H10T1_Parental_H3K36me2_input.CVRG",sep=""), header = FALSE, col.names = c("Chr","Start","End","markdepth","input"),sep = "\t",colClasses = c("character","integer","integer","integer","integer"))
#To try a new way, add the line below otherwise, comment it out.
DIR=paste(DIR,"Segmentation_",numberOfWindows,"/",sep="")
system(paste("mkdir -p ",DIR,sep=""))


moment<-c()
for (CHR in c("chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"))
{
  print (CHR)
  input_mark0<-subset(input_mark00,(input_mark00$Chr==CHR))
  
  input_mark0<-subset(input_mark0,(input_mark0$input>0))
  
  input_mark0$input<-input_mark0$input+1
  input_mark0$markdepth<-input_mark0$markdepth+1
  input_mark0$normalized<-(input_mark0$markdepth/sum(input_mark0$markdepth))/(input_mark0$input/sum(input_mark0$input))
  
  input_mark0$scaled<-input_mark0$normalized
  #(input_mark0$normalized-min(input_mark0$normalized))/(max(input_mark0$normalized)-min(input_mark0$normalized))
  input_mark<-input_mark0
  
  moment<-rbind(moment,input_mark0)
  
  
  cpts.cvrg=cpt.meanvar(input_mark$scaled,minseglen = numberOfWindows,penalty="SIC",method="PELT",test.stat="Normal")
  
  
  
  sink(file = paste(DIR,paste(intermediatename,mark,CHR,numberOfWindows,"peaks.bedgraph",sep = "-"),sep=""), append = FALSE)
  cat(paste(CHR,1,input_mark[cpts(cpts.cvrg)[1],"End"],(mean(input_mark$scaled[0:cpts(cpts.cvrg)[1]])), sep = "\t"),"\n",sep="")
  numberOfpeaks=length(cpts(cpts.cvrg))-1
  for (i in 1:numberOfpeaks)
  {
    cat(paste(CHR,input_mark[cpts(cpts.cvrg)[i],"Start"],input_mark[cpts(cpts.cvrg)[i+1],"Start"],(mean(input_mark$scaled[cpts(cpts.cvrg)[i]:cpts(cpts.cvrg)[i+1]])), sep = "\t"),"\n",sep="")
  }
  
  cat(paste(CHR,input_mark[cpts(cpts.cvrg)[length(cpts(cpts.cvrg))],"Start"],input_mark[cpts(cpts.cvrg)[length(cpts(cpts.cvrg))],"End"],(mean(input_mark$scaled[cpts(cpts.cvrg)[length(cpts(cpts.cvrg))]:nrow(input_mark)])), sep = "\t"),"\n",sep="")
  sink()

}

system(paste("cat ",DIR,intermediatename,"-",mark,"-chr*",numberOfWindows,"-peaks.bedgraph  | sed 's/^chr//' | sed 's/X/23/' | sort -nk1,1 -nk2,2 | sed 's/^/chr/' | sed 's/chr23\\s/chrX\t/' | grep -v -e '^chrY' | awk '{print $1\"\t\"$2\"\t\"$3-1\"\t\"$4}' > ",DIR,samplename,"-",numberOfWindows,"kb-",mark,"-1stComplete_Genome-peaks.bedgraph", sep=""))



system(paste("grep -v Inf ",DIR,samplename,"-",numberOfWindows,"kb-",mark,"-1stComplete_Genome-peaks.bedgraph > ",DIR,samplename,"-",numberOfWindows,"kb-",mark,"-1stComplete_Genome-peaks.corrected.bedgraph",sep=""))


to_fix<-read.csv(paste(DIR,samplename,"-",numberOfWindows,"kb-",mark,"-1stComplete_Genome-peaks.corrected.bedgraph",sep=""),sep = "\t",header = FALSE,col.names = c("chr","start","end","seg_value"))
fixed<-to_fix
till=nrow(fixed)-1
for (i in 1:till)
{
  endz=fixed[i,3]
  startz=fixed[(i+1),2]
  if (endz != startz && startz > 1)
  { 
  fixed[i,3] = fixed[(i+1),2]
  }
}

finalfile<-paste(DIR,samplename,"-",numberOfWindows,"kb-",mark,"-1stComplete_Genome-peaks.fixed.bedgraph",sep="")
write.table(fixed,finalfile,sep="\t",col.names = FALSE,row.names = FALSE)

system(paste("sed -i 's/\"//g' ", finalfile,sep = ""))

