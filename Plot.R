library("ggplot2")
library("hexbin")
options(scipen=999)
par(mfrow=c(3,3))
for (mark in c("K36me2","K27me3","K9me3")){
for (domain in c("LARGE","GENICLARGE","NONGENICLARGE")){


WTFILE<-paste("/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_3_Cells/Copy_of_Final_Results_",domain,"/Parental_C3H10T_BS_",domain,"-Parental-",mark,"-domain_CpG_counts_sum_average.dat",sep="")

WT<-read.csv(WTFILE, header = FALSE, col.names = c("Chr","Start","End","Mark","seg_length","seg_num","CpG","methylation","average_methylation"),sep = "\t",colClasses = c("character","integer","integer","numeric","integer","integer","integer","numeric","numeric"))

K36MFILE<-paste("/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_3_Cells/Copy_of_Final_Results_",domain,"/K36M_C3H10T_BS_",domain,"-Parental-",mark,"-domain_CpG_counts_sum_average.dat",sep="")
K36M<-read.csv(K36MFILE, header = FALSE, col.names = c("Chr","Start","End","Mark","seg_length","seg_num","CpG","methylation","average_methylation"),sep = "\t",colClasses = c("character","integer","integer","numeric","integer","integer","integer","numeric","numeric"))

DKOFILE<-paste("/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_3_Cells/Copy_of_Final_Results_",domain,"/C1_C3H10T_BS_",domain,"-Parental-",mark,"-domain_CpG_counts_sum_average.dat",sep="")
DKO<-read.csv(DKOFILE, header = FALSE, col.names = c("Chr","Start","End","Mark","seg_length","seg_num","CpG","methylation","average_methylation"),sep = "\t",colClasses = c("character","integer","integer","numeric","integer","integer","integer","numeric","numeric"))


WTAVG=sum(WT$methylation)/sum(WT$CpG)
K36MAVG=sum(K36M$methylation)/sum(K36M$CpG)
DKOAVG=sum(DKO$methylation)/sum(DKO$CpG)
xx<-barplot(main=paste(gsub("GENIC","GENIC ",domain)," domains of ",mark," (WT)"), c(WTAVG,K36MAVG,DKOAVG), ylim=c(0,10+max(WTAVG,K36MAVG,DKOAVG)),xlab = "Cell lines", ylab = "Average methylation", names.arg=c("WT","K36M","DKO"))
text(x = xx, y = c(WTAVG,K36MAVG,DKOAVG), label = c(round(WTAVG,digits = 2),round(K36MAVG,digits = 2),round(DKOAVG,digits = 2)), pos = 3, cex = 1.7, col = "red")


}
}
