MRC_B1_Rx_H3K27ac=29508831
MRC_B1_Rx_INPUTK27ac=23469486

MRC_B1_Rx_H3K4me1=51885134
MRC_B1_Rx_Input=55979718

MRC_C1_Rx_H3K27ac=35968139
MRC_C1_Rx_INPUTK27ac=30139533

MRC_C1_Rx_H3K4me1=70786886
MRC_C1_Rx_Input=36517535

MRC_K36M_Rx_H3K27ac=46194343
MRC_K36M_Rx_INPUTK27ac=19674991

MRC_K36M_Rx_H3K4me1=56295648
MRC_K36M_Rx_Input=27005358

MRC_Pa_III_Rx_H3K4me1=45739496
MRC_Pa_III_Rx_Input=32147992

MRC_Pa_Rx_H3K27ac=31530289
MRC_Pa_Rx_INPUTK27ac=6766509


pdf("/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/H3K4me1_H3K27ac_in_WT_Domains/K27ac_in_Large_Domains.pdf")
par(mfrow=c(3,3))
for (mark in c("K36me2","K27me3","K9me3")){
for (domain in c("LARGE","GENIC_LARGE","NONGENIC_LARGE")){
inputfile=paste("/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/H3K4me1_H3K27ac_in_WT_Domains/",domain,"_",mark,".k4me1k27acmarks.tsv",sep="")
input<-read.csv(inputfile, header=TRUE,sep = "\t",colClasses = c("character","integer","integer","numeric","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer"))

normalized_B1ac<-(sum(input$B1K27ac)/MRC_B1_Rx_H3K27ac)/(sum(input$B1K27acinput)/MRC_B1_Rx_INPUTK27ac)

normalized_C1ac<-(sum(input$C1K27ac)/MRC_C1_Rx_H3K27ac)/(sum(input$C1K27acinput)/MRC_C1_Rx_INPUTK27ac)

normalized_K36Mac<-(sum(input$K36MK27ac)/MRC_K36M_Rx_H3K27ac)/(sum(input$K36MK27acinput)/MRC_K36M_Rx_INPUTK27ac)

normalized_Paac<-(sum(input$PaK27ac)/MRC_Pa_Rx_H3K27ac)/(sum(input$PaK27acinput)/MRC_Pa_Rx_INPUTK27ac)

colnamesac<-c(7,11,15,19)
xx<-barplot(c(normalized_B1ac,normalized_C1ac,normalized_K36Mac,normalized_Paac),ylim=c(0,2.5),ylab = "Normalized K27ac",names.arg=gsub("K27ac","",colnames(input)[colnamesac]) ,las=2,main=paste(domain," ",mark," Domains",sep=""),cex.main = 1)

text(x = xx, y = c(normalized_B1ac,normalized_C1ac,normalized_K36Mac,normalized_Paac), label = c(round(normalized_B1ac,digits=2),round(normalized_C1ac,digits=2), round(normalized_K36Mac,digits=2),round(normalized_Paac,digits=2)), pos = 3, cex = 1, col = "red")

}
}
dev.off()

pdf("/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/H3K4me1_H3K27ac_in_WT_Domains/K4me1_in_Large_Domains.pdf")
par(mfrow=c(3,3))

for (mark in c("K36me2","K27me3","K9me3")){
  for (domain in c("LARGE","GENIC_LARGE","NONGENIC_LARGE")){
    inputfile=paste("/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/H3K4me1_H3K27ac_in_WT_Domains/",domain,"_",mark,".k4me1k27acmarks.tsv",sep="")
    input<-read.csv(inputfile, header=TRUE,sep = "\t",colClasses = c("character","integer","integer","numeric","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer"))
    
    normalized_B1k4me1<-(sum(input$B1K4me1)/MRC_B1_Rx_H3K4me1)/(sum(input$B1K4me1.1)/MRC_B1_Rx_Input)
    
    normalized_C2_k4me1<-(sum(input$C1K4me1)/MRC_C1_Rx_H3K4me1)/(sum(input$C1K4me1.1)/MRC_C1_Rx_Input)
    
    normalized_K36MK4me1<-(sum(input$K36MK4me1)/MRC_K36M_Rx_H3K4me1)/(sum(input$K36MK4me1.1)/MRC_K36M_Rx_Input)
    
    normalized_Pak4me1<-(sum(input$PaK4me1)/MRC_Pa_III_Rx_H3K4me1)/(sum(input$PaK4me1.1)/MRC_Pa_III_Rx_Input)
    

    colnamesk4me1<-c(9,13,17,21)
    xx<-barplot(c(normalized_B1k4me1,normalized_C2_k4me1,normalized_K36MK4me1,normalized_Pak4me1),ylim=c(0,2.5),ylab = "Normalized K4me1",names.arg=gsub("K4me1","",colnames(input)[colnamesk4me1]) ,las=2,main=paste(domain," ",mark," Domains",sep=""),cex.main = 1)
    
    text(x = xx, y = c(normalized_B1k4me1,normalized_C2_k4me1,normalized_K36MK4me1,normalized_Pak4me1), label = c(round(normalized_B1k4me1,digits=2), round(normalized_C2_k4me1,digits=2),round(normalized_K36MK4me1,digits=2),round(normalized_Pak4me1,digits=2)), pos = 3, cex = 1, col = "red")
    
  }
}
dev.off()
