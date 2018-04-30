## Creating the 1kb window bed files:
if [ -e "mm10.chrom.sizes" ]
then
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes &&
cat mm10.chrom.sizes  | grep -v -e '_\|chrM' | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/chr//' | sort -k1,1n -k2,2n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > mm10_21chr.genome
fi
if [ -e "mm10_1kb_intervals.bed" ]
then
bedtools makewindows -g mm10_21chr.genome -w 1000 | awk '{print $1"\t"$2+1"\t"$3}' > mm10_1kb_intervals.bed
fi
## Creating 1kb windows from the bam files (adjust the path accordingly) using the bed file created above:
bedtools multicov -q 50 -bed mm10_1kb_intervals.bed -bams 4252_4295_C3H10T1_K36M_H3K36me2.sorted.bam 4252_4295_C3H10T1_K36M_input.sorted.bam > K36M_1kb_K36me2_input.cvrg

# This file will have two additional columns (in addition to position columns), one is the read count for hsitone mark and one the read count for input.
head -200000 /media/behnam/Black_Seagate2/Mouse/Analysis3/K36M/K36me2_200b/K36M_1kb_K36me2_input.cvrg | tail
#chr1	39998401	39998600	20	7
#chr1	39998601	39998800	21	10
#chr1	39998801	39999000	19	2
#chr1	39999001	39999200	22	3
#chr1	39999201	39999400	22	3
#chr1	39999401	39999600	26	4
#chr1	39999601	39999800	32	8
#chr1	39999801	40000000	28	3
#chr1	40000001	40000200	32	3
#chr1	40000201	40000400	23	3



#The result of this preparation will be fed to the segmentation script (R script).



