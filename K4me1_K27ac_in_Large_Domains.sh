# To calculate the RPKM in large domains in the 4 cell lines, first we need to know the number of mapped reads per sample.
# go to the directory where you store rnaseq bam files with their index files and then:
for file in *.bam;do echo $file;samtools flagstat $file | grep mapped | grep -e '(' | grep -e 'N/A' | awk '{print $1}';done

#Here are the results:
#B2-1.sorted.bam
81126550
#B2-2.sorted.bam
54292876
#B2-3.sorted.bam
49819820
#C22-1.sorted.bam
72154994
#C22-2.sorted.bam
62585132
#C22-3.sorted.bam
54084870
#K36M-1.sorted.bam
74441204
#K36M-2.sorted.bam
55633502
#K36M-3.sorted.bam
44524276
#P2-1.sorted.bam
69928354
#P2-2.sorted.bam
51730034
#P2-3.sorted.bam
55686982
#D1_rep1.sorted.bam
99382598
#D1_rep2.sorted.bam
102524850
#D1_rep3.sorted.bam
102674864
#K36R-1.sorted.bam
47665668
#K36R-2.sorted.bam
44597754
#K36R-3.sorted.bam
42799344


# Next, running bedtools multicov on all the bam files at the same time using the domains and subdomains (genic/nongenic) as the input bed file. This needs to be done for all WT-domain-marks

# After each multicov run, add a line with column names:

sed -i '1s/^/Chr\tStart\tEnd\tMark\tSeg_legnth\tSeg_num\tB2-1\tB2-2\tB2-3\tC22-1\tC22-2\tC22-3\tD1_rep1\tD1_rep2\tD1_rep3\tK36M-1\tK36M-2\tK36M-3\tK36R-1\tK36R-2\tK36R-3\tP2-1\tP2-2\tP2-3\n/' file

## The bedtools are here:
###### K36me2

# LARGE K36me2
/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K36me2_Domains_vs_Methylation_in_4_Cells/K36me2_Domains/H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.ALL.LARGE.bedgraph

# NONGENIC LARGE K36me2
/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K36me2_Domains_vs_Methylation_in_4_Cells/K36me2_Domains/H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.NONGENIC.LARGE.bedgraph

# GENIC LARGE K36me2
/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K36me2_Domains_vs_Methylation_in_4_Cells/K36me2_Domains/H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.GENIC.LARGE.bedgraph


###### K9me3
# LARGE K9me3
/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K9me3_Domains_vs_Methylation_in_4_Cells/K9me3_Domains/H3K9me3-50kb-H3K9me3-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.ALL.LARGE.bedgraph

# NONGENIC LARGE K9me3
/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K9me3_Domains_vs_Methylation_in_4_Cells/K9me3_Domains/H3K9me3-50kb-H3K9me3-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.NONGENIC.LARGE.bedgraph

# GENIC LARGE K9me3
/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K9me3_Domains_vs_Methylation_in_4_Cells/K9me3_Domains/H3K9me3-50kb-H3K9me3-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.GENIC.LARGE.bedgraph

###### K27me3

# LARGE K27me3
/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K27me3_Domains_vs_Methylation_in_4_Cells/K27me3_Domains/H3K27me3-50kb-H3K27me3-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.ALL.LARGE.bedgraph

# NONGENIC LARGE K27me3
/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K27me3_Domains_vs_Methylation_in_4_Cells/K27me3_Domains/H3K27me3-50kb-H3K27me3-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.NONGENIC.LARGE.bedgraph

# GENIC LARGE K27me3
/media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K27me3_Domains_vs_Methylation_in_4_Cells/K27me3_Domains/H3K27me3-50kb-H3K27me3-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.GENIC.LARGE.bedgraph

# I put all these locations in a file called: bedfiles_run.sh and then run this script (Stored in complete_run.sh)
# To created the complete_run.sh I used this command below:

cat bedfiles_run.sh | awk '{if ($0 !~ "#" && $0 != "") {print "bedtools multicov -q 50 -bed "$2" -bams B2-1.sorted.bam B2-2.sorted.bam B2-3.sorted.bam C22-1.sorted.bam C22-2.sorted.bam C22-3.sorted.bam D1_rep1.sorted.bam D1_rep2.sorted.bam D1_rep3.sorted.bam K36M-1.sorted.bam K36M-2.sorted.bam K36M-3.sorted.bam K36R-1.sorted.bam K36R-2.sorted.bam K36R-3.sorted.bam P2-1.sorted.bam P2-2.sorted.bam P2-3.sorted.bam > "$1".0.tsv"} else print $0}' >> complete_run.sh

# and added this line to the end of complete_run.sh to add headers:
for file in *.0.tsv
do
sed '1s/^/Chr\tStart\tEnd\tMark\tSeg_legnth\tSeg_num\tB2-1\tB2-2\tB2-3\tC22-1\tC22-2\tC22-3\tD1_rep1\tD1_rep2\tD1_rep3\tK36M-1\tK36M-2\tK36M-3\tK36R-1\tK36R-2\tK36R-3\tP2-1\tP2-2\tP2-3\n/' $file > ${file/.0.tsv/.tsv}
done


### Next repeat the same for the merged RNASeq bam files (repeats 1,2 and 3)
## after merging, sorting and indexing
## Create the running shell script:

cat ../bedfiles_run.sh | awk '{if ($0 !~ "#" && $0 != "") {print "bedtools multicov -q 50 -bed "$2" -bams B2.Sorted.bam  C22.Sorted.bam  D1.Sorted.bam  K36M.Sorted.bam  K36R.Sorted.bam  P2.Sorted.bam > "$1".merged.0.tsv"} else print $0}' >> complete_run_mergedBams.sh

# After each multicov run (complete_run_mergedBams.sh), add a line with column names:
for file in *.0.tsv
do
sed '1s/^/Chr\tStart\tEnd\tMark\tSeg_legnth\tSeg_num\tB2\tC22\tD1\tK36M\tK36R\tP2\n/' $file > 
${file/.0.tsv/.tsv}
done                                                                                                                              

# For plotting, the mapped reads for rpkm
#B2.Sorted.bam
185239246
#C22.Sorted.bam
188824996
#D1.Sorted.bam
304582312
#K36M.Sorted.bam
174598982
#K36R.Sorted.bam
135062766
#P2.Sorted.bam
177345370



## To look into the H3K4me1 and H3K27ac marks within large domains (mentioned above) in different cell lines:
# count the number of mapped read counts in each bam file (H3K4me1 and H3K27ac in different cell lines):
for file in *.bam;do echo $file;samtools flagstat $file | grep mapped | grep -e '(' | grep -e 'N/A' | awk '{print $1}';done
# Here are the results


# Next count the number of reads in each mentioned marks and their input samples within the large domains:
# first created the run script:

cat ../../Data/bedfiles_run.sh | awk '{if ($0 !~ "#" && $0 != "") {print "echo \""$1"\"\n""bedtools multicov -q 50 -bed "$2" -bams C3H10T1-2-B1-Rx_H3K27ac.sorted.bam C3H10T1-2-B1-Rx_INPUTK27ac.sorted.bam C3H10T1-2-B1-Rx_H3K4me1.sorted.bam C3H10T1-2-B1-Rx_Input.sorted.bam C3H10T1-2-C1-Rx_H3K27ac.sorted.bam C3H10T1-2-C1-Rx_INPUTK27ac.sorted.bam C3H10T1-2-C1-Rx_H3K4me1.sorted.bam C3H10T1-2-C1-Rx_Input.sorted.bam C3H10T1-2-K36M-Rx_H3K27ac.sorted.bam C3H10T1-2-K36M-Rx_INPUTK27ac.sorted.bam C3H10T1-2-K36M-Rx_H3K4me1.sorted.bam C3H10T1-2-K36M-Rx_Input.sorted.bam C3H10T1-2-Pa-Rx_H3K27ac.sorted.bam C3H10T1-2-Pa-Rx_INPUTK27ac.sorted.bam C3H10T1-2-Pa-III-Rx_H3K4me1.sorted.bam C3H10T1-2-Pa-III-Rx_Input.sorted.bam > "$1".k4me1k27acmarks.0.tsv"} else print $0}' >> complete_run_marks.sh

# then add this line to the run scrrip:
for file in *.0.tsv
do
sed '1s/^/Chr\tStart\tEnd\tMark\tSeg_legnth\tSeg_num\tB1K27ac\tB1K27acinput\tB1K4me1\tB1K4me1\tC1K27ac\tC1K27acinput\tC1K4me1\tC1K4me1\tK36MK27ac\tK36MK27acinput\tK36MK4me1\tK36MK4me1\tPaK27ac\tPaK27acinput\tPaK4me1\tPaK4me1\n/' $file > ${file/.0.tsv/.tsv}
done   
# next take the created table to R for plotting:
K27ac_and_K4me1_in_Large_Domains.R



##### For K27me3, K36me3 and K9m3 marks in large domainas (similar to what you did for K4me1 and K27ac
# First, get the total mapped reads:
for file in *.bam;do echo $file;samtools flagstat $file  | grep mapped | grep -e '(' | grep -e 'N/A' | awk '{print $1}';done

#4252_4295_C3H10T1_Parental_H3K36me2.sorted.bam
54341330
#4252_4295_C3H10T1_Parental_input.sorted.bam
39822015
#4252_4295_C3H10T1_K36M_H3K36me2.sorted.bam
47646999
#4252_4295_C3H10T1_K36M_input.sorted.bam
43905805
#4252_4295_C3H10T1_DKO_H3K36me2.sorted.bam
47612786
#4252_4295_C3H10T1_DKO_input.sorted.bam
45335522
#4252_4295_C3H10T1_SETD2-KO_H3K36me2.sorted.bam
51066735
#4252_4295_C3H10T1_SETD2-KO_input.sorted.bam
40044470
#4252_4295_C3H10T1_TKO_H3K36me2.sorted.bam
12305874
#4252_4295_C3H10T1_TKO_input.sorted.bam
30205361
#3607_C3H10T1_Parental_H3K27me3.sorted.bam
42089195
#3641_C3H10T1_Parental_input.sorted.bam
45470506
#3607_C3H10T1_K36M_H3K27me3.sorted.bam
54398433
#3641_C3H10T1_K36M_input.sorted.bam
31466631
#3607_C3H10T1_DKO_H3K27me3.sorted.bam
58079645
#3641_C3H10T1_DKO_input.sorted.bam
44031740
#3607_C3H10T1_SETD2-KO_H3K27me3.sorted.bam
68184802
#3641_C3H10T1_SETD2-KO_input.sorted.bam
53346926
#4290_C3H10T1_Parental_H3K9me3.sorted.bam
69694859
#4290_C3H10T1_Parental_input.sorted.bam
46857429
#4290_C3H10T1_K36M_H3K9me3.sorted.bam
50394716
#4290_C3H10T1_K36M_input.sorted.bam
62088889
#4290_C3H10T1_DKO_H3K9me3.sorted.bam
44335876
#4290_C3H10T1-2-C1-Rx_Input_H3K9me3.sorted.bam
39163222


# Next, count the reads mapped in each interval:
# first created the run script:

cat ../../Data/bedfiles_run.sh | awk '{if ($0 !~ "#" && $0 != "") {print "echo \""$1"\"\n""bedtools multicov -q 50 -bed "$2" -bams 4252_4295_C3H10T1_Parental_H3K36me2.sorted.bam 4252_4295_C3H10T1_Parental_input.sorted.bam 4252_4295_C3H10T1_K36M_H3K36me2.sorted.bam 4252_4295_C3H10T1_K36M_input.sorted.bam 4252_4295_C3H10T1_DKO_H3K36me2.sorted.bam 4252_4295_C3H10T1_DKO_input.sorted.bam 4252_4295_C3H10T1_SETD2-KO_H3K36me2.sorted.bam 4252_4295_C3H10T1_SETD2-KO_input.sorted.bam 4252_4295_C3H10T1_TKO_H3K36me2.sorted.bam 4252_4295_C3H10T1_TKO_input.sorted.bam 3607_C3H10T1_Parental_H3K27me3.sorted.bam 3641_C3H10T1_Parental_input.sorted.bam 3607_C3H10T1_K36M_H3K27me3.sorted.bam 3641_C3H10T1_K36M_input.sorted.bam 3607_C3H10T1_DKO_H3K27me3.sorted.bam 3641_C3H10T1_DKO_input.sorted.bam 3607_C3H10T1_SETD2-KO_H3K27me3.sorted.bam 3641_C3H10T1_SETD2-KO_input.sorted.bam 4290_C3H10T1_Parental_H3K9me3.sorted.bam 4290_C3H10T1_Parental_input.sorted.bam 4290_C3H10T1_K36M_H3K9me3.sorted.bam 4290_C3H10T1_K36M_input.sorted.bam 4290_C3H10T1_DKO_H3K9me3.sorted.bam 4290_C3H10T1-2-C1-Rx_Input_H3K9me3.sorted.bam > "$1".K36me2K27me3K9me3Marks.0.tsv"} else print $0}' >> complete_run_marks.sh

# then add this line to the run scrrip:
for file in *.0.tsv
do
sed '1s/^/Chr\tStart\tEnd\tMark\tSeg_legnth\tSeg_num\tParental_H3K36me2\tParental_input\tK36M_H3K36me2\tK36M_input\tDKO_H3K36me2\tDKO_input\tSETD2-KO_H3K36me2\tSETD2-KO_input\tTKO_H3K36me2\tTKO_input\tParental_H3K27me3\tParental_input\tK36M_H3K27me3\tK36M_input\tDKO_H3K27me3\tDKO_input\tSETD2-KO_H3K27me3\tSETD2-KO_input\tParental_H3K9me3\tParental_input\tK36M_H3K9me3\tK36M_input\tDKO_H3K9me3\tC1-Rx_Input_H3K9me3\n/' $file > ${file/.0.tsv/.tsv}
done  
# next take the created table to R for plotting:
K36me2K27me3K4me1_in_Large_Domains.R



