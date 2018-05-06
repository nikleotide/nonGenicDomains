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
for file in *.bam;do echo $file;samtools flagstat $file;done

# Next count the number of reads in each mentioned marks and their input samples within the large domains:
# first created the run script:


# then add this line to the run scrrip:

# next take the created table to R for plotting:

