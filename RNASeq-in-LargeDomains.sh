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


# To do the counting
# bedtools multicov -q 50 -bed one_of_the_bedgraphs_above -bams list_of_bams


