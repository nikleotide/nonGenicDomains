#!/bin/bash


# For the sake of figure 1c, we don't need the trailing domains. Only domains that are very obvious. Thus I only take the > 1 values.

####################################################################################################
## A) Segmentation and merging the histone marks, filtering for 500kb size and then intersecting them with genic/nongenic 
####################################################################################################
## 1kbWindow > Segment50w > Filter > Merge > Intersect
#NOTE: In order to call the trailing streteches of high peaks as separate entities, I use to different threshold for filtering the segments (based on their values). One of them, I used threshold of 1.0 and the other one 0.5-1.0 and then merge each subset separately. Next I will concatenate the resulting segments and intersect them with genic/nongenic
for file in */Segmentation_50/*50kb-*fixed.bedgraph
do
echo $file 
awk '{if ($4 >= 1.0) print $0"\t"($3-$2)/1000}' $file > ${file/.fixed/.fixed.filtered_VALUE1.0.bedgraph}
echo
bedtools merge -d 10 -i ${file/.fixed/.fixed.filtered_VALUE1.0.bedgraph} -c 4,5,1 -o mean,sum,count > ${file/fixed.bedgraph/-FILTERED_VALUE1.0-MERGED-ALL.bedgraph}
echo
cat ${file/fixed.bedgraph/-FILTERED_VALUE1.0-MERGED-ALL.bedgraph} | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed -e 's/^chr//' | sort -k1,1n -k2,2n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.ALL.bedgraph}
echo
awk '{if ($3-$2 > 500000) print $0}' ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.ALL.bedgraph} > ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.LARGE.bedgraph}
echo

bedtools intersect -b /media/behnam/Black_Seagate2/Mouse/Data/mm10hgTables_NonGenic.bed -a ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.LARGE.bedgraph} | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/chr//' | sort -k1,1n -k2,2n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.NONGENIC.LARGE.bedgraph}
echo
bedtools intersect -b /media/behnam/Black_Seagate2/Mouse/Data/mm10hgTables_Genic.bed -a ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.LARGE.bedgraph} | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/chr//' | sort -k1,1n -k2,2n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.GENIC.LARGE.bedgraph}

rm ${file/fixed.bedgraph/-FILTERED_VALUE1.0-MERGED-ALL.bedgraph}
rm ${file/.fixed/.fixed.filtered_VALUE1.0.bedgraph}
done


####################################################################################################
## B) Calculating the average methylation in the large, large-nongenic and large-genic domain 
####################################################################################################
# This is only for K9me3, for other marks only replace all the K9me3s with the name of the replacement mark
## ALL Large

bedtools intersect -a /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K9me3_Domains_vs_Methylation_in_4_Cells/K9me3_Domains/H3K9me3-50kb-H3K9me3-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.ALL.LARGE.bedgraph -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Raw_Methylation_bedgraphs_3_Cells/B1_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -c > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K9me3_Domains_vs_Methylation_in_4_Cells/CpG-Counts/B1_C3H10T_BS_ALL.LARGE-K9me3_CpG_counts.dat



bedtools intersect -a /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K9me3_Domains_vs_Methylation_in_4_Cells/CpG-Counts/B1_C3H10T_BS_ALL.LARGE-K9me3_CpG_counts.dat -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Raw_Methylation_bedgraphs_3_Cells/B1_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -wa -wb | groupBy -g 1,2,3,4,5,6,7 -c 11 -o sum | awk '{print $0"\t"$8/$7}' > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K9me3_Domains_vs_Methylation_in_4_Cells/CpG_counts_sum_average/B1_C3H10T_BS_LARGE-Parental-K9me3-domain_CpG_counts_sum_average.dat


# NONGENIC
bedtools intersect -a /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K9me3_Domains_vs_Methylation_in_4_Cells/K9me3_Domains/H3K9me3-50kb-H3K9me3-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.NONGENIC.LARGE.bedgraph -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Raw_Methylation_bedgraphs_3_Cells/B1_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -c > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K9me3_Domains_vs_Methylation_in_4_Cells/CpG-Counts/B1_C3H10T_BS_NONGENIC.LARGE-K9me3_CpG_counts.dat



bedtools intersect -a /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K9me3_Domains_vs_Methylation_in_4_Cells/CpG-Counts/B1_C3H10T_BS_NONGENIC.LARGE-K9me3_CpG_counts.dat -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Raw_Methylation_bedgraphs_3_Cells/B1_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -wa -wb | groupBy -g 1,2,3,4,5,6,7 -c 11 -o sum | awk '{print $0"\t"$8/$7}' > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K9me3_Domains_vs_Methylation_in_4_Cells/CpG_counts_sum_average/B1_C3H10T_BS_NONGENICLARGE-Parental-K9me3-domain_CpG_counts_sum_average.dat


# GENIC
bedtools intersect -a /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K9me3_Domains_vs_Methylation_in_4_Cells/K9me3_Domains/H3K9me3-50kb-H3K9me3-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.GENIC.LARGE.bedgraph -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Raw_Methylation_bedgraphs_3_Cells/B1_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -c > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K9me3_Domains_vs_Methylation_in_4_Cells/CpG-Counts/B1_C3H10T_BS_GENIC.LARGE-K9me3_CpG_counts.dat



bedtools intersect -a /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K9me3_Domains_vs_Methylation_in_4_Cells/CpG-Counts/B1_C3H10T_BS_GENIC.LARGE-K9me3_CpG_counts.dat -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Raw_Methylation_bedgraphs_3_Cells/B1_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -wa -wb | groupBy -g 1,2,3,4,5,6,7 -c 11 -o sum | awk '{print $0"\t"$8/$7}' > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Parental_K9me3_Domains_vs_Methylation_in_4_Cells/CpG_counts_sum_average/B1_C3H10T_BS_GENICLARGE-Parental-K9me3-domain_CpG_counts_sum_average.dat


####################################################################################################
## B-2) Looping over all cells and marks (MAKE SURE ABOUT THE MARK NAMES)
####################################################################################################
# the directory strucyture looks like this:
# main= /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Methylation_in_3_Parental_Domains_3_Cells/Parental_K27me3_Domains_vs_Methylation_in_3_Cells

# first creating the cpg count files
#main/CpG-Counts
#go to the methylation directory and run these three loops

for file in *.profile.cg_strand_combined.bedgraph; do bedtools intersect -a ../Parental_K9me3_Domains_vs_Methylation_in_3_Cells/K9me3_Domains/H3K9me3-50kb-H39me3-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.LARGE.bedgraph -b $file -c > ../Parental_K9me3_Domains_vs_Methylation_in_3_Cells/CpG-Counts/${file/1.profile.cg_strand_combined.bedgraph/LARGE-K36me2_CpG_counts.dat};done


for file in *.profile.cg_strand_combined.bedgraph; do bedtools intersect -a ../Parental_K9me3_Domains_vs_Methylation_in_3_Cells/K9me3_Domains/H3K9me3-50kb-H39me3-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.GENIC.LARGE.bedgraph -b $file -c > ../Parental_K9me3_Domains_vs_Methylation_in_3_Cells/CpG-Counts/${file/1.profile.cg_strand_combined.bedgraph/GENIC.LARGE-K36me2_CpG_counts.dat};done


for file in *.profile.cg_strand_combined.bedgraph; do bedtools intersect -a ../Parental_K9me3_Domains_vs_Methylation_in_3_Cells/K9me3_Domains/H3K9me3-50kb-H39me3-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.NONGENIC.LARGE.bedgraph -b $file -c > ../Parental_K9me3_Domains_vs_Methylation_in_3_Cells/CpG-Counts/${file/1.profile.cg_strand_combined.bedgraph/NONGENIC.LARGE-K36me2_CpG_counts.dat};done



####### NEXT ###################3
# once done, go to the cpg count directory and run this
for file in *;do echo "bedtools intersect -a $file -b ../../Raw_Methylation_bedgraphs_3_Cells/${file/_BS_*/_BS_1.profile.cg_strand_combined.bedgraph} -wa -wb | groupBy -g 1,2,3,4,5,6,7 -c 11 -o sum | awk '{print \$0\"\t\"\$8/\$7}' > ../CpG_counts_sum_average/${file/-K36me2_CpG_counts.dat/-Parental-K36me2-domain_CpG_counts_sum_average.dat}" | bash;done

#The results will be here:
#main/CpG_counts_sum_average








####################################################################################################
## C) PLOTTING
####################################################################################################





####################################################################################################
## EXTRA) NOT NECESSARY BUT FOR FUTURE REFERENCE
## To segment methylation
####################################################################################################

# First, before running the segmentation script:
intersectBed -a /media/behnam/Black_Seagate2/Mouse/Data/mm10_1kb_intervals.bed -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -c > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1kb_CpG_counts.dat
intersectBed -a /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1kb_CpG_counts.dat -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -wa -wb | groupBy -g 1,2,3,4 -c 8 -o sum | awk '{print $0"\t"$5/$4}' > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_Methylation_CpGcount_sum_average.bedgraph

#Then run the Segmentation_Methylation.R on the output


