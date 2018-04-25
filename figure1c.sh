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
awk '{if ($3-$2 > 500000) print $0}' ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.ALL.bedgraph} > ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.ALL.LARGE.bedgraph}
echo

bedtools intersect -b /media/behnam/Black_Seagate2/Mouse/Data/mm10hgTables_NonGenic.bed -a ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.ALL.LARGE.bedgraph} |  cut -f 1-4 | awk '{print $0"\t"($3-$2)/1000}' | awk '{if ($4 > 0.5) print $0}' | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/chr//' | sort -k1,1n -k2,2n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.NONGENIC.LARGE.bedgraph}
echo
bedtools intersect -b /media/behnam/Black_Seagate2/Mouse/Data/mm10hgTables_Genic.bed -a ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.ALL.LARGE.bedgraph} |  cut -f 1-4 | awk '{print $0"\t"($3-$2)/1000}' | awk '{if ($4 > 0.5) print $0}' | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/chr//' | sort -k1,1n -k2,2n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.GENIC.LARGE.bedgraph}

rm ${file/fixed.bedgraph/-FILTERED_VALUE1.0-MERGED-ALL.bedgraph}
rm ${file/.fixed/.fixed.filtered_VALUE1.0.bedgraph}
done


####################################################################################################
## B) Calculating the average methylation in the large, large-nongenic and large-genic domain 
####################################################################################################

## ALL Large
bedtools intersect -a H3K36me2/Segmentation_50/H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.ALL.LARGE.bedgraph -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -c > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_ALL.LARGE-K36me2_CpG_counts.dat
intersectBed -a /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_ALL.LARGE-K36me2_CpG_counts.dat -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -wa -wb | groupBy -g 1,2,3,4,5,6,7 -c 11 -o sum | awk '{print $0"\t"$8/$7}' > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_ALL.LARGE-K36me2_CpG_counts_sum_average.dat

## Large-nongenic
bedtools intersect -a H3K36me2/Segmentation_50/H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.NONGENIC.LARGE.bedgraph -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -c > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_NONGENIC.LARGE-K36me2_CpG_counts.dat
intersectBed -a /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_NONGENIC.LARGE-K36me2_CpG_counts.dat -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -wa -wb | groupBy -g 1,2,3,4,5,6,7 -c 11 -o sum | awk '{print $0"\t"$8/$7}' > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_NONGENIC.LARGE-K36me2_CpG_counts_sum_average.dat

# Large-genic
bedtools intersect -a H3K36me2/Segmentation_50/H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.-FILTERED_fixed.filtered_VALUES1.GENIC.LARGE.bedgraph -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -c > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_GENIC.LARGE-K36me2_CpG_counts.dat
intersectBed -a /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_GENIC.LARGE-K36me2_CpG_counts.dat -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -wa -wb | groupBy -g 1,2,3,4,5,6,7 -c 11 -o sum | awk '{print $0"\t"$8/$7}' > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_GENIC.LARGE-K36me2_CpG_counts_sum_average.dat


####################################################################################################
## C) NOT NECESSARY BUT FOR FUTURE REFERENCE
## To segment methylation
####################################################################################################

# First, before running the segmentation script:
intersectBed -a /media/behnam/Black_Seagate2/Mouse/Data/mm10_1kb_intervals.bed -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -c > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1kb_CpG_counts.dat
intersectBed -a /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1kb_CpG_counts.dat -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -wa -wb | groupBy -g 1,2,3,4 -c 8 -o sum | awk '{print $0"\t"$5/$4}' > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/Data/Only_Parental/Methylation/Parental_C3H10T_BS_Methylation_CpGcount_sum_average.bedgraph

#Then run the Segmentation_Methylation.R on the output


