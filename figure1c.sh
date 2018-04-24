# For the sake of figure 1c, we don't need the trailing domains. Only domains that are very obvious. Thus I only take the > 1 values.
####################################################################################################
## 1kbWindow > Segment50w > Filter > Merge > Intersect
#NOTE: In order to call the trailing streteches of high peaks as separate entities, I use to different threshold for filtering the segments (based on their values). One of them, I used threshold of 1.0 and the other one 0.5-1.0 and then merge each subset separately. Next I will concatenate the resulting segments and intersect them with genic/nongenic

awk '{if ($4 >= 1.0) print $0"\t"($3-$2)/1000}' 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/4252_4295_C3H10T1_Parental_H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.fixed.bedgraph > 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/4252_4295_C3H10T1_Parental_H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.fixed.filtered_VALUE1.0.bedgraph
#   segments

# Next is merging them
bedtools merge -d 10 -i 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/4252_4295_C3H10T1_Parental_H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.fixed.filtered_VALUE1.0.bedgraph -c 4,5,1 -o mean,sum,count > 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/Parental_H3K36me2-50kb-H3K36me2-FILTERED_VALUE1.0-MERGED-ALL.bedgraph
#  merged segments including both genic and nongenic

awk '{if ($4 >= 0.5 && $4 < 1) print $0"\t"($3-$2)/1000}' 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/4252_4295_C3H10T1_Parental_H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.fixed.bedgraph > 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/4252_4295_C3H10T1_Parental_H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.fixed.filtered_VALUE0.5-1.bedgraph
#   segments

# Next is merging them
bedtools merge -d 10 -i 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/4252_4295_C3H10T1_Parental_H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.fixed.filtered_VALUE0.5-1.bedgraph -c 4,5,1 -o mean,sum,count > 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/Parental_H3K36me2-50kb-H3K36me2-FILTERED_fixed.filtered_VALUE0.5-1.ALL.bedgraph

#Now,to concatenate and sort them
cat 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/Parental_H3K36me2-50kb-H3K36me2-FILTERED_VALUE1.0-MERGED-ALL.bedgraph 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/Parental_H3K36me2-50kb-H3K36me2-FILTERED_fixed.filtered_VALUE0.5-1.ALL.bedgraph | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed -e 's/^chr//' | sort -k1,1n -k2,2n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/Parental_H3K36me2-50kb-H3K36me2-FILTERED_fixed.filtered_VALUES0.5AND0.5-1.ALL.bedgraph


#Now intersecting them:
bedtools intersect -b NonGenic/mm10hgTables_NonGenic.bed -a 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/Parental_H3K36me2-50kb-H3K36me2-FILTERED_fixed.filtered_VALUES0.5AND0.5-1.ALL.bedgraph |  cut -f 1-4 | awk '{print $0"\t"($3-$2)/1000}' | awk '{if ($4 > 0.5) print $0}' | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/chr//' | sort -k1,1n -k2,2n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/Parental_H3K36me2-50kb-H3K36me2-FILTERED_fixed.filtered_VALUES0.5AND0.5-1.NONGENIC.bedgraph
#and for genic
bedtools intersect -b NonGenic/mm10hgTables_Genic.bed -a 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/Parental_H3K36me2-50kb-H3K36me2-FILTERED_fixed.filtered_VALUES0.5AND0.5-1.ALL.bedgraph |  cut -f 1-4 | awk '{print $0"\t"($3-$2)/1000}' | awk '{if ($4 > 0.5) print $0}' | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/chr//' | sort -k1,1n -k2,2n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/Parental_H3K36me2-50kb-H3K36me2-FILTERED_fixed.filtered_VALUES0.5AND0.5-1.GENIC.bedgraph

#

