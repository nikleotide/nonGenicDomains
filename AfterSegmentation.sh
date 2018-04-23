## 1kbWindow > Segment50w > Intersect > Filter > Merge 
# results from 50kb minimum
# to merge neighboring segments of high values, first we have to filter the genic regions out and also filter out those that are smaller than 100kb and have low values of mark
bedtools intersect -a mm10hgTables_NonGenic.bed -b 4252_4295_C3H10T1_Parental_H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.fixed.bedgraph -wb | cut -f 1-3,7 | awk '{print $0"\t"($3-$2)/1000}' | awk '{if ($4 > 0.5 && $5 > 50) print $0}' | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/chr//' | sort -k1,1n -k2,2n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > Parental_H3K36me2-50kb-H3K36me2-NonGenic-NOT-MERGED.bedgraph 
# next using bedtools merge, we merge the neighboring segments (-d 10 is okay since we don't have any segment smaller than 50k anyway)
bedtools merge -d 10 -i Parental_H3K36me2-50kb-H3K36me2-NonGenic-NOT-MERGED.bedgraph -c 4,5,1 -o mean,sum,count > Parental_H3K36me2-50kb-H3K36me2-NonGenic.bedgraph



####################################################################################################
## 1kbWindow > Segment50w > Filter > Merge > Intersect

## In case of merging before intersectioning with genic/intergenic bed files
## First adding the length of each segment and filtering out those below 0.5 and smaller than 50 (trailing ones)
awk '{if (($3-$2)/1000 > 50 && $4 > 0.5) print $0"\t"($3-$2)/1000}' ../4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/4252_4295_C3H10T1_Parental_H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.fixed.bedgraph > ../4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/4252_4295_C3H10T1_Parental_H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.fixed-filtered-LENGTH.bedgraph
#  9647 segments

# Next is merging them
bedtools merge -d 10 -i ../4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/4252_4295_C3H10T1_Parental_H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.fixed-filtered-LENGTH.bedgraph -c 4,5,1 -o mean,sum,count > Parental_H3K36me2-50kb-H3K36me2-MERGED-ALL.bedgraph
# 2483 merged segments including both genic and nongenic

# Now time for intersecting with Genic/NonGenic
# For NonGenic
behnam:Figure1C$ bedtools intersect -b NonGenic/mm10hgTables_NonGenic.bed -a 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/Parental_H3K36me2-50kb-H3K36me2-FILTERED_VALUE-MERGED-ALL.bedgraph |  cut -f 1-4 | awk '{print $0"\t"($3-$2)/1000}' | awk '{if ($4 > 0.5) print $0}' | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/chr//' | sort -k1,1n -k2,2n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/Parental_H3K36me2-50kb-H3K36me2-FILTERED_VALUE-MERGED-NONGENIC.bedgraph
# And for Genic
bedtools intersect -b NonGenic/mm10hgTables_Genic.bed -a 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/Parental_H3K36me2-50kb-H3K36me2-FILTERED_VALUE-MERGED-ALL.bedgraph |  cut -f 1-4 | awk '{print $0"\t"($3-$2)/1000}' | awk '{if ($4 > 0.5) print $0}' | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/chr//' | sort -k1,1n -k2,2n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > 4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/Parental_H3K36me2-50kb-H3K36me2-FILTERED_VALUE-MERGED-GENIC.bedgraph

#

