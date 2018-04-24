#!/bin/bash
# For the sake of figure 1c, we don't need the trailing domains. Only domains that are very obvious. Thus I only take the > 1 values.
####################################################################################################
## 1kbWindow > Segment50w > Filter > Merge > Intersect
#NOTE: In order to call the trailing streteches of high peaks as separate entities, I use to different threshold for filtering the segments (based on their values). One of them, I used threshold of 1.0 and the other one 0.5-1.0 and then merge each subset separately. Next I will concatenate the resulting segments and intersect them with genic/nongenic


for file in test/Segmentation_50/*50kb-*fixed.bedgraph
do
echo 
awk '{if ($4 >= 1.0) print $0"\t"($3-$2)/1000}' $file > ${file/.fixed/.fixed.filtered_VALUE1.0.bedgraph}
echo
bedtools merge -d 10 -i ${file/.fixed/.fixed.filtered_VALUE1.0.bedgraph} -c 4,5,1 -o mean,sum,count > ${file/fixed.bedgraph/-FILTERED_VALUE1.0-MERGED-ALL.bedgraph}
echo
cat ${file/fixed.bedgraph/-FILTERED_VALUE1.0-MERGED-ALL.bedgraph} | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed -e 's/^chr//' | sort -k1,1n -k2,2n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.ALL.bedgraph}
echo
bedtools intersect -b NonGenic/mm10hgTables_NonGenic.bed -a ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.ALL.bedgraph} |  cut -f 1-4 | awk '{print $0"\t"($3-$2)/1000}' | awk '{if ($4 > 0.5) print $0}' | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/chr//' | sort -k1,1n -k2,2n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.NONGENIC.bedgraph}
echo
bedtools intersect -b NonGenic/mm10hgTables_Genic.bed -a ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.ALL.bedgraph} |  cut -f 1-4 | awk '{print $0"\t"($3-$2)/1000}' | awk '{if ($4 > 0.5) print $0}' | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/chr//' | sort -k1,1n -k2,2n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > ${file/fixed.bedgraph/-FILTERED_fixed.filtered_VALUES1.GENIC.bedgraph}
done
