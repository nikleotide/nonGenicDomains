# results from 50kb minimum
# to merge neighboring segments of high values, first we have to filter the genic regions out and also filter out those that are smaller than 100kb and have low values of mark
 bedtools intersect -a mm10hgTables_NonGenic.bed -b ../4252_4295_C3H10T1_Parental_H3K36me2/Segmentation_50/4252_4295_C3H10T1_Parental_H3K36me2-50kb-H3K36me2-1stComplete_Genome-peaks.fixed.bedgraph -wb | cut -f 1-3,7 | awk '{print $0"\t"($3-$2)/1000}' | awk '{if ($4 > 0.5 && $5 > 50) print $0}' | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/chr//' | sort -k1,1n -k2,2n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > temp0.bedgraph 
 # next using bedtools merge, we merge the neighboring segments (-d 10 is okay since we don't have any segment smaller than 50k anyway)
bedtools merge -d 10 -i temp0.bedgraph -c 4,5,1 -o mean,sum,count > Parental_H3K36me2-50kb-H3K36me2-NonGenic.bedgraph



#
