# Here where I kept the bed files created from merging all the transcripts (genic) and subtracting it from the entire genome to get the non-genic regions
# /media/behnam/Black_Seagate2/Mouse/Data/mm10hgTables_Genic.bed
# /media/behnam/Black_Seagate2/Mouse/Data/mm10hgTables_NonGenic.bed

# We use these files to intersect with methylation and measure their methylation in each gene/non-gene region and compare their changes from cell to cell
### Repeat these steps below for all the cell lines


## Genic
# 1) Count the number of CpGs
bedtools intersect -a /media/behnam/Black_Seagate2/Mouse/Data/mm10hgTables_Genic.bed -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Raw_Methylation_bedgraphs_4_Cells/Parental_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -c > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_genic-nongenic-regions_in_4_Cells/CpG-Counts/Parental_C3H10T_BS_GENIC_CpG_counts.dat

# 2) Calculate the average methylation
bedtools intersect -a /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_genic-nongenic-regions_in_4_Cells/CpG-Counts/Parental_C3H10T_BS_GENIC_CpG_counts.dat -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Raw_Methylation_bedgraphs_4_Cells/Parental_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -wa -wb | groupBy -g 1,2,3,4 -c 8 -o sum | awk '{print $0"\t"$5/$4}' > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_genic-nongenic-regions_in_4_Cells/CpG_counts_sum_average/Parental_C3H10T_BS_GENIC_CpG_counts_sum_average.dat


## Non-Genic
# 1) Count the number of CpGs
bedtools intersect -a /media/behnam/Black_Seagate2/Mouse/Data/mm10hgTables_NonGenic.bed -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Raw_Methylation_bedgraphs_4_Cells/Parental_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -c > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_genic-nongenic-regions_in_4_Cells/CpG-Counts/Parental_C3H10T_BS_NONGENIC_CpG_counts.dat

# 2) Calculate the average methylation
bedtools intersect -a /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_genic-nongenic-regions_in_4_Cells/CpG-Counts/Parental_C3H10T_BS_NONGENIC_CpG_counts.dat -b /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_3_Parental_Domains_in_4_Cells/Raw_Methylation_bedgraphs_4_Cells/Parental_C3H10T_BS_1.profile.cg_strand_combined.bedgraph -wa -wb | groupBy -g 1,2,3,4 -c 8 -o sum | awk '{print $0"\t"$5/$4}' > /media/behnam/Black_Seagate2/Mouse/Final_Figure1C/ANALYSIS/Methylation_in_genic-nongenic-regions_in_4_Cells/CpG_counts_sum_average/Parental_C3H10T_BS_NONGENIC_CpG_counts_sum_average.dat
