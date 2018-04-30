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

# Next, running bedtools multicov on all the bam files at the same time using the domains and subdomains (genic/nongenic) as the input bed file

