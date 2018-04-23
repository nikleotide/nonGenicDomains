# nonGenicDomains
Creating domains of broad histone marks in nonGenic regions

# 1kbWindow > Segment50w > Intersect > Filter > Merge
## Merging after intersecing
First I tried different scaled for segmentation. From segments of haboring minimum two windows of 1kb to 200 windows. Below is a screenshot of the results on IGV.
(You can click on the image for a higher resolution picture)

You can also download the bedgraph files from <a href="http://nikleotide.com/wp-content/uploads//2018/04/different_scales.zip">here</a> and the bigwig file for this samples and its input can be found on graham, here: <b>  /project/6007495/nikbakht/nonGenic/ </b>. <br><br>

![alt text](http://nikleotide.com/wp-content/uploads//2018/04/igv_snapshot.png)

Out of different scales, it seems 50kb to be a suitable threshold for large domains. I used the 50kb segmentation and intersect it with the nongenic regions actoss mm10 genome. The results were filtered for segments (now all intergenic) that their segmentation value is > 0.5 (i.e. average of normalized read counts over the 1kb windows included in each segment). Below is a screenshot of how its bedgraph looks like. 

![alt text](http://nikleotide.com/wp-content/uploads//2018/04/igv_snapshot2.png)


However, it seems it still exist some segments neighboring each other which can be put together as a bigger segment since their values are very close. What I did was to merge the bed file containing the segment data and use the average of the value of segments merging together as their final value. Below is the same screenshot with the new merged bed file added. The smoothing (merging of neighboring segments) in this example is more clear on the far left of the screenshot.

![alt text](http://nikleotide.com/wp-content/uploads//2018/04/igv_snapshot3.png)

#### Out of 3698 segments, after merging there are 2696 larger segments left. After filtering for minimum 500kb size of the intergenic segments, there are 114 segments left.

You can download the segments and merged segments <a href="http://nikleotide.com/wp-content/uploads//2018/04/Parental_H3K36me2-50kb-H3K36me2-NonGenic-NOT-MERGED.bedgraph.zip">before merging</a> and <a href="http://nikleotide.com/wp-content/uploads//2018/04/Parental_H3K36me2-50kb-H3K36me2-NonGenic.bedgraph.zip">after merging</a>.

-----------------------------------------------------------------------------------------------------------------------
To replicate the results, please consult the BeforeSegmentation.sh, Segmentation.R and AfterSegmentation.sh scripts in that order.
<b>Note: The necessary bed files for mm10 can be downloaded from <a href="http://nikleotide.com/wp-content/uploads//2018/04/mm10.bed_.files_.zip">here</a>.
<br>
 The bam files also can be found here:
 /project/6007495/hchen009/Data/K36M/Simon/K36M/chip_seq/C3H10T1/
</b>

-----------------------------------------------------------------------------------------------------------------------
# 1kbWindow > Segment50w > Filter > Merge > Intersect
## Merging the segments before intersecting them with Genic/NonGenic regions
<br>
The bedgraph files for this igv session can be downloaded from <a href="http://nikleotide.com/wp-content/uploads//2018/04/Merged_then_Intersected.zip">here</a>
Below is a screenshot of the results. The tracks are ordered as:

1. Actual bw from histone mark 
2. Input bw track 
3. Only segmentation (no filtering, not merging, no intersection) 
4. Only segmentation plus filtering (if I don't don filtering all segments will merge together; i.e. we will have one big segment per chromosome at the end) 
5. Segmentation + filtering + merging 
6. Segmentation + filtering + merging + intersected with NonGenic regions 
7. Segmentation + filtering + merging + intersected with Genic regions 
<br>

![alt text](http://nikleotide.com/wp-content/uploads//2018/04/igv_snapshot-merged-intersected.png)

-----------------------------------------------------------------------------------------------------------------------
# 1kbWindow > Segment50w > Filter and subset for value  > Merge within subsets > Concatanation > Intersect
## In this approach, I will use different hight threshold values for merging. Then, concatante the merged results and sort them before intersecting them with Genic/NonGenic regions

This is to account for lower level trailing stretches of marks which are higher than zero but not as high as their neighboring mark domains.

Below is an example of how they look like:
1. Mark
2. Input
3. SICER results to compare (garbage results from sicer for this mark even though it worked okay with K27me3)
4. Resulting segments before intersecting with genic/nongenic 
5. Genic segments
6. Genic annotation from UCSC (mm10)
7. Intergenic segments
8. Intergenic annotation from UCSC (mm10)
<br>

![alt text](http://nikleotide.com/wp-content/uploads//2018/04/igv_snapshot4.png)







