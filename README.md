# nonGenicDomains
Creating domains of broad histone marks in nonGenic regions

#### First I tried different scaled for segmentation. From segments of haboring minimum two windows of 1kb to 200 windows. Below is a screenshot of the results on IGV.
(You can click on the image for a higher resolution picture)

You can also download the bedgraph files from <a href="http://nikleotide.com/wp-content/uploads//2018/04/different_scales.zip">here</a> and the bigwig file for this samples and its input can be found on graham, here: <b> /project/6007495/nikbakht/nonGenic/ </b>. <p>

![alt text](http://nikleotide.com/wp-content/uploads//2018/04/igv_snapshot.png)

Out of different scales, it seems 50kb to be a suitable threshold for large domains. I used the 50kb segmentation and intersect it with the nongenic regions actoss mm10 genome. The results were filtered for segments (now all intergenic) that their segmentation value is > 0.5 (i.e. average of normalized read counts over the 1kb windows included in each segment). Below is a screenshot of how its bedgraph looks like. 

![alt text](http://nikleotide.com/wp-content/uploads//2018/04/igv_snapshot2.png)

However, it seems it still
