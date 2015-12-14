Working With Coverage
========================================================
author: MRC Clinical Sciences Centre
width: 1440
height: 1100
autosize: true
font-import: <link href='http://fonts.googleapis.com/css?family=Slabo+27px' rel='stylesheet' type='text/css'>
font-family: 'Slabo 27px', serif;
css:style.css




Coverting BAM files to Coverage
========================================================
Coverage provides a useful genome wide picture of signal.


Bam files can easily be converted to coverage like objects with the coverage() function.
This by default produces a Run Length Encoding of signal over the genome.

========================================================

```r
library(rtracklayer)
```

So we use rtracklayer package to work with coverage.






```r
covExample <- coverage("/home/ubuntu//chipseqcourseData/sortedbams/wgEncodeSydhTfbsCh12CmycIggrabRawDataRep1.bam")
```

The RLElist object
========================================================
 The RLElist stores the runs of scores across the genomes

```r
length(covExample)
```

```
[1] 22
```

```r
covExample[[1]]
```

```
integer-Rle of length 129993255 with 1465520 runs
  Lengths: 3000000      36      25      36 ...       3       3     813
  Values :       0       1       0       2 ...       3       2       0
```

Arithmetric Operations on RLElist objects
========================================================
Arithmetric operations and running window operations are easily performed

```r
cov10 <- covExample+10
smoothedcov <- runsum(covExample[[10]],50)

cov10[1:10]
RleList of length 10
$`10`
numeric-Rle of length 129993255 with 1465520 runs
  Lengths: 3000000      36      25      36 ...       3       3     813
  Values :      10      11      10      12 ...      13      12      10

$`11`
numeric-Rle of length 121843856 with 2266855 runs
  Lengths: 3000063       1       2       1 ...       1      20     158
  Values :      10      11      12      13 ...      12      11      10

$`12`
numeric-Rle of length 121257530 with 1219909 runs
  Lengths: 3000020      36      33      36 ...      11      36      37
  Values :      10      11      10      11 ...      10      11      10

$`13`
numeric-Rle of length 120284312 with 1352336 runs
  Lengths: 3000325      36      55       1 ...       6      30      42
  Values :      10      11      10      11 ...      12      11      10

$`14`
numeric-Rle of length 125194864 with 1251007 runs
  Lengths: 3001566      27     726       8 ...      17      36    2724
  Values :      10      11      10      11 ...      10      11      10

...
<5 more elements>

smoothedcov[1:10]
integer-Rle of length 10 with 1 run
  Lengths: 10
  Values :  0
```



RLElist functions
========================================================
 The RLElist stores the runs of scores across the genomes. Here calculate base pairs above a certain height.

```r
tableOfDepths <- table(covExample[[1]])
sum(tableOfDepths[names(tableOfDepths) > 50])
```

```
[1] 142192
```

RLElist functions
========================================================
We can also get summary statistics

```r
sdCov <- sd(covExample[[1]])
meanCov <- mean(covExample[[1]])

sdCov
```

```
[1] 0.7877756
```

```r
meanCov
```

```
[1] 0.2263528
```


RLElist functions
========================================================
We can use slice() function to select islands above mean + 1SD

```r
myViews <- slice(covExample[[1]],lower=meanCov+sdCov)
```

Views functions
========================================================
Views objects hold the position of slice and scores along slice

```r
myViews[1:4]
```

```
Views on a 129993255-length Rle subject

views:
      start     end width
[1] 3000062 3000097    36 [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 ...]
[2] 3000601 3000636    36 [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 ...]
[3] 3001849 3001884    36 [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 ...]
[4] 3004224 3004236    13 [2 2 2 2 2 2 2 2 2 2 2 2 2]
```
Creating views from IRanges
========================================================

```r
exampleIR <- IRanges(start=c(3004224,8004225),end=c(3004236,8004237))
myNewViews <- Views(covExample[[2]],exampleIR)
```
========================================================

```r
myNewViews
```

```
Views on a 121843856-length Rle subject

views:
      start     end width
[1] 3004224 3004236    13 [5 5 5 5 5 5 5 5 5 5 5 5 5]
[2] 8004225 8004237    13 [0 0 0 0 0 0 0 0 0 0 0 0 0]
```

Views can be very useful in summarising score across a region.
========================================================

```r
viewMeans(myViews)[1:4]
```

```
[1] 2 2 2 2
```

```r
viewSums(myViews)[1:4]
```

```
[1] 72 72 72 26
```

Views also have functions to extract the value and position of max or min signal in a region.

```r
unique(viewMins(myViews))
```

```
[1]  2  3  4  6  5  9  8  7 11
```

```r
unique(viewMaxs(myViews))
```

```
 [1]   2   3   4   5   6  10   7   8  13   9  11  16  20  22  18  41  12
[18]  24  61  15  19  17  23  48  21  14  44  27  25  35  39  95 292  56
[35]  30 807  33  32  31 142  43  28  29 148  37  89  26  34 149  71  51
[52] 103  64  90  58  38  36
```

Exporting bigWig
========================================================
Finally we can export our coverage RLElist as a bigWig for visualising in IGV.

```r
export.bw(covExample, "ch12Myc.bw")
```

Typically we would want to extend fragments first. 

You could do this with resize on the GRanges converted GAlignments and passing the GRanges to export.bw.

Example extension of fragments
=========================================================

We know from ChIP-qc earlier the predicted fragment length



Out of bounds reads should trimmed automatically



```r
covExample <- readGAlignments("/home/ubuntu//chipseqcourseData/sortedbams/wgEncodeSydhTfbsCh12CmycIggrabRawDataRep1.bam")
covExample <- GRanges(covExample)
extended <- resize(covExample,165,fix="start")
export.bw(extended,con="Extended_ch12Myc.bw")
```

sessioninfo
============

```r
sessionInfo()
```

```
R version 3.2.2 (2015-08-14)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.11 (El Capitan)

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
[1] rtracklayer_1.30.1   GenomicRanges_1.22.1 GenomeInfoDb_1.6.1  
[4] IRanges_2.4.4        S4Vectors_0.8.3      BiocGenerics_0.16.1 
[7] knitr_1.11          

loaded via a namespace (and not attached):
 [1] XVector_0.10.0             magrittr_1.5              
 [3] zlibbioc_1.16.0            GenomicAlignments_1.6.1   
 [5] BiocParallel_1.4.0         stringr_1.0.0             
 [7] tools_3.2.2                SummarizedExperiment_1.0.1
 [9] Biobase_2.30.0             lambda.r_1.1.7            
[11] futile.logger_1.4.1        digest_0.6.8              
[13] formatR_1.2.1              futile.options_1.0.0      
[15] bitops_1.0-6               RCurl_1.95-4.7            
[17] evaluate_0.8               stringi_1.0-1             
[19] Biostrings_2.38.2          Rsamtools_1.22.0          
[21] XML_3.98-1.3              
```
