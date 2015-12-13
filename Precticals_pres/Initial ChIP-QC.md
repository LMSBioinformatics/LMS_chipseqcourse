Initial ChIP-QC
========================================================
author: MRC Clinical Sciences Centre
date:http://mrccsc.github.io/r_course/introToR_Session1.html
width: 1440
height: 1100
autosize: true
font-import: <link href='http://fonts.googleapis.com/css?family=Slabo+27px' rel='stylesheet' type='text/css'>
font-family: 'Slabo 27px', serif;
css:style.css

ChIP-QC
========================================================

ChIP-seq data is very noisy and can be variable in specificity and sensitivity depending on the ChIP antibody used.

Several tools exist in R to allow us to evaluate our ChIPQC prior to peak calling and post visualisation in a browser.

ChIPQC package
========================================================

The ChIPQC package was developed to wrap up most of the useful metrics and place them in the right context.

Given unfiltered data it processes the data and gathers the QC metrics at the appropriate filtering steps.

Loading packages
===============
First lets load the packages we need


```r
library(DiffBind)
library(ChIPQC)
library(GenomeInfoDb)
library(BiocParallel)
```



ChIPQCsample
========================================================

The most basic function in ChIPQC is ChIPQCsample() which will gather the information on a single sample.

Just a BAM file name, genome for annotation and any areas of blacklisted signal.



ChIPQC accecpts a BAM file (with or without index). 
It is also recommended to provide a genome to test for enrichment in genomic features and a Blacklist.


```r
indexedBams[1]
```

```
[1] "/Users/tcarroll/chipseqcourse/chipseqDataFQ//wgEncodeSydhTfbsCh12CmycIggrabRawDataRep1sorted.bam.bam"
```

```r
myQC <- ChIPQCsample(indexedBams[1],
             annotation="mm9",
             blacklist="/Users/tcarroll/chipseqcourse//referenceData/mm9-blacklist.bed",chromosomes=paste0(1:2))
```

```


[1] 1
[1] 1
```


Strange Results for Reads in Blacklist and features.
========================================================

```r
myQC
```

```
                        ProportionOfCounts
BlackList                                0
LongPromoter20000to2000                  0
Promoters2000to500                       0
Promoters500                             0
All5utrs                                 0
Alltranscripts                           0
Allcds                                   0
Allintrons                               0
All3utrs                                 0
GRanges object with 0 ranges and 0 metadata columns:
   seqnames    ranges strand
      <Rle> <IRanges>  <Rle>
  -------
  seqinfo: no sequences
```

Fixing differences in Genome Contig Annotation.
========================================================
Lets check names of contigs in our BAMS


```r
myBam <- BamFile(indexedBams[1])
names(scanBamHeader(myBam)$targets)
```

```
 [1] "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "1"  "2"  "3"  "4" 
[15] "5"  "6"  "7"  "8"  "9"  "MT" "X"  "Y" 
```

Fixing differences in Genome Contig Annotation
========================================================
The names in our files start with a chr so we will have to fix these using the **renameSeqlevels** function


```r
myBL <- ChIPQC:::GetGRanges("/Users/tcarroll/chipseqcourse//referenceData/mm9-blacklist.bed")
myBL <- renameSeqlevels(myBL,gsub("chr","",seqlevels(myBL)))
```

We will also need to the same for gene annotation
========================================================
We can extract,alter and provide custom feature annotation for ChIPQC


```r
mm9Anno <- ChIPQC:::getAnnotation("mm9",AllChr=NULL)
mm9AnnoNew <- lapply(mm9Anno[-1],
                  function(x)
                  renameSeqlevels(x,gsub("chr","",seqlevels(x))
                       )
              )
```

Lets now provide fixed annotation
========================================================
We can now try with fixed annotation




```r
myQC <- ChIPQCsample(indexedBams[1],
             annotation=mm9AnnoNew,
             blacklist = myBL,
             chromosomes=paste0(1:10)
             )
```

These results look better
========================================================

```r
load("../robjects/singleSampleQC.RData")
```


```r
myQC
```

```
                        ProportionOfCounts
BlackList                      0.138962494
LongPromoter20000to2000        0.164999877
Promoters2000to500             0.027274214
Promoters500                   0.021972120
All5utrs                       0.009702207
Alltranscripts                 0.421248126
Allcds                         0.017216495
Allintrons                     0.393648791
All3utrs                       0.013905895
GRanges object with 0 ranges and 0 metadata columns:
   seqnames    ranges strand
      <Rle> <IRanges>  <Rle>
  -------
  seqinfo: no sequences
```

Some useful plots
=========================================================
In the CC plot we see fragment length and artefact peak. Note in this samnple the artefact peak is higher then fragment length peak


```r
p <- plotCC(myQC)
p$layers[[2]] <- NULL
p
```

![plot of chunk unnamed-chunk-12](Initial ChIP-QC-figure/unnamed-chunk-12-1.png) 

Some more useful plots
=========================================================
By comparing SSD before and after Blacklisting we get an understanding of artefact and remaining signal.


```r
plotSSD(myQC)+xlim(0,10)
```

![plot of chunk unnamed-chunk-13](Initial ChIP-QC-figure/unnamed-chunk-13-1.png) 

Some more useful plots
=========================================================
By using Regi plot we can see where the distribution of signal in our features.


```r
plotRegi(myQC)
```

![plot of chunk unnamed-chunk-14](Initial ChIP-QC-figure/unnamed-chunk-14-1.png) 

Running across an experiment
=========================================================
Lets set up a ChIPQC evaluation across all the experiments




Plotting across an entire experiment
=========================================================
You should have an object loaded called ss. This Samplesheet is required to merge individual ChIPQCsample results.




```r
ss
```

```
    SampleID Tissue Factor Treatment Replicate Condition
1 myc_ch12_1   ch12    myc        NA         1        NA
2 myc_ch12_2   ch12    myc        NA         2        NA
3 input_ch12   ch12  input        NA         1        NA
4  myc_Mel_1    mel    myc        NA         1        NA
5  myc_Mel_2    mel    myc        NA         2        NA
6  input_Mel    mel  input        NA         1        NA
                                                                                              bamReads
1 /Users/tcarroll/chipseqcourse/chipseqDataFQ//wgEncodeSydhTfbsCh12CmycIggrabRawDataRep1sorted.bam.bam
2 /Users/tcarroll/chipseqcourse/chipseqDataFQ//wgEncodeSydhTfbsCh12CmycIggrabRawDataRep2sorted.bam.bam
3    /Users/tcarroll/chipseqcourse/chipseqDataFQ//wgEncodeSydhTfbsCh12InputIggmusRawDatasorted.bam.bam
4  /Users/tcarroll/chipseqcourse/chipseqDataFQ//wgEncodeSydhTfbsMelCmycIggrabRawDataRep1sorted.bam.bam
5  /Users/tcarroll/chipseqcourse/chipseqDataFQ//wgEncodeSydhTfbsMelCmycIggrabRawDataRep2sorted.bam.bam
6     /Users/tcarroll/chipseqcourse/chipseqDataFQ//wgEncodeSydhTfbsMelInputIggmusRawDatasorted.bam.bam
  bamControl ControlID Peaks
1       <NA>        NA  <NA>
2       <NA>        NA  <NA>
3       <NA>        NA  <NA>
4       <NA>        NA  <NA>
5       <NA>        NA  <NA>
6       <NA>        NA  <NA>
```

Plotting across an entire experiment
=========================================================
You should have an object loaded called ss


```r
myRes2 <- ChIPQC(ss,samples=myRes)
```


```r
p <- plotCC(myRes2)
p$layers[[2]] <- NULL
p
```

![plot of chunk unnamed-chunk-20](Initial ChIP-QC-figure/unnamed-chunk-20-1.png) 
     
We will come back to QC when we have some peaks.