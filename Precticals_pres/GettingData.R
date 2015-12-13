GettingData
========================================================
  
  As we saw there are many places to grab data from.
========================================================
  
  Some of the most popular are GEO, SRA and ENA.

Interfacing with Repositories in R.
========================================================
  Lets load the packages we need, the repo should be obvious.

```{r}
library(GEOquery)
library(SRAdb)
```

GEOquery
========================================================
  So lets go look for some samples we have already downloaded
We are using Encode data so we find some GEOinforma tion [here](http://hgdownload.cse.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeSydhTfbs/files.txt)

One we are interest in is GSM912906

```{r, echo=T,cahce=T}
gds <- getGEO("GSM11805")
```

Have a look at gds object

GEOquery
========================================================
  So lets go look for some samples we have already downloaded
We are using Encode data so we find some GEOinforma tion [here](http://hgdownload.cse.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeSydhTfbs/files.txt)

One we are interest in is GSM912906

```{r, echo=T,cahce=T}
gds <- getGEO("GSM912906")
```

Have a look at gds object


GEOquery
========================================================
  GEO has much of the meta data available try meta(gds)

```{r, echo=T,cahce=T}
Meta(gds)$supplementary_file_1
Meta(gds)$supplementary_file_2
Meta(gds)$supplementary_file_3
```
With the third file pointing to SRA

SRAdb
========================================================
  SRAdb acts an interface between SRA and ENA

```{r, echo=T,cache=T}
fileToGen <- basename(Meta(gds)$supplementary_file_3) 
fileToGen
```
SRAdb
========================================================
  The first step in querying SRAdb is to download their schema.
```{r, echo=T,cahce=T}
sqlfile <<- getSRAdbFile()

Running across an experiment
=========================================================
  Lets set up a ChIPQC evaluation across all the experiments post-peak calling

```{r,echo=FALSE,eval=F,cache=T}
SampleID <- c("myc_ch12_1","myc_ch12_2","myc_Mel_1","myc_Mel_2")
Tissue <- c(rep("ch12",2),rep("mel",2))
Treatment <- rep(NA,4)
Condition <- rep(NA,4)
Factor <- c(rep("myc",2),rep("myc",2))
Replicate <- c(1,2,1,2)
bamReads  <- indexedBams[-c(3,4)]
bamControl= c("/Users/tcarroll/chipseqcourse/chipseqDataFQ//wgEncodeSydhTfbsCh12InputIggmusRawDatasorted.bam.bam",
              "/Users/tcarroll/chipseqcourse/chipseqDataFQ//wgEncodeSydhTfbsCh12InputIggmusRawDatasorted.bam.bam",
              "/Users/tcarroll/chipseqcourse/chipseqDataFQ//wgEncodeSydhTfbsMelInputIggmusRawDatasorted.bam.bam",
              "/Users/tcarroll/chipseqcourse/chipseqDataFQ//wgEncodeSydhTfbsMelInputIggmusRawDatasorted.bam.bam")
ControlID = c("ch12","ch12","ch12","ch12")
Peaks <- dir("/Users/tcarroll/chipseqcourse/MacsPeaks/",full.names=T)
PeakCaller <- rep("macs",4)
ss <- data.frame(SampleID,Tissue,Factor,Treatment,Replicate,Condition,bamReads,bamControl,ControlID,Peaks,PeakCaller)
ss2 <- merge(ss,temp,by=1,all.x=T,all.y=F)
```

```