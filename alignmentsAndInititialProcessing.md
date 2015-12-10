Alignments and initial QC
========================================================
author: MRC Bioinformatics Core
date: 

First Steps in Alignment
========================================================

We have already looked at FastQ files. [FASTQ Description](http://mrccsc.github.io/genomicFormats.html#/8)

Fastq is one of the more basic datatypes representing sequence information from your fragment and the quality associated to it

You can find some examples in the directory.
"/Users/tcarroll/chipseqcourse/chipseqData"

First lets load the libraries we will need for this practical


```r
library(BiocParallel)
library(ShortRead)
library(Rsubread)
library(Rsamtools)
library(QuasR)
```



A quick peak at FastQ
========================================================


```r
pathToFQ <- "/Users/tcarroll/chipseqcourse/chipseqDataFQ/"
fastQ <- "wgEncodeSydhTfbsCh12CmycIggrabRawDataRep1.fastq.gz"
sampleOfFastQ <- FastqSampler(file.path(pathToFQ,fastQ), 10)
fq <- yield(sampleOfFastQ)
fq
```

```
class: ShortReadQ
length: 10 reads; width: 36 cycles
```

ShortReadQ class
========================================================
The ShortReadQ class hold information on each FASTQ as well as general statistics.


```
  A DNAStringSet instance of length 1
    width seq
[1]    36 TTGCAAGCTTTCTGTTGTTTTTGTTGTTGTTGTTGA
```

```
class: SFastqQuality
quality:
  A BStringSet instance of length 1
    width seq
[1]    36 a_GUHU_WaUa_T[UQ_T`_aQ`PXXLLFT`[][BB
```


ShortReadQ Quality
========================================================
The ShortReadQ class also helps in predicting encoding


```
 ;  <  =  >  ?  @  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S 
-5 -4 -3 -2 -1  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 
 T  U  V  W  X  Y  Z  [ \\  ]  ^  _  `  a  b  c  d  e  f  g  h  i 
20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 
```

```
[1] "Quality is phred64 Scale"
```
Phred 64 information will be used with the aligner.

Starting with alignments.
=========================================================

For this course we will be using rsubread but first we need a FASTA file.

We saw FASTA previously in another [course](http://mrccsc.github.io/genomicFormats.html#/6).

An aligner specific index must be created for this FASTA file

```r
ref = "mm9_Fasta.fa"
buildindex(basename="PATH_TO_INDEX/mm9_Index",
           reference=ref)
```

Starting with alignments. (The aligning)
=========================================================

Here we simply need to point the subread aligner at the reference genome and fastq files to align too. It doesnt matter if they are zipped or not.
Lets just do one for brevity.

We will use type 1 (for DNA) and stop the search for complex indels

```r
fastqToAlign <-dir("/Users/tcarroll/chipseqcourse/chipseqDataFQ/",
                   pattern="*.fastq.gz$",full.names=T)
index = "/Users/tcarroll/chipseqcourse/referenceData//mm9_index"
# Get our number of cores
BiocParallel::register(BiocParallel::MulticoreParam(20))


for(i in 1:length(fastqToAlign)){
align(index,fastqToAlign[i],nthreads=2,type=1,complexIndels=FALSE)
message("file", i, "aligned")
}  
```
Working with alignments. (Sorting)
=========================================================

Bam files need to be sorted in order to be used by most programs and properly indexed

```r
BiocParallel::register(BiocParallel::MulticoreParam(8))
toSort <- dir("/Users/tcarroll/chipseqcourse/chipseqDataBAM/",pattern="*.BAM$",full.names=T)
bplapply(toSort,function(x) sortBam(file=x,destination=gsub("\\.BAM","sorted\\.BAM",x),maxMemory=1024))
```



Working with alignments. (Indexing)
=========================================================

Indexing a file allows for random access of only the parts of the files of interest

```r
BiocParallel::register(BiocParallel::MulticoreParam(8))
toIndex <- dir("/Users/tcarroll/chipseqcourse/chipseqDataBAM/",pattern="*subreadsorted\\.BAM\\.bam$",full.names=T)
bplapply(toIndex,function(x) indexBam(file=x))
```


Working with alignments. (Mapping Rates)
=========================================================

```r
BiocParallel::register(BiocParallel::MulticoreParam(8))
indexedBams <- dir("/Users/tcarroll/chipseqcourse/chipseqDataBAM/",pattern="*subreadsorted\\.BAM\\.bam$",full.names=T)
stats <- bplapply(indexedBams,function(x) alignmentStats(x))
stats[[5]]
```


Working with BAM files
========================



Inspecting the BAM header.
===========================
The BAM header contains information on contigs aligned to and potentially programs used.









