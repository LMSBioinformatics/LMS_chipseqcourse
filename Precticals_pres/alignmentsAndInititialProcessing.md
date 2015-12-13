Alignments and initial QC
========================================================
author: MRC Clinical Sciences Centre
date:http://mrccsc.github.io/r_course/introToR_Session1.html
width: 1440
height: 1100
autosize: true
font-import: <link href='http://fonts.googleapis.com/css?family=Slabo+27px' rel='stylesheet' type='text/css'>
font-family: 'Slabo 27px', serif;
css:style.css
First Steps in Alignment
========================================================

We have already looked at FastQ files. [FASTQ Description](http://mrccsc.github.io/genomicFormats.html#/8)

Fastq is one of the more basic datatypes representing sequence information from your fragment and the quality associated to it

You can find some examples in the directory.
"/chipseqData/FQs"

First lets load the libraries we will need for this practical


```r
library(BiocParallel)
library(ShortRead)
library(Rsubread)
library(Rsamtools)
library(QuasR)
library(GenomicAlignments)
library(GenomicRanges)
```



A quick peak at FastQ
========================================================
Here we use **FastqSampler()** function to select 10 reads at random.
This can be useful for QC when you may want 1,000,000 reads instead of all reads.


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

Hre we can see the sequence and quality of the read using the **sread()** and **quality()** functions respectively


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
The ShortReadQ class also helps in predicting encoding.
You cqan always compare to [chart](http://en.wikipedia.org/wiki/FASTQ_format) to see which quality scores yours are.

R also offers a means to compare your encoding to that of a predefined set.


```r
encoding(quality(fq[1]))
```

```
 ;  <  =  >  ?  @  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S 
-5 -4 -3 -2 -1  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 
 T  U  V  W  X  Y  Z  [ \\  ]  ^  _  `  a  b  c  d  e  f  g  h  i 
20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 
```

```r
if(
all(encoding(quality(fq[1])) == encoding(SFastqQuality()))){
  print("Quality is phred64 Scale")
}
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

But this takes a long time so we wont bother here for time constraints.

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
```

Lets start with one sample
=========================================================

```r
fastqToAlign[1]
```

```
[1] "/Users/tcarroll/chipseqcourse/chipseqDataFQ//wgEncodeSydhTfbsCh12CmycIggrabRawDataRep1.fastq.gz"
```

```r
align(index,x,output_file = gsub("\\.fastq\\.gz$","\\.bam",x),nthreads=4,unique=F,type=1,complexIndels=FALSE,phredOffset=64)))
```

Important parameters
=========================================================
Take a look at **?align** for all the parameters available

Some of the most important here are:
- nthreads == _Number of threads/cpus to use in aligning_
- unique == _When TRUE will only report unique reads_
- phredOffset == *It is important to supply correct Phred scale (although you are warned by subread if it looks wrong)*

Now for all aligning all files in the experiment 
=========================================================

```r
BiocParallel::register(MulticoreParam(workers=3))

bplapply(fastqToAlign,function(x)
  capture.output(align(index,x,output_file = gsub("\\.fastq\\.gz$","\\.bam",x),
  nthreads=4,unique=F,type=1,complexIndels=FALSE,phredOffset=64)))
```
Working with alignments. (Sorting)
=========================================================

Bam files need to be sorted in order to be used by most programs and properly indexed
The **sortBam()** function will perform this with a user supplied max memory supplied in the **maxMemory** argument.


```r
toSort <- dir("/Users/tcarroll/chipseqcourse/chipseqDataFQ/",pattern="*.bam$",full.names=T)
sortBam(file=x,destination=gsub("\\.bam","sorted\\.bam",x),maxMemory=1024))
```
Now for all sorting all files in the experiment 
=========================================================

```r
BiocParallel::register(BiocParallel::MulticoreParam(8))
toSort <- dir("/Users/tcarroll/chipseqcourse/chipseqDataFQ/",pattern="*.bam$",full.names=T)
bplapply(toSort,function(x) sortBam(file=x,destination=gsub("\\.bam","sorted\\.bam",x),maxMemory=1024))
```


Working with alignments. (Indexing)
=========================================================

Indexing a file allows for random access of only the parts of the files of interest.
In R we can use **indexBam()** function. 

```r
toIndex <- dir("/Users/tcarroll/chipseqcourse/chipseqDataFQ/",pattern="*sorted\\.bam",full.names=T)
indexBam(file=toIndex))
```

Note that most programs need indexed data. **ChIPQC** will create it if it doesnt exists using same method

Indexing all files
=========================================================


```r
BiocParallel::register(BiocParallel::MulticoreParam(8))
toIndex <- dir("/Users/tcarroll/chipseqcourse/chipseqDataFQ/",pattern="*sorted\\.bam",full.names=T)
bplapply(toIndex,function(x) indexBam(file=x))
```


Working with alignments. (Mapping Rates)
=========================================================

The **alignmentStats()** function will provide simple statistics on mapping.


```r
indexedBams <- dir("/Users/tcarroll/chipseqcourse/chipseqDataFQ/",pattern="*sorted\\..*bam$",full.names=T)
stats <- alignmentStats(indexedBams[1])
stats
```

```
                                                         seqlength
wgEncodeSydhTfbsCh12CmycIggrabRawDataRep1sorted.bam.bam 2654911517
                                                          mapped unmapped
wgEncodeSydhTfbsCh12CmycIggrabRawDataRep1sorted.bam.bam 18144863  9471721
```

Working with alignments. (Mapping Rates)
=========================================================

```r
BiocParallel::register(BiocParallel::MulticoreParam(8))
indexedBams <- dir("/Users/tcarroll/chipseqcourse/chipseqDataFQ/",pattern="*sorted\\..*bam$",full.names=T)
stats <- bplapply(indexedBams,function(x) alignmentStats(x))
stats[[5]]
```

```
                                                        seqlength   mapped
wgEncodeSydhTfbsMelCmycIggrabRawDataRep2sorted.bam.bam 2654911517 17840276
                                                       unmapped
wgEncodeSydhTfbsMelCmycIggrabRawDataRep2sorted.bam.bam 10249532
```


Working with BAM files
========================



Retrieving global information with RSamtools.
===========================

The Rsamtools library is the main library to deal with BAM files.
Creating a **BamFile** object is easy and will allow finer control in the future over simply supplying strings of BAM file names.




```r
bamToRead <- BamFile(indexedBams[1])
bamToRead
```

```
class: BamFile 
path: /Users/tca.../wgEncodeSydhTfbsCh12CmycIggrabRawDataRep1sorted.bam.bam
index: /User.../wgEncodeSydhTfbsCh12CmycIggrabRawDataRep1sorted.bam.bam.bai
isOpen: FALSE 
yieldSize: NA 
obeyQname: FALSE 
asMates: FALSE 
qnamePrefixEnd: NA 
qnameSuffixStart: NA 
```


Bam Headers 
===========================
</br>

The BAM header contains information on contigs aligned to and potentially programs used.




```r
bamToRead <- BamFile(indexedBams[1])
targets <- scanBamHeader(bamToRead)$targets
```
***

|Contigs |   Lengths|
|:-------|---------:|
|10      | 129993255|
|11      | 121843856|
|12      | 121257530|
|13      | 120284312|
|14      | 125194864|
|15      | 103494974|
|16      |  98319150|
|17      |  95272651|
|18      |  90772031|
|19      |  61342430|
|1       | 197195432|
|2       | 181748087|
|3       | 159599783|
|4       | 155630120|
|5       | 152537259|
|6       | 149517037|
|7       | 152524553|
|8       | 131738871|
|9       | 124076172|
|MT      |     16299|
|X       | 166650296|
|Y       |  15902555|

Retrieving global information with RSamtools.
===========================
Here we have sort order description, version and programs version used to align




```r
scanBamHeader(bamToRead)$text["@HD"]
```

```
$`@HD`
[1] "VN:1.0"        "SO:coordinate"
```

```r
scanBamHeader(bamToRead)$text["@PG"]
```

```
$`@PG`
[1] "ID:subread"         "PN:subread"         "VN:Rsubread 1.20.2"
```


Getting the alignments
===========================

We can read in BAM file to R thorugh many methods. 
One of the most common is **readGAlignment()** function


```r
temp <- readGAlignments(bamToRead)
length(temp)
```

```
[1] 18144863
```

```r
temp[1]
```

```
GAlignments object with 1 alignment and 0 metadata columns:
      seqnames strand       cigar    qwidth     start       end     width
         <Rle>  <Rle> <character> <integer> <integer> <integer> <integer>
  [1]       10      -         36M        36   3000001   3000036        36
          njunc
      <integer>
  [1]         0
  -------
  seqinfo: 22 sequences from an unspecified genome
```

Converting to GRanges
===========================
GenomicAlignments can be easily converted to GRanges to remove information we are not interested in for ChIP-seq.


```r
GRanges(temp[1])
```

```
GRanges object with 1 range and 0 metadata columns:
      seqnames             ranges strand
         <Rle>          <IRanges>  <Rle>
  [1]       10 [3000001, 3000036]      -
  -------
  seqinfo: 22 sequences from an unspecified genome
```

Finer control with ScanBamParam
===========
The ScanBamParam() and it's related functions scanBamFlag() and scanBamWhich() allow much finer control over reading from Bam files. 

Here we select reads on negative strand


```r
temp <- readGAlignments(bamToRead,
        param=ScanBamParam(flag=scanBamFlag(isMinusStrand = T)))
```

====


```r
temp[1]
```

```
GAlignments object with 1 alignment and 0 metadata columns:
      seqnames strand       cigar    qwidth     start       end     width
         <Rle>  <Rle> <character> <integer> <integer> <integer> <integer>
  [1]       10      -         36M        36   3000001   3000036        36
          njunc
      <integer>
  [1]         0
  -------
  seqinfo: 22 sequences from an unspecified genome
```
