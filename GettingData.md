GettingData
========================================================

As we saw there are many places to grab data from.
========================================================

Some of the most popular are GEO, SRA and ENA.

Interfacing with Repositories in R.
========================================================
Lets load the packages we need, the repo should be obvious.


```r
library(GEOquery)
library(SRAdb)
```

GEOquery
========================================================
So lets go look for some samples we have already downloaded
We are using Encode data so we find some GEOinforma tion [here](http://hgdownload.cse.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeSydhTfbs/files.txt)

One we are interest in is GSM912906


```r
gds <- getGEO("GSM11805")
```

Have a look at gds object

GEOquery
========================================================
So lets go look for some samples we have already downloaded
We are using Encode data so we find some GEOinforma tion [here](http://hgdownload.cse.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeSydhTfbs/files.txt)

One we are interest in is GSM912906


```r
gds <- getGEO("GSM912906")
```

Have a look at gds object


GEOquery
========================================================
GEO has much of the meta data available try meta(gds)


```r
Meta(gds)$supplementary_file_1
```

```
[1] "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM912nnn/GSM912906/suppl/GSM912906_mm9_wgEncodeSydhTfbsCh12CmycIggrabPk.narrowPeak.gz"
```

```r
Meta(gds)$supplementary_file_2
```

```
[1] "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM912nnn/GSM912906/suppl/GSM912906_mm9_wgEncodeSydhTfbsCh12CmycIggrabSig.bigWig"
```

```r
Meta(gds)$supplementary_file_3
```

```
[1] "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX140/SRX140369"
```
With the third file pointing to SRA

SRAdb
========================================================
SRAdb acts an interface between SRA and ENA


```r
fileToGen <- basename(Meta(gds)$supplementary_file_3) 
fileToGen
```

```
[1] "SRX140369"
```
 SRAdb
========================================================
 The first step in querying SRAdb is to download their schema.


