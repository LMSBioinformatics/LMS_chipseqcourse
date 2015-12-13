GettingData
========================================================
author: MRC Clinical Sciences Centre
date:http://mrccsc.github.io/r_course/introToR_Session1.html
width: 1440
height: 1100
autosize: true
font-import: <link href='http://fonts.googleapis.com/css?family=Slabo+27px' rel='stylesheet' type='text/css'>
font-family: 'Slabo 27px', serif;
css:style.css

As we saw there are many places to grab data from.
========================================================

Some of the most popular are :-
- [GEO](http://www.ncbi.nlm.nih.gov/geo/)
- [SRA](http://www.ncbi.nlm.nih.gov/sra)
- [ENA.](http://www.ebi.ac.uk/ena)


Interfacing with Repositories in R.
========================================================
Lets load the packages we need, which repositories are associated with which tool should be obvious.


```r
library(GEOquery)
library(SRAdb)
```

GEOquery
========================================================
GEOquery allows you to find out information about known GEO id.
Here is an example of an expression array sample


```r
gdsArray <- getGEO("GSM11805")
```


```r
show(gdsArray)
```

Have a look at gds object by typing **show(gds)** and we can see it contains alot of information



GEOquery
========================================================
So lets go look for some samples we have already downloaded.
We are using Encode data so we find some GEOinformation  [here](http://hgdownload.cse.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeSydhTfbs/files.txt)

One we are interest in is GSM912906


```r
gdsChIP <- getGEO("GSM912906")
```

Have again look at gds object


GEOquery
========================================================
GEO has much of the meta data available try **meta(gds)**


```r
names(Meta(gdsChIP))
```

```
 [1] "biomaterial_provider_ch1" "channel_count"           
 [3] "characteristics_ch1"      "contact_address"         
 [5] "contact_city"             "contact_country"         
 [7] "contact_email"            "contact_institute"       
 [9] "contact_name"             "contact_state"           
[11] "contact_zip/postal_code"  "data_processing"         
[13] "data_row_count"           "extract_protocol_ch1"    
[15] "geo_accession"            "instrument_model"        
[17] "last_update_date"         "library_selection"       
[19] "library_source"           "library_strategy"        
[21] "molecule_ch1"             "organism_ch1"            
[23] "platform_id"              "relation"                
[25] "series_id"                "source_name_ch1"         
[27] "status"                   "submission_date"         
[29] "supplementary_file_1"     "supplementary_file_2"    
[31] "supplementary_file_3"     "taxid_ch1"               
[33] "title"                    "treatment_protocol_ch1"  
[35] "type"                    
```

GEOquery
========================================================
You can also just see a subsection of information by using **head(Meta(object))**
Here we show subsection bits of information

```r
head(Meta(gdsChIP),3)
```

```
$biomaterial_provider_ch1
[1] "Weissman lab"

$channel_count
[1] "1"

$characteristics_ch1
 [1] "lab: Stanford-m"                                                                                                                               
 [2] "lab description: Snyder - Stanford University"                                                                                                 
 [3] "datatype: ChipSeq"                                                                                                                             
 [4] "datatype description: Chromatin IP Sequencing"                                                                                                 
 [5] "cell: CH12"                                                                                                                                    
 [6] "cell organism: mouse"                                                                                                                          
 [7] "cell description: B-cell lymphoma (GM12878 analog)"                                                                                            
 [8] "cell sex: F"                                                                                                                                   
 [9] "antibody: c-Myc"                                                                                                                               
[10] "antibody antibodydescription: rabbit polyclonal to amino acids 1-262 of c-Myc human origin. Antibody Target: c-Myc"                            
[11] "antibody targetdescription: transcription factor; c-Myc-encoded proteins function in cell proliferation,differentiation and neoplastic disease"
[12] "antibody vendorname: Santa Cruz Biotech"                                                                                                       
[13] "antibody vendorid: sc-764"                                                                                                                     
[14] "treatment: None"                                                                                                                               
[15] "treatment description: No special treatment or protocol applies"                                                                               
[16] "control: IgG-rab"                                                                                                                              
[17] "control description: Input signal from Normal Rabbit IgG ChIP-seq."                                                                            
[18] "age: immortalized"                                                                                                                             
[19] "age description: Immortal cells"                                                                                                               
[20] "control: IgG-rab"                                                                                                                              
[21] "control description: Input signal from Normal Rabbit IgG ChIP-seq."                                                                            
[22] "controlid: CH12/None/Input/IgG-rab"                                                                                                            
[23] "replicate: 1"                                                                                                                                  
[24] "strain: B10.H-2aH-4bp/Wts"                                                                                                                     
[25] "strain description: Derived by inbreeding from selected F2 progeny of B10.A X B10.129"                                                         
```

GEOquery
========================================================
We can see now that this file has 3 supplementary files attached.


```r
Meta(gdsChIP)$supplementary_file_1
[1] "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM912nnn/GSM912906/suppl/GSM912906_mm9_wgEncodeSydhTfbsCh12CmycIggrabPk.narrowPeak.gz"
Meta(gdsChIP)$supplementary_file_2
[1] "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM912nnn/GSM912906/suppl/GSM912906_mm9_wgEncodeSydhTfbsCh12CmycIggrabSig.bigWig"
Meta(gdsChIP)$supplementary_file_3
[1] "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX140/SRX140369"
```
With the third file pointing to SRA

SRAdb
========================================================
SRAdb acts an interface between SRA and ENA


```r
fileToGen <- basename(Meta(gdsChIP)$supplementary_file_3) 
fileToGen
```

```
[1] "SRX140369"
```
 SRAdb
========================================================
 The first step in querying SRAdb is to download their schema.
 
 ```r
 if(!file.exists("/Users/tcarroll/chipseqcourseNew/SRAmetadb.sqlite")){
  sqlfile <<- getSRAdbFile()
 }else{
  sqlfile <- "/Users/tcarroll/chipseqcourseNew/SRAmetadb.sqlite"
 }
 sra_con <- dbConnect(SQLite(),sqlfile)
 ```
You can try this but it takes a few minutes
 SRAdb
========================================================
Now we can retrieve the location of file we want in fastq format
 
 ```r
 rs = listSRAfile( c(fileToGen), sra_con, fileType = "fastq") 
 rs$ftp
 ```
 
 ```
 [1] "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR478/SRR478478/SRR478478.fastq.gz"
 [2] "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR478/SRR478479/SRR478479.fastq.gz"
 ```
  SRAdb - more info
========================================================
Now we can retrieve a little more information about the sample with 
 
 ```r
 sra_con <- dbConnect(SQLite(),sqlfile) 
 info <- getSRAinfo("SRR478478",sra_con=sra_con)
 info
 ```
 
AnnotationHub
========================================================
AnnotationHub provides a nice interface to retrieve data into R
Try the code below for yourself.

```r
ah = AnnotationHub()
display(ah)
```
 
 
