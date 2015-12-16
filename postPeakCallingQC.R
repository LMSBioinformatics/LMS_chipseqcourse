## ---- echo=T-------------------------------------------------------------
library(ChIPQC)
library(GenomeInfoDb)

## ---- echo=T,eval=T------------------------------------------------------
myBL <- ChIPQC:::GetGRanges("/home/ubuntu/chipseqcourseData/referencedata/mm9-blacklist.bed")
myBL <- renameSeqlevels(myBL,gsub("chr","",seqlevels(myBL)))

## ---- echo=T-------------------------------------------------------------
mm9Anno <- ChIPQC:::getAnnotation("mm9",AllChr=NULL)
mm9AnnoNew <- lapply(mm9Anno[-1],
                  function(x)
                  renameSeqlevels(x,gsub("chr","",seqlevels(x))
                       )
              )

## ----echo=T,eval=T,cache=T-----------------------------------------------
indexedBams <- dir("/home/ubuntu//chipseqcourseData/sortedbams/",
                   pattern="*sorted\\..*bam$",full.names=T)

SampleID <- c("myc_ch12_1","myc_ch12_2","myc_Mel_1","myc_Mel_2")
Tissue <- c(rep("ch12",2),rep("mel",2))
Treatment <- rep(NA,4)
Condition <- rep(NA,4)
Factor <- c(rep("myc",2),rep("myc",2))
Replicate <- c(1,2,1,2)

## ----echo=T,eval=T,cache=T-----------------------------------------------

bamReads  <- indexedBams[-c(3,6)]
bamControl= c("/home/ubuntu//chipseqcourseData/sortedbams//wgEncodeSydhTfbsCh12InputIggmusRawDatasorted.bam.bam",
              "/home/ubuntu//chipseqcourseData/sortedbams//wgEncodeSydhTfbsCh12InputIggmusRawDatasorted.bam.bam",
              "/home/ubuntu//chipseqcourseData/sortedbams//wgEncodeSydhTfbsMelInputIggmusRawDatasorted.bam.bam",
              "/home/ubuntu//chipseqcourseData/sortedbams//wgEncodeSydhTfbsMelInputIggmusRawDatasorted.bam.bam")
ControlID = c("ch12","ch12","MEL","MEL")

## ----echo=T,eval=T,cache=T-----------------------------------------------


Peaks <- dir("/home/ubuntu/chipseqcourseNew/MacsPeaks/",full.names=T)
PeakCaller <- rep("macs",4)


ss <- data.frame(SampleID,Tissue,Factor,Treatment,Replicate,
                 Condition,bamReads,
                 bamControl,ControlID,
                 Peaks,
                 PeakCaller)



## ----echo=F,eval=T,cache=F-----------------------------------------------
load("/home/ubuntu//chipseqcourseNew/robjects/ChIPQCwithPeaks.RData")

## ----echo=T,eval=T,cache=T-----------------------------------------------
ss[1,]

## ----echo=F,eval=T,cache=F-----------------------------------------------
load("../robjects/ChIPQCwithPeaks.RData")

## ----echo=T,eval=F,cache=T-----------------------------------------------
## register(SerialParam(), default=TRUE)
## #p <- MulticoreParam(workers = 4)
## #register(p)
## 
## res <- ChIPQC(ss,annotation=mm9AnnoNew,
##               chromosomes=paste0(1:10),
##               blacklist=myBL,consensus=F,bCount=F)

## ----echo=T,eval=T,cache=T-----------------------------------------------
QCmetrics(res)

## ----echo=T,eval=T,cache=T,fig.width=20, fig.height=10-------------------
plotCC(res)

## ----echo=T,eval=T,cache=T,fig.width=20, fig.height=10-------------------
plotSSD(res)

## ----echo=T,eval=T,cache=T,fig.width=20, fig.height=10-------------------
plotRegi(res)

## ----echo=T,eval=T,cache=T,fig.width=20, fig.height=10-------------------
plotFribl(res)

## ----echo=T,eval=T,cache=T,fig.width=20, fig.height=10-------------------
plotFrip(res)

## ----echo=T,eval=T,cache=T,fig.width=20, fig.height=10-------------------
plotPeakProfile(res)

## ----echo=T,eval=T,cache=T-----------------------------------------------
sessionInfo()


