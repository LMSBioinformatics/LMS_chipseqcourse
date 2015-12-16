## ------------------------------------------------------------------------
library(ChIPQC)
library(DESeq2)
library(GenomicRanges)
library(Rsubread)

## ----eval=T,echo=T-------------------------------------------------------
macsPeaks <- dir("/home/ubuntu/chipseqcourseNew/MacsPeaks/",full.names=T)
singlePeakSet <- ChIPQC:::GetGRanges(macsPeaks[1],sep="\t",simplify=T)
singlePeakSet

## ---- echo=TRUE,cache=F,eval=T-------------------------------------------

listOfPeaks <- GRangesList(lapply(macsPeaks,function(x)ChIPQC:::GetGRanges(x,sep="\t",simplify=T)))
flattenedPeaks <- reduce(unlist(listOfPeaks))

## ---- echo=TRUE----------------------------------------------------------
matOfOverlaps <- sapply(listOfPeaks,function(x)
(flattenedPeaks %over% x)+0
)

colnames(matOfOverlaps) <- basename(gsub("_peaks\\.xls","",macsPeaks))


elementMetadata(flattenedPeaks) <- matOfOverlaps

flattenedPeaks[1:2,]

## ---- echo=TRUE----------------------------------------------------------
limma:::vennCounts(as.data.frame(elementMetadata(flattenedPeaks)))

## ---- echo=T,fig.width=20, fig.height=20---------------------------------
limma:::vennDiagram(as.data.frame(elementMetadata(flattenedPeaks)))

## ---- echo=TRUE----------------------------------------------------------

mych12Peaks <- flattenedPeaks[elementMetadata(flattenedPeaks)$mycch12rep1 + 
                 elementMetadata(flattenedPeaks)$mycch12rep2 == 2]
mycMelPeaks <- flattenedPeaks[elementMetadata(flattenedPeaks)$mycmelrep1 +                  elementMetadata(flattenedPeaks)$mycmelrep2 == 2]


## ---- echo=TRUE----------------------------------------------------------

mycMelPeaks_Only <- flattenedPeaks[elementMetadata(flattenedPeaks)$mycmelrep1 +                  elementMetadata(flattenedPeaks)$mycmelrep2 == 2 &
elementMetadata(flattenedPeaks)$mycch12rep1 + 
                 elementMetadata(flattenedPeaks)$mycch12rep2 == 0]

mycMelPeaks_Only[1,]

## ---- echo=TRUE----------------------------------------------------------

highConfidence_Only <- flattenedPeaks[elementMetadata(flattenedPeaks)$mycmelrep1 +                  elementMetadata(flattenedPeaks)$mycmelrep2 == 2 |
elementMetadata(flattenedPeaks)$mycch12rep1 + 
                 elementMetadata(flattenedPeaks)$mycch12rep2 == 2]


## ---- echo=TRUE----------------------------------------------------------
boxplot(width(highConfidence_Only))
abline(h=400,col="red")

## ---- echo=TRUE----------------------------------------------------------
PeaksToCount <- resize(highConfidence_Only,width = 400,fix = "center")
PeaksToCount[1:2,]

## ----eval=T,echo=TRUE,message=T,warning=F--------------------------------
Bams <- dir("/home/ubuntu/chipseqcourseData/sortedbams/",pattern="*sorted\\..*bam$",full.names=T)

simplePeaksToCount <- ChIPQC:::GetGRanges(PeaksToCount,simple=T,simplify=T)


toCount <- data.frame(GeneID = paste0("ID",seq(1,length(simplePeaksToCount))),
                      Chr=seqnames(simplePeaksToCount),Start=start(simplePeaksToCount),
                      End=end(simplePeaksToCount),Strand="-")

out  <- capture.output(myCountTable <- featureCounts(Bams[1],annot.ext = toCount,nthreads=1))

myCountTable$counts[1:10]

## ---- echo=FALSE---------------------------------------------------------
if(!file.exists("/home/ubuntu/chipseqcourseNew/robjects/MycCounts.RData")){
out <- capture.output(
  myCountTableList <- lapply(Bams,function(x)
                               featureCounts(x,annot.ext = toCount,nthreads=4))
)
}

## ---- echo=TRUE,eval=F---------------------------------------------------
## 
## out <- capture.output(
##   myCountTableList <- lapply(Bams,function(x)
##                                featureCounts(x,annot.ext = toCount,nthreads=4))
## )
## 

## ---- echo=TRUE,eval=FALSE-----------------------------------------------
## 
## countTable <- sapply(myCountTableList,function(x)x$counts)
## rownames(countTable) <- paste0(toCount[,1],"-",toCount[,2],";",toCount[,3],"-",toCount[,4])
## colnames(countTable) <- c("ch12myc","ch12myc","ch12input","melmyc","melmyc","meinput")

## ---- echo=FALSE,eval=FALSE----------------------------------------------
## if(!file.exists("/home/ubuntu/chipseqcourseNew/robjects/MycCounts.RData")){
## save(countTable,file="/home/ubuntu/chipseqcourseNew/robjects/MycCounts.RData")
## }

## ---- echo=FALSE,eval=T--------------------------------------------------

load("/home/ubuntu/chipseqcourseNew/robjects/MycCounts.RData")

## ---- echo=TRUE----------------------------------------------------------
library("DESeq2")

colData <- data.frame(SampleName=colnames(countTable[,-c(3,6)]),CellLine=c("ch12","ch12","mel","mel"))
dds <- DESeqDataSetFromMatrix(countData = countTable[,-c(3,6)],
                              colData = colData,
                              design = ~ CellLine)

dds <- DESeq(dds)
testcellline <- results(dds, contrast=c("CellLine","ch12","mel"))

## ---- echo=T,eval=T------------------------------------------------------

testcellline <- testcellline[order(testcellline$pvalue),]
testcellline <- as.data.frame(testcellline)
testcellline[1:3,]

## ---- echo=T,eval=T------------------------------------------------------
resultsMat <- matrix(unlist(strsplit(rownames(testcellline),"-")),byrow=T,ncol=3)
resultMatPart2 <- matrix(unlist(strsplit(resultsMat[,2],";")),byrow=T,ncol=2)
DEasGRanges <- GRanges(seqnames=resultMatPart2[,1],IRanges(as.numeric(as.vector(resultMatPart2[,2])),as.numeric(as.vector(resultsMat[,3]))),elementMetadata=testcellline)
colnames(elementMetadata(DEasGRanges)) <- gsub("elementMetadata","",colnames(elementMetadata(DEasGRanges)))
DEasGRanges[1:5,]

## ---- echo=F,eval=F------------------------------------------------------
## 
## save(DEasGRanges,file="/home/ubuntu/chipseqcourseNew/robjects/DEasGRanges.RData")

## ---- echo=T,eval=T------------------------------------------------------
sessionInfo()

