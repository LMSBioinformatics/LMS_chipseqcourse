## ----echo=T,eval=T-------------------------------------------------------
library(GenomicRanges)
library(DESeq2)
library(goseq)

load("/home/ubuntu/chipseqcourseNew/robjects/DEasGRanges.RData")
DEasGRanges[1,]

## ---- echo=TRUE,collapse=T-----------------------------------------------
library(GenomicRanges)
library(GenomeInfoDb)
library(DESeq2)

UpInMel <- DEasGRanges[DEasGRanges$.padj < 0.05 
                       & !is.na(DEasGRanges$.padj) 
                       & DEasGRanges$.log2FoldChange > 0]

DownInMel <- DEasGRanges[DEasGRanges$.padj < 0.05 
                         & !is.na(DEasGRanges$.padj) 
                         & DEasGRanges$.log2FoldChange < 0]

length(UpInMel)

length(DownInMel)

## ---- echo=TRUE,collapse=F-----------------------------------------------

mm9Genes <- read.delim("/home/ubuntu/chipseqcourseNew/robjects/mm9Genes_May2012.txt",sep="\t",h=T)
mm9GeneRanges <- GRanges(seqnames=paste0("chr",mm9Genes[,3]),ranges=IRanges(start=mm9Genes[,1],end=mm9Genes[,2]),strand=mm9Genes[,4],name=mm9Genes[,5],biotype=mm9Genes[,6])

JustChrOfInterest <- unique(as.vector(seqnames(mm9GeneRanges)))[grep("\\d.|\\d|X|Y",unique(as.vector(seqnames(mm9GeneRanges))))]

mm9PC <- mm9GeneRanges[mm9GeneRanges$biotype == "protein_coding"]
mm9PC <- mm9PC[order(width(mm9PC),decreasing=T)]
mm9PC <- mm9PC[match(unique(mm9PC$name),mm9PC$name)]
mm9PC <- mm9PC[!mm9PC$name == ""]
mm9PC <- mm9PC[seqnames(mm9PC) %in% JustChrOfInterest]
mm9PC[1:3]


## ---- echo=TRUE,collapse=T-----------------------------------------------

mm9Promoters <- promoters(mm9PC,1000,1000)
mm9PC[1:2,]
mm9Promoters[1:2,]
mm9Promoters <- renameSeqlevels(mm9Promoters,gsub("chr","",seqlevels(mm9Promoters)))
mm9Promoters

## ---- echo=TRUE,collapse=F-----------------------------------------------

GeneswithUpInMel <- mm9Promoters %over% UpInMel + 0
GeneswithDownInMel <- mm9Promoters %over% DownInMel + 0

names(GeneswithUpInMel) <- names(GeneswithDownInMel) <- mm9Promoters$name

## ------------------------------------------------------------------------


library(KEGG.db)
library(goseq)
xx <- as.list(KEGGPATHID2NAME)
temp <- cbind(names(xx),unlist(xx))
addKeggTogoseq <- function(JX,temp){
  for(l in 1:nrow(JX)){
    if(JX[l,1] %in% temp[,1]){
      JX[l,"term"] <- temp[temp[,1] %in% JX[l,1],2]
      JX[l,"ontology"] <- "KEGG"
    }
    
  }
  return(JX)
}

## ---- echo=TRUE,collapse=F-----------------------------------------------

pwf=nullp(GeneswithUpInMel,"mm9","geneSymbol")
GeneswithUpInMelEnrich <- goseq(pwf,"mm9","geneSymbol",test.cats=c("GO:BP","GO:MF","KEGG"),method="Hypergeometric")
pwf=nullp(GeneswithDownInMel,"mm9","geneSymbol")
GeneswithDownInMelEnrich <- goseq(pwf,"mm9","geneSymbol",test.cats=c("GO:BP","GO:MF","KEGG"),method="Hypergeometric")

GeneswithUpInMelEnrich <- addKeggTogoseq(GeneswithUpInMelEnrich,temp)
GeneswithDownInMelEnrich <- addKeggTogoseq(GeneswithDownInMelEnrich,temp)


## ---- echo=TRUE,collapse=F-----------------------------------------------
GeneswithUpInMelEnrich[1:10,]



## ---- echo=TRUE,collapse=F-----------------------------------------------
GeneswithDownInMelEnrich[1:10,]

## ---- echo=TRUE,collapse=F-----------------------------------------------

library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm9)
genome <- BSgenome.Mmusculus.UCSC.mm9
UpInMel <- renameSeqlevels(UpInMel,paste0("chr",seqlevels(UpInMel)))
UpInMel <- UpInMel[seqnames(UpInMel) != "chrMT"]

DownInMel <- renameSeqlevels(DownInMel,paste0("chr",seqlevels(DownInMel)))
DownInMel <- DownInMel[seqnames(DownInMel) != "chrMT"]

UpinMelSequences <- getSeq(genome,GRanges(UpInMel))
DowninMelSequences <- getSeq(genome,GRanges(DownInMel))

DowninMelSequences[1:2,]


## ---- echo=TRUE,collapse=F-----------------------------------------------
writeXStringSet(DowninMelSequences,file="DowninMel.fa")
writeXStringSet(UpinMelSequences,file="UpinMel.fa")

## ---- eval=F,echo=TRUE,collapse=F----------------------------------------
## source("https://bioconductor.org/biocLite.R")
## biocLite("MotifDB")

## ---- eval=T,echo=TRUE,collapse=F----------------------------------------
library(seqLogo)
library(MotifDb)
PotentialMyc <- as.list(subset (MotifDb, tolower (geneSymbol) == "myc"))
seqLogo(PotentialMyc[[9]])

## ---- eval=T,echo=TRUE,collapse=F----------------------------------------
UpinMelCount <- lapply(UpinMelSequences,function(x)countPWM(PotentialMyc[[9]],x))
DowninMelCount <- lapply(DowninMelSequences,function(x)countPWM(PotentialMyc[[9]],x))
RevComUpinMelCount <- lapply(reverseComplement(UpinMelSequences),function(x)countPWM(PotentialMyc[[9]],x))
RevComDowninMelCount <- lapply(reverseComplement(DowninMelSequences),function(x)countPWM(PotentialMyc[[9]],x))


UpinMelFreq <- (sum(unlist(c(UpinMelCount,RevComUpinMelCount))))/(sum(width(UpinMelSequences))*2)
DowninMelFreq <- (sum(unlist(c(DowninMelCount,RevComDowninMelCount))))/(sum(width(DowninMelSequences))*2)



## ---- echo=TRUE,collapse=F-----------------------------------------------

chr1MycCounts <- countPWM(PotentialMyc[[9]],genome$chr1)
RevComchr1MycCounts <- countPWM(PotentialMyc[[9]],reverseComplement(genome$chr1))
chr1Freq <- (sum(c(chr1MycCounts,RevComchr1MycCounts)))/(length(genome$chr1)*2)

## ---- echo=TRUE,collapse=T-----------------------------------------------

UpinMelFreq
DowninMelFreq
chr1Freq

## ---- echo=TRUE,collapse=F-----------------------------------------------
sessionInfo()

