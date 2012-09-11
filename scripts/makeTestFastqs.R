library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ShortRead)
eg <- org.Hs.egSYMBOL2EG[["TP53"]]
tx <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene,
                  vals = list(gene_id = eg))
newRange <- GRanges(ranges=IRanges(start(range(tx)), end(range(tx))), seqnames="17")
param <- ScanBamParam(what = c("seq", "qual"),
                      which = newRange + 1e6)
writeFastqForReads <- function(x, file) {
  seqs <- mcols(x)$seq
  quals <- mcols(x)$qual
  negStranded <- which(strand(x) == "-")
  seqs[negStranded] <- reverseComplement(seqs[negStranded])
  quals[negStranded] <- reverse(quals[negStranded])
  
  writeFastq(ShortReadQ(seqs,
                        BStringSet(quals),
                        BStringSet(names(x))),
             file)
}

bamDirStem <- "/gnet/is3/research/data/bioinfo/ngs_analysis/RNA-seq/R116_pipeline3.4"
bamDirs <- file.path(bamDirStem, c("H1993_SAM633816", "H2073_SAM633819"), "bams")
if (!all(file.exists(bamDirs))) {
  stop("Could not find all BAM-file directories")
}

##take only the concordant unique reads
concUniqs <- dir(bamDirs, pattern="\\.concordant_uniq\\.bam$", full.names=TRUE)
names(concUniqs) <- c("H1993", "H2073")

filePath <- file.path(bamDirs[1], "H1993_SAM633816.concordant_uniq.bam")
if(!file.exists(filePath)) stop("Could not find the file at ", filePath)
gap <- readGappedAlignmentPairs(filePath, param = param, use.names = TRUE)
##subsampling because alignment took too long for testing and examples
set.seed(1)
sampledGap <- gap[sample(seq_len(length(gap)), 2500L)]
writeFastqForReads(first(sampledGap), "H1993.reads-first.fastq")
writeFastqForReads(last(sampledGap), "H1993.reads-last.fastq")
system("gzip H1993.reads-first.fastq")
system("gzip H1993.reads-last.fastq")

filePath <- file.path(bamDirs[2], "H2073_SAM633819.concordant_uniq.bam")
gap <- readGappedAlignmentPairs(filePath, param = param, use.names = TRUE)
set.seed(2)
sampledGap <- gap[sample(seq_len(length(gap)), 2500L)]
writeFastqForReads(first(sampledGap), "H2073.reads-first.fastq")
writeFastqForReads(last(sampledGap), "H2073.reads-last.fastq")
system("gzip H2073.reads-first.fastq")
system("gzip H2073.reads-last.fastq")
