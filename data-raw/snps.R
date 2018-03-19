# Download snp150Common.txt.gz from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
# Common SNPs(150) - SNPs with >= 1% minor allele frequency (MAF), mapping only once to reference assembly.
# Double click file to unzip

#------
# bash
#------------------------------------------------------------------------------
awk -F "\t" '$12 == "single" && \
$7 == "+" && \
$9 ~ /^A$|^T$|^C$|^G$/ && \
$10 ~ /^A\/[TCG]$|^T\/[ACG]$|^C\/[TAG]$|^G\/[TCA]$/ && \
$14 !~ /^-/  && \
$2 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ \
{ print $2"\t"$3"\t"$4"\t"$5"\t"$9"\t"$10"\t"$14 }' snp150Common.txt | \
sort -k 7 -r | \
head -n 1100000 > \
snp150Common.prefiltered.txt
#------------------------------------------------------------------------------


#---
# R
#-----------------------------------------------------------------------------------------------------
filtered <- read.table("/dcl01/scharpf1/data/dbruhm/svpipeline/data/snp150Common.prefiltered.txt",
                       header = FALSE,
                       stringsAsFactors = FALSE)
colnames(filtered) <- c("chrom", "chromStart", "chromEnd", "name", "refUCSC", "observed", "avHet")

# Extracting altAllele
filtered$observed <- gsub("/", "", filtered$observed)
allele1 <- sapply(strsplit(filtered$observed, ""), `[`, 1)
allele2 <- sapply(strsplit(filtered$observed, ""), `[`, 2)

a1match <- filtered$refUCSC == allele1
a2match <- filtered$refUCSC == allele2

filtered$altAllele[a1match] <- allele2[a1match]
filtered$altAllele[a2match] <- allele1[a2match]

# Removing rows where refUCSC is not present in observed
filtered <- filtered[-which(is.na(filtered$altAllele)),]

# Convert to GRanges

library(GenomicRanges)

gr <- GRanges(seqnames = filtered$chrom,
              ranges = IRanges(start = filtered$chromEnd,
                               end = filtered$chromEnd),
              refUCSC = filtered$refUCSC,
              altAllele = filtered$altAllele)

# Making seqinfo match the seqinfo from svfilters.hg19
data(bins1kb, package = "svfilters.hg19")

gr <- keepSeqlevels(gr, seqlevels(bins1kb), pruning.mode = "coarse")
gr <- sortSeqlevels(gr)
genome(gr) <- "hg19"
seqlengths(gr) <- seqlengths(svfilters.hg19::bins1kb)
strand(gr) <- "+" ## in dnSNP the refUCSC column is always the + strand

# liftOver to hg38

library(rtracklayer)

chain <- import.chain("/dcl01/scharpf1/data/dbruhm/svpipeline/data/hg19ToHg38.over.chain")
gr38 <- unlist(liftOver(gr, chain))
gr38 <- sortSeqlevels(gr38)

# Remove any overlapping positions
gr38 <- gr38[which(countOverlaps(gr38, gr38) == 1)]

# Filling in seqinfo
genome(gr38) <- "hg38"

library(BSgenome.Hsapiens.UCSC.hg38)
seqlengths(gr38) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[1:24]

# Random sample 1,000,000 SNPs from the hg38 GRanges
gr38 <- gr38[sample(1:length(gr38), size = 1000000, replace = FALSE)]

# Sort
gr38 <- sort(gr38)

# Save object under the name snps

snps <- gr38
save(snps, file = "/dcl01/scharpf1/data/dbruhm/svpipeline/data/hg38/snps.rda")
#-----------------------------------------------------------------------------------------------------
