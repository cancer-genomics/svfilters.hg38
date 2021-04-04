

# the library TxDb.Hsapiens.UCSC.hg38.refGene does not exist and has been replaced with
# TxDb.Hsapiens.UCSC.hg38.knownGene
#
# to keep things consistent, I (Noushin) Downloaded the hg38 refGene label from the
# UCSC table browser and placed it at /dcl01/scharpf/data/pipeline-hub/public-resources/ucsc/hg38/hgTables_refGene_hg38.txt
#
#
#
#

.libPaths('/dcl01/scharpf/data/pipeline-hub/pipeline-lib/R-libs/3.10-bioc-release')

library(GenomicRanges)
library(svfilters.hg19)

data(transcripts)

cancer_genes <- unique(subset(transcripts, cancer_connection == TRUE)$gene_name)
biol_genes <-  unique(subset(transcripts, biol_sign == TRUE)$gene_name)

# a couple of gene names have been updated
# fix these manually
cancer_genes <- c(cancer_genes, 'H3-3A', 'H3-3B')

# even more for biol_sign
biol_genes <- c(biol_genes, 'H3-3A', 'H3-3B', 'COQ8A', 'NSD2', 'GUCY1A1', 'ERBIN', 'H4C9',
                'CRYBG1', 'CEP43', 'AFDN', 'CILK1', 'NSD3', 'JCAD', 'GRK2', 'CARS1', 'MRE11',
                'LHFPL6', 'KNL1', 'HASPIN', 'SEPTIN9', 'H3-3B', 'COQ8B', 'PAK5', 'SEPTIN5',
                'GRK3', 'MRTFA', 'NEXMIF', 'SEPTIN6', 'PRAG1', 'RSKR')



tx <- data.table::fread('/dcl01/scharpf/data/pipeline-hub/public-resources/ucsc/hg38/hgTables_refGene_hg38.txt')

tx <- tx[,c('chrom', 'txStart', 'txEnd', 'strand', 'name', 'name2'), with = FALSE]
tx$cancer_connection <- tx$name2 %in% cancer_genes
tx$biol_sign <- tx$name2 %in% biol_genes

data.table::setnames(tx, c('name2', 'name'), c('gene_name', 'tx_name'))


tx = GRanges(subset(tx, chrom %in% levels(seqnames(transcripts))))
tx <- sort(tx)

transcripts <- tx
save(transcripts, file="~/Software/svpackages/svfilters.hg38/data/transcripts.rda")
