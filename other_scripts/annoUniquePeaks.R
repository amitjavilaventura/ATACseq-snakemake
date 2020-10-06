# Script that calls the unique peaks in a certain sample after overlaps between two different samples. 

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(dplyr)


# ----- Define the command line arguments of this script ----- #
args = commandArgs(trailingOnly=TRUE)

inputfile <- args[1]
out1 <- args[2]
out2 <- args[3]
out3 <- args[4]
out4 <- args[5]
out5 <- args[6]
out6 <- args[7]
out7 <- args[8]

peakfile <- readPeakFile(inputfile)
peakanno <- annotatePeak(peak = peakfile, tssRegion = c(-2500, 2500), TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb = "org.Mm.eg.db")
peakanno@anno$annotation <- sub(" \\(.*\\)", "", peakanno@anno$annotation)
peakanno <- as.GRanges(peakanno)


peakanno <- as.data.frame(peakanno)
# Filter annotared bed file to obtain promoter target genes, bed files from peaks overlaping with promoters/distal regions and the coordinates of the promoter target genes.
#  promo and distal peaks: coordinates of peaks that are linked to a gene (promo = promoter peaks; distal = all non promoter peaks) + the name of the (nearest) gene. 
distal.peaks <- peakanno %>% subset(annotation != "Promoter") %>% dplyr::select(c("seqnames", "start", "end", "V4", "V5", "ENSEMBL", "SYMBOL")) 
promo.peaks <- peakanno %>% subset(annotation == "Promoter") %>% dplyr::select(c("seqnames", "start", "end", "V4", "V5", "ENSEMBL", "SYMBOL")) 
#Genes (and their cordinates) with a peak in the promoter. 
promo.targets.bed <-  peakanno %>% subset(annotation == "Promoter") %>% 
  dplyr::select(c("seqnames", "geneStart", "geneEnd", "ENSEMBL", "SYMBOL", "geneStrand")) %>% 
  mutate(geneStrand = replace(geneStrand, geneStrand == 1, "+")) %>%
  mutate(geneStrand = replace(geneStrand, geneStrand == 2, "-")) %>% unique()
promo.targets <- peakanno %>% subset(annotation == "Promoter") %>% dplyr::select("SYMBOL") %>% unique()
#Genes (and their coordinates) with a peak in a distal (all that aren't promoters) reagions. 
distal.targets.bed <- peakanno %>% subset(annotation == "Promoter") %>%
  dplyr::select(c("seqnames", "geneStart", "geneEnd", "ENSEMBL", "SYMBOL", "geneStrand")) %>% 
  mutate(geneStrand = replace(geneStrand, geneStrand == 1, "+")) %>%
  mutate(geneStrand = replace(geneStrand, geneStrand == 2, "-")) %>% unique()
distal.targets <- peakanno %>% subset(annotation != "Promoter") %>% dplyr::select("SYMBOL") %>% unique()

write.table(peakanno, file = out1, sep = "\t", quote = F, row.names = F)
write.table(promo.targets.bed, file = out2, sep = "\t", quote = F, row.names = F, col.names = F)
write.table(promo.targets, file = out3, sep = "\t", quote = F, row.names = F, col.names = F)
write.table(promo.peaks, file = out4, sep = "\t", quote = F, row.names = F, col.names = F)
write.table(distal.peaks, file = out5, sep = "\t", quote = F, row.names = F, col.names = F)
write.table(distal.targets.bed, file = out6, sep = "\t", quote = F, row.names = F, col.names = F)
write.table(distal.targets, file = out7, sep = "\t", quote = F, row.names = F, col.names = F)