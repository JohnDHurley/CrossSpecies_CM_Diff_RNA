library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readr)
library(tidyverse)

### Loading in txt file


Peaks_file <- read.delim("~/ASHG_Data/FeatureCounts/ASHG_h_samples_counts_2pSTP_noD4.txt", comment.char="#")
# #or use rtracklayer
# gr <- rtracklayer::import("~/ASHG_Data/FeatureCounts/ASHG_h_samples_counts_2pSTP_noD4.txt", format = "BED")

head(Peaks_file)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
### check if these exist:
seqnames(gr)
ranges(gr)
mcols(gr)

Peaks_file %>%
  column_to_rownames()

gr <- GRanges(
  seqnames = Peaks_file$Chr,
  ranges   = IRanges(start = Peaks_file$Start + 1, end = Peaks_file$End),
  strand   = "*")

## keep extra columns if present
mcols(gr)$name <- Peaks_file$Geneid


# peakannotation ----------------------------------------------------------

my_anno_peaks <- annotatePeak(gr, tssRegion =c(-2000,2000), TxDb=txdb)

###Viola!   save as RDS
plotAnnoBar(my_anno_peaks)
plotAnnoPie(my_anno_peaks)

####  now pick only autosomes
autosomes <- paste0("chr", 1:22)
gr <- gr[seqnames(gr) %in% autosomes]
Collapse

TSS_Regions <- my_anno_peaks@anno %>% 
  as.data.frame() %>% 
  dplyr::filter(str_starts(annotation,"Promoter")) %>% 
  distinct(name)














