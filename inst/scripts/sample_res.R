library(annoLinker)
library(rtracklayer)
library(TxDb.Drerio.UCSC.danRer10.refGene)
library(org.Dr.eg.db)
txdb <- TxDb.Drerio.UCSC.danRer10.refGene
org <- org.Dr.eg.db
extPath <- system.file('extdata', package='annoLinker')
## load peaks
peaks <- rtracklayer::import(file.path(extPath, 'peaks.bed'))
## load interactions
interactions <- rtracklayer::import(file.path(extPath, 'interaction.bedpe'))
## load annotation data
annoData <- genes(txdb)
anno <- annoLinker(peaks, annoData, interactions, verbose=TRUE)
saveRDS(anno, 'sample_res.rds')
