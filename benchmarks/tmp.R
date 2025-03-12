main_path <- "/home/grocamora/RytenLab-Research/37-UTRome_pipeline" # here::here()

### Input Paths
input_gtf <- file.path(main_path, "Results/Sqanti3_Rescue/sq3.annotated_rescued.gtf")
txdb <- GenomicFeatures::makeTxDbFromGFF(input_gtf)

set.seed(0)
grlExons <- GenomicFeatures::exonsBy(txdb, use.names = T)
mcols(grlExons) <- NULL

BPPARAM = BiocParallel::MulticoreParam(workers = 16, progressbar = T)

grlExons <- sample(grlExons, 1000)
dfTxGene <- AnnotationDbi::select(txdb, keys = names(grlExons), keytype = "TXNAME", columns = "GENEID")
mapTxToGene <- setNames(dfTxGene$GENEID, dfTxGene$TXNAME)



a %>% sort()
b %>% sort()




a[1120:1130]
b[1120:1130]

gr <- grlExons[["ENST00000522546.1"]]

gr <- GRanges(
  seqnames = c("chr8", "chr8"),
  ranges = IRanges(start = c(138779509, 138777907), end = c(138779562, 138778405)),
  strand = c("-", "-")
)



library(txcutr)
library(gtxcutr)
library(tidyverse)

main_path <- "/home/grocamora/RytenLab-Research/37-UTRome_pipeline" # here::here()

### Input Paths
input_gtf <- file.path(main_path, "Results/Sqanti3_Rescue/sq3.annotated_rescued.gtf")
gtf <- rtracklayer::import.gff(input_gtf)

set.seed(0)
filter_transcripts <- GenomicRanges::mcols(gtf)[["transcript_id"]] %>% unique %>% sample(1000)
gtf_filter <- gtf %>% plyranges::filter(transcript_id %in% filter_transcripts)
txdb <- GenomicFeatures::makeTxDbFromGRanges(gtf_filter)

system.time({txdb_def <- txcutr::truncateTxome(txdb, maxTxLength = 500)})
system.time({txdb_mod <- gtxcutr::truncateTxome(txdb, maxTxLength = 500)})

benchmark_results <- microbenchmark::microbenchmark(
  default_method = txcutr::truncateTxome(txdb, maxTxLength = 500, BPPARAM = BiocParallel::MulticoreParam(workers = 20)),
  modded_method = gtxcutr::truncateTxome(txdb, maxTxLength = 500, BPPARAM = BiocParallel::MulticoreParam(workers = 20)),
  times = 5
)


library(GenomicFeatures)
all.equal(genes(txdb_def) %>% tibble::as_tibble(), genes(txdb_mod) %>% tibble::as_tibble())
all.equal(transcripts(txdb_def) %>% tibble::as_tibble(), transcripts(txdb_mod) %>% tibble::as_tibble())
all.equal(exons(txdb_def) %>% tibble::as_tibble(), exons(txdb_mod) %>% tibble::as_tibble())