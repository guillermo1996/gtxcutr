library(GenomicRanges)
library(GenomicFeatures)
BiocParallel::register(BiocParallel::SerialParam())

############
## Mock Data
############

## Single Exon Gene - 3'
gr_contig_3p <- GRanges(
  seqnames=rep("chr1", 5),
  strand="+",
  ranges=IRanges(end=5000,
                 width=c(1000, 1000, 1000, 800, 800)),
  type=c("gene", "transcript", "exon", "transcript", "exon"),
  ID=c("gene_1", "tx_1", "exon_1", "tx_2", "exon_2"),
  Parent=c(NA, "gene_1", "tx_1", "gene_1", "tx_2"),
  gene_id="gene_1",
  tx_id=c(NA, "tx_1", "tx_1", "tx_2", "tx_2"),
  exon_id=c(NA, NA, "exon_1", NA, "exon_2"))

txdb_contig_3p <- makeTxDbSafe(gr_contig_3p)
txdb_contig_3p_inv <- makeTxDbSafe(invertStrand(gr_contig_3p))

## Negative Strand - 3'
gr_contig_neg_3p <- GRanges(
  seqnames=rep("chr1", 5),
  strand="-",
  ranges=IRanges(start=5000,
                 width=c(1000, 1000, 1000, 800, 800)),
  type=c("gene", "transcript", "exon", "transcript", "exon"),
  ID=c("gene_1", "tx_1", "exon_1", "tx_2", "exon_2"),
  Parent=c(NA, "gene_1", "tx_1", "gene_1", "tx_2"),
  gene_id="gene_1",
  tx_id=c(NA, "tx_1", "tx_1", "tx_2", "tx_2"),
  exon_id=c(NA, NA, "exon_1", NA, "exon_2"))

txdb_contig_neg_3p <- makeTxDbSafe(gr_contig_neg_3p)
txdb_contig_neg_3p_inv <- makeTxDbSafe(invertStrand(gr_contig_neg_3p))

## Overlapping Genes - 3'
gr_multigene_3p <- GRanges(
  seqnames=rep("chr1", 6),
  strand="+",
  ranges=IRanges(end=5000,
                 width=c(1000, 1000, 1000, 1200, 1200, 1200)),
  type=c("gene", "transcript", "exon",
         "gene", "transcript", "exon"),
  ID=c("gene_1", "tx_1", "exon_1",
       "gene_2", "tx_2", "exon_2"),
  Parent=c(NA, "gene_1", "tx_1",
           NA, "gene_2", "tx_2"),
  gene_id=c("gene_1", "gene_1", "gene_1",
            "gene_2", "gene_2", "gene_2"),
  tx_id=c(NA, "tx_1", "tx_1",
          NA, "tx_2", "tx_2"),
  exon_id=c(NA, NA, "exon_1",
            NA, NA, "exon_2"))

txdb_multigene_3p <- makeTxDbSafe(gr_multigene_3p)

## Single Exon Gene - 5'
gr_contig_5p <- GRanges(
  seqnames=rep("chr1", 5),
  strand="+",
  ranges=IRanges(start=5000,
                 width=c(1000, 1000, 1000, 800, 800)),
  type=c("gene", "transcript", "exon", "transcript", "exon"),
  ID=c("gene_1", "tx_1", "exon_1", "tx_2", "exon_2"),
  Parent=c(NA, "gene_1", "tx_1", "gene_1", "tx_2"),
  gene_id="gene_1",
  tx_id=c(NA, "tx_1", "tx_1", "tx_2", "tx_2"),
  exon_id=c(NA, NA, "exon_1", NA, "exon_2"))

txdb_contig_5p <- makeTxDbSafe(gr_contig_5p)
txdb_contig_5p_inv <- makeTxDbSafe(invertStrand(gr_contig_5p))

## Negative Strand - 5'
gr_contig_neg_5p <- GRanges(
  seqnames=rep("chr1", 5),
  strand="-",
  ranges=IRanges(end=5000,
                 width=c(1000, 1000, 1000, 800, 800)),
  type=c("gene", "transcript", "exon", "transcript", "exon"),
  ID=c("gene_1", "tx_1", "exon_1", "tx_2", "exon_2"),
  Parent=c(NA, "gene_1", "tx_1", "gene_1", "tx_2"),
  gene_id="gene_1",
  tx_id=c(NA, "tx_1", "tx_1", "tx_2", "tx_2"),
  exon_id=c(NA, NA, "exon_1", NA, "exon_2"))

txdb_contig_neg_5p <- makeTxDbSafe(gr_contig_neg_5p)
txdb_contig_neg_5p_inv <- makeTxDbSafe(invertStrand(gr_contig_neg_5p))


## Overlapping Genes - 5'
gr_multigene_5p <- GRanges(
  seqnames=rep("chr1", 6),
  strand="+",
  ranges=IRanges(start=5000,
                 width=c(1000, 1000, 1000, 1200, 1200, 1200)),
  type=c("gene", "transcript", "exon",
         "gene", "transcript", "exon"),
  ID=c("gene_1", "tx_1", "exon_1",
       "gene_2", "tx_2", "exon_2"),
  Parent=c(NA, "gene_1", "tx_1",
           NA, "gene_2", "tx_2"),
  gene_id=c("gene_1", "gene_1", "gene_1",
            "gene_2", "gene_2", "gene_2"),
  tx_id=c(NA, "tx_1", "tx_1",
          NA, "tx_2", "tx_2"),
  exon_id=c(NA, NA, "exon_1",
            NA, NA, "exon_2"))

txdb_multigene_5p <- makeTxDbSafe(gr_multigene_5p)

########
## Tests
########

test_that("identical transcripts are merged, positive strand, 3'", {
  LENGTHS_TO_TEST <- c(100, 500)

  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate3primeTxome(txdb_contig_3p, maxTxLength=n, quiet = T)

    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 1)
    expect_equal(n_exons, 1)
  }
})

test_that("identical transcripts are merged, negative strand, 3'", {
  LENGTHS_TO_TEST <- c(100, 500)

  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate3primeTxome(txdb_contig_neg_3p, maxTxLength=n, quiet = T)

    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 1)
    expect_equal(n_exons, 1)
  }
})

test_that("non-identical transcripts are retained, positive strand, 3'", {
  LENGTHS_TO_TEST <- c(900)

  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate3primeTxome(txdb_contig_3p, maxTxLength=n, quiet = T)

    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 2)
    expect_equal(n_exons, 2)
  }
})

test_that("non-identical transcripts are retained, negative strand, 3'", {
  LENGTHS_TO_TEST <- c(900)

  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate3primeTxome(txdb_contig_neg_3p, maxTxLength=n, quiet = T)

    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 2)
    expect_equal(n_exons, 2)
  }
})

test_that("APA transcripts are retained, positive strand, 3'", {
  LENGTHS_TO_TEST <- c(100, 500)

  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate3primeTxome(txdb_contig_neg_3p_inv, maxTxLength=n, quiet = T)

    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 2)
    expect_equal(n_exons, 2)
  }
})

test_that("APA transcripts are retained, negative strand, 3'", {
  LENGTHS_TO_TEST <- c(100, 500)

  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate3primeTxome(txdb_contig_3p_inv, maxTxLength=n, quiet = T)

    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 2)
    expect_equal(n_exons, 2)
  }
})

test_that("identical txs from different genes are retained, positive strand, 3'", {
  LENGTHS_TO_TEST <- c(500)

  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate3primeTxome(txdb_multigene_3p, maxTxLength=n, quiet = T)

    n_genes <- length(transcripts(txdb_res))
    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_genes, 2)
    expect_equal(n_txs, 2)
    expect_equal(n_exons, 2)
  }
})

########
## Tests 5'
########

test_that("identical transcripts are merged, positive strand, 5'", {
  LENGTHS_TO_TEST <- c(100, 500)
  
  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate5primeTxome(txdb_contig_5p, maxTxLength=n, quiet = T)
    
    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 1)
    expect_equal(n_exons, 1)
  }
})

test_that("identical transcripts are merged, negative strand, 5'", {
  LENGTHS_TO_TEST <- c(100, 500)
  
  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate5primeTxome(txdb_contig_neg_5p, maxTxLength=n, quiet = T)
    
    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 1)
    expect_equal(n_exons, 1)
  }
})

test_that("non-identical transcripts are retained, positive strand, 5'", {
  LENGTHS_TO_TEST <- c(900)
  
  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate5primeTxome(txdb_contig_5p, maxTxLength=n, quiet = T)
    
    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 2)
    expect_equal(n_exons, 2)
  }
})

test_that("non-identical transcripts are retained, negative strand, 5'", {
  LENGTHS_TO_TEST <- c(900)
  
  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate5primeTxome(txdb_contig_neg_5p, maxTxLength=n, quiet = T)
    
    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 2)
    expect_equal(n_exons, 2)
  }
})

test_that("APA transcripts are retained, positive strand, 5'", {
  LENGTHS_TO_TEST <- c(100, 500)
  
  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate5primeTxome(txdb_contig_neg_5p_inv, maxTxLength=n, quiet = T)
    
    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 2)
    expect_equal(n_exons, 2)
  }
})

test_that("APA transcripts are retained, negative strand, 5'", {
  LENGTHS_TO_TEST <- c(100, 500)
  
  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate5primeTxome(txdb_contig_5p_inv, maxTxLength=n, quiet = T)
    
    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_txs, 2)
    expect_equal(n_exons, 2)
  }
})

test_that("identical txs from different genes are retained, positive strand, 5'", {
  LENGTHS_TO_TEST <- c(500)
  
  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate5primeTxome(txdb_multigene_5p, maxTxLength=n, quiet = T)
    
    n_genes <- length(transcripts(txdb_res))
    n_txs <- length(transcripts(txdb_res))
    n_exons <- length(exons(txdb_res))
    expect_equal(n_genes, 2)
    expect_equal(n_txs, 2)
    expect_equal(n_exons, 2)
  }
})


