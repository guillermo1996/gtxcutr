library(GenomicRanges)
library(GenomicFeatures)
BiocParallel::register(BiocParallel::SerialParam())

############
## Mock Data
############

## Single Exon Gene
gr_single_contig_w1000_5p <- GRanges(
  seqnames=rep("chr1", 3),
  strand="+",
  ranges=IRanges(start=5000,
                 width=1000),
  type=c("gene", "transcript", "exon"),
  ID=c("gene_1", "tx_1", "exon_1"),
  Parent=c(NA, "gene_1", "tx_1"),
  gene_id="gene_1",
  tx_id=c(NA, "tx_1", "tx_1"),
  exon_id=c(NA, NA, "exon_1"))

txdb_single_contig_w1000_5p <- makeTxDbFromGRanges(gr_single_contig_w1000_5p)

## Negative Strand
txdb_single_contig_w1000_5p_neg <- makeTxDbFromGRanges(invertStrand(gr_single_contig_w1000_5p))

## Multi Exon
gr_multi_exon_5p <- GRanges(
  seqnames=rep("chr1", 5),
  strand="+",
  ranges=IRanges(start=c(1000, 1000, 1000, 5000, 9000),
                 end=c(9299, 9299, 1299, 5299, 9299)),
  type=c("gene", "transcript", "exon", "exon", "exon"),
  ID=c("gene_1", "tx_1", "exon_1", "exon_2", "exon_3"),
  Parent=c(NA, "gene_1", "tx_1", "tx_1", "tx_1"),
  gene_id="gene_1",
  tx_id=c(NA, "tx_1", "tx_1", "tx_1", "tx_1"),
  exon_id=c(NA, NA, "exon_1", "exon_2", "exon_3"))

txdb_multi_exon_5p <- makeTxDbFromGRanges(gr_multi_exon_5p)

## Negative Strand
txdb_multi_exon_5p_neg <- makeTxDbFromGRanges(invertStrand(gr_multi_exon_5p))

########
## Tests
########

test_that("simple truncation works, positive strand, 5'", {
  LENGTHS_TO_TEST <- c(100, 500)
  
  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate5primeTxome(txdb_single_contig_w1000_5p, maxTxLength=n, quiet = T)
    
    ## correct lengths
    expect_equal(width(genes(txdb_res)), n)
    expect_equal(width(transcripts(txdb_res)), n)
    expect_equal(width(exons(txdb_res)), n)
    
    ## correct 5' ends
    expect_equal_applied(txdb_res, txdb_single_contig_w1000_5p, fns=list(
      function (x) { start(genes(x)) },
      function (x) { start(transcripts(x)) },
      function (x) { start(exons(x)) }))
  }
})

test_that("simple truncation works, negative strand, 5'", {
  LENGTHS_TO_TEST <- c(100, 500)
  
  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate5primeTxome(txdb_single_contig_w1000_5p_neg, maxTxLength=n, quiet = T)
    
    ## correct lengths
    expect_equal(width(genes(txdb_res)), n)
    expect_equal(width(transcripts(txdb_res)), n)
    expect_equal(width(exons(txdb_res)), n)
    
    ## correct 3' ends
    expect_equal_applied(txdb_res, txdb_single_contig_w1000_5p_neg, fns=list(
      function (x) { end(genes(x)) },
      function (x) { end(transcripts(x)) },
      function (x) { end(exons(x)) }))
  }
})

test_that("idempotent, 5'", {
  txdb_res_w500_1 <- truncate5primeTxome(txdb_single_contig_w1000_5p, maxTxLength=500, quiet = T)
  txdb_res_w500_2 <- truncate5primeTxome(txdb_res_w500_1, maxTxLength=500, quiet = T)
  
  expect_equal_applied(txdb_res_w500_1, txdb_res_w500_2, fns=list(
    function (x) { width(genes(x)) },
    function (x) { width(transcripts(x)) },
    function (x) { width(exons(x)) },
    function (x) { start(genes(x)) },
    function (x) { start(transcripts(x)) },
    function (x) { start(exons(x)) },
    function (x) { end(genes(x)) },
    function (x) { end(transcripts(x)) },
    function (x) { end(exons(x)) }))
})

test_that("spliced truncation works, positive strand, 5'", {
  LENGTHS_TO_TEST <- c(100, 500)
  
  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate5primeTxome(txdb_multi_exon_5p, maxTxLength=n, quiet = T)
    
    ## correct total transcript length
    tx_widths <- unname(sum(width(exonsBy(txdb_res))))
    expect_equal(tx_widths, n)
    
    ## correct 3' ends
    expect_equal_applied(txdb_res, txdb_multi_exon_5p, fns=list(
      function (x) { start(genes(x)) },
      function (x) { start(transcripts(x)) },
      function (x) { min(start(exons(x))) }))
  }
})

test_that("spliced truncation works, negative strand, 5'", {
  LENGTHS_TO_TEST <- c(100, 500)
  
  for (n in LENGTHS_TO_TEST) {
    txdb_res <- truncate5primeTxome(txdb_multi_exon_5p_neg, maxTxLength=n, quiet = T)
    
    ## correct total transcript length
    tx_widths <- unname(sum(width(exonsBy(txdb_res))))
    expect_equal(tx_widths, n)
    
    ## correct 3' ends
    expect_equal_applied(txdb_res, txdb_multi_exon_5p, fns=list(
      function (x) { end(genes(x)) },
      function (x) { end(transcripts(x)) },
      function (x) { max(end(exons(x))) }))
  }
})

