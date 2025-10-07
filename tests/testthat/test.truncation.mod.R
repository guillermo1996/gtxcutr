library(GenomicRanges)
library(GenomicFeatures)
BiocParallel::register(BiocParallel::SerialParam())

############
## Mock Data
############

## Single Exon Gene
single_exon_gr <- GRanges(
  seqnames=rep("chr1", 3),
  strand="+",
  ranges=IRanges(start=5000,
                 width=200),
  type=c("gene", "transcript", "exon"),
  ID=c("gene_1", "tx_1", "exon_1"),
  Parent=c(NA, "gene_1", "tx_1"),
  gene_id="gene_1",
  tx_id=c(NA, "tx_1", "tx_1"),
  exon_id=c(NA, NA, "exon_1"))

txdb_single_exons <- makeTxDbSafe(single_exon_gr)
txdb_single_exons_neg <- makeTxDbSafe(invertStrand(single_exon_gr))

## Triple Exon Gene
triple_exon_gr <- GRanges(
  seqnames=rep("chr1", 5),
  strand="+",
  ranges=IRanges(start = c(5000, 5000, 5000, 5200, 5400),
                 end = c(5499, 5499, 5099, 5299, 5499)),
  type=c("gene", "transcript", "exon", "exon", "exon"),
  ID=c("gene_1", "tx_1", "exon_1", "exon_2", "exon_3"),
  Parent=c(NA, "gene_1", "tx_1", "tx_1", "tx_1"),
  gene_id="gene_1",
  tx_id=c(NA, "tx_1", "tx_1", "tx_1", "tx_1"),
  exon_id=c(NA, NA, "exon_1", "exon_2", "exon_3"))

txdb_triple_exons <- makeTxDbSafe(triple_exon_gr)
txdb_triple_exons_neg <- makeTxDbSafe(invertStrand(triple_exon_gr))


########
## Tests
########

test_that("single exon truncation, positive/negative strand, 3'", {
  # Positive strand
  txdb_res_w150 <- truncate3primeTxome(txdb_single_exons, maxTxLength = 150, quiet = T)
  txdb_res_w200 <- truncate3primeTxome(txdb_single_exons, maxTxLength = 200, quiet = T)
  txdb_res_w250 <- truncate3primeTxome(txdb_single_exons, maxTxLength = 250, quiet = T)
  
  ## Both w200 and w250 should be the same
  expect_equal_applied(txdb_res_w200, txdb_res_w250, fns=list(
    function (x) { start(genes(x)) },
    function (x) { start(transcripts(x)) },
    function (x) { start(exons(x)) },
    function (x) { end(genes(x)) },
    function (x) { end(transcripts(x)) },
    function (x) { end(exons(x)) })
  )
  
  ## Ending position of w150 and w200 should be the same
  expect_equal_applied(txdb_res_w150, txdb_res_w200, fns=list(
    function (x) { end(genes(x)) },
    function (x) { end(transcripts(x)) },
    function (x) { end(exons(x)) })
  )
  
  ## Ending position against reference
  expect_equal(end(exons(txdb_res_w150)), end(exons(txdb_single_exons)))
  expect_equal(end(exons(txdb_res_w200)), end(exons(txdb_single_exons)))
  expect_equal(end(exons(txdb_res_w250)), end(exons(txdb_single_exons)))
  
  ## Starting position against reference
  expect_equal(start(exons(txdb_res_w200)), start(exons(txdb_single_exons)))
  expect_equal(start(exons(txdb_res_w250)), start(exons(txdb_single_exons)))
  
  ## Exons widths
  expect_equal(unname(sum(width(exonsBy(txdb_res_w150)))), 150)
  expect_equal(unname(sum(width(exonsBy(txdb_res_w200)))), 200)
  expect_equal(unname(sum(width(exonsBy(txdb_res_w250)))), 200)
  
  
  # Negative strand
  txdb_res_w150 <- truncate3primeTxome(txdb_single_exons_neg, maxTxLength = 150, quiet = T)
  txdb_res_w200 <- truncate3primeTxome(txdb_single_exons_neg, maxTxLength = 200, quiet = T)
  txdb_res_w250 <- truncate3primeTxome(txdb_single_exons_neg, maxTxLength = 250, quiet = T)
  
  ## Both w200 and w250 should be the same
  expect_equal_applied(txdb_res_w200, txdb_res_w250, fns=list(
    function (x) { start(genes(x)) },
    function (x) { start(transcripts(x)) },
    function (x) { start(exons(x)) },
    function (x) { end(genes(x)) },
    function (x) { end(transcripts(x)) },
    function (x) { end(exons(x)) })
  )
  
  ## Starting position of w150 and w200 should be the same
  expect_equal_applied(txdb_res_w150, txdb_res_w200, fns=list(
    function (x) { start(genes(x)) },
    function (x) { start(transcripts(x)) },
    function (x) { start(exons(x)) })
  )
  
  ## Start position against reference
  expect_equal(start(exons(txdb_res_w150)), start(exons(txdb_single_exons_neg)))
  expect_equal(start(exons(txdb_res_w200)), start(exons(txdb_single_exons_neg)))
  expect_equal(start(exons(txdb_res_w250)), start(exons(txdb_single_exons_neg)))
  
  ## Starting position against reference
  expect_equal(end(exons(txdb_res_w200)), end(exons(txdb_single_exons_neg)))
  expect_equal(end(exons(txdb_res_w250)), end(exons(txdb_single_exons_neg)))
})


test_that("single exon truncation, positive/negative strand, 5'", {
  # Positive strand
  txdb_res_w150 <- truncate5primeTxome(txdb_single_exons, maxTxLength = 150, quiet = T)
  txdb_res_w200 <- truncate5primeTxome(txdb_single_exons, maxTxLength = 200, quiet = T)
  txdb_res_w250 <- truncate5primeTxome(txdb_single_exons, maxTxLength = 250, quiet = T)
  
  ## Both w200 and w250 should be the same
  expect_equal_applied(txdb_res_w200, txdb_res_w250, fns=list(
    function (x) { start(genes(x)) },
    function (x) { start(transcripts(x)) },
    function (x) { start(exons(x)) },
    function (x) { end(genes(x)) },
    function (x) { end(transcripts(x)) },
    function (x) { end(exons(x)) })
  )
  
  ## Starting position of w150 and w200 should be the same
  expect_equal_applied(txdb_res_w150, txdb_res_w200, fns=list(
    function (x) { start(genes(x)) },
    function (x) { start(transcripts(x)) },
    function (x) { start(exons(x)) })
  )
  
  ## Starting position against reference
  expect_equal(start(exons(txdb_res_w150)), start(exons(txdb_single_exons)))
  expect_equal(start(exons(txdb_res_w200)), start(exons(txdb_single_exons)))
  expect_equal(start(exons(txdb_res_w250)), start(exons(txdb_single_exons)))
  
  ## Ending position against reference
  expect_equal(end(exons(txdb_res_w200)), end(exons(txdb_single_exons)))
  expect_equal(end(exons(txdb_res_w250)), end(exons(txdb_single_exons)))
  
  ## Exons widths
  expect_equal(unname(sum(width(exonsBy(txdb_res_w150)))), 150)
  expect_equal(unname(sum(width(exonsBy(txdb_res_w200)))), 200)
  expect_equal(unname(sum(width(exonsBy(txdb_res_w250)))), 200)
  
  
  # Negative strand
  txdb_res_w150 <- truncate5primeTxome(txdb_single_exons_neg, maxTxLength = 150, quiet = T)
  txdb_res_w200 <- truncate5primeTxome(txdb_single_exons_neg, maxTxLength = 200, quiet = T)
  txdb_res_w250 <- truncate5primeTxome(txdb_single_exons_neg, maxTxLength = 250, quiet = T)
  
  ## Both w200 and w250 should be the same
  expect_equal_applied(txdb_res_w200, txdb_res_w250, fns=list(
    function (x) { start(genes(x)) },
    function (x) { start(transcripts(x)) },
    function (x) { start(exons(x)) },
    function (x) { end(genes(x)) },
    function (x) { end(transcripts(x)) },
    function (x) { end(exons(x)) })
  )
  
  ## Ending position of w150 and w200 should be the same
  expect_equal_applied(txdb_res_w150, txdb_res_w200, fns=list(
    function (x) { end(genes(x)) },
    function (x) { end(transcripts(x)) },
    function (x) { end(exons(x)) })
  )
  
  ## Ending position against reference
  expect_equal(end(exons(txdb_res_w150)), end(exons(txdb_single_exons_neg)))
  expect_equal(end(exons(txdb_res_w200)), end(exons(txdb_single_exons_neg)))
  expect_equal(end(exons(txdb_res_w250)), end(exons(txdb_single_exons_neg)))
  
  ## Starting position against reference
  expect_equal(start(exons(txdb_res_w200)), start(exons(txdb_single_exons_neg)))
  expect_equal(start(exons(txdb_res_w250)), start(exons(txdb_single_exons_neg)))
})


test_that("triple exon truncation, positive/negative strand, 3'", {
  # Positive strand
  txdb_res_w150 <- truncate3primeTxome(txdb_triple_exons, maxTxLength = 150, quiet = T)
  txdb_res_w200 <- truncate3primeTxome(txdb_triple_exons, maxTxLength = 200, quiet = T)
  txdb_res_w250 <- truncate3primeTxome(txdb_triple_exons, maxTxLength = 250, quiet = T)
  
  ## All ending positions are equal
  expect_equal_applied(txdb_res_w150, txdb_res_w200, fns=list(
    function (x) { end(genes(x)) },
    function (x) { end(transcripts(x)) },
    function (x) { max(end(exons(x))) })
  )
  expect_equal_applied(txdb_res_w200, txdb_res_w250, fns=list(
    function (x) { end(genes(x)) },
    function (x) { end(transcripts(x)) },
    function (x) { max(end(exons(x))) })
  )
  
  ## Starting positions of exons are pre-defined
  expect_equal(min(start(exons(txdb_res_w150))), 5250)
  expect_equal(min(start(exons(txdb_res_w200))), 5200)
  expect_equal(min(start(exons(txdb_res_w250))), 5050)
  
  ## Ensure that transcripts and genes are also truncated
  expect_equal(start(transcripts(txdb_res_w150)), 5250)
  expect_equal(start(transcripts(txdb_res_w200)), 5200)
  expect_equal(start(transcripts(txdb_res_w250)), 5050)
  
  expect_equal(start(genes(txdb_res_w150)), 5250)
  expect_equal(start(genes(txdb_res_w200)), 5200)
  expect_equal(start(genes(txdb_res_w250)), 5050)
  
  ## Ensure total exon width
  expect_equal(unname(sum(width(exonsBy(txdb_res_w150)))), 150)
  expect_equal(unname(sum(width(exonsBy(txdb_res_w200)))), 200)
  expect_equal(unname(sum(width(exonsBy(txdb_res_w250)))), 250)
  
  
  # Negative strand
  txdb_res_w150 <- truncate3primeTxome(txdb_triple_exons_neg, maxTxLength = 150, quiet = T)
  txdb_res_w200 <- truncate3primeTxome(txdb_triple_exons_neg, maxTxLength = 200, quiet = T)
  txdb_res_w250 <- truncate3primeTxome(txdb_triple_exons_neg, maxTxLength = 250, quiet = T)
  
  ## All starting positions are equal
  expect_equal_applied(txdb_res_w150, txdb_res_w200, fns=list(
    function (x) { start(genes(x)) },
    function (x) { start(transcripts(x)) },
    function (x) { min(start(exons(x))) })
  )
  expect_equal_applied(txdb_res_w200, txdb_res_w250, fns=list(
    function (x) { start(genes(x)) },
    function (x) { start(transcripts(x)) },
    function (x) { min(start(exons(x))) })
  )
  
  ## Ending positions of exons are pre-defined
  expect_equal(max(end(exons(txdb_res_w150))), 5249)
  expect_equal(max(end(exons(txdb_res_w200))), 5299)
  expect_equal(max(end(exons(txdb_res_w250))), 5449)
  
  ## Ensure that transcripts and genes are also truncated
  expect_equal(end(transcripts(txdb_res_w150)), 5249)
  expect_equal(end(transcripts(txdb_res_w200)), 5299)
  expect_equal(end(transcripts(txdb_res_w250)), 5449)
  
  expect_equal(end(genes(txdb_res_w150)), 5249)
  expect_equal(end(genes(txdb_res_w200)), 5299)
  expect_equal(end(genes(txdb_res_w250)), 5449)
  
  ## Ensure total exon width
  expect_equal(unname(sum(width(exonsBy(txdb_res_w150)))), 150)
  expect_equal(unname(sum(width(exonsBy(txdb_res_w200)))), 200)
  expect_equal(unname(sum(width(exonsBy(txdb_res_w250)))), 250)
})


test_that("triple exon truncation, positive/negative strand, 5'", {
  # Positive strand
  txdb_res_w150 <- truncate5primeTxome(txdb_triple_exons, maxTxLength = 150, quiet = T)
  txdb_res_w200 <- truncate5primeTxome(txdb_triple_exons, maxTxLength = 200, quiet = T)
  txdb_res_w250 <- truncate5primeTxome(txdb_triple_exons, maxTxLength = 250, quiet = T)
  
  ## All starting positions are equal
  expect_equal_applied(txdb_res_w150, txdb_res_w200, fns=list(
    function (x) { start(genes(x)) },
    function (x) { start(transcripts(x)) },
    function (x) { min(start(exons(x))) })
  )
  expect_equal_applied(txdb_res_w200, txdb_res_w250, fns=list(
    function (x) { start(genes(x)) },
    function (x) { start(transcripts(x)) },
    function (x) { min(start(exons(x))) })
  )
  
  ## Ending positions of exons are pre-defined
  expect_equal(max(end(exons(txdb_res_w150))), 5249)
  expect_equal(max(end(exons(txdb_res_w200))), 5299)
  expect_equal(max(end(exons(txdb_res_w250))), 5449)
  
  ## Ensure that transcripts and genes are also truncated
  expect_equal(end(transcripts(txdb_res_w150)), 5249)
  expect_equal(end(transcripts(txdb_res_w200)), 5299)
  expect_equal(end(transcripts(txdb_res_w250)), 5449)
  
  expect_equal(end(genes(txdb_res_w150)), 5249)
  expect_equal(end(genes(txdb_res_w200)), 5299)
  expect_equal(end(genes(txdb_res_w250)), 5449)
  
  ## Ensure total exon width
  expect_equal(unname(sum(width(exonsBy(txdb_res_w150)))), 150)
  expect_equal(unname(sum(width(exonsBy(txdb_res_w200)))), 200)
  expect_equal(unname(sum(width(exonsBy(txdb_res_w250)))), 250)
  
  
  # Negative strand
  txdb_res_w150 <- truncate5primeTxome(txdb_triple_exons_neg, maxTxLength = 150, quiet = T)
  txdb_res_w200 <- truncate5primeTxome(txdb_triple_exons_neg, maxTxLength = 200, quiet = T)
  txdb_res_w250 <- truncate5primeTxome(txdb_triple_exons_neg, maxTxLength = 250, quiet = T)
  
  ## All ending positions are equal
  expect_equal_applied(txdb_res_w150, txdb_res_w200, fns=list(
    function (x) { end(genes(x)) },
    function (x) { end(transcripts(x)) },
    function (x) { max(end(exons(x))) })
  )
  expect_equal_applied(txdb_res_w200, txdb_res_w250, fns=list(
    function (x) { end(genes(x)) },
    function (x) { end(transcripts(x)) },
    function (x) { max(end(exons(x))) })
  )
  
  ## Starting positions of exons are pre-defined
  expect_equal(min(start(exons(txdb_res_w150))), 5250)
  expect_equal(min(start(exons(txdb_res_w200))), 5200)
  expect_equal(min(start(exons(txdb_res_w250))), 5050)
  
  ## Ensure that transcripts and genes are also truncated
  expect_equal(start(transcripts(txdb_res_w150)), 5250)
  expect_equal(start(transcripts(txdb_res_w200)), 5200)
  expect_equal(start(transcripts(txdb_res_w250)), 5050)
  
  expect_equal(start(genes(txdb_res_w150)), 5250)
  expect_equal(start(genes(txdb_res_w200)), 5200)
  expect_equal(start(genes(txdb_res_w250)), 5050)
  
  ## Ensure total exon width
  expect_equal(unname(sum(width(exonsBy(txdb_res_w150)))), 150)
  expect_equal(unname(sum(width(exonsBy(txdb_res_w200)))), 200)
  expect_equal(unname(sum(width(exonsBy(txdb_res_w250)))), 250)
})

