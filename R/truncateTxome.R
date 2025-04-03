#' @rdname truncateTxome
#' @export
#'
#' @importFrom methods setGeneric
#' @importFrom BiocParallel bpparam
setGeneric("truncateTxome", signature=c("txdb", "maxTxLength"),
           function(txdb, maxTxLength=500, txEnd="3prime", BPPARAM = bpparam()) standardGeneric("truncateTxome")
)


#' @rdname truncateTxome
#' @param quiet silent output messages
#' 
#' @export
#' @importFrom methods setGeneric
#' @importFrom BiocParallel bpparam
setGeneric("truncate3primeTxome", signature=c("txdb", "maxTxLength"),
           function(txdb, maxTxLength=500, BPPARAM = bpparam(), quiet = F) standardGeneric("truncate3primeTxome")
)


#' @rdname truncateTxome
#' @param quiet silent output messages
#' 
#' @export
#' @importFrom methods setGeneric
#' @importFrom BiocParallel bpparam
setGeneric("truncate5primeTxome", signature=c("txdb", "maxTxLength"),
           function(txdb, maxTxLength=300, BPPARAM = bpparam(), quiet = F) standardGeneric("truncate5primeTxome")
)

#' @rdname truncateTxome
#' @export
#' @importFrom methods setGeneric
#' @importFrom BiocParallel bpparam
setMethod("truncate3primeTxome", "TxDb", function(txdb,
                                                  maxTxLength = 500,
                                                  BPPARAM = bpparam(),
                                                  quiet = F) {
  if(quiet){
    suppressMessages(truncateTxome(txdb, maxTxLength = maxTxLength, txEnd = "3prime", BPPARAM = BPPARAM))
  }else{
    truncateTxome(txdb, maxTxLength = maxTxLength, txEnd = "3prime", BPPARAM = BPPARAM)
  }
})

#' @rdname truncateTxome
#' @export
#' @importFrom methods setGeneric
#' @importFrom BiocParallel bpparam
setMethod("truncate5primeTxome", "TxDb", function(txdb,
                                                  maxTxLength = 300,
                                                  BPPARAM = bpparam(),
                                                  quiet = F) {
  if(quiet){
    suppressMessages(truncateTxome(txdb, maxTxLength = maxTxLength, txEnd = "5prime", BPPARAM = BPPARAM))
  }else{
    truncateTxome(txdb, maxTxLength = maxTxLength, txEnd = "5prime", BPPARAM = BPPARAM)
  }
})

#' Truncate Transcriptome
#' 
#' Truncates the provided transcriptome
#'
#' @rdname truncateTxome
#'
#' @param txdb a \code{TxDb} object
#' @param maxTxLength the maximum length of transcripts
#' @param txEnd transcript end to truncate
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether
#'   and how the method should be parallelized.
#' @return a \code{TxDb} object
#'
#' @examples
#' library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
#'
#' ## load annotation
#' txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
#'
#' ## restrict to 'chrI' transcripts
#' seqlevels(txdb) <- c("chrI")
#'
#' ## last 500 nts per tx
#' txdb_w500 <- truncate3primeTxome(txdb)
#' txdb_w500
#'
#' ## last 100 nts per tx
#' txdb_w100 <- truncate3primeTxome(txdb, maxTxLength=100)
#' txdb_w100
#' @export

#' @importFrom BiocParallel bpparam
setMethod("truncateTxome", "TxDb", function(txdb,
                                            maxTxLength = 500,
                                            txEnd = "3prime",
                                            BPPARAM = bpparam()) {
  ############################################################################
  # Ensure correct values of `txEnd`
  valid_3prime = c("3", "3'", "3prime", "3_prime")
  valid_5prime = c("5", "5'", "5prime", "5_prime")
  if(txEnd %in% valid_3prime) txEnd <- "3prime"
  if(txEnd %in% valid_5prime) txEnd <- "5prime"
  
  if(!txEnd %in% c("3prime", "5prime")) stop("txEnd parameter not valid - only '3prime' or '5prime'.")
  
  ############################################################################
  # Split exons by transcripts and create a mapping dictionary from transcript_id
  # to gene_id
  grlExons <- exonsBy(txdb, use.names=TRUE)
  dfTxGene <- AnnotationDbi::select(txdb, keys=names(grlExons), keytype="TXNAME", columns="GENEID")
  mapTxToGene <- setNames(dfTxGene$GENEID, dfTxGene$TXNAME)
  
  ############################################################################
  # Transcript truncation
  message("Truncating transcripts...")
  clipped <- .clipTranscript(grlExons, maxTxLength = maxTxLength, txEnd = txEnd, BPPARAM = BPPARAM)
  seqinfo(clipped) <- seqinfo(grlExons)
  
  ## GR: Split the exons by transcript to facilitate finding overlaps
  grlC_clipped <- S4Vectors::split(clipped, mcols(clipped)["transcript_id"])
  
  message("Done.")
  
  ############################################################################
  # Remove overlapping transcripts
  message("Checking for duplicate transcripts...")
  overlaps <- findOverlaps(grlC_clipped, minoverlap=maxTxLength,
                           ignore.strand=FALSE,
                           drop.self=TRUE, drop.redundant=TRUE)
  grlC_clipped_names <- names(grlC_clipped)
  
  ## GR: Ensure that overlaps are from the same gene
  matched_overlaps <- tibble::as_tibble(overlaps) %>%
    dplyr::mutate(queryGene = mapTxToGene[grlC_clipped_names[queryHits]],
                  subjectGene = mapTxToGene[grlC_clipped_names[subjectHits]]) %>% 
    dplyr::filter(queryGene == subjectGene) %>% 
    dplyr::select(-queryGene, -subjectGene)
  
  duplicates <- unique(matched_overlaps$queryHits)
  if (length(duplicates) > 0) {
    grlC_clipped <- grlC_clipped[-duplicates]
  }
  message(sprintf("Removed %d duplicates.", length(duplicates)))
  
  ############################################################################
  # Create the final exon ranges
  message("Creating exon ranges...")
  
  ## flatten with tx_id in metadata
  grExons <- slot(grlC_clipped, "unlistData") 
  grExons$transcript_id <- as.character(grExons$transcript_id)
  grExons$type <- "exon"
  
  ## add gene id
  grExons$gene_id <- mapTxToGene[grExons$transcript_id]
  
  ## reindex exon info
  grExons <- sort(grExons)
  mcols(grExons)["exon_id"] <- seq_along(grExons)
  mcols(grExons)["exon_name"] <- NULL
  ## TODO: include `exon_rank`
  
  message("Done.")
  
  ############################################################################
  # Create the final transcript ranges
  message("Creating tx ranges...")
  
  grTxs <- .generateTranscriptRanges(grExons)
  mcols(grTxs)["type"] <- "transcript"
  mcols(grTxs)["gene_id"] <- mapTxToGene[grTxs$transcript_id]
  
  message("Done.")
  
  ############################################################################
  # Create the final gene ranges
  message("Creating gene ranges...")
  
  grGenes <- .generateGeneRanges(grExons)
  mcols(grGenes)["type"] <- "gene"
  message("Done.")
  
  ############################################################################
  # Generate the final TxDb object
  dfMetadata <- data.frame(
    name=c("Truncated by", "Maximum Transcript Length"),
    value=c("gtxcutr", maxTxLength)
  )
  
  makeTxDbFromGRanges(c(grGenes, grTxs, grExons),
                      taxonomyId = taxonomyId(txdb),
                      metadata = dfMetadata)
})


#' Clip Transcript to given length
#'
#' Internal function for operating on a \code{GRangesList}, where ranges
#' represent exons in a transcript. This is designed to be used with a
#' \code{GRangesList} object.
#'
#' @param grlExons a \code{CompressedGRangesList} object
#' @param maxTxLength a positive integer
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether
#'   and how the method should be parallelized.
#'
#' @return the clipped \code{GRanges} object
#'
#' @importFrom methods slot
#' @importFrom magrittr %>% %<>%
#' 
.clipTranscript <- function(grlExons, maxTxLength, txEnd, BPPARAM){
  # GR: Merge all exons into a single GRanges while adding the transcript ID
  exons <- slot(.mutateEach(grlExons, transcript_id = names(grlExons)), "unlistData")
  exons$transcript_id <- forcats::fct_inorder(exons$transcript_id)
  # transcript_order <- names(grlExons)
  
  # GR: Convert the GRanges into a tibble and remove transcripts with inconsistent strands
  exons_df <- tibble::as_tibble(exons[, "transcript_id"])
  exons_df <- .pruneInconsistentStrand(exons_df)
  
  # GR: Split the exon tibble into groups by the number of threads specified
  # with BiocParallel and by strand.
  n_cores <- BPPARAM$workers
  
  tx_to_threads <- exons_df %>% 
    dplyr::distinct(transcript_id) %>% 
    dplyr::mutate(thread_id = as.numeric(dplyr::row_number() %% n_cores))
  
  exons_by_strand_threads <- exons_df %>% 
    dplyr::left_join(tx_to_threads, by = "transcript_id", relationship = "many-to-one") %>% 
    dplyr::group_by(thread_id, strand) %>% 
    dplyr::group_split()
  
  # GR: Apply the truncation pipeline to each group individually. The resulting
  # exons are merged into a tibble again.
  exons_truncated <- BiocParallel::bplapply(exons_by_strand_threads, BPPARAM = BPPARAM, function(exons_split){
    # GR: The pipeline to truncate the exons can be the same to both strands if
    # we make the following modifications to the reverse strand:
    #
    # 1. Swap the start and ending positions.
    # 2. Set the genomic coordinates as negatives.
    # 3. Apply the truncation process.
    # 4. Revert the sign of the genomic coordinates and the start/end positions.
    
    split_strand = unique(exons_split$strand)
    if(split_strand == "-") exons_split %<>% dplyr::mutate(start = -start, end = -end) %>% dplyr::rename(start = end, end = start)
    
    # GR: The main idea behind truncation is to extract the exons with a
    # cumulative width less than the desired maxTxLength and then prune the
    # following exon so that the total width per transcript is exactly
    # maxTxLength. No truncation is performed if overall width per transcript is
    # less than maxTxLength.
    #
    # Changes between 3' and 5' truncation are:
    #
    # + Cumulative width starts counting from the transcript start in 5' and
    # from the transcript end in 3'.    #
    # + If truncation happens, it does so in the transcript in 5' and in the
    # transcript beginning in 3'.

    if(txEnd == "3prime"){
      exons_truncated_split <- exons_split %>% 
        dplyr::group_by(transcript_id) %>% 
        dplyr::arrange(-end, .by_group = TRUE) %>% 
        dplyr::mutate(cumLength = cumsum(width)) %>% 
        # Filter the first N + 1 exons, where N is the number of exons with cumulative width less than maxTxLength
        dplyr::filter(dplyr::row_number() <= sum(cumLength < maxTxLength) + 1) %>%  
        # Modify only the exon with a cumulative width higher than maxTxLength
        dplyr::mutate(start = dplyr::if_else(cumLength > maxTxLength, start + (cumLength - maxTxLength), start)) %>%  
        dplyr::arrange(end, .by_group = TRUE) %>% 
        dplyr::select(-cumLength, -thread_id) %>% 
        dplyr::ungroup()
      
    }else if(txEnd == "5prime"){
      exons_truncated_split <- exons_split %>% 
        dplyr::group_by(transcript_id) %>% 
        dplyr::arrange(start, .by_group = TRUE) %>% 
        dplyr::mutate(cumLength = cumsum(width)) %>% 
        # Filter the first N + 1 exons, where N is the number of exons with cumulative width less than maxTxLength
        dplyr::filter(dplyr::row_number() <= sum(cumLength < maxTxLength) + 1) %>% 
        dplyr::mutate(end = dplyr::if_else(cumLength > maxTxLength, end - (cumLength - maxTxLength), end)) %>% 
        dplyr::arrange(start, .by_group = TRUE) %>% 
        dplyr::select(-cumLength, -thread_id) %>% 
        dplyr::ungroup()
    }
    
    # GR: Revert start/end modifications to the reverse strand transcripts.
    if(split_strand == "-") exons_truncated_split %<>% dplyr::mutate(start = -start, end = -end) %>% dplyr::rename(start = end, end = start)
    
    return(exons_truncated_split)
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::arrange(transcript_id)
  
  return(GenomicRanges::GRanges(exons_truncated))
}

#' Prune GRanges of transcripts with inconsistent strand
#'
#' @param exons_df a \code{tibble} of exons with \code{strand} and
#'   \code{transcript_id} fields.
#'
#' @return \code{tibble} with inconsistent transcripts removed
#'
.pruneInconsistentStrand <- function(exons_df){
  # GR: Group by transcripts and filter based on the number of distinct strands.
  multistrand_tx <- exons_df %>% 
    dplyr::group_by(transcript_id) %>% 
    dplyr::summarise(n = dplyr::n_distinct(strand)) %>% 
    dplyr::filter(n > 1)
  
  if(nrow(multistrand_tx) > 0){
    warning("Some transcripts have inconsistend strand annotation! These will be ignored")
    remove_tx <- unique(multistrand_tx$transcript_id)
    exons_df <- exons_df %>% dplyr::filter(!transcript_id %in% remove_tx)
  }
  
  return(exons_df)
}

#' Generates the transcript ranges from a \code{GRanges}
#'
#' @param grlC_clipped a \code{GRanges} of truncated exons.
#'
#' @return \code{GRanges} with transcripts defined from the truncated exons
#'
#' @importFrom methods slot
#' @importFrom GenomicRanges GRanges
.generateTranscriptRanges <- function(grExons){
  grTxs <- grExons %>% 
    tibble::as_tibble() %>% 
    dplyr::select(seqnames, start, end, strand, transcript_id) %>% 
    dplyr::group_by(transcript_id) %>% 
    dplyr::mutate(start = min(start), end = max(end)) %>% 
    dplyr::ungroup() %>% 
    dplyr::distinct() %>% 
    GRanges() %>% 
    sort()
  
  seqinfo(grTxs) <- seqinfo(grExons)
  return(grTxs)
}

#' Generates the genes ranges from a \code{GRanges}
#'
#' @param grlC_clipped a \code{GRanges} of truncated exons.
#'
#' @return \code{GRanges} with genes defined from the truncated exons
#'
#' @importFrom methods slot
#' @importFrom GenomicRanges GRanges
.generateGeneRanges <- function(grExons){
  grGenes <- grExons %>% 
    tibble::as_tibble() %>% 
    dplyr::select(seqnames, start, end, strand, gene_id) %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::mutate(start = min(start), end = max(end)) %>% 
    dplyr::ungroup() %>% 
    dplyr::distinct() %>% 
    GRanges() %>% 
    sort()
  
  seqinfo(grGenes) <- seqinfo(grExons)
  return(grGenes)
}