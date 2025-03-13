#' @rdname truncateTxome
#' @param txdb an object representing a transcriptome
#' @param maxTxLength the maximum length of resulting transcripts
#' @param ... additional arguments
#'
#' @return a \code{TxDb} object
#' @export
#'
#' @importFrom methods setGeneric
setGeneric("truncateTxome", signature=c("txdb", "maxTxLength"),
           function(txdb, maxTxLength=500, ...) standardGeneric("truncateTxome")
)

#' Truncate Transcriptome
#'
#' @rdname truncateTxome
#'
#' @param txdb a \code{TxDb} object
#' @param maxTxLength the maximum length of transcripts
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
#' txdb_w500 <- truncateTxome(txdb)
#' txdb_w500
#'
#' ## last 100 nts per tx
#' txdb_w100 <- truncateTxome(txdb, maxTxLength=100)
#' txdb_w100
#'
#' @importFrom GenomicRanges GRangesList mcols
#' @importFrom GenomicFeatures exonsBy makeTxDbFromGRanges
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom AnnotationDbi select taxonomyId
#' @importFrom S4Vectors queryHits subjectHits split
#' @importFrom magrittr %>%
#' @importFrom methods setMethod
#' @importFrom plyranges mutate
#' @importFrom tibble as_tibble
#' @importFrom dplyr group_by mutate filter slice ungroup
#' @export
setMethod("truncateTxome", "TxDb", function(txdb,
                                            maxTxLength=500,
                                            BPPARAM=bpparam()) {
    grlExons <- exonsBy(txdb, use.names=TRUE)
    mcols(grlExons) <- NULL  # Remove metadata
    dfTxGene <- AnnotationDbi::select(txdb, keys=names(grlExons), keytype="TXNAME", columns="GENEID")
    mapTxToGene <- setNames(dfTxGene$GENEID, dfTxGene$TXNAME)

    message("Truncating transcripts...")
    clipped <- .clipTranscript_modded2(grlExons, maxTxLength = maxTxLength, BPPARAM = BPPARAM)
    seqinfo(clipped) <- seqinfo(grlExons)
    grlC_clipped <- S4Vectors::split(clipped, mcols(clipped)["transcript_id"])
    
    # clipped2 <- bplapply(grlExons, .clipTranscript, maxTxLength=maxTxLength, BPPARAM=BPPARAM)
    # grlC_clipped2 <- GRangesList(clipped2, compress = T)
    message("Done.")
    
    message("Checking for duplicate transcripts...")
    overlaps <- findOverlaps(grlC_clipped, minoverlap=maxTxLength,
                             ignore.strand=FALSE,
                             drop.self=TRUE, drop.redundant=TRUE)
    
    grlC_clipped_names <- names(grlC_clipped)
    ## ensure genes match
    if (length(overlaps) > 0) {
      idx_genes_match <- mapply(function (idx1, idx2) {
        mapTxToGene[grlC_clipped_names[idx1]] == mapTxToGene[grlC_clipped_names[idx2]]
      }, idx=queryHits(overlaps), idx2=subjectHits(overlaps))
      overlaps <- overlaps[idx_genes_match]
    }
    
    ## get duplicate indices
    duplicates <- unique(queryHits(overlaps))
    if (length(duplicates) > 0) {
      grlC_clipped <- grlC_clipped[-duplicates]
    }
    message(sprintf("Removed %d duplicates.", length(duplicates)))
    
    message("Creating exon ranges...")
    ## flatten with tx_id in metadata
    grExons <- unlist(.mutateEach(grlC_clipped, transcript_id=names(grlC_clipped)))
    names(grExons) <- NULL
    mcols(grExons)["type"] <- "exon"
    
    ## add gene id
    mcols(grExons)["gene_id"] <- mapTxToGene[mcols(grExons)$transcript_id]
    
    ## reindex exon info
    grExons <- sort(grExons)
    mcols(grExons)["exon_id"] <- seq_along(grExons)
    mcols(grExons)["exon_name"] <- NULL
    ## TODO: include `exon_rank`
    
    message("Done.")
    
    message("Creating tx ranges...")
    ## generate transcripts GRanges with clipped bounds
    # grTxs <- unlist(GRangesList(bplapply(clipped, .fillReduce, BPPARAM=BPPARAM)))
    
    grTxs <- unlist(grlC_clipped) %>% 
      tibble::as_tibble() %>% 
      dplyr::group_by(transcript_id) %>% 
      dplyr::mutate(end = max(end)) %>% 
      dplyr::filter(start == min(start)) %>% 
      dplyr::slice(1) %>% 
      dplyr::ungroup() %>%
      dplyr::select(seqnames, start, end, strand, transcript_id) %>% 
      GRanges()
    
    mcols(grTxs)["type"] <- "transcript"
    mcols(grTxs)["gene_id"] <- mapTxToGene[grTxs$transcript_id]
    
    message("Done.")
    
    message("Creating gene ranges...")
    # grGenes <- unlist(GRangesList(bplapply(split(grTxs, grTxs$gene_id), .fillReduce, BPPARAM=BPPARAM)))
    
    grGenes <- grTxs %>% 
      tibble::as_tibble() %>% 
      dplyr::group_by(gene_id) %>% 
      dplyr::mutate(end = max(end)) %>% 
      dplyr::filter(start == min(start)) %>% 
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(seqnames, start, end, strand, gene_id) %>% 
      GRanges()
    
    mcols(grGenes)["type"] <- "gene"
    message("Done.")
    
    dfMetadata <- data.frame(
      name=c("Truncated by", "Maximum Transcript Length"),
      value=c("gtxcutr", maxTxLength)
    )
    
    makeTxDbFromGRanges(c(grGenes, grTxs, grExons),
                        taxonomyId=taxonomyId(txdb),
                        metadata=dfMetadata)
}
)

#' Clip Transcript to Given Length
#'
#' Internal function for operating on individual \code{GRanges}, where ranges
#' represent exons in a transcript. This is designed to be used in an
#' \code{*apply} function over a \code{GRangesList} object.
#'
#' @param gr a \code{GRanges} object
#' @param maxTxLength a positive integer
#'
#' @return the clipped \code{GRanges} object
#'
#' @importFrom GenomicRanges GRanges width strand start end intersect
#' @importFrom IRanges IRanges
#'
.clipTranscript <- function (gr, maxTxLength) {
    if (sum(width(gr)) <= maxTxLength) { ## already short enough
      gr
    } else { ## need to adjust
        ## adjustment is directed
        txStrand <- strand(gr)
        if (all(txStrand == "+")) {
          ## order by 3' ends
          idx <- order(-end(gr))

          ## compute cumulative lengths
          cumLength <- cumsum(width(gr[idx]))

          ## index of exon that exceeds maximum length
          idxLast <- min(which(cumLength > maxTxLength))

          ## compute cutoff (genomic position)
          startNew <- start(gr[idx[idxLast]]) + (cumLength[idxLast] - maxTxLength)

          ## new transcript interval
          grMask <- GRanges(seqnames(gr[1]),
                            IRanges(startNew, max(end(gr))),
                            strand="+")

          ## clip exons with interval
          GenomicRanges::intersect(gr, grMask)
      } else if (all(txStrand == "-")) {
          ## order by 3' ends
          idx <- order(start(gr))

          ## compute cumulative lengths
          cumLength <- cumsum(width(gr[idx]))

          ## index of exon that exceeds maximum length
          idxLast <- min(which(cumLength > maxTxLength))

          ## compute cutoff (genomic position)
          endNew <- end(gr[idx[idxLast]]) - (cumLength[idxLast] - maxTxLength)

          ## new transcript interval
          grMask <- GRanges(seqnames(gr[1]),
                            IRanges(min(start(gr)), endNew),
                            strand="-")

          ## clip exons with interval
          GenomicRanges::intersect(gr, grMask)
      } else {
          warning("Skipping Transcript: Encountered inconsistent strand annotation!", gr)
          gr
      }
    }
}


#' Clip Transcript to Given Length - Modded by Guillermo R.
#'
#' Internal function for operating on individual \code{GRanges}, where ranges
#' represent exons in a transcript. This is designed to be used in an
#' \code{*apply} function over a \code{GRangesList} object.
#'
#' @param gr a \code{GRanges} object
#' @param maxTxLength a positive integer
#'
#' @return the clipped \code{GRanges} object
#'
#' @importFrom GenomicRanges GRanges width strand start end intersect
#' @importFrom IRanges IRanges
#'
.clipTranscript_modded <- function (gr, maxTxLength) {
  if (sum(width(gr)) <= maxTxLength) { ## already short enough
    gr
  } else { ## need to adjust
    ## adjustment is directed
    txStrand <- strand(gr)
    if (all(txStrand == "+")) {
      ## order by 3' ends
      idx <- order(-end(gr))
      
      ## compute cumulative lengths
      cumLength <- cumsum(width(gr[idx]))
      
      ## index of exon that exceeds maximum length
      idxLast <- min(which(cumLength >= maxTxLength))
      
      ## compute cutoff (genomic position)
      start(gr[idx[idxLast]]) <- start(gr[idx[idxLast]]) + (cumLength[idxLast] - maxTxLength)
      
      ## Return object
      gr[idx[idxLast:1]]
    } else if (all(txStrand == "-")) {
      ## order by 3' ends
      idx <- order(start(gr))
      
      ## compute cumulative lengths
      cumLength <- cumsum(width(gr[idx]))
      
      ## index of exon that exceeds maximum length
      idxLast <- min(which(cumLength >= maxTxLength))
      
      ## compute cutoff (genomic position)
      end(gr[idx[idxLast]]) <- end(gr[idx[idxLast]]) - (cumLength[idxLast] - maxTxLength)
      
      ## Return object
      gr[idx[1:idxLast]]
    } else {
      warning("Skipping Transcript: Encountered inconsistent strand annotation!", gr)
      gr
    }
  }
}


#' Clip Transcript to Given Length - Modded 2 by Guillermo R.
#'
#' Internal function for operating on individual \code{GRanges}, where ranges
#' represent exons in a transcript. This is designed to be used in an
#' \code{*apply} function over a \code{GRangesList} object.
#'
#' @param gr a \code{GRanges} object
#' @param maxTxLength a positive integer
#'
#' @return the clipped \code{GRanges} object
#'
#' @importFrom GenomicRanges GRanges width strand start end intersect
#' @importFrom IRanges IRanges
#' @importFrom tibble as_tibble enframe
#' @importFrom dplyr mutate left_join group_by arrange filter lag row_number ungroup bind_rows summarize n_distinct
#' @importFrom forcats fct_inorder
#'
.clipTranscript_modded2 <- function (grlExons, maxTxLength, BPPARAM) {
  exons <- unlist(grlExons)
  if(!"transcript_id" %in% colnames(mcols(exons))){
    mcols(exons) <- NULL
    mcols(exons)["transcript_id"] <- names(exons)
  }
    
  truncateSplits <- function(exons_strand, maxTxLength, BPPARAM){
    ## Split by cores
    transcript_id <- tibble::enframe(unique(exons_strand$transcript_id) %>% 
                               setNames(., seq_along(.)), 
                             name = "core_id", value = "transcript_id") %>% 
      dplyr::mutate(core_id = as.numeric(core_id) %% BPPARAM$workers)
    split_exons <- exons_strand %>% 
      dplyr::left_join(transcript_id, by = "transcript_id") %>% 
      split(., .$core_id)
    
    truncated_split_exons <- BiocParallel::bplapply(
      split_exons,
      BPPARAM = BPPARAM,
      function(split_exons_i) {
        strandTx = split_exons_i[[1, "strand"]]
        
        if(strandTx == "+"){
          split_exons_i %>% 
            dplyr::select(-core_id) %>% 
            dplyr::group_by(transcript_id) %>% 
            dplyr::arrange(-end, .by_group = TRUE) %>% 
            dplyr::mutate(cumLength = cumsum(width)) %>% 
            dplyr::filter(cumLength < maxTxLength | dplyr::lag(cumLength < maxTxLength, default = F) | dplyr::row_number() == 1) %>% 
            dplyr::mutate(start = ifelse(cumLength > maxTxLength, start + (cumLength - maxTxLength), start)) %>% 
            dplyr::select(-cumLength) %>% 
            dplyr::arrange(end, .by_group = TRUE) %>% 
            dplyr::ungroup()
        }else if(strandTx == "-"){
          split_exons_i %>% 
            dplyr::select(-core_id) %>% 
            dplyr::group_by(transcript_id) %>% 
            dplyr::arrange(start, .by_group = TRUE) %>% 
            dplyr::mutate(cumLength = cumsum(width)) %>% 
            dplyr::filter(cumLength < maxTxLength | dplyr::lag(cumLength < maxTxLength, default = F) | dplyr::row_number() == 1) %>% 
            dplyr::mutate(end = ifelse(cumLength > maxTxLength, end - (cumLength - maxTxLength), end)) %>% 
            dplyr::select(-cumLength) %>% 
            dplyr::arrange(start, .by_group = TRUE) %>% 
            dplyr::ungroup()
        }else{
          return(NULL)
        }
      }
    ) %>% dplyr::bind_rows()
  }
  
  exons_df <- tibble::as_tibble(exons) %>% dplyr::mutate(transcript_id = forcats::fct_inorder(as.character(transcript_id)))
  
  ## Validate strand consistency
  multistrand_tx <- exons_df %>% 
    dplyr::group_by(transcript_id) %>% 
    dplyr::summarise(n = dplyr::n_distinct(strand)) %>% 
    dplyr::filter(n > 1)
  
  if(nrow(multistrand_tx) > 0){
    warning("Some transcripts have inconsistend strand annotation! These will be ignored")
    remove_tx <- unique(multistrand_tx$transcript_id)
    exons_df <- exons_df %>% dplyr::filter(!transcript_id %in% remove_tx)
    rm(remove_tx)
  }
  rm(multistrand_tx)
  
  ## Split by strand
  exons_split_df <- split(exons_df, exons_df$strand)
  exons_pos <- exons_split_df[["+"]]
  exons_neg <- exons_split_df[["-"]]
  
  rm(exons_split_df)
  
  truncated_pos <- truncateSplits(exons_pos, maxTxLength, BPPARAM)
  truncated_neg <- truncateSplits(exons_neg, maxTxLength, BPPARAM)
  
  truncated_exons <- dplyr::bind_rows(truncated_pos, truncated_neg) %>% GenomicRanges::GRanges()
  return(truncated_exons)
}


#' Convert GRanges to Single Range
#'
#' @param gr a \code{GRanges} with ranges to be merged.
#' @param validate logical determining whether entries should be checked for compatible
#' seqnames and strands.
#'
#' @return \code{GRanges} with single interval
#'
#' @details The validation assumes seqnames and strand are \code{Rle} objects.
#'
#' @importFrom GenomicRanges seqnames start end strand reduce start<- end<-
#' @importFrom S4Vectors nrun
.fillReduce <- function (gr, validate=TRUE) {
  if (validate) {
    stopifnot(nrun(seqnames(gr)) ==  1,
              nrun(strand(gr)) ==  1)
  }

  ## TODO: Check if faster to construct new GRanges
  ## Current implementation makes retention of seqinfo simple.
  start(gr) <- min(start(gr))
  end(gr) <- max(end(gr))
  reduce(gr)
}
