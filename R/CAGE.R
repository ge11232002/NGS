### -----------------------------------------------------------------
### promotersFromGFF: build promotersFromGFF from GFF
### Exported!
promotersFromGFF <- function(gffFn, fastaFn, ...,
                             upstream=500L, downstream=500L,
                             output="data.frame"){
  seqlengths <- fasta.seqlengths(fastaFn)
  txdb <- makeTxDbFromGFF(gffFn, chrominfo=data.frame(chrom=names(seqlengths),
                                                      length=seqlengths),
                          ...)
  promotersGranges <- promoters(txdb, upstream=upstream, downstream=downstream,
                                columns=c("tx_id", "tx_name", "gene_id"))
  promotersGranges <- trim(promotersGranges)
  promotersGranges <- sort(promotersGranges)
  promoters_df <- as.data.frame(promotersGranges)
  promoters_df$gene_id <- as.character(promoters_df$gene_id)
  # select only unique promoters for each gene - there can be multiple promoters
  # from different transcripts, but they should not have the same start and end
  promoters_df <- promoters_df[order(promoters_df$seqnames,
                                      promoters_df$start,
                                      promoters_df$strand,
                                      promoters_df$gene_id),
                                ]
  promoters_df <- promoters_df[!duplicated(
                     promoters_df[c("seqnames", "start", "strand")]), ]
  return(promoters_df)
}

### -----------------------------------------------------------------
### assignTCsToPromoters:
### 
assignTCsToPromoters <- function(tc_df, promoters_df, 
                                 keep_all_promoters = FALSE){
  tc_gr <- makeGRangesFromDataFrame(tc_df, keep.extra.columns = TRUE)
  promoters_gr <- makeGRangesFromDataFrame(promoters_df, 
                                           keep.extra.columns = TRUE)
  tc_prom_overlap <- findOverlaps(tc_gr, promoters_gr)
  tc_prom_overlap_df <- as.data.frame(tc_prom_overlap)
  tc_prom_overlap_df$tpm <- tc_gr$tpm[tc_prom_overlap_df$queryHits]
  tc_prom_overlap_df$gene_id <-
    as.character(promoters_gr$gene_id)[tc_prom_overlap_df$subjectHits]
  
  #retain only one transcript per hit:
  tc_prom_overlap_df <-
    tc_prom_overlap_df[!duplicated(
      tc_prom_overlap_df[,c("queryHits", "gene_id")]), ]
  
  # if multiple TCs hit the same gene, retain only the TC with most tags:
  tc_prom_overlap_df <- tc_prom_overlap_df[order(tc_prom_overlap_df$gene_id, - tc_prom_overlap_df$tpm), ]
  tc_prom_overlap_df <- tc_prom_overlap_df[!duplicated(tc_prom_overlap_df[,c("gene_id")]), ]
  tc_df$gene_id <- NA
  tc_df[tc_prom_overlap_df$queryHits, ]$gene_id <-
    promoters_df[tc_prom_overlap_df$subjectHits, ]$gene_id
  genes_df <- promoters_df[ ,c("gene_id"), drop=FALSE]
  genes_df <- genes_df[!duplicated(genes_df), , drop=FALSE]
  genes_with_tcs <- merge(genes_df, tc_df[tc_prom_overlap_df$queryHits,],
                          by = c("gene_id"),
                          all=keep_all_promoters)
  return(genes_with_tcs)
}

















