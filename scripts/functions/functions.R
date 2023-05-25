map_ID <- function(metadata, files){
  mapped_meta <- transform(
    metadata,
    filename = unlist(files)[adist(metadata$sample_id, files, fixed = FALSE) %>% 
                               sweep(1,0,"==") %>% 
                               apply(1,function(x) which(x)[1])]
  )
  return(mapped_meta)
}

"%noin%" <- Negate("%in%")

block.impute.row_mean <- function(block){
  row_means <- DelayedMatrixStats::rowMeans2(block, na.rm = TRUE)
  for(i in 1:nrow(block)){
    nas <- is.na(block[i,])
    if(any(nas)){
      block[i,nas] <- row_means[i]
    }
  }
  return(block)
}

block.row_naSum <- function(block){
  sumnas = 0
  for(i in 1:nrow(block)){
    nas <- is.na(block[i,])
    if(any(nas)){
      sumnas = sumnas + 1
    }
  }
  return(sumnas)
}

missingness.delayed <- function(delayedmat){
  sum(is.na(delayedmat))/length(delayedmat) * 100
}

complete.cases.sum.delayed <- function(delayedmat, dim, progressbar = FALSE){
  if(dim == "row"){
    rowNaSum <- blockApply(delayedmat, grid = rowAutoGrid(delayedmat, nrow = 10000),
               function(block){block.row_naSum(block)},
               BPPARAM = BiocParallel::MulticoreParam(workers = 10, progressbar = progressbar))
    rowNaSum <- sum(unlist(rowNaSum))
    
    return(nrow(delayedmat) - rowNaSum)
  }
}

block.missingness <- function(delayedmat, progressbar = FALSE){
  nas <- blockApply(delayedmat, grid = rowAutoGrid(delayedmat, nrow = 10000),
                    function(block){sum(is.na(block))},
                    BPPARAM = BiocParallel::MulticoreParam(workers = 10, progressbar = progressbar))
  nas <- sum(unlist(nas))
  return(nas/length(meth_raw_filter_75))
}

block.rowSum <- function(delayedmat, progressbar = FALSE, grid, na.rm = TRUE){
  rowsumblock <- blockApply(delayedmat, grid = grid,
                            function(block){rowSums(block, na.rm = na.rm)},
                            BPPARAM = BiocParallel::MulticoreParam(workers = 10, progressbar = progressbar))
  return(do.call(c, rowsumblock))
}

rowSums_Parallel <- function(DelayedMatrix, cores = 1, nblocks = 1, na.rm = FALSE, progressbar = FALSE){
  rowSums <- blockApply(DelayedMatrix, grid = rowAutoGrid(DelayedMatrix, nrow = nrow(DelayedMatrix)/nblocks), function(block){DelayedMatrixStats::rowSums2(block, na.rm = na.rm)}, 
                        BPPARAM = BiocParallel::MulticoreParam(workers = cores, progressbar = progressbar))
  return(unlist(rowSums))
}

rowMeans_Parallel <- function(DelayedMatrix, cores = 1, nblocks = 1, na.rm = FALSE, progressbar = FALSE){
  rowSums <- blockApply(DelayedMatrix, grid = rowAutoGrid(DelayedMatrix, nrow = nrow(DelayedMatrix)/nblocks), 
                        function(block){DelayedMatrixStats::rowMeans2(block, na.rm = na.rm)}, 
                        BPPARAM = BiocParallel::MulticoreParam(workers = cores, progressbar = progressbar))
  return(unlist(rowSums))
}

rowVars_Parallel <- function(DelayedMatrix, cores = 1, nblocks = 1, na.rm = FALSE, progressbar = FALSE){
  rowVars <- blockApply(DelayedMatrix, grid = rowAutoGrid(DelayedMatrix, nrow = nrow(DelayedMatrix)/nblocks), 
                        function(block){DelayedMatrixStats::rowVars(block, na.rm = na.rm)}, 
                        BPPARAM = BiocParallel::MulticoreParam(workers = cores, progressbar = progressbar))
  return(unlist(rowVars))
}

rowAnyNA_Parallel <- function(DelayedMatrix, cores = 1, nblocks = 1, progressbar = FALSE){
  naIndices <- blockApply(DelayedMatrix, grid = rowAutoGrid(DelayedMatrix, nrow = nrow(DelayedMatrix)/nblocks),
                          function(block){DelayedMatrixStats::rowAnyNAs(block, useNames = TRUE)},
                          BPPARAM = BiocParallel::MulticoreParam(workers = cores, progressbar = progressbar))
  return(unlist(naIndices))
}

realize_Parallel <- function(DelayedMatrix, cores = 1, nblocks = 1, progressbar = FALSE){
  chunks <- blockApply(DelayedMatrix, grid = rowAutoGrid(DelayedMatrix, nrow = nrow(DelayedMatrix)/nblocks),
                       function(block){as.matrix(block)},
                       BPPARAM = BiocParallel::MulticoreParam(workers = cores, progressbar = progressbar))
  
  return(do.call(rbind, chunks))
}

filter_rrbs <- function(se, percent_coverage = .66, min_coverage = 10, 
                        cores = 5, nblocks = 100, 
                        progressbar = FALSE, ignore_samples = "all()",
                        cpg_select = NULL){
  
  
  #Select "All CpGs", "CpG Islands", or "Non CpG Islands"
  if(cpg_select == "cgi"){
      session <- rtracklayer::browserSession("UCSC",url = 'http://genome-euro.ucsc.edu/cgi-bin/')
      GenomeInfoDb::genome(session) <- "mm10"
      query <- rtracklayer::ucscTableQuery(session, track="CpG Islands",table="cpgIslandExt",
                                           range=rtracklayer::GRangesForUCSCGenome("mm10"))
      track <- rtracklayer::track(query)
      se <- IRanges::subsetByOverlaps(se, track)
  }
  if(cpg_select == "noncgi"){
      session <- rtracklayer::browserSession("UCSC",url = 'http://genome-euro.ucsc.edu/cgi-bin/')
      GenomeInfoDb::genome(session) <- "mm10"
      query <- rtracklayer::ucscTableQuery(session, track="CpG Islands",table="cpgIslandExt",
                                           range=rtracklayer::GRangesForUCSCGenome("mm10"))
      track <- rtracklayer::track(query)
      se <- IRanges::subsetByOverlaps(se, track, invert = TRUE)
  }
  if(cpg_select == "allcg"){
    se <- se
  }
  
  message("filtering out x chromosome...")
  xchr <- chrSelectBSseq(se, seqnames = "chrX") %>% granges()
  se <- subsetByOverlaps(se, xchr, invert = TRUE)
  
  cov <- bsseq::getCoverage(se[,eval(parse(text=ignore_samples))], type = "Cov")
  
  message("filtering loci for coverage...")
  keepLoci <- which(rowSums_Parallel(cov >= min_coverage, cores = cores, nblocks = nblocks, na.rm = TRUE, progressbar = progressbar) >= ncol(cov)*percent_coverage)
  
  se <- se[keepLoci,]
  message("done")
  
  meth <- bsseq::getMeth(se[,eval(parse(text=ignore_samples))], type = "raw")
  
  message("filtering loci for variance...")
  var <- rowVars_Parallel(meth, cores = cores, nblocks = nblocks, na.rm = TRUE, progressbar = progressbar)
  o <- order(var, decreasing = TRUE)
  o <- o <= (nrow(meth) * 0.5)
  
  se <- se[o,]
  message("done")
  
  return(se)
}
