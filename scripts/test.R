a <- array(runif(1500000), dim=c(10000, 30, 5, 1))
A <- DelayedArray(a)
b <- array(runif(1500000), dim=c(10000, 30, 5, 1))
B <- DelayedArray(b)
C <- aperm(arbind(aperm(A, 4:1), aperm(B, 4:1)), 4:1)


sessionInfo()
#test combineList function bsseq
#Noticed it still loads assays into memory after combination
M <- matrix(0:8, 3, 3)
Cov <- matrix(1:9, 3, 3)
hdf5_M <- writeHDF5Array(M)
hdf5_Cov <- writeHDF5Array(Cov)
hdf5_BS1 <- BSseq(chr = c("chr1", "chr2", "chr1"),
                  pos = c(1, 2, 3),
                  M = hdf5_M,
                  Cov = hdf5_Cov,
                  sampleNames = c("A", "B", "C"))

hdf5_BS1

hdf5_BS2 <- BSseq(chr = c("chr1", "chr1", "chr1"),
                  pos = c(3, 4, 5),
                  M = hdf5_M,
                  Cov = hdf5_Cov,
                  sampleNames = c("D", "E", "F"))

hdf5_BS2

x <- combineList(list(hdf5_BS1, hdf5_BS2), BACKEND = "HDF5Array")
x

HDF5Array::saveHDF5SummarizedExperiment(hdf5_BS2, 
                                        dir = here("results","h5","test2"))
                                        
test1 <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","test1"))
test2 <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","test2"))
test1;test2
x <- bsseq::combineList(test1, test2)
x
getCoverage(x)


system.time(bsseq::combine(test1, test2))
system.time(bsseq::combineList(test1, test2))                                        


library(rtracklayer)
mcols(rrbs_org)$index <- 1:length(rrbs_org)

mm9tomm10 <- import.chain(here("data", "chainfiles","mm9ToMm10.over.chain"))

mm9 <- rowRanges(rrbs_org)
mm9$index <- 1:length(mm9)

mm10 <- mm9 %>% liftOver(mm9tomm10) %>% unlist()

# Subset based on regions that lifted over
rrbs_org.mm10 <- rrbs_org[which(mm10$index %in% mm9$index)]

rowRanges(rrbs_org.mm10) <- mm10

print(glue::glue("{length(rrbs_org.mm10)} out of {length(rrbs_org)} were lifted over"))
print(glue::glue("{length(rrbs_org) - length(rrbs_org.mm10)} did not liftOver"))

genome(rrbs_org.mm10) <- "mm10"

HDF5Array::setHDF5DumpDir(here("results","h5","HDF5DumpDir"))
DelayedArray::setAutoRealizationBackend(BACKEND="HDF5Array")
library(ExperimentHub)
hub <- ExperimentHub()
query(hub, "TENxBrainData")
fname <- hub[["EH1040"]]
tenx <- HDF5Array(filepath = fname, name = "counts")
tenx_subset <- tenx[, 1:1000]
tenx_subset_hdf5 <- writeHDF5Array(tenx_subset)
DelayedArray:::set_verbose_block_processing(verbose = TRUE)
rowAutoGrid(x =tenx_subset, nrow = 10000)
DelayedMatrixStats::colSums2(tenx_subset)
DelayedArray::colSums(tenx_subset)

source(here::here("scripts", "01_data_preprocessing.R"))
se2 <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_compile"), prefix = "filter66")
set.seed(123)
se <- se[sample(1:nrow(se), size = 100000),1:10]
md_rrbs <- read.csv(here("data", "metadata", "rrbs_metadata.csv"), row.names = 1)
gfwt_samples <- md_rrbs %>% filter(batch != "Organoid" & sample_id != "liv_4402") %>% rownames()

se2 <- se2[,gfwt_samples]
se <- se[,1:10]
m <- getMeth(se, type = "raw")
cov <- getCoverage(se, type = c("Cov"))

setAutoGridMaker(function(x) rowAutoGrid(x, nrow = nrow(x)/10))
DelayedArray:::set_verbose_block_processing(verbose = TRUE)
var <- DelayedMatrixStats::rowVars(m, na.rm = TRUE)

keepLoci <- which(rowSums_blockApply_Parallel(cov >= 10) >= ncol(cov)*.66)
keepLoci2 <- which(DelayedMatrixStats::rowSums2(cov >= 10) >= ncol(cov)*.66)

blockApply(m, function(block){DelayedMatrixStats::rowSums2(block)}, 
           BPPARAM = BiocParallel::MulticoreParam(workers = 10, progressbar = TRUE))


filter_rrbs <- function(se, percent_coverage = .66, min_coverage = 10, cores = 5){
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  cov <- bsseq::getCoverage(se, type = "Cov")
  
  message("filtering loci for coverage...")
  keepLoci <- which(rowSums_Parallel(cov >= min_coverage, na.rm = TRUE, cores = cores) >= ncol(cov)*percent_coverage)
  
  se <- se[keepLoci,]
  
  meth <- bsseq::getMeth(se, type = "raw")
  
  message("filtering loci for variance...")
  var <- rowVars_Parallel(meth, na.rm = TRUE, cores = cores)
  o <- order(var, decreasing = TRUE)
  o <- o <= (nrow(meth) * 0.5)
  
  se <- se[o,]
  message("done")
  
  return(se)
}

test <- filter_rrbs(se)

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
DelayedArray::setAutoRealizationBackend(BACKEND="HDF5Array")
answer <- blockApply(meth, grid = rowAutoGrid(meth, block.length=NULL, nrow = 1000),
                     function(block) {block.impute.row_mean(block)}, 
                     BPPARAM = BiocParallel::MulticoreParam(workers = 5, progressbar = TRUE))

answer_compile <- Reduce(rbind, answer)

HDF5Array::writeHDF5Array(answer_compile, 
                          here("results","h5","test.h5"), 
                          with.dimnames = TRUE, name = "meth")


newfiles <- list.files(here("data","covfiles"),pattern="*.cov.gz", full.names = TRUE)
oldfiles <- sample_table_rrbs$filename

sum(newfiles %noin% oldfiles)
newfiles[newfiles %noin% oldfiles] %>% view()

job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  #optimize filtering step
  rrbs_compile <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_compile"))
  
  test_filter <- filter_rrbs(rrbs_compile, percent_coverage = .75, cores = 10, nblocks = 20, progressbar = TRUE, 
                             ignore_samples = "se$batch != 'Organoid' & se$sample_id != 'liv_4402'")
  job::export(c(test_filter))
}, title = "Filtering methylation matrix 75% coverage")

set.seed(123)
rand <- sample(1:nrow(meth), size = 10000)
meth_rand <- meth[rand, 1:20]
output2 <- PCAtools::pca(meth_rand_naomit, BSPARAM = FastAutoParam(), rank = min(c(ncol(meth_rand_naomit), 30)))


ngenes <- 1000
means <- 2^runif(ngenes, 6, 10)
dispersions <- 10/means + 0.2
nsamples <- 50
counts <- matrix(rnbinom(ngenes*nsamples, mu=means, 
                         size=1/dispersions), ncol=nsamples)

# Choosing the number of PCs
lcounts <- log2(counts + 1)
output <- parallelPCA(meth_rand_naomit)
output$n

#### Distance matrix on samples to identify duplicates ####
library(rdist)

job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  set.seed(123)
  
  #Configure delayed array options
  HDF5Array::setHDF5DumpDir(here("results","h5","HDF5DumpDir"))
  DelayedArray::setAutoRealizationBackend(BACKEND="HDF5Array")
  DelayedArray:::set_verbose_block_processing(verbose = TRUE)
  
  #Load in SE and filter
  covfiles <- paste(here("data/covfiles/"), md_rrbs %>% filter(sampletype != "Organoid") %>% pull(filename), sep = "")
  se <- filter_rrbs(se, cores = 40, nblocks = 20, progressbar = TRUE)
  message("measuring dist")
  d <- mia::calculateJSD(se, exprs_values = "M",
                         BPPARAM = BiocParallel::MulticoreParam(workers = 40, progressbar = TRUE))
  job::export(c(d))
})


#### Test Colon Filtering ####
se <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
se <- se[,se$organ %in% c("DCOL", "PCOL", "COL") &
           se$microbiome %in% c("gf", "wt")]
se <- filter_rrbs(se, cores = 40, nblocks = 20, cpg_select = "allcg")

#Filter for aging sites
se <- subsetByOverlaps(se, aging_sites); se

#Get meth
meth <- getMeth(se, type = "raw") %>% as.matrix()
age <- colData(se)[,c("age"), drop = FALSE] %>% as.matrix()
cbind(age, t(meth)) -> x



