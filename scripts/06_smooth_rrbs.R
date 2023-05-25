job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  HDF5Array::setHDF5DumpDir(here("results","h5","HDF5DumpDir"))
  message("loading data")
  rrbs_org <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_org"))
  message("smoothing data")
  rrbs_org.s <- bsseq::BSmooth(BSseq = rrbs_org,
          BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
          verbose = TRUE)
  message("saving data")
  HDF5Array::saveHDF5SummarizedExperiment(x = rrbs_org.s, 
                                          dir = here("results","h5","rrbs_org"), 
                                          replace = FALSE,
                                          prefix = "smooth")
  
  job::export(c())
}, title = "Smooth ORG")

#SMOOTH RRBS GFWT
job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  HDF5Array::setHDF5DumpDir(here("results","h5","HDF5DumpDir"))
  message("loading data")
  rrbs_gfwt <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_gfwt"), prefix = "update")
  message("smoothing data")
  rrbs_gfwt.s <- bsseq::BSmooth(BSseq = rrbs_gfwt,
                               BPPARAM = BiocParallel::MulticoreParam(workers = 40, progressbar = TRUE),
                               verbose = FALSE)
  message("saving data")
  HDF5Array::saveHDF5SummarizedExperiment(x = rrbs_gfwt.s, 
                                          dir = here("results","h5","rrbs_gfwt"), 
                                          replace = FALSE,
                                          prefix = "smooth")
  
  job::export(c())
}, title = "Smooth GFWT")

job::empty({
  source(here::here("scripts", "01_data_preprocessing.R"))
  HDF5Array::setHDF5DumpDir(here("results","h5","HDF5DumpDir"))
  message("loading data")
  rrbs_combine <- HDF5Array::loadHDF5SummarizedExperiment(dir = here("results","h5","rrbs_combine"))
  message("smoothing data")
  rrbs_combine.s <- bsseq::BSmooth(BSseq = rrbs_combine,
                                BPPARAM = BiocParallel::MulticoreParam(workers = 8, progressbar = TRUE),
                                verbose = TRUE)
  message("saving data")
  HDF5Array::saveHDF5SummarizedExperiment(x = rrbs_combine.s, 
                                          dir = here("results","h5","rrbs_combine"), 
                                          replace = FALSE,
                                          prefix = "smooth")
  
  job::export(c())
}, title = "Smooth Combined")
