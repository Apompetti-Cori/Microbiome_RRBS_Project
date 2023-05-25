message("preprocessing")
library(bsseq) |> suppressPackageStartupMessages()
library(methylSig) |> suppressPackageStartupMessages()
library(MethylResolver) |> suppressPackageStartupMessages()
library(tidyverse) |> suppressPackageStartupMessages()
library(PCAtools) |> suppressPackageStartupMessages()
library(ggplot2) |> suppressPackageStartupMessages()
library(plotly) |> suppressPackageStartupMessages()
library(here) |> suppressPackageStartupMessages()
library(ggeasy) |> suppressPackageStartupMessages()
library(rhdf5) |> suppressPackageStartupMessages()
library(HDF5Array) |> suppressPackageStartupMessages()
library(DelayedArray) |> suppressPackageStartupMessages()
library(DelayedMatrixStats) |> suppressPackageStartupMessages()
source(here::here("scripts", "functions", "functions.R"))

bedfiles <- list.files(here("data", "bedfiles"), full.names = TRUE, pattern = "*.gz")

names(bedfiles) <- basename(tools::file_path_sans_ext(tools::file_path_sans_ext(bedfiles)))

covfiles <- list.files(here("data", "covfiles"),
                             pattern = "*merged|bismark.cov.gz",
                             full.names = FALSE,
                             include.dirs = TRUE)

md_rrbs <- read.csv(here("data", "metadata", "rrbs_metadata.csv"), row.names = 1)

md_gfwt <- read.csv(here("data", "metadata", "rrbs_metadata.csv"), row.names = 1) %>% filter(batch != "Organoid")

md_encode <- read.table(here("data","metadata","encode_metadata.tsv"), sep = "\t", header = T)
