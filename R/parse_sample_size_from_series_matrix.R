#' ---
#'  title: Estimate sample size from GEO series matrix files
#'  author: Taavi Päll
#' ---
#'
#' Loading libraries.
#+ libs
library(readr)
library(dplyr)
library(purrr)
library(Biobase)
library(stringr)
library(tidyr)
library(ggplot2)
library(here)

#' Import series matrix dataset.
#+ data
gsem <- read_rds(here("output/gsem.rds"))
gsem <- gsem %>%
    mutate(gse = map(gse, "result"))
gsem <- filter(gsem, !map_lgl(gse, is.null))
gsem <- mutate(gsem, pdata = map(gse, pData))
gsem <- select(gsem, -gse)

#' Helper function to find informative columns in experiment metadata.
#+ infofun
is_informative <- function(x) {
    not_all_same <- n_distinct(x) > 1
    not_all_unique <- n_distinct(x) < length(x)
    all(not_all_same, not_all_unique)
}

#' Helper function to parse sample sizes.
#+ getter
get_treatments <- function(pdata) {
    select_if(pdata, is_informative) %>%
        group_by_all() %>%
        count() %>%
        ungroup()
}

get_treatments_safely <- safely(get_treatments)

#' Compute sample sizes.
#+ compute
gsem <- gsem %>%
    mutate(N = map(pdata, get_treatments_safely),
           N = map(N, "result"))

#+ munge
sample_n <- gsem %>%
    mutate(sample_size = map(N, "n")) %>%
    select(Accession, series_matrix_file, sample_size) %>%
    unnest()

#' Plot __mean number__ of replicates in treatment arm(s).
#+
sample_n %>%
    group_by(Accession, series_matrix_file) %>%
    summarise(sample_size = mean(sample_size)) %>%
    ggplot() +
    geom_histogram(mapping = aes(x = sample_size), bins = 30) +
    scale_x_log10() +
    labs(x = "Mean sample size", caption = "All sets!")
