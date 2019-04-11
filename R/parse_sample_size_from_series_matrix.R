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
is_informative <- function(df) {
    not_all_same <- n_distinct(df) > 1
    not_all_different <- n_distinct(df) < length(df)
    all(not_all_same, not_all_different)
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

#' Unnest sample size data frame and arrange sample size with experiment ids as first three columns. 
treatments <- gsem %>% 
    select(Accession, series_matrix_file, N) %>% 
    unnest() %>%
    select(Accession, series_matrix_file, n, everything())

#' Variable names defining experiment groups.
treatments_vars <- colnames(treatments)
treatments_vars <- x_vars[4:length(x_vars)]

#' Paste together experiment groups whereas removing NAs.
treatments <- treatments %>% 
    unite(col = "group", treatments_vars) %>% 
    mutate_at("group", str_replace_all, "_?NA_?", "")
write_csv(treatments, here("output/sample_size.csv"))

#' Plot __mean number__ of replicates in treatment arm(s).
#+
treatments_summary <- treatments %>%
    group_by(Accession, series_matrix_file) %>%
    summarise_at("n", list("mean", "min", "max", n = "length"))

treatments_summary %>%
    mutate(n = case_when(
        n >= 10 ~ "10+",
        TRUE ~ as.character(n)
    ),
    n = factor(n, , levels = c(as.character(1:9), "10+"))) %>% 
    ggplot() +
    geom_histogram(mapping = aes(x = mean), bins = 30) +
    facet_wrap(~ n, scales = "free") +
    scale_x_log10() +
    labs(x = "Mean sample size", caption = "All sets! Facet label is number of experimental groups.")

#' Sets with p values.
# Import P value stats
pvals_pub <- read_csv("output/pvalues_pool_pub.csv",
                      col_types = "cccdddddd")

pvals_pub %>% 
    select(Accession, suppdata_id, Type, pi0) %>% 
    left_join(treatments_summary) %>% 
    ggplot(aes(Type, mean)) +
    geom_jitter(height = 0.1, size = 0.5) +
    labs(x = "Mean sample size", caption = "Sets with P values!")
    

