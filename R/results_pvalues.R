
#' ## Load libs
source("R/_common.R")

#' ## Import tabular supplementary files

## ---- loadst -----
#' Import parsed supplementary files.
#+
st <- read_rds(here("output/suppdata.rds"))

#' Imported files object.
#+
imported_geos <- select(st, Accession) %>%  n_distinct()

st_unnested <- unnest(st, result)

#' Unnest xls file sheet names.
#+
st_unnested <- unnest(st_unnested, sheets)

#' Append sheet names to xls file names for unique table id.
#+
st_unnested <- st_unnested %>% 
  mutate(suppdata_id = case_when(
    str_length(sheets) > 0 ~ str_c(suppfiles, "-sheet-", sheets),
    TRUE ~ suppfiles
  ))

#' Import gsem table
#+
gsem <- read_rds(here("output/gsem.rds"))

#' Pull out GEO series matrixes.
#+
gsem <- mutate(gsem, series_matrix = map(gse, "result"))

#' Remove errored GEO series matrixes.
#+
gsem_error <- filter(gsem, map_lgl(series_matrix, is.null))

#+
gsem <- gsem %>%
  filter(!map_lgl(series_matrix, is.null)) %>% 
  mutate(samples = map_int(series_matrix, ~ ncol(exprs(.x))),
         annot = map_chr(series_matrix, annotation)) %>% 
  select(Accession, annot, series_matrix, samples, everything())

#' Get taxon info.
#+
get_tax <- function(data) {
  as_tibble(data) %>% 
    select(matches("organism|taxid")) %>%
    distinct() %>%
    rename_all(str_replace, "_ch1", "") %>% 
    mutate_all(as.character)
}

taxons <- gsem %>% 
  mutate(pd = map(series_matrix, pData),
         tax = map(pd, get_tax)) %>% 
  select(Accession, annot, tax) %>% 
  unnest()

multiple_taxons <- taxons %>% 
  select(organism_ch2:`organism part:ch1`) %>% 
  apply(., 1, function(x) !all(is.na(x)))

multitaxa <- taxons %>% 
  filter(multiple_taxons) %>% 
  group_by(Accession) %>% 
  transmute(organism, taxid,
         organism2 = pmap(list(organism_ch2, organism.1, organism.2, organism.3), ~unique(na.omit(c(..1, ..2, ..3, ..4)))),
         organism2 = map2(organism2, organism, setdiff),
         organism2 = map(organism2, str_c, collapse = ","),
         taxid2 = pmap(list(taxid_ch2, taxid.1, taxid.2, taxid.3), ~unique(na.omit(c(..1, ..2, ..3, ..4)))),
         taxid2 = map2(taxid2, taxid, setdiff),
         taxid2 = map(taxid2, str_c, collapse = ",")) %>% 
  unnest() %>% 
  ungroup() %>% 
  distinct()
  
taxons <- left_join(taxons, multitaxa) %>% 
  select(Accession, annot, organism, taxid, organism2, taxid2) %>% 
  distinct()

#' Save taxons to csv file.
write_csv(taxons, here("output/taxons.csv"))

#' Match samples to correct table/assay. 
#+
dims <- st_unnested %>% 
  left_join(gsem) %>%
  group_by(suppdata_id) %>%
  mutate(idcols = columns - samples) %>%
  filter(idcols >= 0, min_rank(idcols) == 1) %>%
  ungroup() %>%
  select(-series_matrix_file)

#+
dims %>% 
  group_by(samples) %>% 
  summarise(N = n()) %>% 
  write_csv(here("output/number_of_samples_in_geo.csv"))

#+
n_samples <- dims %>% 
  summarise_at( 
    vars(samples, features), 
    list(mean = mean, 
         median = median, 
         Mode = Mode, 
         min = min, 
         max = max)
  )

#' Filter datasets with p values. Keep only numeric P value sets.
#+
p_value_dims <- dims %>% 
  filter(
    !map_lgl(pvalues, is.null), 
    !map_lgl(pvalues, inherits, "try-error")
  )

## ---- plotdims -----

#' ## Plot features versus samples
dims_tabp <- ggplot() + 
  geom_hex(aes(log10(samples), log10(features), fill = ..count.. / sum(..count..)), data = p_value_dims) +
  geom_hex(aes(log10(samples), log10(features), fill = ..count.. / sum(..count..)), data = dims, alpha = 0.5) +
  labs(x = bquote(N~samples~(log[10])),
       y = bquote(N~features~(log[10]))) +
  geom_hline(yintercept = log10(nrowthreshold), linetype = 2) +
  scale_fill_viridis(name = "Fraction")

legend <- g_legend(dims_tabp)
dims_tabp <- dims_tabp + theme(legend.position = "none")

n_imported_sets <- n_distinct(dims$Accession)

dims_featuresp <- dims %>%
  ggplot(aes(log10(features))) + 
  geom_histogram(bins = 60, alpha = 0.5) +
  geom_histogram(bins = 60, aes(log10(features)), data = p_value_dims) +
  labs(x = bquote(N~features~(log[10])),
       y = "Count") +
  geom_vline(xintercept = log10(nrowthreshold), linetype = 2)

dims_samplesp <- dims %>%
  ggplot(aes(log10(samples))) + 
  geom_histogram(bins = 60, alpha = 0.5) +
  geom_histogram(bins = 60, aes(log10(samples)), data = p_value_dims) +
  labs(x = bquote(N~samples~(log[10])),
       y = "Count")

#' Compose plot.
pg <- lapply(list(dims_featuresp, dims_samplesp, dims_tabp), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = length(pg), widths = c(3, 3, 5))
lwidth <- sum(legend$widths)
pga <- arrangeGrob(pga, legend, 
                   ncol = 2, 
                   widths = unit.c(unit(1, "npc") - lwidth, lwidth))

#' Draw plot.
grid.draw(pga)

## ---- pvalues -----

#' ## Collapse rows from GEO series with multiple annotations.
#' Some of the tables have multiple columns with p-values, we need to assign 
#' unique ids to these p-value sets and respectively update suppdata_ids.
#' add id counts
p_value_dims_suppdata_id <- p_value_dims %>% 
  group_by(suppdata_id) %>% 
  add_count()

#' First, collapse annotations in GEO series with multiple annotations.
multi_annot <- p_value_dims_suppdata_id %>% 
  filter(n > 1) %>%
  mutate(annot = str_c(annot, collapse = ",")) %>% 
  slice(1)

#' Then, merge dataset with collapsed annotations back.
p_value_dims <- p_value_dims_suppdata_id %>% 
  filter(n == 1) %>% 
  bind_rows(multi_annot) %>% 
  select(-n) %>% 
  ungroup()

#'
#' Create P values dataset. Filter doubles.
#+ 
p_values <- p_value_dims %>% 
  mutate(pvalues = map(pvalues, make_unique_colnames),
         pvalues = map(pvalues, select, matches("pvalue")),
         pvalues = map(pvalues, as.list)) %>% 
  select(Accession, suppfiles, suppdata_id, annot, features, columns, samples, pvalues)

#' Subset of pvalues with basemean. 
#' Ad hoc filter to remove uninformative features with less than 100 counts 
#' per 1M reads.
p_values_bm <- p_value_dims %>% 
  mutate(basemean = map(pvalues, "basemean"),
         basemean = map_lgl(basemean, ~ !all(is.na(.x)))) %>% 
  filter(basemean) %>%
  mutate(pvalues = map(pvalues, make_unique_colnames),
         pvalues = map(pvalues, filter, basemean / (sum(basemean) / 1e6) > 100),
         pvalues = map(pvalues, select, matches("pvalue")),
         pvalues = map(pvalues, as.list)) %>% 
  select(Accession, suppfiles, suppdata_id, annot, features, columns, samples, pvalues)

#' Proportion of pvalue GSE sets with basemean.
prop_pvalues_with_basemean <- percent(nrow(p_values_bm) / nrow(p_values))

#' Unnest GSE series with multiple pvalue sets and make unique unique ids, omit NA rows from pvalue sets 
#' for each set.
add_set_to_suppdata_id <- function(pvalue_dataset) {
  pvalue_dataset %>% 
    mutate(pvalues = map(pvalues, ~ na.omit(tibble(set = names(.x), pvalues = .x)))) %>% 
    unnest(pvalues) %>% 
    add_count(suppdata_id) %>% 
    mutate(set = str_extract(set, "\\d+"),
           set = case_when(
             is.na(set) ~ "0",
             TRUE ~ set
           ),
           suppdata_id = if_else(n > 1, str_c(suppdata_id, "-set-", set), suppdata_id)) %>% 
    select(-set, -n)
}

#' Append set id to suppdata_ids
p_values <- p_values %>% 
  add_set_to_suppdata_id %>% 
  filter(map_lgl(pvalues, ~typeof(.x) == "double"))

p_values_bm <- p_values_bm %>% 
  add_set_to_suppdata_id %>% 
  filter(map_lgl(pvalues, ~typeof(.x) == "double"))

#' Curated p value histogram classes ---------------------------------------
#' ## Read in classifications
nnet_classes <- read_csv(here("output/pvalue_histogram_nnet_classification.csv")) 
nnet_classes <- nnet_classes %>% 
  mutate(Type = case_when(
    is.na(human) ~ nnet,
    TRUE ~ human
  ))


#' Pool raw and basemean corrected data.
#+ pvalues-pool
pvalues_pool <- bind_rows(
  raw = p_values, 
  basemean = p_values_bm, 
  .id = "Filter"
)

#' ## Calculate ecdf and cluster histograms
#+ add-ecdf
pvalues_pool <- pvalues_pool %>% 
  inner_join(nnet_classes) %>% 
  mutate(eCDF = map(pvalues, ecdf),
         values = map(eCDF, function(Fn) Fn(seq(0, 1, 3 / nrowthreshold))))

## ---- pi0histends -----

#' Rearrange pvalue probs to columns.
#+ pvalues-bins
pvalues_bins <- pvalues_pool %>% 
  mutate(values = map(values, matrix, nrow = 1),
         values = map(values, as.tibble)) %>% 
  unnest(values) %>% 
  select(Filter, suppdata_id, starts_with("v"))

#' Create formula for pca and run pca.
vars <- paste(grep("V", colnames(pvalues_bins), value = TRUE), collapse = "+")
fo <- formula(paste("~", vars))
pca <- prcomp(fo, data = pvalues_bins)
pca_sum <- summary(pca)
pca_screeplot <- pca_sum$importance %>% 
  as_tibble(rownames = "var") %>% 
  filter(str_detect(var, "Prop")) %>% 
  select(1:5) %>% 
  gather(key, value, PC1:PC4) %>% 
  ggplot(aes(key, value)) +
  geom_point() +
  geom_line(aes(group = 1)) +
  facet_wrap(~ var) +
  labs(x = "",
       y = "Proportion")

#' Plot pca results 2D. Use histogram classes, when human assigned class is 
#' missing use nnet.
#+ merge-pca-nnet
pca_tb <- as_tibble(pca$x)
pc12 <- pvalues_pool %>% 
  bind_cols(pca_tb) %>% 
  select(Filter, suppdata_id, PC1, PC2) %>% 
  left_join(nnet_classes)

barcolors <- viridis::viridis(6)
pca_plot <- ggplot(pc12) +
  geom_point(aes(x = PC1, y = PC2, color = Type)) +
  scale_color_manual(values = barcolors) +
  coord_fixed() +
  theme(legend.position = "bottom",
        legend.key = element_blank(),
        legend.title = element_blank())

#' Cluster histograms based on PCA first three components. 
#' Use default euclidean distance.
#+
hc <- pca$x[,1:3] %>% dist() %>% hclust(method = "complete")
hc_phylo <- ape::as.phylo(hc)
types <- as.numeric(factor(pc12$Type))

ggt <- hc_phylo %>%
  ggtree::ggtree(linetype = 2, 
                 color = "steelblue",
                 layout = "slanted") + 
  ggtree::geom_tippoint(color = barcolors[types]) +
  coord_flip()

#' Compose clusters plot. Needs further tweaking!
#+
pg <- lapply(list(pca_screeplot, pca_plot, ggt), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
lay <- rbind(c(1,1,2,2),
             c(1,1,2,2),
             c(3,3,3,3),
             c(3,3,3,3))
pga <- arrangeGrob(grobs = pg, layout_matrix = lay)

#' Draw pi0 plot.
grid.draw(pga)

## ---- pi0hist -----
#' ## Calculate retrospective power (SRP - shitty retrospective power)
#' 
#' Merge pvalues_pool with histogra classes and calculate pi0 and SRP only for
#' anti-conservative sets.
#' 
#' Filter out one dataset with pvalue threshold info column? (logical)
#' and datasets where number of features is less than nrowthreshold
# calculate pi0, the proportion of true nulls
# devtools::install_github("tpall/SRP")
# alternatively, use remotes::install_github()
library(SRP)
safe_srp <- safely(srp)
calculate_pi0_and_srp <- function(pvalue_dataset) {
  pvalue_dataset %>% 
    mutate(pi0 = map_dbl(pvalues, propTrueNull),
           srp = map(pvalues, safe_srp),
           srp = map(srp, "result"))
}

#' Calculate pi0 and srp only for anti-conservative sets.
#+ pvalues-antic
p_values_antic <- pvalues_pool %>% 
  filter(Type == "anti-conservative") %>% 
  calculate_pi0_and_srp

#' Bind anti-conservative sets back to pooled dataset.
#+ merge-antic
pvalues_pool <- pvalues_pool %>% 
  filter(Type != "anti-conservative") %>% 
  bind_rows(p_values_antic)

#' Calculate bins for table histograms. 
pvalues_pool <- pvalues_pool %>%
  mutate(bins = map(pvalues, ntile, 60))

# Save set of vars from pvalues_pool for use in results_publication.R
pvalues_pool_pub <- pvalues_pool %>% 
  filter(Filter == "raw") %>% 
  select(Accession, suppdata_id, Type, pi0, srp) %>% 
  mutate(has_srp = map_lgl(srp, ~!is.null(.x[[1]]))) %>% 
  split(.$has_srp) %>% 
  map_if(~all(.$has_srp), unnest) %>% 
  bind_rows() %>% 
  select(-srp, -has_srp)

write_csv(pvalues_pool_pub, "output/pvalues_pool_pub.csv")

# Histogram of pi0 distribution.
pi0hist <- ggplot(data = p_values_antic) +
  geom_histogram(aes(x = pi0, y = ..count.. / sum(..count..)), bins = 30) +
  labs(x = bquote(Proportion~of~true~nulls~(pi*0)),
       y = "Fraction of anti-conservative P value sets")

pi0_features <- p_values_antic %>% 
  mutate(samples = case_when(
    samples < 4 ~ "2 to 3",
    between(samples, 4, 6) ~ "4 to 6",
    between(samples, 7, 10) ~ "7 to 10",
    TRUE ~ "11+"
  )) %>%
  ggplot(aes(x = pi0, y = log10(features))) +
  geom_point(aes(color = samples)) +
  geom_smooth(method = 'loess', se = FALSE) +
  scale_color_viridis(discrete = TRUE,
                      name = "N samples", 
                      limits = c("2 to 3", "4 to 6", "7 to 10", "11+")) +
  labs(x = bquote(Proportion~of~true~nulls~(pi*0)),
       y = bquote(N~features~(log[10])))

#' Compose pi0 plot.
pg <- lapply(list(pi0hist, pi0_features), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = length(pg), widths = c(1, 1))

#' Draw pi0 plot.
grid.draw(pga)

## ---- sparklines -----

#' Rearrange nnet classes.
# nnet_classes_spread <- nnet_classes %>% 
#   select(suppdata_id, Filter, Type) %>% 
#   spread(Filter, Type) %>% 
#   select(suppdata_id, Type = raw, `Type filtered` = basemean)
make_spark <- . %>% 
  mutate(values = map(pvalues, ~ hist(.x, breaks = seq(0, 1, 1/40), plot = FALSE)$counts),
         pi0 = digits(pi0, 2)) %>%
  mutate(Histogram = map_chr(values, ~ spk_chr(.x,
                                chartRangeMin = 0,
                                type = "bar"))) %>% 
  select(Accession, suppdata_id, Histogram, Type, pi0, srp)

spark_table <- pvalues_pool %>%
  filter(Filter == "raw") %>% 
  make_spark()

spark_table_bm <- pvalues_pool %>%
  filter(Filter == "basemean") %>% 
  make_spark()

#' Original histogram types.
hist_types <- spark_table %>% 
  ungroup() %>% 
  count(Type) %>%
  mutate(`%` = percent(n / sum(n), digits = 1))

#' Histogram types in basemean subset
hist_types_bm <- spark_table_bm %>%
  ungroup() %>% 
  count(Type) %>%
  mutate(`%` = percent(n / sum(n), digits = 1))

#' Recalculated frequencies after basemean filtering.
type_conversion_rate <- select(spark_table, suppdata_id, Type) %>% 
  inner_join(select(spark_table_bm, suppdata_id, `Type filtered`  = Type)) %>% 
  mutate(`Type conversion` = case_when(
    `Type filtered` == Type & Type == "anti-conservative" ~ "same good",
    `Type filtered` == Type & Type != "anti-conservative" ~ "same bad",
    `Type filtered` != Type & `Type filtered` == "anti-conservative" ~ "improvement",
    (`Type filtered` != "uniform" || `Type filtered` != "anti-conservative") & Type == "anti-conservative" ~ "worse",
    `Type filtered` == "uniform" & Type == "anti-conservative" ~ "effects were lost",
    TRUE ~ "no improvement"
  )) %>% 
  count(`Type conversion`) %>% 
  mutate(`%` = percent(n / sum(n), digits = 1))

## ---- histtypes-tab ----
hist_types_caption <- "Summary of histogram types in supplementary files of GEO HT-seq submissions."
hist_types %>% 
  knitr::kable("html", escape = FALSE, caption = hist_types_caption) %>%
  kable_styling(full_width = FALSE)

## ---- histtypesbm-tab ----
hist_types_bm_caption <- "Histogram types after filtering out non-informative features."
hist_types_bm %>% 
  knitr::kable("html", escape = FALSE, caption = hist_types_bm_caption) %>%
  kable_styling(full_width = FALSE)

## ---- typeconversion-tab ----
type_conversion_cap <- "Proportion of histograms that could be fixed by filtering out non-informative features."
type_conversion_rate %>% 
  knitr::kable("html", escape = FALSE, 
               caption = type_conversion_cap) %>%
  kable_styling(full_width = FALSE)

## ---- pvaluespark-tab ----
#' Merge basemean filtered results with spark_table   
spark_table_bm_renamed <- spark_table_bm %>% 
  select(Accession:pi0) %>% 
  rename_at(vars(Histogram, Type, pi0), str_c, "\nafter filter")

pvalue_spark <- spark_table %>% 
  select(Accession:pi0) %>% 
  left_join(spark_table_bm_renamed) %>% 
  rename(`Supplementary file name\n(with p-value set id)` = suppdata_id) 

#' Add taxon info to table
taxons <- rename_all(taxons, str_to_title)
pvalue_spark <- pvalue_spark %>% left_join(taxons)

pv_spark_caption <- glue::glue("P value histograms and proportion of true nulls. 
                               Histograms are colored according to clustering of
                               their empirical cumulative distribution function 
                               outputs. Supplementary file names for tables 
                               from xls(x) files might be appended with sheet 
                               name. This table contains {nrow(spark_table)} 
                               unique P value histograms from 
                               {n_distinct(spark_table$Accession)} GEO 
                               Accessions. pi0 was calculated using limma::propTrueNull() function.") %>% 
  str_replace_all("\n", "")

pvalue_spark %>%
  knitr::kable("html", escape = FALSE, caption = pv_spark_caption) %>%
  kable_styling(full_width = FALSE)

#' Write table for manual reclassification.
pvalue_spark %>% 
  select(1, 2, 4, 7) %>% 
  write_csv(here("output/pvalue_sets_classes.csv"))

## ---- srp-stats ----

spark_table_split <- spark_table %>% 
  split(!map_lgl(.$srp, is.null))

spark_table_bm_split <- spark_table_bm %>% 
  split(!map_lgl(.$srp, is.null))

anti_conservatives <- which(as.logical(names(spark_table_split)))
spark_table_srp <- spark_table_split %>% 
  pluck(anti_conservatives) %>% 
  unnest() %>% 
  arrange(Accession)

anti_conservatives <- which(as.logical(names(spark_table_bm_split)))
spark_table_bm_srp <- spark_table_bm_split %>% 
  pluck(anti_conservatives) %>% 
  unnest() %>% 
  arrange(Accession) %>%
  rename_at(vars(Histogram, Type, pi0, SRP, pi01, fp, rs, ud), str_c, "\nafter filter")

srp_stats <- full_join(spark_table_srp, spark_table_bm_srp)
# Save for further analysis, fix-rename cols for importing
rename_all(srp_stats, str_replace_all, "\\n| ", "_") %>% 
  write_csv(here("output/srp_stats.csv"))

srp_stats_caption <- "SRP and related stats for anti-conservative P value histograms. pi0 was calculated using limma::propTrueNull() function."
srp_stats %>%
  select(-pi0, -`pi0\nafter filter`) %>% 
  rename_all(str_replace, "1", "") %>% 
  mutate_if(is.double, digits, 2) %>% 
  mutate_at(vars(matches("fp|rs|ud")), digits, 0) %>% 
  knitr::kable("html", escape = FALSE, caption = srp_stats_caption) %>%
  kable_styling(full_width = FALSE)

