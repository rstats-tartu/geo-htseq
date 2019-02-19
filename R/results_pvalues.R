
#' ## Load libs
source("R/_common.R")

#' ## Import tabular supplementary files

## ---- loadst -----
#' Import parsed supplementary files.
st <- read_rds(here("output/suppdata.rds"))

#' Imported files object.
imported_geos <- select(st, Accession) %>%  n_distinct()

st_unnested <- unnest(st, result)

#' Unnest xls file sheet names.
st_unnested <- unnest(st_unnested, sheets)

#' Append sheet names to xls file names for unique table id.
st_unnested <- st_unnested %>% 
  mutate(suppdata_id = case_when(
    str_length(sheets) > 0 ~ str_c(suppfiles, "-sheet-", sheets),
    TRUE ~ suppfiles
  ))

#' Import gsem table
gsem <- read_rds(here("output/gsem.rds"))

#' Pull out GEO series matrixes.
gsem <- mutate(gsem, series_matrix = map(gse, "result"))

#' Remove errored GEO series matrixes.
gsem_error <- filter(gsem, map_lgl(series_matrix, is.null))

gsem <- gsem %>%
  filter(!map_lgl(series_matrix, is.null)) %>% 
  mutate(samples = map_int(series_matrix, ~ ncol(exprs(.x))),
         annot = map_chr(series_matrix, annotation)) %>% 
  select(Accession, annot, series_matrix, samples, everything())

#' Match samples to correct table/assay. 
dims <- left_join(st_unnested, gsem) %>%
  group_by(suppdata_id) %>%
  mutate(idcols = columns - samples) %>%
  filter(idcols >= 0, min_rank(idcols) == 1) %>%
  ungroup() %>%
  select(-series_matrix_file)

group_by(dims, samples) %>% 
  summarise(N = n()) %>% 
  write_csv(here("output/number_of_samples_in_geo.csv"))

n_samples <- summarise_at(dims, 
                          vars(samples, features), 
                          list(mean = mean, median = median, Mode = Mode, min = min, max = max))

#' Filter datasets with p values.
p_value_dims <- filter(dims, 
                       !map_lgl(pvalues, is.null), 
                       !map_lgl(pvalues, inherits, "try-error"))

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
#' Plot geo series with p values size distribution over all series.
#' Ad hoc filter to remove uninformative features with less than 100 counts 
#' per 1M reads.
p_values <- p_value_dims %>% 
  mutate(pvalues = map(pvalues, make_unique_colnames),
         pvalues = map(pvalues, select, matches("pvalue")),
         pvalues = map(pvalues, as.list)) %>% 
  select(Accession, suppfiles, suppdata_id, annot, features, columns, samples, pvalues)

#' Subset of pvalues with basemean.
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

#' Unnest GSE series with multiple pvalue sets and make unique unique ids 
#' for each set.
add_set_to_suppdata_id <- function(pvalue_dataset) {
  mutate(pvalue_dataset, pvalues = map(pvalues, ~tibble(set = names(.x), pvalues = .x))) %>% 
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
p_values <- p_values %>% add_set_to_suppdata_id
p_values_bm <- p_values_bm %>% add_set_to_suppdata_id

#' ## Calculate retrospective power (SRP - shitty retrospective power)
#' Filter out one dataset with pvalue threshold info column? (logical)
#' and datasets where number of features is less than nrowthreshold
# calculate pi0, the proportion of true nulls
# devtools::install_github("tpall/SRP")
# alternatively, use remotes::install_github()
library(SRP)
safe_srp <- safely(srp)
calculate_pi0_and_srp <- function(pvalue_dataset) {
  filter(pvalue_dataset, map_lgl(pvalues, is.numeric)) %>% 
    mutate(pi0 = map_dbl(pvalues, propTrueNull),
           srp = map(pvalues, safe_srp),
           srp = map(srp, "result"))
}

p_values <- p_values %>% calculate_pi0_and_srp
p_values_bm <- p_values_bm %>% calculate_pi0_and_srp

#' Filter out too small sets of P values.
p_values <- p_values %>% 
  filter(features > nrowthreshold)

p_values_bm <- p_values_bm %>% 
  filter(features > nrowthreshold)

## ---- pi0hist -----

#' Calculate bins for building a histogram. 
p_values <- p_values %>%
  mutate(bins = map(pvalues, ntile, 60))

p_values_bm <- p_values_bm %>%
  mutate(bins = map(pvalues, ntile, 60))

# Histogram of pi0 distribution
pi0hist <- ggplot(data = p_values) +
  geom_histogram(aes(x = pi0, y = ..count.. / sum(..count..)), bins = 30) +
  labs(x = bquote(Proportion~of~true~nulls~(pi*0)),
       y = "Fraction of P value sets")

pi0hist_bm <- pi0hist %+% p_values_bm
recode_pi0_features <- function(pvalue_dataset) {
  mutate(pvalue_dataset, 
         samples = case_when(
           samples < 4 ~ "2 to 3",
           samples > 3 & samples < 7 ~ "4 to 6",
           samples > 6 & samples < 11 ~ "7 to 10",
           samples > 10 ~ "10+"
         ))
}

pi0_features_recoded <- p_values %>% recode_pi0_features()
pi0_features_recoded_bm <- p_values_bm %>% recode_pi0_features()

pi0_features <- ggplot(pi0_features_recoded, aes(x = pi0, y = log10(features))) +
  geom_point(aes(color = samples)) +
  geom_smooth(method = 'loess', se = FALSE) +
  scale_color_viridis(discrete = TRUE,
                      name = "N samples", 
                      limits = c("2 to 3", "4 to 6", "7 to 10", "10+")) +
  labs(x = bquote(Proportion~of~true~nulls~(pi*0)),
       y = bquote(N~features~(log[10])))

pi0_features_bm <- pi0_features %+% pi0_features_recoded_bm
pi0_features_bm <- pi0_features_bm + theme(legend.position = "none")

legend <- g_legend(pi0_features)
pi0_features <- pi0_features + theme(legend.position = "none")

#' Compose pi0 plot.
pg <- lapply(list(pi0hist, pi0_features), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = length(pg), widths = c(1, 1))
lwidth <- sum(legend$widths)
pga <- arrangeGrob(pga, legend, 
                   ncol = 2, 
                   widths = unit.c(unit(1, "npc") - lwidth, lwidth))

#' Draw pi0 plot.
grid.draw(pga)

## ---- pi0histends -----

#' Save pvalue datasets for classification using machine learning
write_rds(p_values, here("output/pvalues.rds"))
write_rds(p_values_bm, here("output/pvalues_bm.rds"))

#' ## Calculate ecdf and cluster histograms
p_values <- p_values %>%
  mutate(eCDF = map(pvalues, ecdf))

p_values <- p_values %>% 
  mutate(probs = map(eCDF, function(Fn) Fn(seq(0, 1, 1/40))))
probs_mtrx <- p_values$probs %>% 
  unlist %>% 
  matrix(nrow = nrow(p_values), byrow = T)
rownames(probs_mtrx) <- p_values$Accession

#' Use default euclidean distance 
# barcolors <- c("#070d0d", "#feb308", "#9b5fc0", "#6ecb3c", "#1d5dec", "#fe4b03")
hc <- probs_mtrx %>% dist(method = "euclidean") %>% hclust() 
treecut <- hc %>% cutree(k = 6)
hc_phylo <- hc %>% ape::as.phylo()

barcolors <- viridis::viridis(6)
ggt <- hc_phylo %>%
  ggtree::ggtree(linetype = 2, 
                 color = "steelblue",
                 layout = "circular") + 
  ggtree::geom_tippoint(color = barcolors[treecut], size = 1) +
  ggtree::geom_tiplab(aes(angle = angle), size = 1.1, hjust = -0.2)

ggt

## ---- sparklines -----

#' ## Merge clusters to p value dataframe and create sparklines
treecut <- tibble(histclus = map_chr(treecut, ~ barcolors[.x]))
p_values <- p_values %>% bind_cols(treecut)

spark_table <- p_values %>%
  mutate(values = map(pvalues, ~ hist(.x, breaks = seq(0, 1, 1/40), plot = FALSE)$counts),
         pi0 = digits(pi0, 2)) %>%
  unnest(values) %>%
  group_by(Accession, suppdata_id, pi0, histclus) %>%
  summarise(Histogram = spk_chr(values,
                                chartRangeMin = 0,
                                barColor = histclus,
                                type = "bar"))

spark_table_bm <- p_values_bm %>%
  mutate(values = map(pvalues, ~ hist(.x, breaks = seq(0, 1, 1/40), plot = FALSE)$counts),
         pi0 = digits(pi0, 2)) %>%
  unnest(values) %>%
  group_by(Accession, suppdata_id, pi0) %>%
  summarise(Histogram = spk_chr(values,
                                chartRangeMin = 0,
                                type = "bar"))

# Curated p value histogram classes ---------------------------------------


# Histogram types summary table
#types_legend <- read_csv("data/pvalue_hist_types.csv")

# Load manually assigned classes
his_all <- read_delim(here("data/pvalue_hist_UM_190121.csv"), 
                      delim = ",", 
                      locale = locale(decimal_mark = ","))
# Code used to identify rows for fixing in data/pvalue_hist_UM_190121.csv
# his_all_to_be_updated <- bind_rows(p_values %>% 
#             filter(suppfiles %in% his_all$`Supplementary file name`, suppfiles != suppdata_id) %>% 
#             select(suppfiles, suppdata_id, pi0),
#           p_values_bm %>% 
#             filter(suppfiles %in% his_all$`Supplementary file name`, suppfiles != suppdata_id) %>% 
#             select(suppfiles, suppdata_id, pi0), .id = "id")


# Split manualy assigned classes to filtered and non-filtered
his_bm <- his_all %>%
  select(Accession,
         suppdata_id = `Supplementary file name`,
         pi0 = `True nulls proportion filtered`,
         Type = `Type filtered`) %>%
  filter(!is.na(Type))

his <- his_all %>%
  select(Accession,
         suppdata_id = `Supplementary file name`,
         pi0 = `True nulls proportion`,
         Type) %>%
  filter(!is.na(Type))

# Read in classifications from random forest.
his_pre <- read_csv(here("data/histogram_classification_prefiltration.csv")) 

his_post <- read_csv(here("data/histogram_classification_postfiltration.csv")) 

# combine classifications and rename low pi0 to all_effects
#his <- full_join(distinct(his),distinct(his_pre)) %>%
#  mutate(Type = replace(Type,as.numeric(pi0) <= pi0threshold, "all_effects"))
#his_bm <- full_join(distinct(his_bm),distinct(his_post)) %>%
#  mutate(Type = replace(Type,as.numeric(pi0) <= pi0threshold, "all_effects"))

# Module to use only classification from random forest.
his <- his_pre
his_bm <- his_post

# Merge with classes and generate output table
spark_table <- spark_table %>%
  ungroup() %>% 
  left_join(select(his, -pi0)) %>%
  arrange(Accession) %>% 
  distinct() %>%
  rename('P value histogram' = Histogram,
         'True nulls proportion' = pi0,
         'Supplementary file name' = suppdata_id) %>% 
  arrange(Type) %>%
  filter(!map_int(Type,is.na)) #if you want to show only classified stuff

# Save empty table for manual classification. 
if (!any(str_detect(list.files(here("output")), "pvalue_histogram_classes.csv"))) {
  spark_table %>% 
    select(-`P value histogram`) %>% 
    write_excel_csv(here("output/pvalue_histogram_classes.csv"))
}

hist_types <- spark_table %>% 
  group_by(Type) %>% 
  summarise(N = n()) %>%
  ungroup() %>% 
  mutate(`%` = percent(N / sum(N), digits = 1))

# Recalculated frequencies after basemean filtering
# Does nto work properly as for the same suppdata_id can be several sets which will get combined multiple times,
# therefore increasing number of rows used for comparision. It means that all combinations from the same set will be 
# compared.
type_conversion_rate <- 
  left_join(
  select(spark_table, Accession, `Supplementary file name`, Type), 
  select(his_bm, Accession, `Supplementary file name` = suppdata_id, type_after_filter = Type)
  ) %>% 
  filter(!is.na(type_after_filter)) %>% 
  mutate(type_conversion = case_when(
    type_after_filter == Type ~ "same",
    type_after_filter != Type & type_after_filter == "anti-conservative" ~ "improvement",
    (type_after_filter != "uniform" || type_after_filter != "anti-conservative") & Type == "anti-conservative" ~ "worse",
    type_after_filter == "uniform" & Type == "anti-conservative" ~ "effects were lost",
    TRUE ~ "no improvement"
  )) %>% 
  group_by(type_conversion) %>% 
  summarise(N = n()) %>%
  ungroup() %>% 
  mutate(`%` = percent(N / sum(N), digits = 1))

# Merge with classes and generate output table
spark_table_bm_renamed <- spark_table_bm %>%
  ungroup() %>% 
  left_join(select(his_bm, -pi0)) %>%
  rename('P value histogram,\nfiltered' = Histogram,
         'True nulls proportion,\nfiltered' = pi0,
         'Supplementary file name' = suppdata_id) %>%
  rename('Type,\nfiltered' = Type) %>%
  arrange(Accession) %>% 
  filter(!map_lgl(`Type,\nfiltered`, is.na)) # used to show only manually classified sets

# This part was not correct
# instead of his_mb spark_table_bm_renamed was used
# this is important when reference list is bigger than actual list of sets.
# By using anti_join we keep only sets with Type conversion after basemean filtering.
his_bm_converted <- 
  anti_join(
    select(spark_table_bm_renamed, Accession, `Supplementary file name`, Type = 'Type,\nfiltered'),
    select(spark_table, Accession, `Supplementary file name`, Type)
    )

# This part is not very correct: sets are cross checked at the level 
# of supplementary file names (actually suppdata_id), but often there are more than 1 set per
# suppdata_id
hist_types_bm <- select(spark_table, Accession, `Supplementary file name`, Type) %>% 
  filter(!(`Supplementary file name` %in% his_bm_converted$`Supplementary file name`)) %>% 
  bind_rows(his_bm_converted)

write_csv(hist_types_bm, here("output/hist_types_after_bm_filter.csv"))

hist_types_bm <- hist_types_bm %>% 
  group_by(Type) %>% 
  summarise(`N after filter` = n()) %>%
  ungroup() %>% 
  mutate(`% after filter` = percent(`N after filter` / sum(`N after filter`), digits = 1))

# Merge summary tables with pre and post filtration summaries.
hist_types <- full_join(hist_types, hist_types_bm) 

# replace NA-s (missing types) with zeros.
hist_types[is.na(hist_types)] <- 0

hist_types_caption <- "Summary of histogram types in supplementary files of GEO HT-seq submissions."

hist_types %>% 
  knitr::kable("html", escape = FALSE, caption = hist_types_caption) %>%
  kable_styling(full_width = FALSE)

type_conversion_rate %>% 
  knitr::kable("html", escape = FALSE, caption = "Fraction of histograms that could be fixed by filtering out non-informative features.") %>%
  kable_styling(full_width = FALSE)

pv_hist_caption <- glue::glue("P value histograms and proportion of true nulls. Histograms are colored according to clustering of their empirical cumulative distribution function outputs. Supplementary file names for tables from xls(x) files might be appended with sheet name. This table contains {nrow(spark_table)} unique P value histograms from {length(unique(spark_table$'Supplementary file name'))} supplementary tables related to {length(unique(spark_table$Accession))} GEO Accessions.")


## type and comment has to be removed to avoid conflict between spark_table and spark_table_bm_renamed  
pvalue_spark <- left_join(spark_table, spark_table_bm_renamed) %>%
  select(Accession, 
         `Supplementary file name`, 
         `P value histogram`, 
         Type, 
         `True nulls proportion`,
         `P value histogram,\nfiltered`,
         `Type,\nfiltered`,
         `True nulls proportion,\nfiltered`) 

pvalue_spark %>% 
  select(-matches("histogram")) %>%
  write_csv(here("output/pvalue_spark_table.csv"))

# to group histograms by the type
value_spark <- arrange(pvalue_spark, Type) 

(pvalue_spark <- pvalue_spark %>%
  knitr::kable("html", escape = FALSE, caption = pv_hist_caption) %>%
  kable_styling(full_width = FALSE))

write_rds(pvalue_spark, here("output/pvalue_spark_table.rds"))

# srp stats ---------------------------------------------------------------

p_values_srp <- select(p_values, Accession, suppfiles, suppdata_id, annot, features, samples, pi0, srp) %>% 
  filter(map_lgl(srp, is.data.frame)) %>% 
  unnest(srp) %>% 
  rename(pi0_raw = pi0, SRP_raw = SRP, pi01_raw = pi01, fp_raw = fp, rs_raw = rs, ud_raw = ud)

p_values_srp_bm <- select(p_values_bm, Accession, suppfiles, suppdata_id, annot, features, samples, pi0, srp) %>% 
  filter(map_lgl(srp, is.data.frame)) %>% 
  unnest(srp) %>% 
  rename(pi0_bm = pi0, SRP_bm = SRP, pi01_bm = pi01, fp_bm = fp, rs_bm = rs, ud_bm = ud)

srp_stats <- full_join(p_values_srp, p_values_srp_bm)

write_csv(srp_stats, here("output/srp_stats.csv"))
