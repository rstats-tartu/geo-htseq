
# Load libs
source("R/_common.R")

# Import of tabular supplementary files -----------------------------------

## ---- loadst -----
st <- readRDS("output/suppdata.rds")

# imported files
imported_geos <- select(st, Accession) %>%  n_distinct()

st_unnested <- unnest(st, result)

# Unnest xls file sheet names
st_unnested <- unnest(st_unnested, sheets)

# Append sheet names to xls file names for unique table id
st_unnested <- st_unnested %>% 
  mutate(suppdata_id = case_when(
    str_length(sheets) > 0 ~ str_c(suppfiles, "-sheet-", sheets),
    str_length(sheets) == 0 ~ suppfiles
  ))

## Let's use gsem table
gsem <- readRDS("output/gsem.rds")

# Pull out series matrixes
gsem <- mutate(gsem, series_matrix = map(gse, "result"))

## Remove errored matrixes
gsem_error <- filter(gsem, map_lgl(series_matrix, is.null))

gsem <- gsem %>%
  filter(!map_lgl(series_matrix, is.null)) %>% 
  mutate(samples = map_int(series_matrix, ~ ncol(exprs(.x))),
         annot = map_chr(series_matrix, annotation)) %>% 
  select(Accession, annot, series_matrix, samples, everything())

## Match samples to right table/assay 
dims <- left_join(st_unnested, gsem) %>%
  group_by(suppdata_id) %>%
  mutate(idcols = columns - samples) %>%
  filter(idcols >= 0, idcols == min(idcols)) %>%
  ungroup %>%
  select(-series_matrix_file)

group_by(dims, samples) %>% 
  summarise(N = n()) %>% 
  write_csv("output/number_of_samples_in_geo.csv")

n_samples <- summarise_at(dims, 
                          vars(samples, features), 
                          funs(mean, median, Mode, min, max))

#' Datasets with p values
p_value_dims <- filter(dims, !map_lgl(pvalues, is.null), 
                       !map_lgl(pvalues, inherits, "try-error"))

## ---- plotdims -----

## Plot features versus samples
dims_tabp <- ggplot() + 
  geom_hex(aes(log10(samples), log10(features), fill = ..count.. / sum(..count..)), data = p_value_dims) +
  geom_hex(aes(log10(samples), log10(features), fill = ..count.. / sum(..count..)), data = dims, alpha = 0.5) +
  labs(x = bquote(N~samples~(log[10])),
       y = bquote(N~features~(log[10]))) +
  geom_hline(yintercept = log10(nrowthreshold), linetype = 2) +
  scale_fill_viridis(name = "Fraction")

legend <- g_legend(dims_tabp)
dims_tabp <- dims_tabp + theme(legend.position = "none")

n_imported_sets <- length(unique(dims$Accession))

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

# Compose plot
pg <- lapply(list(dims_featuresp, dims_samplesp, dims_tabp), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = length(pg), widths = c(3, 3, 5))
lwidth <- sum(legend$widths)
pga <- arrangeGrob(pga, legend, 
                   ncol = 2, 
                   widths = unit.c(unit(1, "npc") - lwidth, lwidth))

# draw plot
grid.draw(pga)

## ---- pvalues -----
## 
# plot geo series with p values size distribution over all series

# ad hoc filter to remove uninformative features with less than 100 counts per 1M reads on average
p_values_bm <- mutate(p_value_dims, 
                     basemean = map(pvalues, "basemean"),
                     basemean = map_lgl(basemean, ~ !all(is.na(.x)))) %>% 
  filter(basemean) %>%
  mutate(pvalues = map(pvalues, make_unique_colnames),
         pvalues = map(pvalues, filter, basemean / (sum(basemean) / 1e6) > 100),
         pvalues = map(pvalues, select, matches("pvalue")),
         pvalues = map(pvalues, as.list)) %>% 
  select(Accession, suppfiles, suppdata_id, annot, features, columns, samples, pvalues)

p_values <- p_value_dims %>% 
  mutate(pvalues = map(pvalues, make_unique_colnames),
         pvalues = map(pvalues, select, matches("pvalue")),
         pvalues = map(pvalues, as.list)) %>% 
  select(Accession, suppfiles, suppdata_id, annot, features, columns, samples, pvalues)

p_values <- unnest_listcol(p_values, pvalues)

p_values_bm <- unnest_listcol(p_values_bm, pvalues)

#' Filter out one dataset with pvalue threshold info (logical)
#' and datasets where number of features is less than nrowthreshold
# calculate pi0, the proportion of true nulls
# devtools::install_github("tpall/SRP")
# alternatively, use pacman to (install and) load SRP package
library(SRP)
safe_srp <- safely(srp)
p_values <- p_values %>% 
  filter(map_lgl(pvalues, is.numeric)) %>% 
  mutate(pi0 = map_dbl(pvalues, propTrueNull),
         srp = map(pvalues, safe_srp),
         srp = map(srp, "result"))

p_values_bm <- p_values_bm %>% 
  filter(map_lgl(pvalues, is.numeric)) %>% 
  mutate(pi0 = map_dbl(pvalues, propTrueNull),
         srp = map(pvalues, safe_srp),
         srp = map(srp, "result"))

# Identify pvalue sets with low pi0
p_values_low <- p_values %>% 
  filter(features > nrowthreshold,
         pi0 <= pi0threshold)
p_values_bm_low <- p_values_bm %>% 
  filter(features > nrowthreshold,
         pi0 <= pi0threshold)

p_values_low <- select(p_values,Accession, suppdata_id)
p_values_low <- p_values_low %>%
  mutate(type = 6)
write_csv(p_values_low, "output/histogram_classes_all_efects.csv")

p_values_bm_low <- select(p_values_bm,Accession, suppdata_id)
p_values_bm_low <- p_values_bm_low %>%
  mutate(type = 6)
write_csv(p_values_bm_low, "output/histogram_classes_all_efects_filtered.csv")

# Filter out sets with supposedly bad sets of P values

#p_values <- p_values %>% 
#  filter(features > nrowthreshold,
#         pi0 > pi0threshold)
#
#p_values_bm <- p_values_bm %>% 
#  filter(features > nrowthreshold,
#         pi0 > pi0threshold)

p_values <- p_values %>% 
  filter(features > nrowthreshold)

p_values_bm <- p_values_bm %>% 
  filter(features > nrowthreshold)


## ---- pi0hist -----

#' Calculate bins 
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

pi0_features_recoded <- p_values %>%
  mutate(samples = case_when(
    samples < 4 ~ "2 to 3",
    samples > 3 & samples < 7 ~ "4 to 6",
    samples > 6 & samples < 11 ~ "7 to 10",
    samples > 10 ~ "10+"
  ))

pi0_features_recoded_bm <- p_values_bm %>%
  mutate(samples = case_when(
    samples < 4 ~ "2 to 3",
    samples > 3 & samples < 7 ~ "4 to 6",
    samples > 6 & samples < 11 ~ "7 to 10",
    samples > 10 ~ "10+"
  ))

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

# Compose plot
pg <- lapply(list(pi0hist, pi0_features), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = length(pg), widths = c(1, 1))
lwidth <- sum(legend$widths)
pga <- arrangeGrob(pga, legend, 
                   ncol = 2, 
                   widths = unit.c(unit(1, "npc") - lwidth, lwidth))

# draw plot
grid.draw(pga)

## ---- pi0histends -----

# Calculate ecdf and cluster histograms -----------------------------
p_values <- p_values %>%
  mutate(eCDF = map(pvalues, ecdf))

p_values <- p_values %>% 
  mutate(probs = map(eCDF, function(Fn) Fn(seq(0, 1, 1/40))))
probs_mtrx <- p_values$probs %>% 
  unlist %>% 
  matrix(nrow = nrow(p_values), byrow = T)
rownames(probs_mtrx) <- p_values$Accession

# Ok, let's use default eucl + ward.d 
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

# Merge clusters to p value dataframe and create sparklines -------
treecut <- data_frame(histclus = map_chr(treecut, ~ barcolors[.x]))
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

# Save pvalue histograms in 40 bin format
# Used for classification with maschine learning approach

spark_table_write <- p_values %>%
  mutate(values = map(pvalues, ~ hist(.x, breaks = seq(0, 1, 1/40), plot = FALSE)$counts)) %>%
  mutate(values = map(values, ~ .x/sum(.x))) %>%
  select(Accession,suppdata_id,annot,values)
write_rds(spark_table_write,"output/pvalue_spark_bins.rds")

spark_table_write <- p_values_bm %>%
  mutate(values = map(pvalues, ~ hist(.x, breaks = seq(0, 1, 1/40), plot = FALSE)$counts)) %>%
  mutate(values = map(values, ~ .x/sum(.x))) %>%
  select(Accession,suppdata_id,annot,values)
write_rds(spark_table_write,"output/pvalue_bm_spark_bins.rds")

# Curated p value histogram classes ---------------------------------------

# Load manually assigned classes
his <- read_delim("data/pvalue_hist_UM.csv", 
                  delim = ";", 
                  locale = locale(decimal_mark = ","))
colnames(his) <- c("Accession", "suppdata_id", "pi0", "code", "legend")

# remove column 5/legend 
his <- select(his, -"legend")

# Second set of p value histograms
his_added <- read_delim("data/20180411_histogram_classes.csv", 
                        delim = ";", 
                        locale = locale(decimal_mark = ","))

his_added <- select(his_added, 
       Accession, 
       suppdata_id = `Supplementary file name`,
       pi0 = `True nulls proportion`,
       code)
his_added <- mutate(his_added, pi0 = parse_double(pi0))
his <- full_join(his, his_added)
his <- distinct(his)

# parse histogram types from codes
his <- his %>%
  select(Accession, suppdata_id, code) %>% 
  distinct() %>% 
  mutate(type = case_when(
    str_detect(code, "6") ~ 3,
    str_detect(code, "2") & str_detect(code, "1") ~ 5,
    str_detect(code, "2") ~ 2,
    code == 1 ~ 1,
    code == 0 ~ 0,
    suppdata_id == "GSE90615_DifferentialExpression.xlsx-sheet-1dP-MI_vs_Sfrp2_12dP-MI" ~ 4,
    TRUE ~ 4
  ))

# Histogram types summary table
types_legend <- read_csv("data/pvalue_hist_types.csv")
his <- left_join(his, types_legend)

# Histograms after basemean filtering
his_bm <- read_delim("data/pvalue_histogram_codes_after_filtering.csv", delim = ";")
his_bm <- select(his_bm, Accession, `Supplementary file name`, `True nulls proportion filtered`, code_after_filter)
his_bm <- distinct(his_bm) %>% 
  mutate(type = case_when(
    str_detect(code_after_filter, "6") ~ 3,
    str_detect(code_after_filter, "2") & str_detect(code_after_filter, "1") ~ 5,
    str_detect(code_after_filter, "2") ~ 2,
    code_after_filter == 1 ~ 1,
    code_after_filter == 0 ~ 0,
    TRUE ~ 4
  )) %>% 
  left_join(types_legend)

# Merge with classes and generate output table
spark_table <- spark_table %>%
  ungroup() %>% 
  left_join(his) %>%
  rename('P value histogram' = Histogram,
         'True nulls proportion' = pi0,
         'Supplementary file name' = suppdata_id,
         'Type' = typetext) %>% 
  arrange(Accession) %>% 
  distinct() #%>% 
  #filter(!map_int(Type,is.na)) #if you want to show only manually classified stuff

# Save empty table for manual classification. 
if (!any(str_detect(list.files("output"), "pvalue_histogram_classes.csv"))) {
  spark_table %>% 
    select(-`P value histogram`) %>% 
    write_excel_csv("output/pvalue_histogram_classes.csv")
}

hist_types <- spark_table %>% 
  group_by(Type, Comment = comment) %>% 
  summarise(N = n()) %>%
  ungroup() %>% 
  mutate(`%` = percent(N / sum(N), digits = 1))

# Recalculated frequencies after basemean filtering

type_conversion_rate <- left_join(select(spark_table, Accession, `Supplementary file name`, type), 
                                  select(his_bm, Accession, `Supplementary file name`, type_after_filter = type)) %>% 
  filter(!is.na(type_after_filter)) %>% 
  distinct() %>% 
  mutate(type_conversion = case_when(
    type_after_filter == type ~ "same",
    type_after_filter != type & type_after_filter == 1 ~ "improvement",
    type_after_filter > 1 & type == 1 ~ "worse",
    type_after_filter == 0 & type == 1 ~ "effects were lost",
    TRUE ~ "no improvement"
  )) %>% 
  group_by(type_conversion) %>% 
  summarise(N = n()) %>%
  ungroup() %>% 
  mutate(`%` = percent(N / sum(N), digits = 1))
  
his_bm_converted <- anti_join(select(his_bm, Accession, `Supplementary file name`, type, comment),
                              select(spark_table, Accession, `Supplementary file name`, type, comment))

hist_types_bm <- select(spark_table, Accession, `Supplementary file name`, type, comment) %>% 
  filter(!(`Supplementary file name` %in% his_bm_converted$`Supplementary file name`)) %>% 
  bind_rows(his_bm_converted) %>% 
  left_join(types_legend)

write_csv(hist_types_bm, "output/hist_types_after_bm_filter.csv")

hist_types_bm <- hist_types_bm %>% 
  group_by(Type = typetext, Comment = comment) %>% 
  summarise(`N after filter` = n()) %>%
  ungroup() %>% 
  mutate(`% after filter` = percent(`N after filter` / sum(`N after filter`), digits = 1))

hist_types <- left_join(hist_types, hist_types_bm)

hist_types_caption <- "Summary of histogram types in supplementary files of GEO HT-seq submissions."

hist_types %>% 
  knitr::kable("html", escape = FALSE, caption = hist_types_caption) %>%
  kable_styling(full_width = FALSE)

type_conversion_rate %>% 
  knitr::kable("html", escape = FALSE, caption = "Fraction of histograms that could be fixed by filtering out non-informative features.") %>%
  kable_styling(full_width = FALSE)

pv_hist_caption <- glue::glue("P value histograms and proportion of true nulls. Histograms are colored according to clustering of their empirical cumulative distribution function outputs. Supplementary file names for tables from xls(x) files might be appended with sheet name. This table contains {nrow(spark_table)} unique P value histograms from {length(unique(spark_table$'Supplementary file name'))} supplementary tables related to {length(unique(spark_table$Accession))} GEO Accessions.")

spark_table_bm_renamed <- spark_table_bm %>%
  rename('P value histogram,\nfiltered' = Histogram,
         'True nulls proportion,\nfiltered' = pi0,
         'Supplementary file name' = suppdata_id) %>%
  ungroup() %>% 
  left_join(his_bm) %>%
  rename('Type,\nfiltered' = typetext) %>%
  arrange(Accession) %>% 
  distinct() #%>%
  #filter(!map_int('Type,\nfiltered',is.na)) #used to show only manually classified sets

## type and comment has to be removed to avoid conflict between spark_table and spark_table_bm_renamed  
pvalue_spark <- left_join(select(spark_table,-type,-comment),
                          select(spark_table_bm_renamed,-type,-comment)) %>%
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
  write_csv("output/pvalue_spark_table.csv")

#pvalue_spark <- arrange(pvalue_spark, Type)

(pvalue_spark <- pvalue_spark %>%
  knitr::kable("html", escape = FALSE, caption = pv_hist_caption) %>%
  kable_styling(full_width = FALSE))

write_rds(pvalue_spark, "output/pvalue_spark_table.rds")


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
write_csv(srp_stats, "output/srp_stats.csv")
