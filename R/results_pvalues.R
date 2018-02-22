
# Load libs
source("R/_common.R")

# Import of tabular supplementary files -----------------------------------

## ---- loadst -----
st <- readRDS("output/suppdata.rds")

# imported files
imported_geos <- select(st, Accession) %>%  n_distinct()

st_unnested <- st %>% unnest(result)
st_unnested <- st_unnested %>% unnest(sheets)

# Add sheet names to xls files
st_unnested <- st_unnested %>% 
  mutate(suppdata_id = case_when(
    str_length(sheets) > 0 ~ str_c(suppfiles, "-sheet-", sheets),
    str_length(sheets) == 0 ~ suppfiles
  ))

## Let's use gsem table
gsem <- readRDS("output/gsem.rds")

## Remove errored matrixes
library(Biobase)
gsem_error <- gsem %>%
  filter(map_lgl(series_matrix, inherits, "try-error"))

gsem <- gsem %>%
  filter(!map_lgl(series_matrix, inherits, "try-error")) %>% 
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

n_samples <- dims %>% 
  summarise_at(vars(samples, features), funs(mean, median, Mode, min, max))

#' Datasets with p values
p_value_dims <- dims %>% 
  filter(!map_lgl(pvalues, is.null), 
         !map_lgl(pvalues, inherits, "try-error"))

## ---- plotdims -----

## Plot features versus samples
dims_tabp <- ggplot() + 
  geom_hex(aes(log10(samples), log10(features), fill = ..count.. / sum(..count..)), data = p_value_dims) +
  geom_hex(aes(log10(samples), log10(features), fill = ..count.. / sum(..count..)), data = dims, alpha = 0.5) +
  labs(x = bquote(N~samples~(log[10])),
       y = bquote(N~features~(log[10]))) +
  geom_hline(yintercept = log10(nrowthreshold), linetype = 2) +
  scale_fill_continuous(name = "Fraction", type = "viridis")

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

p_values <- p_value_dims %>% 
  mutate(pvalues = map(pvalues, make_unique_colnames),
         pvalues = map(pvalues, select, matches("pvalue")),
         pvalues = map(pvalues, as.list)) %>% 
  select(Accession, suppfiles, suppdata_id, annot, features, columns, samples, pvalues)

p_values <- p_values %>% 
  unnest_listcol(pvalues)

#' Filter out one dataset with pvalue threshold info (logical)
#' and datasets where number of features is less than nrowthreshold
# calculate pi0, the proportion of true nulls
p_values <- p_values %>% 
  mutate(pi0 = map_dbl(pvalues, propTrueNull))

# Filter out sets with supposedly bad sets of P values
p_values <- p_values %>% 
  filter(map_lgl(pvalues, is.numeric), 
         features > nrowthreshold,
         pi0 > pi0threshold)
#' Dataset with pvaluethreshold
# p_values %>% filter(Accession == "GSE83134")

## ---- pi0hist -----

#' Calculate bins 
p_values <- p_values %>%
  mutate(bins = map(pvalues, ntile, 60))

# Histogram of pi0 distribution
pi0hist <- p_values %>% 
  ggplot(aes(pi0, ..count.. / sum(..count..))) +
  geom_histogram(bins = 30) +
  labs(x = bquote(Proportion~of~true~nulls~(pi*0)),
       y = "Fraction of P value sets")

pi0_features <- p_values %>%
  mutate(samples = case_when(
    samples < 4 ~ "2 to 3",
    samples > 3 & samples < 7 ~ "4 to 6",
    samples > 6 & samples < 11 ~ "7 to 10",
    samples > 10 ~ "10+",
  )) %>% 
  ggplot(aes(pi0, log10(features))) +
  geom_point(aes(color = samples)) +
  geom_smooth(method = 'loess', se = FALSE) +
  scale_color_viridis(discrete = TRUE,
                      name = "N samples", 
                      limits = c("2 to 3", "4 to 6", "7 to 10", "10+")) +
  labs(x = bquote(Proportion~of~true~nulls~(pi*0)),
       y = bquote(N~features~(log[10])))

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
  mutate(values = map(pvalues, ~ hist(.x, breaks = seq(0, 1, 1/44), plot = FALSE)$counts),
         pi0 = digits(pi0, 2)) %>%
  unnest(values) %>%
  group_by(Accession, suppdata_id, pi0, histclus) %>%
  summarise(Histogram = spk_chr(values,
                                chartRangeMin = 0,
                                barColor = histclus,
                                type = "bar"))

# Save empty table for manual classification. 
if (!any(str_detect(list.files("output"), "pvalue_histogram_classes.csv"))) {
  spark_table %>% 
    select(-Histogram) %>% 
    write_excel_csv("output/pvalue_histogram_classes.csv")
}

# Curated p value histogram classes ---------------------------------------

# Load manually assigned classes
his <- read_delim("data/pvalue_hist_UM.csv", 
                  delim = ";", 
                  locale = locale(decimal_mark = ","))
colnames(his) <- c("Accession", "suppdata_id", "pi0", "code")

# create legend table
# code_legend <- his[, 5]
# colnames(code_legend) <- "description"
# code_legend <- code_legend %>% 
#   na.omit() %>% 
#   separate(description, c("code", "legend"), sep = "-") %>% 
#   mutate_all(str_trim) %>% 
#   mutate_at("code", parse_integer)

# parse histogram types from codes
his <- his %>%
  select(Accession, suppdata_id, code) %>% 
  distinct() %>% 
  mutate(type = case_when(
    str_detect(code, "2") ~ 2,
    str_detect(code, "6") ~ 3,
    code == 1 ~ 1,
    code == 0 ~ 0,
    suppdata_id == "GSE90615_DifferentialExpression.xlsx-sheet-1dP-MI_vs_Sfrp2_12dP-MI" ~ 4,
    TRUE ~ 4
  ))

# Histogram types summary table
types_legend <- read_csv("data/pvalue_hist_types.csv")
his <- left_join(his, types_legend)

# Merge with classes and generate output table
spark_table <- spark_table %>%
  ungroup() %>% 
  left_join(his) %>%
  rename('P value histogram' = Histogram,
         'True nulls proportion' = pi0,
         'Supplementary file name' = suppdata_id,
         'Type' = typetext) %>% 
  arrange(Accession) %>% 
  distinct()

hist_types <- spark_table %>% 
  group_by(Type, Comment = comment) %>% 
  summarise(N = n()) %>%
  ungroup() %>% 
  mutate(`%` = percent(N / sum(N), digits = 1))

hist_types_caption <- "Summary of histogram types in supplementary files of GEO HT-seq submissions."

hist_types %>% 
  knitr::kable("html", escape = FALSE, caption = hist_types_caption) %>%
  kable_styling(full_width = FALSE)

pv_hist_caption <- glue::glue("P value histograms and proportion of true nulls. Histograms are colored according to clustering of their empirical cumulative distribution function outputs. Supplementary file names for tables from xls(x) files might be appended with sheet name. This table contains {nrow(spark_table)} unique P value histograms from {length(unique(spark_table$'Supplementary file name'))} supplementary tables related to {length(unique(spark_table$Accession))} GEO Accessions.")

spark_table %>% 
  select(Accession, 
         `Supplementary file name`, 
         `P value histogram`, 
         Type, 
         `True nulls proportion`) %>%
  knitr::kable("html", escape = FALSE, caption = pv_hist_caption) %>%
  kable_styling(full_width = FALSE)
