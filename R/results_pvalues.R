
# Load libs
source("R/_common.R")

# Import of tabular supplementary files -----------------------------------

## ---- loadst -----
st <- readRDS("output/suppdata.rds")

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
  mutate(samples = map_int(series_matrix, ~ncol(exprs(.x))),
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
  filter(!map_lgl(pvalues, is.null), !map_lgl(pvalues, inherits, "try-error"))

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
  select(Accession, suppfiles, suppdata_id, annot, features, columns, pvalues)

p_values <- p_values %>% 
  unnest_listcol(pvalues)

#' Filter out one dataset with pvalue threshold info (logical)
#' and datasets where number of features is less than nrowthreshold
p_values <- p_values %>% 
  filter(map_lgl(pvalues, ~ is.numeric(.x)), 
         features > nrowthreshold)
#' Dataset with pvaluethreshold
# p_values %>% filter(Accession == "GSE83134")

## ---- pi0hist -----

#' Calculate bins and pi0 
p_values <- p_values %>%
  mutate(bins = map(pvalues, ntile, 60),
         pi0 = map_dbl(pvalues, propTrueNull))

# Histogram of pi0 distribution
pi0hist <- p_values %>% 
  ggplot(aes(pi0, ..count.. / sum(..count..))) +
  geom_histogram(bins = 30) +
  labs(x = bquote(Proportion~of~true~nulls~(pi*0)),
       y = "Fraction of P value sets")

pi0_features <- p_values %>% 
  ggplot(aes(pi0, log10(features))) +
  geom_point() +
  geom_smooth() +
  labs(x = bquote(Proportion~of~true~nulls~(pi*0)),
       y = bquote(N~features~(log[10])))

# Compose plot
pg <- lapply(list(pi0hist, pi0_features), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = length(pg), widths = c(1, 1))

# draw plot
grid.draw(pga)

## ---- pi0histends -----

# Calculate ecdf and cluster histograms -----------------------------
p_values <- p_values %>%
  filter(map_lgl(pvalues, is.numeric)) %>% 
  mutate(eCDF = map(pvalues, ecdf))

p_values <- p_values %>% 
  mutate(probs = map(eCDF, function(Fn) Fn(seq(0, 1, 1/40))))
probs_mtrx <- p_values$probs %>% 
  unlist %>% 
  matrix(nrow = nrow(p_values), byrow = T)

# Ok, let's use default eucl + ward.d 
library(ape)
barcolors <- c("#070d0d", "#feb308", "#9b5fc0", "#6ecb3c", "#1d5dec", "#fe4b03")
hc <- probs_mtrx %>% dist(method = "euclidean") %>% hclust() 
treecut <- hc %>% cutree(k = 6)
hc_phylo <- hc %>% as.phylo()

library(ggtree)
ggt <- hc_phylo %>%
  ggtree(linetype = 2, 
         color = "steelblue") + 
  geom_tippoint(color = barcolors[treecut], size = 1)

gghist <- function(x) {
  ggplot(x, aes(pvalues, group = tip)) +
    geom_freqpoly(bins = 40, colour = barcolors[unique(x$clus)], size = .1) +
    ggimage::theme_transparent()
}

clus_fp <- data_frame(clus = treecut, 
                      tip = hc_phylo$tip.label, 
                      pvalues = p_values$pvalues) %>% 
  unnest() %>% 
  filter(complete.cases(.)) %>% 
  group_by(clus) %>% 
  do(p = gghist(.))

# Create sparklines, inspect histograms --------------------------------------------------

library(sparklines)
bar_sparkline <- function(x, ...){
  commonopts <- list(barWidth = 2, barSpacing = 0, chartRangeMin = 0)
  hist(x, breaks = seq(0, 1, 1/44), plot = FALSE) %>% 
    .$counts %>% 
    sparklines::sparkline(chart_type = "bar", config = c(commonopts, ...))
}


# Merge clusters to p value dataframe -------
p_values <- p_values %>%
  bind_cols(data_frame(histclus = map_chr(treecut, ~barcolors[.x]))) %>% 
  mutate(spark = map2(pvalues, histclus, ~bar_sparkline(.x, barColor = .y)))

