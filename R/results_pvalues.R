
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

#' dataset with p values
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

# draw plot
pi0hist

## ---- pi0histends -----

p_values_unn <- p_values %>% 
  unnest_listcol(bins) %>% 
  group_by(Accession, suppdata_id, annot, bins, pi0) %>% 
  summarise(count = n())

p_values_unn <- p_values %>% 
  select(Accession, suppdata_id, annot, features, pi0) %>% 
  left_join(p_values_unn) %>% 
  mutate(freq = count / features)



# Calculate ecdf and cluster histograms -----------------------------
ecdf_pv <- dims %>% 
  filter(!map_lgl(pvalues, is.null), features >= nrowthreshold) %>% 
  mutate(pvals = map(pvalues, "pvalue"),
         bins = map(pvals, ntile, 44),
         pi0 = map_dbl(pvals, propTrueNull)) %>% 
  filter(!map_lgl(pvals, is.null)) %>%
  mutate(eCDF = map(pvals, ecdf))

ecdf_pv <- ecdf_pv %>% 
  mutate(probs = map(eCDF, function(Fn) Fn(seq(0, 1, 1/44))))
probs_mtrx <- ecdf_pv$probs %>% 
  unlist %>% 
  matrix(nrow = nrow(ecdf_pv), byrow = T)
tsne_out <- Rtsne(probs_mtrx, check_duplicates = FALSE, 
                  perplexity = 50, dims = 2, theta = 0)
as_data_frame(tsne_out$Y) %>% 
  bind_cols(annot = p_values_unn$annot, 
            pi0 = p_values_unn$pi0, 
            Accession = p_values_unn$Accession) %>%
  ggplot(aes(V1, V2, col = pi0)) +
  geom_point(aes(text = Accession)) +
  scale_color_viridis()
library(plotly)
ggplotly(tooltip = c("Accession", "pi0"))
