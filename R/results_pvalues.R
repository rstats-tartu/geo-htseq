# Import of tabular supplementary files -----------------------------------

## ---- loadst -----
load("data/st.RData")

st_unnested <- st %>% unnest(result)
st_unnested <- st_unnested %>% unnest(sheets)

## Add sheet names to xls files
st_unnested <- st_unnested %>% 
  mutate(suppdata_id = case_when(
    str_length(sheets) > 0 ~ str_c(suppfiles, "-sheet-", sheets),
    str_length(sheets) == 0 ~ suppfiles
  ))


## Let's use gsem table
load("data/gsem.RData")

## Remove errored matrixes
library(Biobase)
gsem <- gsem %>%
  filter(!map_lgl(gsematrix, inherits, "try-error")) %>% 
  mutate(samples = map_int(gsematrix, ~ncol(exprs(.x))),
         annot = map_chr(gsematrix, annotation)) %>% 
  select(Accession, annot, gsematrix, samples, everything())

## Match samples to right table/assay 
dims <- left_join(st_unnested, gsem) %>%
  group_by(suppdata_id) %>%
  mutate(idcols = columns - samples) %>%
  filter(idcols >= 0, idcols == min(idcols)) %>%
  ungroup %>%
  select(-matrixfiles, -suppfiles) %>%
  distinct

n_samples <- dims %>% 
  summarise_at(vars(samples, features), funs(mean, median, Mode, min, max))

## ---- plotdims -----

## Plot features versus samples
dims_tabp <- dims %>%
  ggplot(aes(log10(samples), log10(features))) + 
  geom_hex() +
  labs(x = bquote(N~samples~(log[10])),
       y = bquote(N~features~(log[10]))) +
  geom_hline(yintercept = log10(nrowthreshold), linetype = 2) +
  scale_fill_continuous(name = "Count")

n_imported_sets <- length(unique(dims$Accession))

#' dataset with p values
p_value_dims <- dims %>% 
  filter(!map_lgl(pvalues, is.null))

n_sets_with_pvalues <- length(unique(p_value_dims$Accession))

dims_featuresp <- dims %>%
  ggplot(aes(log10(features))) + 
  geom_histogram(bins = 60, fill = "gray70") +
  geom_histogram(bins = 60, aes(log10(features)), data = p_value_dims) +
  labs(x = bquote(N~features~(log[10])),
       y = "Count") +
  geom_vline(xintercept = log10(nrowthreshold), linetype = 2)

dims_samplesp <- dims %>%
  ggplot(aes(log10(samples))) + 
  geom_histogram(bins = 60, fill = "gray70") +
  geom_histogram(bins = 60, aes(log10(samples)), data = p_value_dims) +
  labs(x = bquote(N~samples~(log[10])),
       y = "Count")

pg <- lapply(list(dims_tabp, dims_featuresp, dims_samplesp), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = length(pg), widths = c(3, 2, 2))
grid.draw(pga)

## ---- pvalues -----
## 
# plot geo series with p values size distribution over all series

p_values <- dims %>% 
  filter(!map_lgl(pvalues, is.null), features >= nrowthreshold) %>% 
  mutate(pvals = map(pvalues, "pvalue"),
         bins = map(pvals, ntile, 60),
         pi0 = map_dbl(pvals, propTrueNull))

p_values_unn <- p_values %>%
  filter(!map_lgl(pvals, is.null)) %>% 
  select(Accession, suppdata_id, annot, pvals, bins, pi0) %>% 
  unnest()

p_values_unn <- p_values_unn %>% 
  group_by(Accession, suppdata_id, annot, bins, pi0) %>% 
  summarise(count = n())

p_values_unn <- p_values %>% 
  select(Accession, suppdata_id, annot, features, pi0) %>% 
  right_join(p_values_unn) %>% 
  mutate(freq = count / features)

library(Rtsne)
library(viridis)
p_values_unn <- p_values_unn %>% 
  select(Accession, suppdata_id, annot, bins, freq, pi0) %>%
  filter(complete.cases(.)) %>% 
  nest(bins, freq) %>% 
  mutate(data = map(data, ~ spread(.x, bins, freq)))

p_matrix <- p_values_unn$data %>% bind_rows() %>% as.matrix()
# p_matrix <- p_matrix[!duplicated(p_matrix),]
tsne_out <- Rtsne(p_matrix, check_duplicates = FALSE, perplexity = 50, dims = 3)
as_data_frame(tsne_out$Y) %>% 
  bind_cols(annot = p_values_unn$annot, pi0 = p_values_unn$pi0) %>%
  mutate(pi0bin = case_when(
    pi0 == 1 ~ "no effects",
    pi0 > 0.8 & pi0 < 1 ~ "pi0 ok",
    pi0 < 0.8 ~ "too many\neffects"
  )) %>% 
  ggplot(aes(V1, V2, col = pi0bin)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE)

library(tsne)
tsne_res <- tsne(p_matrix, perplexity = 50, epoch = 50)
as_data_frame(tsne_res) %>%
  bind_cols(annot = p_values_unn$annot, pi0 = p_values_unn$pi0) %>%
  mutate(pi0bin = case_when(
    pi0 == 1 ~ "no effects",
    pi0 > 0.8 & pi0 < 1 ~ "pi0 ok",
    pi0 < 0.8 ~ "too many\neffects"
  )) %>% 
  ggplot(aes(V1, V2, col = pi0bin)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE)

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
