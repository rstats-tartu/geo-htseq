
# For development only, not run during article compilation

## Load libs
source("R/_common.R")

# Load helper functions
source("lib/helpers.R")

## ---- rna-seq-dynamics ----

# Last date to consider geo series and suppfilenames
last_date <- ymd("2017-06-19")

## R/A01_GEO_query.R
load("data/ds.RData") # mouse and human GEO HT-seq expr datasets

## all HT-seq datasets 
ds_all <- filter(ds, ymd(PDAT) <= last_date)
first_date <- min(ymd(ds$PDAT))

## Human or mouse datasets
ds <- filter(ds_all, str_detect(taxon, "Mus musculus|Homo sapiens"))

## Merge all datasets for plotting
ds_merge <- bind_rows(ds_all, ds, .id = "id") %>% 
  mutate_at("PDAT", ymd)

## Count series with publications
pdat <- ds_merge %>% 
  select(id, PDAT, PubMedIds) %>% 
  mutate(pub = str_length(PubMedIds) != 0) %>% 
  group_by(id, PDAT) %>% 
  summarise(N = n(),
            pub = sum(pub)) %>% 
  mutate_at(vars(N, pub), cumsum)

pdat <- gather(pdat, key, value, -PDAT, -id) 

pdat <- ungroup(pdat) %>% 
  mutate(id = if_else(id == 1, "All taxa", "Human and mouse"))

## Plot submissions
geop <- pdat %>% 
  ggplot(aes(ymd(PDAT), value, linetype = key)) + 
  geom_line() +
  facet_wrap(~id) +
  labs(x = "Publication date", 
       y = "Number of GEO series") +
  scale_linetype_discrete(labels = c("All series","Series with\npublications")) +
  theme(legend.position = c(0.75, 0.77),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank())

## Calculate percent series using human or mouse 
## formattable::percent()
perc_mmhs <- percent(round(nrow(ds) / nrow(ds_all), 1), digits = 0)

# number of publications
ppub_n <- group_by(pdat, id, key) %>% 
  summarise_at("value", max)

# percent with publications
ppub <- ppub_n %>% 
  group_by(id) %>% 
  summarise(ppub = value[key=="pub"] / value[key=="N"]) %>% 
  ungroup %>% 
  filter(str_detect(id, "Human")) %>% 
  .$ppub %>% 
  round(digits = 2) %>% 
  percent(digits = 0)

# ---- missingsuppfiles ----
# In this section we download and analyse supplementary file names
# R/A02_download_suppfiles.R
load("data/suppfilenames_2017-06-19.RData")

# Let's keep study time frame fixed
suppfilenames <- suppfilenames %>% 
  mutate_at("PDAT", ymd) %>% 
  filter(PDAT <= last_date)

# Timeouts need to be checked!!!
failed_suppfiles <- suppfilenames %>% 
  filter(map_lgl(SuppFileNames, ~inherits(.x, "try-error"))) %>% 
  select(Id, Accession, PDAT, PubMedIds, SuppFileNames) %>%
  unnest(SuppFileNames) %>% 
  mutate(pub = str_length(PubMedIds) != 0,
         SuppFileNames = case_when(
           str_detect(SuppFileNames, "Server denied") ~ "Server denied",
           str_detect(SuppFileNames, "timed out") ~ "Timed out"
         ))

fsupp <- failed_suppfiles %>% 
  filter(SuppFileNames == "Server denied") %>% 
  group_by(PDAT) %>% 
  summarise(N = n(),
            pub = sum(pub)) %>% 
  mutate_at(vars(N, pub), cumsum) %>% 
  gather(key, value, -PDAT) %>% 
  ggplot(aes(ymd(PDAT), value, linetype = key)) +
  geom_line() +
  labs(x = "Publication date", 
       y = "Number of GEO series") +
  scale_linetype_discrete(labels = c("All series", "Series with\npublications")) +
  theme(legend.position = c(0.5, 0.77),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank())

has_suppfile <- suppfilenames %>% 
  mutate(has_suppfile = !map_lgl(SuppFileNames, ~inherits(.x, "try-error"))) %>% 
  select(has_suppfile) %>% 
  table

## in 24 cases we failed to get dir listing because time out
# perc_wsuppfile <- percent(round(1 - (has_suppfile[1] / sum(has_suppfile)), 2), 
#                           digits = 0)
# temporary hack
perc_wsuppfile <- percent(round(1 - (nrow(
  filter(failed_suppfiles, SuppFileNames == "Server denied")) / nrow(suppfilenames)), 2), 
                          digits = 0)

# Queryfig ----------------------------------------------------------------

## ---- queryfig -----

pg <- lapply(list(geop, fsupp), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = 2, widths = c(2, 1))
grid.draw(pga)

# Commonfilenames ---------------------------------------------------------
## ---- commonfilenames -----

# Analyse file extensions
file_ext <- suppfilenames %>% 
  filter(!map_lgl(SuppFileNames, ~inherits(.x, "try-error"))) %>% 
  mutate(filext = map(SuppFileNames, get_filext),
         filext = map(filext, str_replace, "[[:punct:]]", "")) %>% 
  select(Accession, filext) %>% 
  unnest() %>% 
  group_by(filext) %>% 
  summarise(N = n()) %>%
  mutate(perc = (N / sum(N)) * 100) %>% 
  arrange(desc(N))

file_ext_plot <- ggplot(file_ext, aes(reorder(filext, N), N)) +
  geom_point() +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
  coord_flip() + 
  labs(y = "N of files") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(color = "black", size = 6))

# Single most common filename: filelist.txt
most_common_filename <- suppfilenames %>% 
  filter(!map_lgl(SuppFileNames, ~inherits(.x, "try-error"))) %>% 
  unnest(SuppFileNames) %>% 
  group_by(SuppFileNames) %>% 
  summarise(N = n())

# Supplemental file names with more than N=10 occurences
cf <- suppfilenames %>% 
  filter(!map_lgl(SuppFileNames, ~inherits(.x, "try-error"))) %>% 
  unnest(SuppFileNames) %>%
  mutate(common_filenames = str_replace(SuppFileNames, "GSE[0-9]*_", ""),
         common_filenames = str_replace(common_filenames, "\\.gz$", ""), 
         common_filenames = tolower(common_filenames))

cfn <- group_by(cf, common_filenames) %>% 
  summarise(N = n()) %>% 
  arrange(desc(N)) %>% 
  filter(N > 10)

cfp <- ggplot(cfn, aes(reorder(common_filenames, N), N)) +
  geom_point() +
  scale_y_log10(breaks = c(40, 200, 1000, 5000)) +
  coord_flip() + 
  labs(y = "N of files") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(color = "black"))

# plot commonfilenames ggplot
pg <- lapply(list(file_ext_plot, cfp), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = 2, widths = c(2, 3))
grid.draw(pga)

# Raw data as supplemental file plot (decide wether to show)
cfraw <- filter(cf, str_detect(common_filenames, "raw.tar")) %>%
  arrange(PDAT) %>%
  group_by(PDAT) %>%
  summarise(N = n()) %>%
  mutate(N = cumsum(N)) %>%
  ggplot(aes(ymd(PDAT), N, group = 1)) +
  geom_line()

# Percent GEO ids submit raw data as supplemental file
n_raw <- cf %>% 
  filter(str_detect(common_filenames, "raw.tar")) %>% 
  n_distinct()

# Total nr of accessions
n_acc <- suppfilenames %>% 
  select(Accession) %>% 
  n_distinct()
supp_raw_perc <- percent(n_raw / n_acc, 0)

## ---- out-strings ----

out_string1 <- c("filelist","annotation","readme","error","raw.tar","csfasta",
                 "bam","sam","bed","[:punct:]hic","hdf5","bismark","map",
                 "barcode","peaks")
out_string2 <- c("tar","gtf","(big)?bed(\\.txt|12|graph|pk)?","bw","wig",
                 "hic","gct(x)?","tdf","gff(3)?","pdf","png","zip","sif",
                 "narrowpeak","fa", "r$", "rda(ta)?$")

## ---- filesofinterest ----
suppf_total <- suppfilenames %>%
  filter(!map_lgl(SuppFileNames, ~inherits(.x, "try-error"))) %>% 
  unnest(SuppFileNames) %>% 
  nrow()

suppfiles_of_interest <- suppfilenames %>% 
  filter(!map_lgl(SuppFileNames, ~inherits(.x, "try-error"))) %>%
  unnest(SuppFileNames) %>%
  filter(!str_detect(tolower(SuppFileNames), 
                     paste0(out_string1, collapse = "|")),
         !str_detect(tolower(SuppFileNames), 
                     paste0(out_string2, "(\\.gz|\\.bz2)?$", 
                            collapse = "|"))) %>% 
  select(Accession, SuppFileNames, FTPLink, PDAT) %>% 
  mutate(filext = str_extract(tolower(SuppFileNames), "\\.[:alpha:]+([:punct:][bgz2]+)?$")) 

## Number accessions with files of interest
filesofinterest <- suppfiles_of_interest %>% 
  select(Accession) %>% 
  n_distinct()

suppf_oi_perc <- percent(1 - (filesofinterest / n_acc), 0)

# Plot number of suppfiles per accession
suppf_per_acc <- suppfiles_of_interest %>% 
  group_by(Accession) %>% 
  summarise(N = n()) %>% 
  select(N) %>% 
  table() %>% 
  as_data_frame()
colnames(suppf_per_acc) <- c("Files", "Count")

suppf_per_acc %>% 
  mutate_at("Files", as.numeric) %>% 
  ggplot(aes(log2(Files), Count)) +
  geom_bar(stat = "identity") + 
  labs(x = "Files per GEO Accession, log2",
       y = "N of GEO Accessions")

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

#' dataset with p values
p_value_dims <- dims %>% 
  filter(!map_lgl(pvalues, is.null))

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
