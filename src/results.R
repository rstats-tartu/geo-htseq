
# For development only, not run during article compilation

## Load libs
source("src/_common.R")

# Load helper functions
source("lib/helpers.R")

## ---- rna-seq-dynamics ----

# Last date to consider geo series and suppfilenames
last_date <- ymd("2017-06-19")

## src/A01_GEO_query.R
load("data/ds.RData") # mouse and human GEO HT-seq expr datasets
## all HT-seq datasets 
ds_all <- filter(ds, ymd(PDAT) <= last_date)
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
  ylab("Number of GEO series") +
  xlab("Publication date") +
  scale_linetype_discrete(labels = c("All series","Series with\npublications")) +
  theme(legend.position = c(0.75, 0.77),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank())

## Calculate percent series using human or mouse 
## formattable::percent()
perc_mmhs <- percent(round(nrow(ds) / nrow(ds_all), 1), digits = 0)
ppub <- group_by(pdat, id, key) %>% 
  summarise_at("value", max) %>% 
  group_by(id) %>% 
  summarise(ppub = value[key=="pub"] / value[key=="N"]) %>% 
  ungroup %>% 
  summarise(ppub = mean(ppub)) %>% 
  .$ppub %>% 
  round(digits = 1) %>% 
  percent(digits = 0)

# ---- missingsuppfiles ----
# In this section we download and analyse supplementary file names
# src/A02_download_suppfiles.R
load("data/suppfilenames_2017-06-19.RData")

# Let's keep study time frame fixed
suppfilenames <- suppfilenames %>% 
  mutate_at("PDAT", ymd) %>% 
  filter(PDAT <= last_date)

failed_suppfiles <- suppfilenames %>% 
  filter(map_lgl(SuppFileNames, ~inherits(.x, "try-error"))) %>% 
  select(PDAT, PubMedIds) %>% 
  mutate(pub = str_length(PubMedIds) != 0) %>% 
  arrange(PDAT) %>% 
  group_by(PDAT) %>% 
  summarise(N = n(),
            pub = sum(pub)) %>% 
  mutate_at(vars(N, pub), cumsum) %>% 
  gather(key, value, -PDAT) 

failed_suppfiles$id <- "Human or mouse series\nlacking supplementary files"

fsupp <- ggplot(failed_suppfiles, aes(ymd(PDAT), value, linetype = key)) +
  geom_line() +
  facet_wrap(~id) +
  xlab("Publication date") +
  ylab("Number of GEO series") +
  scale_linetype_discrete(labels = c("All series", "Series with\npublications")) +
  theme(legend.position = c(0.5, 0.77),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank())

no_suppfile <- suppfilenames %>% 
  mutate(has_suppfile = !map_lgl(SuppFileNames, ~inherits(.x, "try-error"))) %>% 
  select(has_suppfile) %>% 
  table

perc_wsuppfile <- percent(round(1 - (no_suppfile[1] / sum(no_suppfile)), 2), 
                          digits = 0)

# Queryfig ----------------------------------------------------------------

## ---- queryfig -----

pg <- lapply(list(geop, fsupp), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = 2, widths = c(2, 1))
grid.draw(pga)

# Commonfilenames ---------------------------------------------------------
## ---- commonfilenames -----

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

cfp <- ggplot(cfn, aes(common_filenames, N)) +
  geom_point() +
  scale_x_discrete(limits = rev(cfn$common_filenames)) +
  scale_y_log10() +
  coord_flip() + 
  xlab("Common stubs of SuppFileNames\n(>10 occurences) ") +
  ylab("Number of files")

# plot commonfilenames ggplot
cfp

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
## total nr of accessions
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

suppfiles_of_interest <- suppfilenames %>% 
  unnest(SuppFileNames) %>%
  filter(!str_detect(tolower(SuppFileNames), 
                     paste0(out_string1, collapse = "|")),
         !str_detect(tolower(SuppFileNames), 
                     paste0(out_string2, "(\\.gz|\\.bz2)?$", 
                            collapse = "|"))) %>% 
  select(Accession, SuppFileNames, FTPLink, PDAT) %>% 
  mutate(filext = str_extract(tolower(SuppFileNames), "\\.[:alpha:]+([:punct:][bgz2]+)?$")) 

## number accessions with files of interest
filesofinterest <- suppfiles_of_interest %>% 
  select(Accession) %>% 
  n_distinct()

suppf_oi_perc <- percent(1 - (filesofinterest / n_acc), 0)

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

## Remove errored matixes
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
  labs(x = bquote(Number~of~samples~(log[10])),
       y = bquote(Number~of~features~(log[10]))) +
  geom_hline(yintercept = log10(nrowthreshold), linetype = 2) +
  scale_fill_continuous(name = "Count")

dims_featuresp <- dims %>%
  ggplot(aes(log10(features))) + 
  geom_histogram(bins = 60) +
  labs(x = bquote(Number~of~features~(log[10])),
       y = "Count") +
  geom_vline(xintercept = log10(nrowthreshold), linetype = 2)

pg <- lapply(list(dims_tabp, dims_featuresp), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = 2, widths = c(3, 2))
grid.draw(pga)
