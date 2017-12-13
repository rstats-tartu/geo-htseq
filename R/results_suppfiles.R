
# Load libs
source("R/_common.R")

# ---- missingsuppfiles ----
# In this section we download and analyse supplementary file names
# R/A02_download_suppfiles.R
suppfilenames_imported <- readRDS("data/suppfilenames_2017-11-27.rds")

# Let's keep study time frame fixed
# Extract supplemetary file names
suppfilenames <- suppfilenames_imported %>% 
  mutate_at("PDAT", ymd) %>% 
  filter(PDAT <= last_date) %>% 
  mutate(SuppFileNames = map(suppfiles, "suppl"))

# Datasets with suppfiles
suppfilenames_present <- suppfilenames %>% 
  filter(map_lgl(SuppFileNames, ~length(.x) > 0))

# unnested suppfiles
suppfilenames_present_unnested <- suppfilenames_present %>% 
  unnest(SuppFileNames)

# datasets missing public suppfiles
suppfilenames_not_present <- suppfilenames %>% 
  filter(map_lgl(SuppFileNames, ~length(.x) == 0))

# Timeouts need to be checked!!!
failed_suppfiles <- suppfilenames_not_present %>% 
  mutate(pub = str_length(PubMedIds) != 0)

fsupp <- failed_suppfiles %>% 
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

# in 24 cases we failed to get dir listing because time out
# perc_wsuppfile <- percent(round(1 - (has_suppfile[1] / sum(has_suppfile)), 2), 
#                           digits = 0)
# temporary hack
perc_wsuppfile <- percent(round(1 - (nrow(suppfilenames_not_present) / nrow(suppfilenames)), 2), digits = 0)

ppub_n_ci_nopublicsupps <- suppfilenames_not_present %>% 
  mutate(pub = str_length(PubMedIds) != 0) %>% 
  group_by(PDAT) %>% 
  summarise(geoseries = n(), 
            pub = sum(pub)) %>% 
  mutate_at(vars(geoseries, pub), cumsum) %>% 
  gather(key, value, -PDAT) %>% 
  group_by(key) %>% 
  summarise_at("value", max) %>%
  spread(key, value) %>% 
  mutate(pois = map2(pub, geoseries, binom.test),
         estimate = map_dbl(pois, "estimate"),
         ci = map(pois, "conf.int"),
         ci = map(ci, percent, 1),
         ci = map_chr(ci, ~glue("95%CI, {.x[1]} to {.x[2]}")))

# ---- commonfilenames -----

# Analyse file extensions
file_ext <- suppfilenames_present %>%
  mutate(filext = map(SuppFileNames, get_filext),
         filext = map(filext, str_replace, "[[:punct:]]", "")) %>% 
  select(Accession, filext) %>% 
  unnest() %>% 
  group_by(filext) %>% 
  summarise(N = n()) %>%
  mutate(perc = (N / sum(N)) * 100) %>% 
  arrange(desc(N))

# Plot file extensions
file_ext_plot <- ggplot(file_ext, aes(reorder(filext, N), N)) +
  geom_point() +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
  coord_flip() + 
  labs(y = "N of files") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(color = "black", size = 6))

# Single most common filename: filelist.txt
most_common_filename <- suppfilenames_present_unnested %>% 
  group_by(SuppFileNames) %>% 
  summarise(N = n())

# filelist.txt and raw.tar only
suppfiles_rawtar_only <- suppfilenames_present %>% 
  filter(map_lgl(SuppFileNames, ~ all(str_detect(.x, "filelist.txt|RAW.tar"))))
supp_rawtaronly_perc <- percent(nrow(suppfiles_rawtar_only) / nrow(suppfilenames), 0)

# for diagram: files after removing sets with RAW.tar and filelist.txt only
suppfiles_notonly_rawtar_unnested <- suppfilenames_present %>% 
  select(Id, Accession, SuppFileNames) %>% 
  unnest(SuppFileNames) %>% 
  filter(!map_lgl(SuppFileNames, ~ str_detect(.x, "filelist.txt|RAW.tar")))

# Publications of series with RAW.tar and filelist supplementary files
fsupp_rawtar_only <- fsupp %+% (suppfilenames_present %>% 
  filter(map_lgl(SuppFileNames, ~ all(str_detect(.x, "filelist.txt|RAW.tar")))) %>% 
  mutate(pub = str_length(PubMedIds) != 0) %>% 
  group_by(PDAT) %>% 
  summarise(N = n(),
            pub = sum(pub)) %>% 
  mutate_at(vars(N, pub), cumsum) %>% 
  gather(key, value, -PDAT))

# Supplemental file names with more than N = 10 occurences
cf <- suppfilenames_present_unnested %>%
  mutate(common_filenames = str_replace(SuppFileNames, "GSE[0-9]*_", ""),
         common_filenames = str_replace(common_filenames, "\\.gz$", ""), 
         common_filenames = tolower(common_filenames))

cfn <- group_by(cf, common_filenames) %>% 
  summarise(N = n()) %>% 
  arrange(desc(N)) %>% 
  filter(N > 10)

# Plot common file names
cfp <- ggplot(cfn, aes(reorder(common_filenames, N), N)) +
  geom_point() +
  scale_y_log10(breaks = c(40, 200, 1000, 5000)) +
  coord_flip() + 
  labs(y = "N of files") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(color = "black"))

# Plot commonfilenames ggplot
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
supp_raw_perc <- percent(n_raw / nrow(suppfilenames), 0)

# ---- out-strings ----

out_string1 <- c("filelist","annotation","readme","error","raw.tar","csfasta",
                 "bam","sam","bed","[:punct:]hic","hdf5","bismark","map",
                 "barcode","peaks")
out_string2 <- c("tar","gtf","(big)?bed(\\.txt|12|graph|pk)?","bw","wig",
                 "hic","gct(x)?","tdf","gff(3)?","pdf","png","zip","sif",
                 "narrowpeak","fa", "r$", "rda(ta)?$")

# ---- filesofinterest ----
suppfiles_of_interest <- suppfilenames_present_unnested %>%
  filter(!str_detect(tolower(SuppFileNames), str_c(out_string1, collapse = "|")),
         !str_detect(tolower(SuppFileNames), str_c(out_string2, "(\\.gz|\\.bz2)?$", 
                                                   collapse = "|"))) %>% 
  select(Accession, SuppFileNames, FTPLink, PDAT) %>% 
  mutate(filext = str_extract(tolower(SuppFileNames), "\\.[:alpha:]+([:punct:][bgz2]+)?$")) 

# Number accessions with files of interest
filesofinterest <- suppfiles_of_interest %>% 
  select(Accession) %>% 
  n_distinct()

# Percent of Acc with suppfiles of interest
suppf_oi_perc <- percent(1 - (filesofinterest / nrow(suppfilenames_present)), 0)

# Plot number of suppfiles per accession
# Number of files per accession, all files 
suppfall_per_acc <- suppfilenames_present_unnested %>% 
  group_by(Accession) %>% 
  summarise(N = n()) %>% 
  select(N) %>% 
  table() %>% 
  as_data_frame()
colnames(suppfall_per_acc) <- c("Files", "Count")

# Number of files per accession, files of interest
suppfoi_per_acc <- suppfiles_of_interest %>% 
  group_by(Accession) %>% 
  summarise(N = n()) %>% 
  select(N) %>% 
  table() %>% 
  as_data_frame()
colnames(suppfoi_per_acc) <- c("Files", "Count")

# Percent files after filter
files_oi_perc <- percent(1 - (nrow(suppfiles_of_interest) / nrow(suppfilenames_present_unnested)), digits = 0)

bind_rows(suppfall_per_acc, suppfoi_per_acc, .id = "id") %>% 
  mutate(Files = as.numeric(Files),
         id = case_when(
           id == 1 ~ "All files",
           id == 2 ~ "Files after filter"
         ),
         Filebins = if_else(Files > 10, ">10", as.character(Files))) %>% 
  ggplot(aes(reorder(Filebins, Files), Count)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~ id) +
  annotate("text", x = "6", y = 4500, 
           label = c(glue("N = {nrow(suppfilenames_present_unnested)} "), 
                     glue("N = {nrow(suppfiles_of_interest)}"))) +
  labs(x = "Files per GEO Accession",
       y = "N of GEO Accessions")
