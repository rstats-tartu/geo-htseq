
# For development only, not run during article compilation

# Load libs
source("R/_common.R")

## ---- rna-seq-dynamics ----

# R/A01_GEO_query.R
load("data/ds.RData") # mouse and human GEO HT-seq expr datasets

# all HT-seq datasets 
ds_all <- filter(ds, ymd(PDAT) <= last_date)
first_date <- min(ymd(ds$PDAT))

# Human or mouse datasets
ds <- filter(ds_all, str_detect(taxon, "Mus musculus|Homo sapiens"))

# Merge all datasets for plotting
ds_merge <- bind_rows(ds_all, ds, .id = "id") %>% 
  mutate_at("PDAT", ymd)

# Count series with publications
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

# Plot submissions
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

# Calculate percent series using human or mouse 
# formattable::percent()
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

## ---- queryfig -----
geop

# pg <- lapply(list(geop, fsupp), ggplotGrob)
# pg <- add_labels(pg, case = panel_label_case)
# pga <- arrangeGrob(grobs = pg, ncol = 2, widths = c(2, 1))
# grid.draw(pga)


