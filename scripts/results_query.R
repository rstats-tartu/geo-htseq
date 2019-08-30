
# For development only, not run during article compilation

# Load libs
source("scripts/_common.R")

## ---- rna-seq-dynamics ----
ds_redline <- read_csv("output/ds_redline.csv",
                       col_types = "cccccccccccccccccccccccccc")

first_date <- min(ds_redline$PDAT)

# Count series with publications
pdat <- ds_redline %>% 
  mutate(pub = !is.na(PubMedIds)) %>% 
  group_by(model, PDAT) %>% 
  summarise(geoseries = n(), 
            pub = sum(pub)) %>% 
  mutate_at(vars(geoseries, pub), cumsum) %>% 
  gather(key, value, -PDAT, -model) %>% 
  mutate_at("PDAT", as_date)

# Plot submissions
geop <- pdat %>% 
  ggplot(aes(PDAT, value)) + 
  geom_line(aes(group = key, linetype = key)) +
  facet_wrap(~model) +
  labs(x = "Publication date", 
       y = "Number of GEO series") +
  scale_linetype_discrete(labels = c("All series","Series with\npublications")) +
  theme(legend.position = c(0.75, 0.77),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank())

# Yearly submissions
geop_year <- ds_redline %>% 
  mutate(pub = !is.na(PubMedIds)) %>% 
  group_by(model, year = lubridate::year(PDAT)) %>% 
  summarise(geoseries = n(), 
            pub = sum(pub)) %>% 
  mutate(ppub = pub / geoseries) %>% 
  ggplot(aes(x = year, y = geoseries)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(breaks = seq(2006, 2018, by = 2)) +
  facet_wrap(~model) +
  labs(x = "Publication year", 
       y = "Number of GEO series") +
  theme(legend.position = c(0.75, 0.77),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank())

# Calculate percent series using human or mouse 
# formattable::percent()
perc_mmhs <- percent(table(ds_redline$model)[1] / sum(table(ds_redline$model)), digits = 0)

# Number of publications
ppub_n <- pdat %>%
  group_by(model, key) %>% 
  summarise_at("value", max)

# Percent with publications
ppub <- ppub_n %>% 
  group_by(model) %>% 
  summarise(ppub = value[key=="pub"] / value[key=="geoseries"]) %>% 
  ungroup %>% 
  # filter(str_detect(model, "Human")) %>% 
  .$ppub %>% 
  round(digits = 2) %>% 
  percent(digits = 0)

ppub_n_ci <- ppub_n %>% 
  spread(key, value) %>% 
  mutate(test = map2(pub, geoseries, binom.test),
         ci = map(test, "conf.int"),
         ci = map(ci, percent, 1),
         ci = map_chr(ci, ~ glue("95%CI, {.x[1]} to {.x[2]}")))

pg <- lapply(list(geop_year, geop), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
query_fig <- arrangeGrob(grobs = pg, nrow = length(pg))

## ---- queryfig -----
grid.draw(query_fig)
