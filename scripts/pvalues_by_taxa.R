non_filtered <- select(taxons, Accession, Organism) %>% 
  distinct() %>% 
  left_join(spark_table, .) %>% 
  mutate(group = case_when(
    Organism == "Homo sapiens" ~ "human",
    Organism == "Mus musculus" ~ "mouse",
    TRUE ~ "other"
  )) %>% 
  group_by(group, Type) %>% 
  summarise(n = n())

filtered <- select(taxons, Accession, Organism) %>% 
  distinct() %>% 
  left_join(spark_table_bm, .) %>% 
  mutate(group = case_when(
    Organism == "Homo sapiens" ~ "human",
    Organism == "Mus musculus" ~ "mouse",
    TRUE ~ "other"
  )) %>% 
  group_by(group, Type) %>% 
  summarise(n = n()) 

type_props <- bind_rows(non_filtered, filtered, .id = "id") %>% 
  ungroup() %>% 
  mutate(id = case_when(
    id == 1 ~ "raw",
    TRUE ~ "filtered"),
    id = factor(id, levels = c("raw", "filtered"))) %>% 
  complete(Type, nesting(id, group), fill = list(n = 0)) %>% 
  group_by(id, group) %>% 
  mutate(total = sum(n), 
         `%` = n / total)

type_props %>% 
  ggplot() + 
  geom_bar(aes(x = Type, y = `%`, fill = group), stat = "identity", position = "dodge") +
  facet_wrap(~id) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

write_csv(type_props, "output/type_props.csv")

ggplot(type_props) + 
  geom_bar(aes(x = id, y = `%`, fill = Type), stat = "identity", position = "stack") +
  facet_wrap(~ group, scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete(limits = c("raw", "filtered")) +
  labs(y = "Percent of P value histograms") + 
  theme(axis.title.x = element_blank()) +
  scale_fill_viridis_d()


library(brms)
library(tidybayes)
fit <- brm(n | trials(total) ~ Type + id + group + Type:id + Type:group + group:id, 
           data = type_props, 
           iter = 2000,
           cores = 4,
           family = binomial(),
           control = list(max_treedepth = 12))
summary(fit)
plot(marginal_effects(fit))
fitted_values <- posterior_linpred(fit, transform = TRUE)
head(fitted_values)
apply(fitted_values, 2, mean_hdi, .width = 0.95) %>% 
  bind_rows()
merged_props %>% 
  slice(rep(1:n(), each = 2))
