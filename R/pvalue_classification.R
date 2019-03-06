
#' Load libraries. We are going to use parsnip api.
source("R/_common.R")
library(rsample)
library(recipes)
library(parsnip)
library(yardstick)

#' READ IN CLASSIFICATIONS --------------------------------------------
#' Load manually assigned classes
his_all <- read_csv2(here("data/pvalue_sets_classes.csv"))

#' Split manualy assigned classes to filtered and non filtered
his_all <- his_all %>% 
  gather(key = Filter, value = Type, type, type_filtered) %>%
  select(-starts_with("pi0")) %>% 
  na.omit() %>% 
  mutate(Filter = case_when(
    Filter == "type" ~ "raw",
    TRUE ~ "basemean"
  )) %>% 
  mutate_at("Type", str_replace, "[:punct:]$", "") # Ülo introduced "anti-conservative'" class

# READ IN DATASET WITH PVALUE BINS ----------------------------------

#' Import pvalue data.
pvalues <- read_rds(here("output/pvalues.rds"))
pvalues_bm <- read_rds(here("output/pvalues_bm.rds"))

#' Pool raw and basemean corrected data.
pvalues_pool <- bind_rows(raw = pvalues, basemean = pvalues_bm, .id = "Filter")

#' Number of features to sample from
features <- pvalues_pool$features

#' Prepare variables for classification. Keep only sets with at least 
#' nrowthreshold pvalues.
pvalues_pool <- pvalues_pool %>% 
  select(Filter, suppdata_id, pi0, pvalues) %>% 
  filter(map_lgl(pvalues, ~ length(.x) >= nrowthreshold)) %>% 
  mutate(eCDF = map(pvalues, ecdf),
         values = map(eCDF, function(Fn) Fn(seq(0, 1, 3/nrowthreshold))))

#' Rearrange pvalue probs to columns.
pvalues_bins <- pvalues_pool %>% 
  mutate(values = map(values, matrix, nrow = 1),
         values = map(values, as.tibble)) %>% 
  unnest(values) %>% 
  select(Filter, suppdata_id, pi0, starts_with("V"))

#' ADD MANUAL HISTOGRAM CODING ----------------------------------------------------
pvalues_coded <- pvalues_bins %>%
  left_join(his_all) %>%
  filter(!is.na(Type)) %>% 
  select(-Accession)

#' Generate synthetic p-value distributions
set.seed(11)
uni <- replicate(100, runif(sample(features, 1)), simplify = FALSE)
aco <- replicate(10, 
                 replicate(sample(features, 1), t.test(rnorm(3, sample(c(0, 0.5, 1, 1.5), 1, prob = c(0.8, rep(0.2 / 3, 3))), 0.3))$p.value),
                 simplify = TRUE)

#' create anti-conservative sets with fdr adjustment
con <- map(c(uni, aco), p.adjust, method = "fdr")
bim <- map2(aco, con[101:110], ~ sample(c(.x, .y, sample(features, 1))))
synt <- bind_rows(tibble(Type = "uniform", pvalues = uni),
                  tibble(Type = "anti-conservative", pvalues = aco),
                  tibble(Type = "conservative", pvalues = con),
                  tibble(Type = "bimodal", pvalues = bim),
                  
)

synt <- synt %>% 
  filter(map_lgl(pvalues, ~ length(.x) >= nrowthreshold)) %>% 
  mutate(eCDF = map(pvalues, ecdf),
         values = map(eCDF, function(Fn) Fn(seq(0, 1, 3/nrowthreshold))))

synt <- synt %>% 
  mutate(values = map(values, matrix, nrow = 1),
         values = map(values, as.tibble)) %>% 
  unnest(values) %>% 
  select(starts_with("V"), Type)

# RANDOM FOREST ---- CREATE CLASSIFICATION MODEL -----------------------------

# Splitting the dataset into train and validation set in the ratio 70:30
set.seed(11)
data_split <- pvalues_coded %>% 
  select(Type, starts_with("V")) %>% 
  bind_rows(synt) %>% 
  mutate_at("Type", as.factor) %>% 
  initial_split(p = 0.8)
train <- training(data_split)
test <- testing(data_split)
   
pvalue_train <- recipe(Type ~ ., data = train) %>%
  prep(training = train, retain = TRUE)

#' Fit nnet model
nn_fit <- 
  mlp(mode = "classification", penalty = 0.3, epochs = 4000) %>% 
  set_engine("nnet", MaxNWts = 2000) %>% 
  fit(Type ~ ., data = juice(pvalue_train))
test_results <- test %>%
  select(Type) %>%
  as_tibble() %>% 
  mutate(
    class = predict_class(nn_fit, new_data = test)
  ) %>% 
  bind_cols(predict_classprob(nn_fit, new_data = test))
test_results %>% accuracy(Type, class)
test_results %>% conf_mat(Type, class)

# PREDICT AND ADD CLASSIFICATION -----------------------------------------------------

#' Get set without manual classification and exclude
#' sets with bad pi0.
pvalues_bins_unclass <- pvalues_bins %>% 
  select(Filter, suppdata_id, pi0, starts_with("V")) %>% 
  anti_join(select(pvalues_coded, Filter, suppdata_id, pi0, starts_with("V")), 
            by = c("Filter", "suppdata_id")) %>% 
  filter(pi0 > pi0threshold)

#' Get bins only.
pvalues_bins_only <- pvalues_bins_unclass %>% 
  select(starts_with("V"))

#' Predict classes for all unclassified sets sets.
pred_all <- predict_class(nn_fit, pvalues_bins_only)

#' Add predicted classification
pvalues_bins_classes <- pvalues_bins_unclass %>%
  mutate(Type = pred_all) %>%
  select(Filter, suppdata_id, Type)

#' Merge model and human classified data and write result to file. 
bind_rows(human = select(pvalues_coded, Filter, suppdata_id, Type),
          nnet = pvalues_bins_classes, .id = "Method") %>% 
  arrange(suppdata_id) %>% 
  write_csv(here("output/pvalue_histogram_nnet_classification.csv"))
