
#' Load libraries. We are going to use parsnip api.
source("scripts/_common.R")
library(rsample)
library(recipes)
library(parsnip)
library(yardstick)

#' READ IN CLASSIFICATIONS --------------------------------------------
#' Load manually assigned classes'

his_all <- read_delim(here("data/pvalue_sets_classes.csv"),
                                  delim =  ";",
                                  col_types = "cccdcd",
                                  escape_double = FALSE, 
                                  trim_ws = TRUE)

his_18 <- read_delim("data/corrected histograms_taavile.csv", 
                     delim =  ";",
                     col_types = "cccdcd",
                     escape_double = FALSE, 
                     trim_ws = TRUE)

his_all <- bind_rows(his_all, his_18) %>% 
  select(-starts_with("pi0")) %>% 
  distinct()

#' Split manualy assigned classes to filtered and non filtered
his_all <- his_all %>% 
  gather(key = Filter, value = Type, type, type_filtered) %>%
  na.omit() %>% 
  mutate(Filter = case_when(
    Filter == "type" ~ "raw",
    TRUE ~ "basemean"
  )) %>% 
  mutate_at("Type", str_replace, "[:punct:]$", "") # Ülo introduced "anti-conservative'" class

# READ IN DATASET WITH PVALUE BINS ----------------------------------

#' Import pvalue data.
p_values <- read_rds(here("output/pvalues.rds"))
p_values_bm <- read_rds(here("output/pvalues_bm.rds"))

#' Pool raw and basemean corrected data.
pvalues_pool <- bind_rows(raw = p_values, basemean = p_values_bm, .id = "Filter")

#' Number of features to sample from
features <- pvalues_pool$features

#' Prepare variables for classification. Keep only sets with at least 
#' nrowthreshold pvalues and non-truncated sets identified by any values over 0.9
pvalues_pool_ecdf <- pvalues_pool %>% 
  select(Filter, suppdata_id, pvalues) %>% 
  mutate(eCDF = map(pvalues, ecdf), # Error in .f(.x[[i]], ...) : 'x' must have 1 or more non-missing values
         values = map(eCDF, function(Fn) Fn(seq(0, 1, 3 / nrowthreshold))))

#' Rearrange pvalue probs to columns.
pvalues_bins <- pvalues_pool_ecdf %>% 
  mutate(values = map(values, matrix, nrow = 1),
         values = map(values, as.tibble)) %>% 
  unnest(values) %>% 
  select(Filter, suppdata_id, starts_with("V"))

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

#' Create anti-conservative sets with fdr adjustment
con <- map(c(uni, aco), p.adjust, method = "fdr")
bim <- map2(aco, con[101:110], ~ sample(c(.x, .y, sample(features, 1))))
synt <- bind_rows(tibble(Type = "uniform", pvalues = uni),
                  tibble(Type = "anti-conservative", pvalues = aco),
                  tibble(Type = "conservative", pvalues = con),
                  tibble(Type = "bimodal", pvalues = bim))

synt <- synt %>% 
  filter(map_lgl(pvalues, ~ length(.x) >= nrowthreshold)) %>% 
  mutate(eCDF = map(pvalues, ecdf),
         values = map(eCDF, function(Fn) Fn(seq(0, 1, 3/nrowthreshold))))

synt <- synt %>% 
  mutate(values = map(values, matrix, nrow = 1),
         values = map(values, as.tibble)) %>% 
  unnest(values) %>% 
  select(starts_with("V"), Type)

# nnet ---- CREATE CLASSIFICATION MODEL -----------------------------

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

#' Get sets without manual classification.
# pvalues_coded <- select(pvalues_coded, Filter, suppdata_id, starts_with("V"))
# pvalues_bins_unclass <- pvalues_bins %>% 
#   select(Filter, suppdata_id, starts_with("V")) %>% 
#   anti_join(pvalues_coded, by = c("Filter", "suppdata_id"))

#' Get bins only.
pvalues_bins_only <- pvalues_bins %>% 
  select(starts_with("V"))

#' Predict classes for all unclassified sets sets.
pred_all <- predict_class(nn_fit, pvalues_bins_only)

#' Add predicted classification
pvalues_bins_classes <- pvalues_bins %>%
  mutate(Type = as.character(pred_all)) %>%
  select(Filter, suppdata_id, Type)

#' Merge model and human classified data and write result to a file. 
bind_rows(human = select(his_all, -Accession),
          nnet = pvalues_bins_classes, .id = "Method") %>% 
  arrange(suppdata_id) %>% 
  spread(key = Method, value = Type) %>% 
  mutate(Type = case_when(
    is.na(human) ~ nnet,
    TRUE ~ human
  ),
  classifier = case_when(
    is.na(human) ~ "nnet",
    TRUE ~ "human"
  )) %>% 
  write_csv(here("output/pvalue_histogram_nnet_classification.csv"))

