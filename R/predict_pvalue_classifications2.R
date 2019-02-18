# *Reads in:
# a)manual histogram classifications.
## b)all efects automatic classifications.
# c)pvalue histograms in 40 bin format.
## *Combines manual classification with automatic classifications
# *Matches histogram classifications and histograms in 40 bin format.
# *Uses random forest maschine learning algorithm to create classification model
# based on histograms with earlier manual classification.
# *Afterwards, created model is applied on all histograms.
# 
# Both basemean prefiltrated and postfiltrated data are used.
#
# Output:
# 1) tsv-formatted files with manual (if present) and predicted classification.
# 2) tsv-formatted files with pvalue sets where manual and predicted
# classification do no match.
# PS! basemean prefiltrated and postfiltrated data are outputted separately.

library(readr)
library(dplyr)
library(randomForest)
library(rsample)
library(recipes)
library(parsnip)
library(yardstick)
source("R/_common.R")

# Histogram types summary table
types_legend <- read_csv(here("data/pvalue_hist_types.csv"))

# remove comments
types_legend <- types_legend %>%
  select(type, typetext)

# READ IN CLASSIFICATIONS --------------------------------------------

# Load manually assigned classes
his_all <- read_delim(here("data/pvalue_hist_UM_190121.csv"), 
                  delim = ",", 
                  locale = locale(decimal_mark = ","))

# Split manualy assigned classes to filtered and non filtered
his_bm <- his_all %>%
  select(Accession,
         suppdata_id = `Supplementary file name`,
         Type = `Type filtered`) %>%
  filter(!is.na(Type))

his <- his_all %>%
  select(Accession,
         suppdata_id = `Supplementary file name`,
         Type) %>%
  filter(!is.na(Type))

# Skip all_effects as it does not seam to help to improve prediction
# and it is possible to identify them easily.
#
# # Load automatic classification of all_effects
# his_auto <- read_delim("data/histogram_classes_all_efects.csv", 
#                   delim = ",", 
#                   locale = locale(decimal_mark = ",")) %>%
#   mutate(Type = "all_effects") %>%
#   select(-type)
# 
# his_auto_bm <- read_delim("data/histogram_classes_all_efects_filtered.csv", 
#                        delim = ",", 
#                        locale = locale(decimal_mark = ",")) %>%
#   mutate(Type = "all_effects") %>%
#   select(-type)
# 
# #combine manually assigned classifications with automatic classifications
# his = full_join(his, his_auto)
# his_bm = full_join(his_bm, his_auto_bm)
# #his_full = full_join(his_1,his_2)
# #his_full

# READ IN DATASET WITH PVALUE BINS ----------------------------------

# Pre basemean corrected data
pvalues <- read_rds(here("output/pvalue_spark_bins.rds"))

# Sort by id
pvalues <- pvalues %>%
  arrange(suppdata_id) %>% 
  mutate_at("pi0", as.numeric)

# Rearrange pvalue frequency bins to columns
pvalues_bins <- pvalues %>% 
  mutate(values = map(values, ~as.tibble(matrix(.x, nrow = 1)))) %>% 
  unnest(values) %>% 
  select(Accession, suppdata_id, annot, everything())

# Post basemean corrected data
pvalues_bm <- read_rds(here("output/pvalue_bm_spark_bins.rds"))

# Sort by id
pvalues_bm <- pvalues_bm %>%
  arrange(suppdata_id) %>% 
  mutate_at("pi0", as.numeric)

# Rearrange pvalue frequency bins to columns
pvalues_bm_bins <- pvalues_bm %>% 
  mutate(values = map(values, ~as.tibble(matrix(.x, nrow = 1)))) %>% 
  unnest(values) %>% 
  select(Accession, suppdata_id, annot, everything())

# ADD MANUAL HISTOGRAM CODING ----------------------------------------------------
pvalues_coded <- pvalues_bins %>%
  left_join(his) %>%
  filter(!is.na(Type))

pvalues_bm_coded <- pvalues_bm_bins %>%
  left_join(his_bm) %>% 
  filter(!is.na(Type))

# Combine pvalue datasets and convert type text to factor
pvalues_all_coded <- bind_rows(pvalues_coded, pvalues_bm_coded) %>% 
  mutate_at("Type", as.factor)

# RANDOM FOREST ---- CREATE CLASSIFICATION MODEL -----------------------------

# Splitting the dataset into train and validation set in the ratio 70:30
set.seed(11)
data_split <- select(pvalues_all_coded, pi0:Type) %>%
  select(Type, everything()) %>% 
  initial_split(p = 0.8)
train <- training(data_split)
test <- testing(data_split)
   
pvalue_train <- recipe(Type ~ ., data = train) %>%
  prep(training = train, retain = TRUE)

# Maximum or near maximum mtry seems to work well in most of the times. Sometimes the maximum is at 10-20 mtry,
# and it seams to depend from the random nature of creation of training set.
# Probable reason seams to be the (difficult and therefore non stable) classification of spiky histograms. 
# As there are not many spiky histograms (and many of them are not correctly assigned) they tend to be distributed
# between training and validation sets differently in different runs. 
# It Will continue with near maximum mtry. This migth allowe more better to replace the wrong classification
# of filtered datasets later.
rf_fit <- 
  rand_forest(mode = "classification", mtry = .preds(), trees = 3000) %>%
  set_engine("ranger") %>% 
  fit(Type ~ ., data = juice(pvalue_train))
test_results <- test %>%
  select(Type) %>%
  as_tibble() %>% 
  mutate(
    `rf class` = predict_class(rf_fit, new_data = test)
  ) %>% 
  bind_cols(predict_classprob(rf_fit, new_data = test))
test_results %>% roc_auc(Type, `anti-conservative`:uniform)
test_results %>% accuracy(Type, `rf class`)
test_results %>% conf_mat(Type, `rf class`)

fit <- 
  boost_tree(mode = "classification", mtry = 5, trees = 100, min_n = 8) %>%
  set_engine("xgboost") %>% 
  fit(Type ~ ., data = juice(pvalue_train))
test_results <- test %>%
  select(Type) %>%
  as_tibble() %>% 
  mutate(
    class = predict_class(fit, new_data = test)
  ) %>% 
  bind_cols(predict_classprob(fit, new_data = test))
test_results %>% roc_auc(Type, `anti-conservative`:uniform)
test_results %>% accuracy(Type, class)
test_results %>% conf_mat(Type, class)

# PREDICT AND ADD CLASSIFICATION -----------------------------------------------------

## Get set without manual classification and exclude
## sets with bad pi0
pvalues_bins_unclass <- pvalues_bins %>%
  left_join(distinct(his)) %>%
  filter(is.na(Type), pi0 > pi0threshold)

pvalues_bm_bins_unclass <- pvalues_bm_bins %>%
  left_join(distinct(his_bm)) %>%
  filter(is.na(Type), pi0 > pi0threshold)

# Get bins only
pvalues_bins_only <- pvalues_bins_unclass %>% 
  select(pi0:Type)

pvalues_bm_bins_only <- pvalues_bm_bins_unclass %>%
  select(pi0:Type)

pred_all <- predict_class(fit, pvalues_bins_only)
pred_bm_all <- predict_class(fit, pvalues_bm_bins_only)

# Add predicted classification
pvalues_bins_unclass <- pvalues_bins_unclass %>%
  mutate(Type = pred_all) %>%
  select(Accession, suppdata_id, pi0, annot, Type)

pvalues_bm_bins_unclass <- pvalues_bm_bins_unclass %>%
  mutate(Type = pred_bm_all) %>%
  select(Accession, suppdata_id, pi0, annot, Type)

# write new reference files
write_csv(pvalues_bins_unclass, here("output/histogram_classification_prefiltration.csv"))
write_csv(pvalues_bm_bins_unclass, here("output/histogram_classification_postfiltration.csv"))
