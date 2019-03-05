
#' Load libraries. We are going to use parsnip api.
source("R/_common.R")
library(rsample)
library(recipes)
library(parsnip)
library(yardstick)

#' READ IN CLASSIFICATIONS --------------------------------------------
#' Load manually assigned classes
his_all <- read_csv("data/pvalue_table_YM_with_set.csv")

#' Split manualy assigned classes to filtered and non filtered
#' 
# his_all <- his_all %>% 
#   select(Accession, suppdata_id = `Supplementary file name`, starts_with("Type")) %>% 
#   gather(key = Filter, value = Type, Type, `Type filtered`) %>%
#   na.omit() %>% 
#   mutate(Filter = case_when(
#     Filter == "Type" ~ "raw",
#     TRUE ~ "basemean"
#   ))
his_all <- his_all %>% 
  mutate(suppdata_id = case_when(
    is.na(set) ~ suppdata_id,
    !is.na(set) ~ str_c(suppdata_id, "-set-", set)
  )) %>% 
  select(Accession, suppdata_id, starts_with("type")) %>% 
  rename(Type = type, `Type filtered` = type_filt) %>% 
  gather(key = Filter, value = Type, Type, `Type filtered`) %>%
  na.omit() %>% 
  mutate(Filter = case_when(
    Filter == "Type" ~ "raw",
    TRUE ~ "basemean"
  ))

# READ IN DATASET WITH PVALUE BINS ----------------------------------

#' Import pvalue data.
pvalues <- read_rds(here("output/pvalues.rds"))
pvalues_bm <- read_rds(here("output/pvalues_bm.rds"))

#' Pool raw and basemean corrected data.
pvalues_pool <- bind_rows(raw = pvalues, basemean = pvalues_bm, .id = "Filter")

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
  mutate_at("Type", as.factor) %>% 
  select(-Accession)

# RANDOM FOREST ---- CREATE CLASSIFICATION MODEL -----------------------------

# Splitting the dataset into train and validation set in the ratio 70:30
set.seed(11)
data_split <- pvalues_coded %>% 
  select(Type, starts_with("V")) %>% 
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
nn_fit <- 
  mlp(mode = "classification", penalty = 0.01, epochs = 3000) %>% 
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
