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
source("_common.R")

# Histogram types summary table
types_legend <- read_csv("../data/pvalue_hist_types.csv")
# remove comments
types_legend <- types_legend %>%
  select(type,typetext)

# READ IN CLASSIFICATIONS --------------------------------------------

# Load manually assigned classes
his_all <- read_delim("../data/pvalue_hist_UM_190121.csv", 
                  delim = ",", 
                  locale = locale(decimal_mark = ","))

# Split manualy assigned classes to filtered and non filtered
his_bm <- his_all %>%
  select(Accession,
         suppdata_id = `Supplementary file name`,
         pi0 = `True nulls proportion filtered`,
         Type = `Type filtered`) %>%
  filter(!is.na(Type))

his<- his_all %>%
  select(Accession,
         suppdata_id = `Supplementary file name`,
         pi0 = `True nulls proportion`,
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

pvalues <- readRDS("../output/pvalue_spark_bins.rds")

# Sort
pvalues <- pvalues %>%
  arrange(suppdata_id) %>%
  mutate(pi0 = as.character(digits(pi0, 6)))

pvalues_bins <- pvalues

#add pvalue frequency bins to columns
for (i in 1:40) {
  new_col_name <-paste0('V', i)
  pvalues_bins <- pvalues_bins %>%
    add_column(!!new_col_name:=map_dbl(pvalues[,5][[1]], ~ .x[i]))
}

# Post basemean corrected data

pvalues_bm <- readRDS("../output/pvalue_bm_spark_bins.rds")

# Sort
pvalues_bm <- pvalues_bm %>%
  arrange(suppdata_id) %>%
  mutate(pi0 = as.character(digits(pi0, 6)))

pvalues_bm_bins <- pvalues_bm

#add pvalue frequency bins to columns
for (i in 1:40) {
  new_col_name <-paste0('V', i)
  pvalues_bm_bins <- pvalues_bm_bins %>%
    add_column(!!new_col_name:=map_dbl(pvalues_bm[,5][[1]], ~ .x[i]))
}

# ADD MANUAL HISTOGRAM CODING ----------------------------------------------------

pvalues_coded <- pvalues_bins %>%
  left_join(his) %>%
  filter(!is.na(Type)) %>%
  select(-Accession,
         -suppdata_id, -pi0,-annot, -values)

pvalues_bm_coded <- pvalues_bm_bins %>%
  left_join(his_bm) %>%
  filter(!is.na(Type)) %>%
  select(-Accession,
         -suppdata_id, -pi0,-annot, -values) 

# type text has to be factor
pvalues_coded$Type <- as.factor(pvalues_coded$Type)
pvalues_bm_coded$Type <- as.factor(pvalues_bm_coded$Type)

# combine pvalue datasets
pvalues_all_coded <- rbind(pvalues_coded,pvalues_bm_coded)

# RANDOM FOREST ---- CREATE CLASSIFICATION MODEL -----------------------------

# Splitting the dataset into train and validation set in the ratio 70:30
train <- sample(nrow(pvalues_all_coded), 0.7*nrow(pvalues_all_coded), replace = FALSE)
TrainSet <- pvalues_all_coded[train,]
ValidSet <- pvalues_all_coded[-train,]
summary(TrainSet)
summary(ValidSet)

 #Using For loop to identify the right mtry for model
   a=c()
   for (i in 3:40) {
     model1 <- randomForest(Type ~ ., data = TrainSet, ntree = 500, mtry = i, importance = TRUE)
     predValid <- predict(model1, ValidSet, type = "class")
     a[i-2] = mean(predValid == ValidSet$Type)
   }
   a
   plot(3:40,a)

# Maximum or near maximum mtry seams to work well in most of the times. Sometimes the maximum is at 10-20 mtry,
# and it seams to depend from the random nature of creation of training set.
# Probable reason seams to be the (difficult and therefore non stable) classification of spiky histograms. 
# As there are not many spiky histograms (and many of them are not correctly assigned) they tend to be distributed
# between training and validation sets differently in different runs. 
# It Will continue with near maximum mtry. This migth allowe more better to replace the wrong classification
# of filtered datasets later.

model2 <- randomForest(Type ~ ., data = TrainSet, mtry = 40, ntree = 2000, importance = TRUE)
model2

# Predicting on train set
predTrain <- predict(model2, TrainSet, type = "class")
mean(predTrain == TrainSet$Type)
# Checking classification accuracy
table(predTrain, TrainSet$Type)

# Predicting on Validation set
predValid <- predict(model2, ValidSet, type = "class")
# Checking classification accuracy
mean(predValid == ValidSet$Type)                    
table(predValid,ValidSet$Type)

# Predicting on full set
predValid <- predict(model2, pvalues_all_coded, type = "class")
# Checking classification accuracy
mean(predValid == pvalues_all_coded$Type)                    
table(predValid,pvalues_all_coded$Type)


# PREDICT AND ADD CLASSIFICATION -----------------------------------------------------

## Get set without manual classification and exclude
## sets with bad pi0
pvalues_bins_unclass <- pvalues_bins %>%
  left_join(distinct(his)) %>%
  filter(is.na(Type),pi0 > pi0threshold)


pvalues_bm_bins_unclass <- pvalues_bm_bins %>%
  left_join(distinct(his_bm)) %>%
  filter(is.na(Type),pi0 > pi0threshold)

# get bins only
pvalues_bins_only <- pvalues_bins_unclass %>%
  select(-Accession,
         -suppdata_id, -annot, -values, -pi0)

pvalues_bm_bins_only <- pvalues_bm_bins_unclass %>%
  select(-Accession,
         -suppdata_id, -annot, -values, -pi0) 

pred_all <- predict(model2, pvalues_bins_only, Type = "class")
pred_bm_all <- predict(model2, pvalues_bm_bins_only, Type = "class")

# add predicted classification
pvalues_bins_unclass <- pvalues_bins_unclass %>%
  mutate(Type = pred_all) %>%
  select(Accession, suppdata_id,pi0,Type)

pvalues_bm_bins_unclass <- pvalues_bm_bins_unclass %>%
  mutate(Type = pred_bm_all) %>%
  select(Accession, suppdata_id,pi0,Type)

# write new reference files
write_csv(pvalues_bins_unclass,"../output/histogram_classification_prefiltration.csv")
write_csv(pvalues_bm_bins_unclass,"../output/histogram_classification_postfiltration.csv")


