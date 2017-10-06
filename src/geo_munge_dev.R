## Return s3 object --------------------------------------------------------

## Copy test files to HD
if(FALSE){
  local_suppfile_folder <- "/Volumes/Media/srp_example/data/counts"
  geos <- c("GSE78760","GSE68155","GSE88970","GSE17403","GSE74993","GSE76261",
            "GSE40918","GSE67517")
  test <- list.files(local_suppfile_folder, full.names = T)[str_detect(list.files(local_suppfile_folder), paste(geos, collapse = "|"))]
  paste(test, collapse = " ")
  system(sprintf("cp %s %s", paste(test, collapse = " "), "~/Dropbox/geo-htseq-article/data/test/"))
}

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
source("lib/helpers.R")
nrowthreshold <- 4000
source("lib/munge_geo.R")
source("lib/checkFullRank.R")
source("lib/text_funs.R")
library(Biobase)
load("data/gsem.RData")

local_suppfile_folder <- "/Volumes/Media/srp_example/data/counts" # "data/counts"
supptabs <- geofile_df(local_suppfile_folder, "suppfiles")

# Merge gsem ExpressionSets to supptabs
supptabs <- select(gsem, Accession, gsematrix) %>% 
  nest(gsematrix, .key = "matrixfiles") %>% 
  left_join(supptabs, .)

start <- Sys.time()
st <- mutate(supptabs, result = map(suppfiles, ~ try(munge_geo2(.x, dir = local_suppfile_folder))))
end <- Sys.time()
start-end

st_unnested <- select(st, Accession, suppfiles, result) %>% 
  unnest %>%
  unnest(sheets)

st_pvalues <- filter(st_unnested, !map_lgl(pvalues, is_null), !map_lgl(pvalues, inherits, "try-error")) 
filter(st_pvalues, features>1) %>% 
  ggplot(aes(log10(features))) + 
  geom_histogram(bins = 40) +
  geom_vline(xintercept = log10(nrowthreshold), linetype = 2)

filter(st_pvalues, features>nrowthreshold) %>%
  summarise_at("features", c("mean","median","min","max"))

filter(st_pvalues, features>nrowthreshold) %>% 
  sample_n(10)


