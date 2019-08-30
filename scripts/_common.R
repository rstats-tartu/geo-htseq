
# Load installed libraries
library(ggplot2)
library(purrr)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(readxl)
library(broom)
library(lubridate)
library(glue)
library(stringr)
library(limma)
library(ape)
library(ggtree)
library(data.table)
library(bookdown)
library(viridis)
library(kableExtra)
library(knitr)
library(formattable)
library(sparkline)
library(grid)
library(gridExtra)
library(Biobase)
library(hexbin)
library(here)

#' default pacman and devtools install fails in conda, this should work 
library(devtools)
options(unzip = "internal")

if (!require(SRP)) {
        remotes::install_github("tpall/SRP")
}

if (!require(entrezquery)) {
        remotes::install_github("tpall/entrezquery")
}

library(entrezquery)

# Plot options
opts_chunk$set(echo = FALSE, 
               message = FALSE, 
               warning = FALSE,
               fig.align = 'center',
               fig.show = "hold",
               dev = "svg")

options(htmltools.dir.version = FALSE, 
        formatR.indent = 2, 
        width = 55, 
        digits = 4, 
        warnPartialMatchAttr = FALSE, 
        warnPartialMatchDollar = FALSE)

# Set panel label case
panel_label_case <- "upper"

# Load helper functions
source("scripts/helpers.R")

# Last date to consider geo series and suppfilenames
last_date <- ymd("2018-12-31")
nrowthreshold <- 1000
pi0threshold <- 0.05

# Folder where downloaded supplementary and matrix files live
local_suppfile_folder <- "output/suppl" 
local_matrixfile_folder <- "output/matrix"
