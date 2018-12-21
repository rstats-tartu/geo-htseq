
# Load installed libraries
library(ggplot2)
library(purrr)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
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

#' default pacman and devtools install fails in conda, this should work 
library(devtools)
options(unzip = "internal")
devtools::install_github("tpall/SRP")
devtools::install_github("tpall/entrezquery")
library(entrezquery)

# Plot options
opts_chunk$set(echo = FALSE, 
               message = FALSE, 
               warning = FALSE,
               out.width = "80%",
               fig.align = 'center',
               fig.width = 6,
               fig.asp = 0.5,
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
source("R/helpers.R")

# Last date to consider geo series and suppfilenames
last_date <- ymd("2017-12-31")
nrowthreshold <- 1000
pi0threshold <- 0.05

# Folder where downloaded supplementary and matrix files live
local_suppfile_folder <- "output/suppl" 
local_matrixfile_folder <- "output/matrix"
