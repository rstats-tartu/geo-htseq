
library(tidyverse)
library(lubridate)
library(glue)
library(formattable)
library(Biobase)
library(grid)
library(gridExtra)
library(limma)
library(viridis)
library(sparkline)
library(kableExtra)
library(knitr)
library(ape)
library(ggtree)

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

# Load helper functions
source("R/helpers.R")

# Set panel label case
panel_label_case <- "upper"
nrowthreshold <- 4000
pi0threshold <- 0.05

# Last date to consider geo series and suppfilenames
last_date <- ymd("2017-06-19")

# folder where downloaded supplementary and matrix files live
local_suppfile_folder <- "output/suppl" 
local_matrixfile_folder <- "output/matrix"
