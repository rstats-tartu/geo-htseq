
# Install and load libraries
source("https://bioconductor.org/biocLite.R")
if (!require("pacman")) install.packages("pacman", repos = "https://cloud.r-project.org/")

pacman::p_load(bookdown, tidyverse, lubridate, glue, stringr, formattable, Biobase, grid, 
                 gridExtra, limma, viridis, sparkline, kableExtra, knitr,
                 ape, ggtree, data.table)
pacman::p_load_gh("tpall/entrezquery")

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
last_date <- ymd("2017-06-19")
nrowthreshold <- 4000
pi0threshold <- 0.05

# Folder where downloaded supplementary and matrix files live
local_suppfile_folder <- "output/suppl" 
local_matrixfile_folder <- "output/matrix"
