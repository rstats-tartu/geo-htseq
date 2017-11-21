
knitr::opts_chunk$set(echo = FALSE, 
                      message = FALSE, 
                      warning = FALSE,
                      out.width = "70%",
                      fig.align = 'center',
                      fig.width = 6,
                      fig.asp = 0.5,
                      fig.show = "hold",
                      dev = 'svg')

options(htmltools.dir.version = FALSE, 
        formatR.indent = 2, 
        width = 55, 
        digits = 4, 
        warnPartialMatchAttr = FALSE, 
        warnPartialMatchDollar = FALSE)

library(tidyverse)
library(lubridate)
library(stringr)
library(formattable)
library(Biobase)
library(grid)
library(gridExtra)
library(limma)

# Load helper functions
source("lib/helpers.R")
# Set panel label case
panel_label_case <- "upper"
nrowthreshold <- 4000
