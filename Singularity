BootStrap: shub
From: tpall/singularity-r:3.5.2

%labels
  Maintainer tpall

%help
  This will run geo-rnaseq workflow.

%post
  # Install dependencies
  apt-get update
  apt-get install -y unzip
  
  # Install CRAN packages
  Rscript -e "install.packages(c('bookdown','XML','devtools','ggplot2','purrr','tibble','dplyr','tidyr','stringr','readr','lubridate','glue','formattable','gridExtra','gridBase','viridis','knitr','ape','data.table','kableExtra','sparkline','evaluate','hexbin','broom','readxl','digest','tidyverse','BiocManager'), repos = 'https://cloud.r-project.org', dependencies = TRUE)"
  
  # Install Bioconductor packages
  Rscript -e "BiocManager::install(update=FALSE,ask=FALSE)"
  Rscript -e "BiocManager::install(c('GEOquery','Biobase','limma','ggtree'),update=FALSE,ask=FALSE)"

  # Clean up
  rm -rf /var/lib/apt/lists/*
