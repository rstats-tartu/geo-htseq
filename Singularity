BootStrap: shub
From: tpall/singularity-tidyverse

%labels
  Maintainer tpall

%help
  This container will run geo-rnaseq workflow R scripts. Includes Rstudio server.

%post
  # Install dependencies
  apt-get update
  apt-get install -y \
    unzip \
    xorg \
    libx11-dev \
    libglu1-mesa-dev \
    libfreetype6-dev
  
  # Install CRAN packages
  Rscript -e "install.packages(c('bookdown','XML','formattable','gridExtra','gridBase','viridis','knitr','ape','data.table','kableExtra','sparkline','evaluate','hexbin','broom','digest'), repos = 'https://cloud.r-project.org', dependencies = TRUE)"
  
  # Install Bioconductor packages
  Rscript -e "BiocManager::install(update = FALSE, ask = FALSE)"
  Rscript -e "BiocManager::install(c('GEOquery', 'Biobase', 'limma', 'ggtree'), update = FALSE, ask = FALSE)"

  # Install Github packages
  Rscript -e "devtools::install_github('tpall/SRP')"
  Rscript -e "devtools::install_github('tpall/entrezquery')"
  
  # Install gridDiagram
  Rscript -e  "install.packages('https://www.stat.auckland.ac.nz/~paul/R/Diagram/gridDiagram_0.2-1.tar.gz', repos = NULL, type = 'source')" 

  # Clean up
  rm -rf /var/lib/apt/lists/*
