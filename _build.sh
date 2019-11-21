#!/bin/bash

Rscript -e "bookdown::render_book('index.Rmd', output_format = 'bookdown::html_document2', encoding = 'UTF-8')"
