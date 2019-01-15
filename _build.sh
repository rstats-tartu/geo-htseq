#!/bin/sh

Rscript -e "bookdown::render_book('index.Rmd', output_format = 'bookdown::html_document2')"
