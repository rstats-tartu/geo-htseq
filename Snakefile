
singularity: "shub://tpall/geo-rnaseq"
last_date = "2018-12-31"

rule all:
  input: "output/gsem.rds", "output/suppdata.rds", "_main.html"

localrules: all, geo_query

rule geo_query:
  output: 
    "output/document_summaries.rds"
  script:
    "R/01_geo_query.R"

rule download_suppfilenames:
  input: 
    rules.geo_query.output
  output: 
    "output/suppfilenames.rds"
  params: 
    last_date = last_date
  script:
    "R/02_download_suppfilenames.R"

rule filter_suppfilenames:
  input: 
    rules.download_suppfilenames.output
  output: 
    "output/suppfilenames_filtered.rds"
  script:
    "R/03_filter_suppfilenames.R"

rule download_suppfiles:
  input: 
    rules.filter_suppfilenames.output
  output: 
    touch("output/downloading_suppfiles.done")
  script:
    "R/04_download_suppfiles.R"

rule series_matrixfiles:
  input: 
    rules.download_suppfiles.output
  output: 
    "output/gsem.rds"
  script:
    "R/05_series_matrixfiles.R"

n_files = list(range(1, 101, 1))
rule split_suppfiles:
  input: 
    docsums = rules.geo_query.output, 
    suppfilenames_filtered = rules.filter_suppfilenames.output, 
    gsem = rules.series_matrixfiles.output
  output: 
    temp(expand("output/tmp/supptabs_{n}.rds", n = n_files))
  script:
    "R/06_split_suppfiles.R"

rule import_suppfiles:
  input: 
    "output/tmp/supptabs_{n}.rds"
  output: 
    temp("output/tmp/suppdata_{n}.rds")
  script:
    "R/06_import_suppfiles.R"
  
rule merge_suppdata:
  input: 
    expand("output/tmp/suppdata_{n}.rds", n = n_files)
  output: "output/suppdata.rds"
  script:
    "R/07_merge_suppdata.R"

rule download_publications:
  input: rules.geo_query.output
  output: "output/publications.rds"
  params: 
    last_date = last_date
  script:
    "R/07_download_publications.R"

rule report:
  input: "index.Rmd", "01_introduction.Rmd", "02_methods.Rmd", "03_results.Rmd", "04_discussion.Rmd", "05_references.Rmd", "output/document_summaries.rds", "output/suppfilenames.rds", "output/suppfilenames_filtered.rds", "output/gsem.rds", "output/suppdata.rds", "output/publications.rds"
  output: "_main.html"
  shell:
    """
    chmod +x ./_build.sh
    ./_build.sh
    """
