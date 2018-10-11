
rule all:
  input: "output/gsem.rds", "output/suppdata.rds"

rule geo_query:
  output: "output/document_summaries.rds"
  conda:
    "envs/r.yaml"
  script:
    "R/01_geo_query.R"

rule download_suppfilenames:
  input: "output/document_summaries.rds"
  output: "output/suppfilenames.rds"
  params: 
    last_date = "2017-06-19"
  conda:
    "envs/r.yaml"
  script:
    "R/02_download_suppfilenames.R"

rule filter_suppfilenames:
  input: 
    suppfilenames = "output/suppfilenames.rds"
  output: 
    suppfilenames_filtered = "output/suppfilenames_filtered.rds"
  conda:
    "envs/r.yaml"
  script:
    "R/03_filter_suppfilenames.R"

rule download_suppfiles:
  input: "output/suppfilenames_filtered.rds"
  output: touch("output/downloading_suppfiles.done")
  conda:
    "envs/r.yaml"
  script:
    "R/04_download_suppfiles.R"

rule series_matrixfiles:
  input: "downloading_suppfiles.done"
  output: "output/gsem.rds"
  conda:
    "envs/r.yaml"
  script:
    "R/05_series_matrixfiles.R"

rule split_suppfiles:
  input: 
    docsums = "output/document_summaries.rds", 
    suppfilenames_filtered = "output/suppfilenames_filtered.rds", 
    gsem = "output/gsem.rds"
  output: temp(dynamic("output/tmp/supptabs_{n}.rds"))
  conda:
    "envs/r.yaml"
  script:
    "R/06_split_suppfiles.R"

rule import_suppfiles:
  input: "output/tmp/supptabs_{n}.rds"
  output: temp("output/tmp/suppdata_{n}.rds")
  conda:
    "envs/r.yaml"
  script:
    "R/06_import_suppfiles.R"
  
rule merge_suppdata:
  input: dynamic("output/tmp/suppdata_{n}.rds")
  output: "output/suppdata.rds"
  conda:
    "envs/r.yaml"
  script:
    "R/07_merge_suppdata.R"
