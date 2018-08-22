
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
  input: "output/suppfilenames.rds"
  output: expand("output/suppfilenames_filtered.{ext}", ext = ['rds','txt'])
  conda:
    "envs/r.yaml"
  script:
    "R/03_filter_suppfilenames.R"

rule download_suppfiles:
  input: "output/suppfilenames_filtered.rds"
  output: touch("downloading_suppfiles.done")
  conda:
    "envs/r.yaml"
  script:
    "R/04_download_suppfiles.R"

rule munge_matrixfiles:
  input: "downloading_suppfiles.done"
  output: "output/gsem.rds"
  conda:
    "envs/r.yaml"
  script:
    "R/05_munge_series_matrixfiles.R"

rule munge_suppfiles:
  input: docsums = "output/document_summaries.rds", gsem = "output/gsem.rds"
  output: dynamic("output/tmp/suppdata_{n}.rds")
  conda:
    "envs/r.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 8000
  script:
    "R/06_munge_suppfiles_slurm.R"
  
rule merge_suppdata:
  input: dynamic("output/tmp/suppdata_{n}.rds")
  output: "output/suppdata.rds"
  conda:
    "envs/r.yaml"
  script:
    "R/07_merge_suppdata.R"
