import os
SIMG = "shub://tpall/geo-rnaseq"
LAST_DATE = "2019-12-31"
QUERY = 'expression profiling by high throughput sequencing[DataSet Type] AND ("2000-01-01"[PDAT] : "{}"[PDAT])'.format(LAST_DATE)
EMAIL = "taavi.pall@ut.ee"

localrules: all, filter_suppfilenames, suppfiles_list

K = 10
N = 10
rule all:
  input: 
    "output/document_summaries.csv",
    "output/single-cell.csv",
    "output/parsed_suppfiles.csv", 
    "output/publications.csv",
    "output/scopus_citedbycount.csv",
    expand(["output/tmp/suppfilenames_filtered_{k}.txt", 
            "output/downloading_suppfiles_{k}.done"], 
            k = list(range(0, K, 1)))


# Queries HT-seq expression profiling experiments
# Requires NCBI api_key as NCBI_APIKEY environment variable
rule geo_query:
  output: 
    "output/document_summaries.csv"
  params:
    email = EMAIL,
    api_key = os.environ["NCBI_APIKEY"],
    query = QUERY,
    db = "gds",
    retmax = 40000,
    batch_size = 100
  conda:
    "envs/geo-query.yaml"
  resources:
    runtime = 90
  script:
    "scripts/geo_query.py"


# Single-cell experiment accessions
rule single_cell:
  output: 
    "output/single-cell.csv"
  params:
    email = EMAIL,
    api_key = os.environ["NCBI_APIKEY"],
    query = QUERY + ' AND "single-cell"[All Fields]',
    db = "gds",
    retmax = 5000,
    columns = ["Accession"]
  conda:
    "envs/geo-query.yaml"
  resources:
    runtime = 30
  script:
    "scripts/geo_query.py"


# Split GEO document summaries
rule split_document_summaries:
  input:
    rules.geo_query.output
  output:
    expand("output/tmp/document_summaries_{k}.csv", k = list(range(0, K, 1)))
  params:
    chunks = K
  conda:
    "envs/geo-query.yaml"
  resources:
    runtime = 30
  script:
    "scripts/split_df.py"


# Download supplementary file names
rule download_suppfilenames:
  input: 
    "output/tmp/document_summaries_{k}.csv"
  output: 
    "output/tmp/suppfilenames_{k}.txt"
  params:
    email = EMAIL
  conda:
    "envs/geo-query.yaml"
  resources:
    runtime = lambda wildcards, attempt: 100 + (attempt * 20)
  shell:
    """
    python3 -u scripts/download_suppfilenames.py --list {input} --out {output} --email {params.email}
    """


# Filter supplementary file names by filename extension
rule filter_suppfilenames:
  input: 
    rules.download_suppfilenames.output
  output: 
    "output/tmp/suppfilenames_filtered_{k}.txt"
  conda:
    "envs/geo-query.yaml"
  resources:
    runtime = 30
  script:
    "scripts/filter_suppfilenames.py"


# Download filterd supplementary files
rule download_suppfiles:
  input: 
    rules.filter_suppfilenames.output
  output: 
    touch("output/downloading_suppfiles_{k}.done")
  conda:
    "envs/geo-query.yaml"
  resources:
    runtime = lambda wildcards, attempt: 100 + (attempt * 20)
  script:
    "scripts/download_suppfiles.py"


# Split list of supplementary files
# dir is the location of suppl/ folder
# Drop some offending files 
BLACKLIST_FILE = "output/blacklist.txt"
with open(BLACKLIST_FILE) as h:
    BLACKLIST = [os.path.basename(i.rstrip()) for i in h.readlines()]

rule suppfiles_list:
  input: 
    rules.filter_suppfilenames.output,
    rules.download_suppfiles.output
  output: 
    expand("output/tmp/suppfilenames_filtered_{{k}}_{n}.txt", n = list(range(0, N, 1)))
  params:
    chunks = N,
    dir = "output",
    blacklist = BLACKLIST
  conda: 
    "envs/geo-query.yaml"
  resources:
    runtime = 30
  script:
    "scripts/split_lines.py"


# Import supplementary data
rule import_suppfiles:
  input: 
    "output/tmp/suppfilenames_filtered_{k}_{n}.txt"
  output: 
    "output/tmp/parsed_suppfiles_{k}_{n}.csv"
  params:
    "--var basemean=10 logcpm=1 rpkm=1 fpkm=1 aveexpr=3.32 --bins 40 --fdr 0.05 -v --blacklist {}".format(BLACKLIST_FILE)
  conda: 
    "envs/geo-query.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: 8000 + (attempt * 8000),
    runtime = lambda wildcards, attempt: 40 + (attempt * 20)
  shell:
    """
    python3 -u scripts/import_suppfiles.py --list {input} --out {output} {params}
    """


# Merge chunks
rule merge_parsed_suppfiles:
  input: 
    expand("output/tmp/parsed_suppfiles_{k}_{n}.csv", k = list(range(0, K, 1)), n = list(range(1, N, 1)))
  output: 
    "output/parsed_suppfiles.csv"
  conda: 
    "envs/geo-query.yaml"
  resources:
    mem_mb = 4000,
    runtime = 30
  script:
    "scripts/concat_tabs.py"
    

# Download publication metadata
rule download_publications:
  input: 
    rules.geo_query.output
  output: 
    "output/publications.csv"
  params: 
    email = EMAIL,
    api_key = os.environ["NCBI_APIKEY"]
  conda:
    "envs/geo-query.yaml"
  resources:
    runtime = 360
  script:
    "scripts/download_publications.py"


# Download citations
rule download_citations:
  input: 
    rules.geo_query.output
  output: 
    "output/scopus_citedbycount.csv"
  params:
    api_key = os.environ["ELSEVIER_GEOSEQ"]
  conda:
    "envs/geo-query.yaml"
  resources:
    runtime = 120
  script:
    "scripts/download_scopus_citations.py"
