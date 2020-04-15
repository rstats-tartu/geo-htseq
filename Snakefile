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
    retmax = 50000,
    batch_size = 100
  conda:
    "envs/geo-query.yaml"
  resources:
    runtime = 90
  script:
    "scripts/preprocess/geo_query.py"


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
    "scripts/preprocess/geo_query.py"


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
    "scripts/preprocess/split_df.py"


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
    python3 -u scripts/preprocess/download_suppfilenames.py --list {input} --out {output} --email {params.email}
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
    "scripts/preprocess/filter_suppfilenames.py"


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
    "scripts/preprocess/download_suppfiles.py"


# Split list of supplementary files
# dir is the location of suppl/ folder
# Drop some offending files 
MEM=["GSE84086_T.S.meth.txt.gz", 
"GSE118184_Organoid.combined.dge.txt.gz", 
"GSE115746_cells_exon_counts.csv.gz",
"GSE81750_Mosaic-Seq_15SE_71HS_241sgRNA_K562.full_matrix.final.txt.gz",
"GSE121891_OB_6_runs.raw.dge.csv.gz",
"GSE117131_Zhu_SC-10X_mouse_raw-counts.txt.gz"]

TIME=["GSE116470_F_GRCm38.81.P60Cerebellum_ALT.raw.dge.txt.gz", 
"GSE30839_GSM765295_ribo_mesc_harr150s-profiles.tar.gz"]

rule suppfiles_list:
  input: 
    "output/tmp/suppfilenames_filtered_{k}.txt"
  output: 
    expand("output/tmp/suppfilenames_filtered_{{k}}_{n}.txt", n = list(range(0, N, 1)))
  params:
    chunks = N,
    dir = "output",
    drop = MEM + TIME
  conda: 
    "envs/geo-query.yaml"
  resources:
    runtime = 30
  script:
    "scripts/preprocess/split_lines.py"



# Import supplementary data
rule import_suppfiles:
  input: 
    "output/tmp/suppfilenames_filtered_{k}_{n}.txt"
  output: 
    "output/tmp/parsed_suppfiles_{k}_{n}.csv"
  params:
    "--var basemean=10 logcpm=1 rpkm=1 fpkm=1 aveexpr=3.32 --bins 40 --fdr 0.05 -v"
  conda: 
    "envs/geo-query.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: 8000 + (attempt * 8000),
    runtime = lambda wildcards, attempt: 40 + (attempt * 20)
  shell:
    """
    python3 -u scripts/preprocess/import_suppfiles.py --list {input} --out {output} {params}
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
    "scripts/preprocess/concat_tabs.py"
    


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
    "scripts/preprocess/download_publications.py"


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
    "scripts/preprocess/download_scopus_citations.py"
