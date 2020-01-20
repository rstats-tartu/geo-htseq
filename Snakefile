import os
SIMG = "shub://tpall/geo-rnaseq"
LAST_DATE = "2018-12-31"
QUERY = 'expression profiling by high throughput sequencing[DataSet Type] AND ("2000-01-01"[PDAT] : "{}"[PDAT])'.format(LAST_DATE)
EMAIL = "taavi.pall@ut.ee"

K = 10
N = 10
rule all:
  input: 
    expand("output/gsem_{k}.rds", k = list(range(0, K, 1))), 
    "output/suppdata.rds", 
    "_main.html"


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
    retmax = 200,
    batch_size = 1
  conda:
    "envs/geo-query.yaml"
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
    retmax = 2000,
    columns = ["Accession"]
  conda:
    "envs/geo-query.yaml"
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
  script:
    "scripts/preprocess/split_document_summaries.py"


# Download supplementary file names
rule download_suppfilenames:
  input: 
    "output/tmp/document_summaries_{k}.csv"
  output: 
    "output/tmp/suppfilenames_{k}.txt"
  conda:
    "envs/geo-query.yaml"
  script:
    "scripts/preprocess/download_suppfilenames.py"


# Filter supplementary file names by filename extension
rule filter_suppfilenames:
  input: 
    rules.download_suppfilenames.output
  output: 
    "output/tmp/suppfilenames_filtered_{k}.txt"
  conda:
    "envs/geo-query.yaml"
  script:
    "scripts/preprocess/filter_suppfilenames.py"


# Download filterd supplementary files
rule download_suppfiles:
  input: 
    rules.filter_suppfilenames.output
  output: 
    touch("output/downloading_suppfiles_{k}.done")
  singularity: SIMG
  script:
    "scripts/preprocess/download_suppfiles.R"


# Parse series matrix files
rule series_matrixfiles:
  input: 
    rules.download_suppfiles.output
  output: 
    "output/gsem_{k}.rds"
  singularity: SIMG
  script:
    "scripts/preprocess/series_matrixfiles.R"


# Supplementary files that kill R
BAD = ["GSE93374_Merged_all_020816_DGE.txt.gz", 
       "GSE88931_RNA-seq_MergedReadCounts.tsv.gz",
       "GSE55385_transcripts_GSE.tsv.gz",
       "GSE74549_ChIP_1kbWindows_correctedReadCount.txt.gz",
       "GSE77213_Nguyen_GEO_TN03_tallies.xls",
       "GSE77213_Nguyen_GEO_TN05_tallies_total.xls",
       "GSE60012_100bpTiles_RRBS_Mouse.txt.gz",
       "GSE60415_heatmap-upload.xls",
       "GSE67516_RNA_seq_rep1_diffExp_analysis.xls",
       "GSE53298_processed_data.xlsx",
       "GSE53350_SAGE_KomatsuFurukawa_Processed.xls",
       "GSE60483_AQ_exp.tab.gz",
       "GSE53260_Supplementary_S2.xls",
       "GSE53260_Supplementary_S1.xls",
       "GSE89113_T35_vs_T45.T45.UP.0.05.A3SS.MATS.JunctionCountOnly.txt.gz",
       "GSE89113_T35_vs_T45.T45.UP.0.05.MXE.MATS.JunctionCountOnly.txt.gz",
       "GSE89113_T40_vs_T45.A3SS.MATS.JunctionCountOnly.txt.gz",            
       "GSE89113_T40_vs_T45.T40.UP.0.05.A3SS.MATS.JunctionCountOnly.txt.gz",
       "GSE89113_T40_vs_T45.T45.UP.0.05.A3SS.MATS.JunctionCountOnly.txt.gz",
       "GSE89113_T35_vs_T45.A3SS.MATS.JunctionCountOnly.txt.gz",            
       "GSE89113_T35_vs_T45.MXE.MATS.JunctionCountOnly.txt.gz",             
       "GSE89113_T25_vs_T40.A3SS.MATS.JunctionCountOnly.txt.gz",            
       "GSE99484_Sulfolobus_acidocaldarius_DSM_639_CP000077_.gbk.txt.gz",   
       "GSE121228_gene_expression_anno.txt.gz", 
       "GSE113074_Raw_combined.annotated_counts.tsv.gz", 
       "GSE113074_Corrected_combined.annotated_counts.tsv.gz", 
       "GSE121737_early_and_medium_bud.repGene.txt.gz"]


# Split list of supplementary files to chunks for parsing
rule split_suppfiles:
  input:
    suppfilenames_filtered = rules.filter_suppfilenames.output, 
    gsem = rules.series_matrixfiles.output
  output: 
    temp(expand("output/tmp/supptabs_{{k}}-{n}.rds", n = list(range(0, N, 1))))
  params:
    bad = BAD
  singularity: SIMG
  script:
    "scripts/preprocess/split_suppfiles.R"


# Import supplementary data
rule import_suppfiles:
  input: 
    "output/tmp/supptabs_{k}-{n}.rds"
  output: 
    temp("output/tmp/suppdata_{k}-{n}.rds")
  params:
    bad = BAD
  singularity: SIMG
  script:
    "scripts/preprocess/import_suppfiles.R"


# Merge chunks
rule merge_suppdata:
  input: 
    expand("output/tmp/suppdata_{k}-{n}.rds", k = list(range(0, K, 1)), n = list(range(0, N, 1)))
  output: "output/suppdata.rds"
  singularity: SIMG
  script:
    "scripts/preprocess/merge_suppdata.R"


# Download publication metadata
rule download_publications:
  input: rules.geo_query.output
  output: "output/publications.csv"
  params: 
    email = EMAIL,
    api_key = os.environ["NCBI_APIKEY"]
  conda:
    "envs/geo-query.yaml"
  script:
    "scripts/preprocess/download_publications.py"


# Download citations
rule download_citations:
  input: 
    "output/publications.csv"
  output: 
    "output/scopus_citedbycount.csv"
  params:
    api_key = os.environ["ELSEVIER_GEOSEQ"]
  conda:
    "envs/geo-query.yaml"
  script:
    "scripts/preprocess/download_scopus_citations.py"


# Knit report
rule report:
  input: "index.Rmd", "introduction.Rmd", "methods.Rmd", "results.Rmd", "discussion.Rmd", "references.Rmd", "output/document_summaries.csv", "output/single-cell.csv", "output/suppfilenames.rds", "output/suppfilenames_filtered.rds", "output/gsem.rds", "output/suppdata.rds", "output/publications.csv", "output/scopus_citedbycount.csv"
  output: "_main.html"
  singularity: SIMG
  shell:
    """
    Rscript -e "bookdown::render_site(encoding = 'UTF-8')"
    """
