import os
SIMG = "shub://tpall/geo-rnaseq"
last_date = "2018-12-31"

rule all:
  input: "output/gsem.rds", "output/suppdata.rds", "_main.html"

# Queries HT-seq expression profiling experiments
# Requires NCBI api_key as NCBI_APIKEY environment variable
rule geo_query:
  output: 
    "output/document_summaries.csv"
  params:
    email = "taavi.pall@ut.ee",
    api_key = os.environ["NCBI_APIKEY"]
  conda:
    "envs/geo-query.yaml"
  script:
    "scripts/preprocess/01_geo_query.py"

rule download_suppfilenames:
  input: 
    rules.geo_query.output
  output: 
    "output/suppfilenames.rds"
  params: 
    last_date = last_date
  singularity: SIMG
  script:
    "scripts/preprocess/02_download_suppfilenames.R"

rule filter_suppfilenames:
  input: 
    rules.download_suppfilenames.output
  output: 
    "output/suppfilenames_filtered.rds"
  singularity: SIMG
  script:
    "scripts/preprocess/03_filter_suppfilenames.R"

rule download_suppfiles:
  input: 
    rules.filter_suppfilenames.output
  output: 
    touch("output/downloading_suppfiles.done")
  singularity: SIMG
  script:
    "scripts/preprocess/04_download_suppfiles.R"

rule series_matrixfiles:
  input: 
    rules.download_suppfiles.output
  output: 
    "output/gsem.rds"
  singularity: SIMG
  script:
    "scripts/preprocess/05_series_matrixfiles.R"

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

n_files = list(range(1, 101, 1))
rule split_suppfiles:
  input:
    suppfilenames_filtered = rules.filter_suppfilenames.output, 
    gsem = rules.series_matrixfiles.output
  output: 
    temp(expand("output/tmp/supptabs_{n}.rds", n = n_files))
  params:
    bad = BAD
  singularity: SIMG
  script:
    "scripts/preprocess/06_split_suppfiles.R"

rule import_suppfiles:
  input: 
    "output/tmp/supptabs_{n}.rds"
  output: 
    temp("output/tmp/suppdata_{n}.rds")
  params:
    bad = BAD
  singularity: SIMG
  script:
    "scripts/preprocess/06_import_suppfiles.R"
  
rule merge_suppdata:
  input: 
    expand("output/tmp/suppdata_{n}.rds", n = n_files)
  output: "output/suppdata.rds"
  singularity: SIMG
  script:
    "scripts/preprocess/07_merge_suppdata.R"

rule download_publications:
  input: rules.geo_query.output
  output: "output/publications.csv"
  params: 
    last_date = last_date
  singularity: SIMG
  script:
    "scripts/preprocess/07_download_publications.R"

rule download_citations:
  input: 
    pubs = "output/publications.csv", 
    document_summaries = "output/document_summaries.csv"
  output: "output/scopus_citedbycount.csv"
  params: 
    last_date = last_date,
    api_key = os.environ["ELSEVIER_GEOSEQ"]
  singularity: SIMG
  script:
    "scripts/preprocess/download_scopus_citations.R"

rule report:
  input: "index.Rmd", "01_introduction.Rmd", "02_methods.Rmd", "03_results.Rmd", "04_discussion.Rmd", "05_references.Rmd", "output/document_summaries.csv", "output/suppfilenames.rds", "output/suppfilenames_filtered.rds", "output/gsem.rds", "output/suppdata.rds", "output/publications.csv", "output/scopus_citedbycount.csv"
  output: "_main.html"
  singularity: SIMG
  shell:
    """
    Rscript -e "bookdown::render_site(encoding = 'UTF-8')"
    """
