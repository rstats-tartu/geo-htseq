import os
LAST_DATE = "2020-12-31"
QUERY = 'expression profiling by high throughput sequencing[DataSet Type] AND ("2000-01-01"[PDAT] : "{}"[PDAT])'.format(LAST_DATE)
EMAIL = "taavi.pall@ut.ee"

onsuccess:
    print("Workflow finished, no error")
    shell("mail -s 'Forkflow finished successfully' {EMAIL} < {log}")

onerror:
    print("An error occurred")
    shell("mail -s 'An error occurred' {EMAIL} < {log}")

localrules: all, filter_suppfilenames, suppfiles_list

K = 15
N = 10
rule all:
  input: 
    "output/document_summaries.csv",
    "output/single-cell.csv",
    "output/parsed_suppfiles.csv", 
    "output/publications.csv",
    "output/scopus_citedbycount.csv",
    "output/suppfilenames.txt",
    "output/suppfilenames_filtered.txt",
    "output/spots.csv",
    "output/geo-htseq-until-{}.tar.gz".format(LAST_DATE),
    expand("output/downloading_suppfiles_{k}.done", 
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
    retmax = 100000,
    batch_size = 500
  conda:
    "envs/geo-query.yaml"
  resources:
    runtime = 120
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
    retmax = 25000,
    columns = ["Accession"]
  conda:
    "envs/geo-query.yaml"
  resources:
    runtime = 120
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
    runtime = 120
  script:
    "scripts/split_df.py"


# Download supplementary file names
rule download_suppfilenames:
  input: 
    "output/tmp/document_summaries_{k}.csv"
  output: 
    "output/tmp/suppfilenames_{k}.txt"
  log:
    "log/download_suppfilenames_{k}.log"
  params:
    email = EMAIL,
    dirs = "suppl",
    size = 200,
  conda:
    "envs/geo-query.yaml"
  resources:
    runtime = lambda wildcards, attempt: 90 + (attempt * 30)
  shell:
    """
    python3 -u scripts/download_suppfilenames.py --input {input} --output {output} --email {params.email} --dirs {params.dirs} size = {params.size} 2> {log}
    """

# Merge suppfilenames
rule suppfilenames:
  input:
    expand("output/tmp/suppfilenames_{k}.txt", k = list(range(0, K, 1)))
  output:
    "output/suppfilenames.txt"
  resources:
    runtime = 120
  shell:
    """
    for file in {input}; do grep "^suppl" $file >> {output}; done
    """

# Download read run data
rule download_spots:
  input: 
    "output/tmp/document_summaries_{k}.csv"
  output: 
    "output/tmp/spots_{k}.csv"
  params:
    email = EMAIL,
    api_key = os.environ["NCBI_APIKEY"],
    retmax = 100,
    max_tries = 3
  shadow:
    "minimal"
  conda:
    "envs/geo-query.yaml"
  resources:
    runtime = 2440
  script:
    "scripts/read_runs.py"


rule merge_spots:
    input:
        expand("output/tmp/spots_{k}.csv", k = list(range(0, K, 1)))
    output:
        "output/spots.csv"
    resources:
        runtime = 120
    run:
        import pandas as pd
        with open(output[0], "a") as output_handle:
  	      for file in input:
        	  spots = pd.read_csv(file, sep=",")
          	  spots.to_csv(
                	output_handle,
                	sep=",",
                	mode="a",
                	header=not output_handle.tell(),
                	index=False,
           		)	


# Filter supplementary file names by filename extension
rule filter_suppfilenames:
  input: 
    rules.download_suppfilenames.output
  output: 
    "output/tmp/suppfilenames_filtered_{k}.txt"
  conda:
    "envs/geo-query.yaml"
  resources:
    runtime = 120
  script:
    "scripts/filter_suppfilenames.py"


rule suppfilenames_filtered:
  input:
    expand("output/tmp/suppfilenames_filtered_{k}.txt", k = list(range(0, K, 1)))
  output:
    "output/suppfilenames_filtered.txt"
  resources:
    runtime = 20
  shell:
    """
    cat {input} > {output}
    """

# Download filterd supplementary files
rule download_suppfiles:
  input: 
    rules.filter_suppfilenames.output
  output: 
    touch("output/downloading_suppfiles_{k}.done")
  log:
    "log/download_suppfiles_{k}.log"
  params:
    email = EMAIL,
    size = 200,
    dir = "",
  conda:
    "envs/geo-query.yaml"
  resources:
    runtime = 1440 #lambda wildcards, attempt: 90 + (attempt * 30)
  shell:
    """
    python3 -u scripts/download_suppfiles.py --input {input} --email {params.email} --size {params.size} --dir {params.dir} 2> {log}
    """


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
    runtime = 120
  script:
    "scripts/split_lines.py"


# Import supplementary data
rule import_suppfiles:
  input: 
    "output/tmp/suppfilenames_filtered_{k}_{n}.txt"
  output: 
    "output/tmp/parsed_suppfiles_{k}_{n}.csv"
  log:
    "log/import_suppfiles_{k}_{n}.log"
  params:
    f"--var basemean=10 logcpm=1 rpkm=1 fpkm=1 aveexpr=3.32 --bins 40 --fdr 0.05 --pi0method lfdr -v --blacklist {BLACKLIST_FILE}"
  conda: 
    "envs/geo-query.yaml"
  resources:
    mem_mb = 64000,
    runtime = 480,
  shell:
    """
    python3 -u scripts/import_suppfiles.py --list {input} --out {output} {params} 2> {log}
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
    runtime = 120
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
    api_key = os.environ["NCBI_APIKEY"],
    batch_size = 500
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
    runtime = lambda wildcards, attempt: 120 + (attempt * 60)
  script:
    "scripts/download_scopus_citations.py"


rule archive:
    input:
        "output/document_summaries.csv",
        "output/single-cell.csv",
        "output/parsed_suppfiles.csv", 
        "output/publications.csv",
        "output/scopus_citedbycount.csv",
        "output/suppfilenames.txt",
        "output/suppfilenames_filtered.txt",
        "output/spots.csv"
    output:
        "output/geo-htseq-until-{}.tar.gz".format(LAST_DATE)
    shell:
        "tar -czvf {output[0]} {input}"
