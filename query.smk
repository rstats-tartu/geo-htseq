import os

LAST_DATE = "2020-12-31"
EMAIL = "taavi.pall@ut.ee"
QUERY = 'expression profiling by high throughput sequencing[DataSet Type] AND ("2000-01-01"[PDAT] : "{}"[PDAT])'.format(
    LAST_DATE
)


K = 15


rule all:
    input:
        "output/document_summaries.csv",
        "output/single-cell.csv",
        "output/publications.csv",
        "output/scopus_citedbycount.csv",
        "output/suppfilenames.txt",
        "output/spots.csv",


# Queries HT-seq expression profiling experiments
# Requires NCBI api_key as NCBI_APIKEY environment variable
rule geo_query:
    output:
        "output/document_summaries.csv",
    params:
        email=EMAIL,
        api_key=os.environ["NCBI_APIKEY"],
        query=QUERY,
        db="gds",
        retmax=100000,
        batch_size=500,
    conda:
        "envs/environment.yaml"
    resources:
        runtime=120,
    script:
        "scripts/geo_query.py"


# Single-cell experiment accessions
rule single_cell:
    output:
        "output/single-cell.csv",
    params:
        email=EMAIL,
        api_key=os.environ["NCBI_APIKEY"],
        query=QUERY + ' AND "single-cell"[All Fields]',
        db="gds",
        retmax=25000,
        columns=["Accession"],
    conda:
        "envs/environment.yaml"
    resources:
        runtime=120,
    script:
        "scripts/geo_query.py"


# Split GEO document summaries
rule split_document_summaries:
    input:
        rules.geo_query.output,
    output:
        expand("output/tmp/document_summaries_{k}.csv", k=list(range(0, K, 1))),
    params:
        chunks=K,
    conda:
        "envs/environment.yaml"
    resources:
        runtime=120,
    script:
        "scripts/split_df.py"


# Download supplementary file names
rule download_suppfilenames:
    input:
        "output/tmp/document_summaries_{k}.csv",
    output:
        "output/tmp/suppfilenames_{k}.txt",
    log:
        "log/download_suppfilenames_{k}.log",
    params:
        email=EMAIL,
        dirs="suppl",
        size=200,
    conda:
        "envs/environment.yaml"
    resources:
        runtime=lambda wildcards, attempt: 90 + (attempt * 30),
    shell:
        """
        python3 -u scripts/download_suppfilenames.py --input {input} --output {output} --email {params.email} --dirs {params.dirs} --size {params.size} 2> {log}
        """


# Merge suppfilenames
checkpoint suppfilenames:
    input:
        expand("output/tmp/suppfilenames_{k}.txt", k=list(range(0, K, 1))),
    output:
        "output/suppfilenames.txt",
    resources:
        runtime=120,
    shell:
        """
        for file in {input}; do grep "^suppl" $file >> {output}; done
        """


# Download read run data
rule download_spots:
    input:
        "output/tmp/document_summaries_{k}.csv",
    output:
        "output/tmp/spots_{k}.csv",
    params:
        email=EMAIL,
        api_key=os.environ["NCBI_APIKEY"],
        retmax=100,
        max_tries=3,
    shadow:
        "minimal"
    conda:
        "envs/environment.yaml"
    resources:
        runtime=2440,
    script:
        "scripts/read_runs.py"


rule merge_spots:
    input:
        expand("output/tmp/spots_{k}.csv", k=list(range(0, K, 1))),
    output:
        "output/spots.csv",
    resources:
        runtime=120,
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


# Download publication metadata
rule download_publications:
    input:
        rules.geo_query.output,
    output:
        "output/publications.csv",
    params:
        email=EMAIL,
        api_key=os.environ["NCBI_APIKEY"],
        batch_size=500,
    conda:
        "envs/environment.yaml"
    resources:
        runtime=360,
    script:
        "scripts/download_publications.py"


# Download citations
rule download_citations:
    input:
        rules.geo_query.output,
    output:
        "output/scopus_citedbycount.csv",
    params:
        api_key=os.environ["ELSEVIER_GEOSEQ"],
    conda:
        "envs/environment.yaml"
    resources:
        runtime=lambda wildcards, attempt: 120 + (attempt * 60),
    script:
        "scripts/download_scopus_citations.py"
