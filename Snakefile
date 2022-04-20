import os
import re
from scripts.import_suppfiles import drop
from snakemake.utils import min_version

min_version("6.0")
EMAIL = "taavi.pall@ut.ee"
BLACKLIST_FILE = "blacklist.txt"
p = re.compile("GSE\\d+")


module query:
    snakefile:
        "query.smk"


use rule * from query as query_*


# Suppfilenames
def get_parsed_suppfiles(wildcards):
    SUPPFILENAMES_FILE = checkpoints.query_suppfilenames.get().output[0]
    with open(SUPPFILENAMES_FILE, "r") as f:
        SUPPFILENAMES = [
            os.path.basename(line.rstrip())
            for line in f.readlines()
            if line.strip() != ""
        ]

    with open(BLACKLIST_FILE) as h:
        BLACKLIST = [os.path.basename(line.rstrip()) for line in h.readlines()]

    # Drop alignment files, sequence files, binary files, README, etc
    SUPPFILENAMES = [
        file for file in SUPPFILENAMES if not bool(drop.search(file.lower()))
    ]

    # Drop files in blacklist
    SUPPFILENAMES = [file for file in SUPPFILENAMES if file not in BLACKLIST]
    return expand(
        ["output/tmp/parsed_suppfiles__{suppfilename}__.csv"],
        suppfilename=SUPPFILENAMES,
    )


rule all:
    input:
        rules.query_all.input,
        "output/parsed_suppfiles.csv",
        "output/geo-htseq.tar.gz",
    default_target: True


# Download filterd supplementary files
rule download_suppfiles:
    output:
        temp("suppl/{suppfilename}"),
    log:
        "log/download__{suppfilename}__.log",
    params:
        id=lambda wildcards: p.search(wildcards.suppfilename.upper()).group(0),
        stub=lambda wildcards: p.search(wildcards.suppfilename.upper()).group(0)[0:-3],
    resources:
        runtime=lambda wildcards, attempt: attempt * 1440,
    shell:
        """
        curl -sS -H 'Expect:' -o {output[0]} --user anonymous:{EMAIL} ftp://ftp.ncbi.nlm.nih.gov/geo/series/{params.stub}nnn/{params.id}/{output[0]} 2> {log}
        """


# Import supplementary data
rule import_suppfiles:
    input:
        "suppl/{suppfilename}",
    output:
        "output/tmp/parsed_suppfiles__{suppfilename}__.csv",
    log:
        "log/import__{suppfilename}__.log",
    params:
        f"--var basemean=10 logcpm=1 rpkm=1 fpkm=1 aveexpr=3.32 --bins 40 --fdr 0.05 --pi0method lfdr --blacklist {BLACKLIST_FILE}",
    conda:
        "envs/environment.yaml"
    resources:
        mem_mb=lambda wildcards, input: min(
            200000, max([int(32 * (input.size // 1e6)), 4000])
        ),
        runtime=lambda wildcards, attempt, input: min(
            1440, attempt * (5 + int(360 * input.size // 1e9))
        ),
    shell:
        """
        python3 -u scripts/import_suppfiles.py --file {input} --out {output} {params} 2> {log}
        """


rule merge_parsed_suppfiles:
    input:
        get_parsed_suppfiles,
    output:
        "output/parsed_suppfiles.csv",
    conda:
        "envs/environment.yaml"
    resources:
        mem_mb=4000,
        runtime=120,
    script:
        "python3 -u scripts/concat_tabs.py --tabs {input} --out {output}"


rule archive:
    input:
        rules.query_all.input,
        "output/parsed_suppfiles.csv",
    output:
        "output/geo-htseq.tar.gz",
    shell:
        "tar -czvf {output[0]} {input}"
