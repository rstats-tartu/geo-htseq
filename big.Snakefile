import os
import re
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

with open("output/sample_of_giga_suppfiles.txt", "r") as f:
    SUPPFILENAMES=[os.path.basename(line.rstrip()) for line in f.readlines()]

EMAIL="taavi.pall@ut.ee"

rule all:
    input: expand(["output/tmp/parsed_suppfiles__{suppfilename}__.csv"], suppfilename=SUPPFILENAMES), "output/parsed_suppfiles__giga__.csv"

# Download filterd supplementary files
rule download_suppfiles:
    output: 
        temp("suppl/{suppfilename}")
    log:
        "log/download__{suppfilename}__.log"
    conda: 
        "envs/geo-query.yaml"
    resources:
        runtime = lambda wildcards, attempt: attempt * 120
    shell:
        """
        python3 -u scripts/download_suppfiles.py --file {output[0]} --email {EMAIL} --maxfilesize 10000000000 2> {log}
        """


# Split list of supplementary files
# dir is the location of suppl/ folder
# Drop some offending files 
BLACKLIST_FILE = "output/blacklist.txt"
with open(BLACKLIST_FILE) as h:
    BLACKLIST = [os.path.basename(i.rstrip()) for i in h.readlines()]

# Import supplementary data
rule import_suppfiles:
  input: 
    "suppl/{suppfilename}"
  output: 
    "output/tmp/parsed_suppfiles__{suppfilename}__.csv"
  log:
    "log/import__{suppfilename}__.log"
  params:
    f"--var basemean=10 logcpm=1 rpkm=1 fpkm=1 aveexpr=3.32 --bins 40 --fdr 0.05 --pi0method lfdr -v --blacklist {BLACKLIST_FILE}"
  conda: 
    "envs/geo-query.yaml"
  resources:
    mem_mb = lambda wildcards, input: max([int(4 * (input.size // 1e6)), 4000]),
    runtime = lambda wildcards, input: 120 + int(60 * input.size // 1e9),
  shell:
    """
    python3 -u scripts/import_suppfiles.py --file {input} --out {output} {params} 2> {log}
    """

rule merge_parsed_suppfiles:
  input: 
    expand("output/tmp/parsed_suppfiles__{suppfilename}__.csv", suppfilename=SUPPFILENAMES)
  output: 
    "output/parsed_suppfiles__giga__.csv"
  conda: 
    "envs/geo-query.yaml"
  resources:
    mem_mb = 4000,
    runtime = 120,
  script:
    "scripts/concat_tabs.py"
