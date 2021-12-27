import os
import re

# Drop some offending files 
BLACKLIST_FILE = "output/blacklist.txt"
with open(BLACKLIST_FILE) as h:
    BLACKLIST = [os.path.basename(line.rstrip()) for line in h.readlines()]

SUPPFILENAMES_FILE = "output/sample_of_giga_suppfiles.txt"
with open(SUPPFILENAMES_FILE, "r") as f:
    SUPPFILENAMES=[os.path.basename(line.rstrip()) for line in f.readlines() if os.path.basename(line.rstrip()) not in BLACKLIST]

EMAIL="taavi.pall@ut.ee"
p = re.compile("GSE\\d+")

rule all:
    input: expand(["output/tmp/parsed_suppfiles__{suppfilename}__.csv"], suppfilename=SUPPFILENAMES), "output/parsed_suppfiles__giga__.csv"

# Download filterd supplementary files
rule download_suppfiles:
    output: 
        temp("suppl/{suppfilename}")
    log:
        "log/download__{suppfilename}__.log"
    params:
        id=lambda wildcards: p.search(wildcards.suppfilename).group(0),
        stub=lambda wildcards: p.search(wildcards.suppfilename).group(0)[0:-3],
    resources:
        runtime = lambda wildcards, attempt: attempt * 120
    shell:
        """
        curl -sS -o {output[0]} --user anonymous:{EMAIL} ftp://ftp.ncbi.nlm.nih.gov/geo/series/{params.stub}nnn/{params.id}/{output[0]} 2> {log}
        """


# Import supplementary data
rule import_suppfiles:
  input: 
    "suppl/{suppfilename}"
  output: 
    "output/tmp/parsed_suppfiles__{suppfilename}__.csv"
  log:
    "log/import__{suppfilename}__.log"
  params:
    f"--var basemean=10 logcpm=1 rpkm=1 fpkm=1 aveexpr=3.32 --bins 40 --fdr 0.05 --pi0method lfdr --blacklist {BLACKLIST_FILE}"
  conda: 
    "envs/geo-query.yaml"
  resources:
    mem_mb = lambda wildcards, input: max([int(16 * (input.size // 1e6)), 4000]),
    runtime = lambda wildcards, attempt, input: attempt * (120 + int(120 * input.size // 1e9)),
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