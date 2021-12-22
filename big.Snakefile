import os
import re
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

with open("output/sample_of_giga_suppfiles.txt", "r") as f:
    SUPPFILENAMES=[os.path.basename(line.rstrip()) for line in f.readlines()]

EMAIL="taavi.pall@ut.ee"
FTP = FTPRemoteProvider(username="anonymous", password=EMAIL)
p = re.compile("GSE\\d+")

def get_url(wildcards):
    id = p.search(wildcards.suppfilename).group(0)
    return os.path.join("ftp.ncbi.nlm.nih.gov", "geo/series", id[0:-3] + "nnn", id, "suppl", wildcards.suppfilename)

rule all:
    input: expand(["output/tmp/parsed_suppfiles__{suppfilename}__.csv"], suppfilename=SUPPFILENAMES), "output/parsed_suppfiles__giga__.csv"

# Download filterd supplementary files
rule download_suppfiles:
    input: lambda wildcards: FTP.remote(get_url(wildcards), keep_local=True, immediate_close=True)
    output: 
        temp("suppl/{suppfilename}")
    resources:
        runtime = lambda wildcards: 120 + int(60 * input.size // 1e9)
    shell:
        "mv {input} {output}"


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
    runtime = lambda wildcards: 120 + int(60 * input.size // 1e9),
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
