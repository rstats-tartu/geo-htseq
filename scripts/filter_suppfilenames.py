import re
import os


assert os.stat(snakemake.input[0]).st_size > 0, "File is empty!"
input = snakemake.input[0]
output = snakemake.output[0]

keep = "|".join(
    ["\." + i + "(\.gz)?$" for i in "tab xlsx diff tsv xls csv txt rtf tar".split(" ")]
)
keep = re.compile(keep)
drop = "series_matrix\.txt\.gz|filelist\.txt|_raw\.tar$|readme|\.bam|\.sam|\.csfasta|\.fa(sta)?|(big)?wig|\.bed(graph)?|(broad_)?lincs"
drop = re.compile(drop)

with open(input) as i:
    with open(output, "w") as o:
        for line in i:
            ll = line.lower()
            if not drop.search(ll):
                o.write(line)
