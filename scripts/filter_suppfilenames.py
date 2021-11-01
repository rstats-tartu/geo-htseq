import re
import os


assert os.stat(snakemake.input[0]).st_size > 0, "File is empty!"
input = snakemake.input[0]
output = snakemake.output[0]

drop = "series_matrix\.txt\.gz|filelist\.txt|readme|\.bam|\.sam|\.csfasta|\.fa(sta)?|\.f(a|n)a|(big)?wig|\.bed(graph)?|(broad_)?lincs"
drop = re.compile(drop)

with open(input) as i:
    with open(output, "w") as o:
        for line in i:
            ll = line.lower()
            if not drop.search(ll):
                o.write(line)
