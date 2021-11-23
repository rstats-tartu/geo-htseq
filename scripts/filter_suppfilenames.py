import re
import os


assert os.stat(snakemake.input[0]).st_size > 0, "File is empty!"
input = snakemake.input[0]
output = snakemake.output[0]

drop = "series_matrix\.txt\.gz$|filelist\.txt$|readme|\.bam(\.tdf|$)|\.bai(\.gz|$)|\.sam(\.gz|$)|\.csfasta|\.fa(sta)?(\.gz|\.bz2|\.txt\.gz|$)|\.f(a|n)a(\.gz|$)|(big)?wig|\.bw(\.|$)|\.bed(graph)?(\.tdf|\.gz|\.bz2|\.tar\.gz|\.txt\.gz|$)|(broad_)?lincs"
drop = re.compile(drop)

with open(input) as i:
    with open(output, "w") as o:
        for line in i:
            ll = line.lower()
            if not drop.search(ll):
                o.write(line)
