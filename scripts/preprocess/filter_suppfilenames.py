import re

keep = "|".join(
    ["\." + i + "(.gz)?$" for i in "tab xlsx diff tsv xls csv txt rtf tar".split(" ")]
)
keep = re.compile(keep)
drop = "filelist.txt|raw.tar$"
drop = re.compile(drop)

input = snakemake.input[0]
output = snakemake.output[0]

with open(input) as i:
    with open(output, "w") as o:
        for line in i:
            ll = line.lower()
            if keep.search(ll) and not drop.search(ll):
                o.write(line)
