import numpy as np
import os

input = snakemake.input[0]
output = snakemake.output
chunks = int(snakemake.params.get("chunks", 1))
dir = snakemake.params.get("dir", "")
blacklist = snakemake.params.get("blacklist", [])

assert len(output) == chunks, "The number of chunks and outputs don't match!"
assert isinstance(drop, list), "Drop must be list"

with open(input) as i:
    lines = i.readlines()

# Drop series matrix files
lines = [i for i in lines if "series_matrix.txt.gz" not in i]

# Drop some extra files (e.g. too large)
lines = [i for i in lines if os.path.basename(i.rstrip()) not in blacklist]

# Split lines to chunks
input_chunks = np.array_split(lines, chunks)

# Write chunks to files
for item, file in zip(input_chunks, output):
    with open(file, "w") as f:
        for line in item:
            f.write(os.path.join(dir, line))
