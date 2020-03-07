import numpy as np
import os

input = snakemake.input[0]
output = snakemake.output
chunks = int(snakemake.params.get("chunks", 1))
dir = snakemake.params.get("dir", "")

assert len(output) == chunks, "The number of chunks and outputs don't match!"

with open(input) as i:
    lines = i.readlines()

lines = [i for i in lines]
input_chunks = np.array_split(lines, chunks)

for item, file in zip(input_chunks, output):
    with open(file, "w") as f:
        for line in item:
            f.write(os.path.join(dir, line))
