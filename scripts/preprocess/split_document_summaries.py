import pandas as pd
import numpy as np

input = snakemake.input[0]
output = snakemake.output
chunks = int(snakemake.params.get("chunks", 1))
assert len(output) == chunks, "The number of chunks and outputs don't match!"

ds = pd.read_csv(input, index_col="Accession")
input_chunks = np.array_split(ds, chunks)

for item, file in zip(input_chunks, output):
    item.to_csv(file)
