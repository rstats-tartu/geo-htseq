import os

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = False
            if not entry:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

N = snakemake.params.get("n", 10)
dir = snakemake.params.get("dir", "")
s = []
with open(snakemake.input[0], "r") as h:
    for line in h:
        if "series_matrix.txt.gz" not in line:
            s.append(os.path.join(dir, line.rstrip()))
chunks =  N - 1
parts = os.path.splitext(snakemake.input[0])

for i, batch in enumerate(batch_iterator(iter(s), len(s)//chunks), start=1):
    filename = parts[0] + "_{}" + parts[1]
    with open(filename.format(i), "w") as h:
        h.writelines("{}\n".format(l) for l in batch)

