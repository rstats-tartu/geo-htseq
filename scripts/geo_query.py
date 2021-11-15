from Bio import Entrez
import pandas as pd
import json

# Helper functions
def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i : i + n]


def extract(d, keys):
    return dict((k, d[k]) for k in keys if k in d)


COLUMNS = [
    "Id",
    "Accession",
    "GDS",
    "title",
    "summary",
    "GPL",
    "GSE",
    "taxon",
    "entryType",
    "gdsType",
    "ptechType",
    "valType",
    "SSInfo",
    "subsetInfo",
    "PDAT",
    "suppFile",
    "Samples",
    "n_samples",
    "SeriesTitle",
    "PlatformTitle",
    "PlatformTaxa",
    "SamplesTaxa",
    "PubMedIds",
    "Projects",
    "FTPLink",
]


# Parse parameters
Entrez.email = snakemake.params.get("email", None)
if Entrez.email is None:
    raise ValueError("An email must be provided")

api_key = snakemake.params.get("api_key", None)
if api_key is None:
    print(
        "Personal API key from NCBI. If not set, only 3 queries per second are allowed. 10 queries per seconds otherwise with a valid API key."
    )
else:
    Entrez.api_key = api_key

query = snakemake.params.get("query", None)
if query is None:
    raise ValueError("A query string must be provided")

db = snakemake.params.get("db", "gds")
retmax = snakemake.params.get("retmax", 10)
batch_size = snakemake.params.get("batch_size", 1)
columns = snakemake.params.get("columns", COLUMNS)

# Run query
search = Entrez.esearch(db=db, term=query, retmax=retmax)
ids = Entrez.read(search)
chunked_ids = chunks(ids["IdList"], batch_size)


# Fetch and parse summaries
with open(snakemake.output[0], "a") as f:
    for chunk in chunked_ids:
        summary = Entrez.esummary(
            db=db, id=",".join(chunk), retmode="xml", retmax=batch_size
        )
        records = Entrez.parse(summary)
        docsums = pd.DataFrame(columns=columns)
        for record in records:
            record = extract(record, columns)
            docsums = docsums.append(
                {k: v for k, v in record.items()}, ignore_index=True
            )
        if "PubMedIds" in docsums.columns:
            docsums["PubMedIds"] = docsums.PubMedIds.apply(
                lambda x: ";".join(pd.Series(x, dtype="str"))
            )
        if "Projects" in docsums.columns:
            docsums["Projects"] = docsums.Projects.apply(json.dumps)
        if "Samples" in docsums.columns:
            docsums["Samples"] = docsums.Samples.apply(json.dumps)
        docsums.to_csv(f, mode="a", header=not f.tell(), index=False)
