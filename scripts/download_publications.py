import pandas as pd
from Bio import Entrez

# Helper functions
def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i : i + n]


def extract(d, keys):
    return dict((k, d[k]) for k in keys if k in d)


COLUMNS = [
    "Id",
    "PubDate",
    "EPubDate",
    "Source",
    "AuthorList",
    "LastAuthor",
    "Title",
    "Volume",
    "Issue",
    "Pages",
    "LangList",
    "NlmUniqueID",
    "ISSN",
    "ESSN",
    "PubTypeList",
    "RecordStatus",
    "PubStatus",
    "ArticleIds",
    "DOI",
    "History",
    "References",
    "HasAbstract",
    "PmcRefCount",
    "FullJournalName",
    "ELocationID",
    "SO",
]

# Import document summaries
ds = pd.read_csv(snakemake.input[0], sep=",")

# Get list of pubmed ids
pids = ds["PubMedIds"].tolist()
pids = [i for i in pids if str(i) != "nan"]
pids = [i.split(";") for i in pids]
pids = [i for sublist in pids for i in sublist]
pids = list(set(pids))

# Setup query params
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

db = "pubmed"
batch_size = snakemake.params.get("batch_size", 1)


chunked_ids = chunks(pids, batch_size)

# Fetch and parse summaries
with open(snakemake.output[0], "a") as f:
    for chunk in chunked_ids:
        summary = Entrez.esummary(db=db, id=",".join(chunk), retmode="xml")
        records = Entrez.parse(summary)
        docsums = pd.DataFrame(columns=COLUMNS)
        for record in records:
            record = extract(record, COLUMNS)
            docsums = docsums.append(
                {k: v for k, v in record.items()}, ignore_index=True
            )
        for col in docsums.columns:
            if "List" in col:
                docsums[col] = docsums[col].apply(
                    lambda x: ";".join(pd.Series(x, dtype="str"))
                )
        docsums.to_csv(f, mode="a", header=not f.tell(), index=False)
