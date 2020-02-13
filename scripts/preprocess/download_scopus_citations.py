import pandas as pd
import requests
import itertools
import json


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


api_key = snakemake.params.get("api_key", None)
if api_key is None:
    raise Exception("API key is required for request")

# Import document summaries
df = pd.read_csv(snakemake.input[0], sep=",")

# Get list of dois
ids = df["PubMedIds"].dropna().tolist()
ids = [i.split(";") for i in ids]
ids = list(set(list(itertools.chain.from_iterable(ids))))
chunked = chunks(ids, 25)

# Query parameters
url = "https://api.elsevier.com/content/search/scopus"
headers = {"Accept": "application/json", "Encoding": "UTF-8"}


# Run query
with open(snakemake.output[0], "a") as f:
    citations = {}
    for chunk in chunked:
        pmids_query = " OR ".join(["PMID({})".format(i) for i in chunk])
        params = {
            "query": pmids_query,
            "apiKey": api_key,
            "count": 25,
            "field": "pubmed-id,citedby-count",
        }
        r = requests.get(url, headers=headers, params=params)
        r.raise_for_status()
        entry = r.json()["search-results"]["entry"]
        citedby_count = {i["pubmed-id"]: int(i["citedby-count"]) for i in entry}
        citations.update(citedby_count)
    df = (
        pd.DataFrame.from_dict(citations, orient="index", columns=["citations"])
        .reset_index()
        .rename(columns={"index": "PubMedIds"})
    )
    df.to_csv(f, mode="a", header=not f.tell(), index=False)
