import pandas as pd
import requests

api_key = snakemake.params.get("api_key", None)
if api_key is None:
    raise Exception("API key is required for request")

# Import document summaries
pubs = pd.read_csv(snakemake.input[0], sep=",")

# Get list of dois
dois = pubs["DOI"].tolist()

# Query parameters
url = "https://api.elsevier.com/content/search/scopus"
headers = {"Accept": "application/json", "Encoding": "UTF-8"}


# Run query
with open(snakemake.output[0], "a") as f:
    citations = pd.DataFrame(columns=["DOI", "citations"])
    for doi in dois:
        params = {"query": "doi({})".format(doi), "apiKey": api_key}
        r = requests.get(url, headers=headers, params=params)
        entry = r.json()["search-results"]["entry"][0]
        try:
            entry = {"DOI": doi, "citations": entry["citedby-count"]}
        except KeyError:
            entry = {"DOI": doi, "citations": "NA"}
        citations = citations.append(
            {k: v for k, v in entry.items()}, ignore_index=True
        )
    citations.to_csv(f, mode="a", header=not f.tell(), index=False)
