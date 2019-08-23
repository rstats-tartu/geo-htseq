from Bio import Entrez
import pandas as pd
import json

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]

Entrez.email = snakemake.params["email"]
Entrez.api_key = snakemake.params["api_key"]


query = 'expression profiling by high throughput sequencing[DataSet Type]'
db = "gds"

search = Entrez.esearch(db = db, term = query, RetMax = 30000)
ids = Entrez.read(search)
chunked_ids = chunks(ids["IdList"], 200)

columns = ["Id", "Accession", "GDS", "title", "summary", "GPL", "GSE", "taxon", "entryType", "gdsType", 
"ptechType", "valType", "SSInfo", "subsetInfo", "PDAT", "suppFile", "Samples", "n_samples","SeriesTitle","PlatformTitle","PlatformTaxa",
"SamplesTaxa", "PubMedIds", "Projects", "FTPLink"]

with open(snakemake.output[0], "a") as f:
  for chunk in chunked_ids:
    summary = Entrez.esummary(db = db, id = ",".join(chunk), retmode = "xml")
    records = Entrez.parse(summary)
    docsums = pd.DataFrame(columns = columns)
    for record in records:
      docsums = docsums.append({"Id": record["Id"], 
                        "Accession" : record["Accession"], 
                               "GDS": record["GDS"], 
                             "title": record["title"], 
                           "summary": record["summary"], 
                               "GPL": record["GPL"], 
                               "GSE": record["GSE"],
                             "taxon": record["taxon"],
                         "entryType": record["entryType"],
                           "gdsType": record["gdsType"],
                         "ptechType": record["ptechType"],
                           "valType": record["valType"],
                            "SSInfo": record["SSInfo"],
                        "subsetInfo": record["subsetInfo"],
                              "PDAT": record["PDAT"],
                          "suppFile": record["suppFile"],
                           "Samples": record["Samples"],
                         "n_samples": record["n_samples"],
                       "SeriesTitle": record["SeriesTitle"],
                     "PlatformTitle": record["PlatformTitle"],
                      "PlatformTaxa": record["PlatformTaxa"],
                       "SamplesTaxa": record["SamplesTaxa"],
                         "PubMedIds": record["PubMedIds"],
                          "Projects": record["Projects"],
                           "FTPLink": record["FTPLink"]}, 
                                ignore_index=True)
    docsums["PubMedIds"] = docsums.PubMedIds.apply(lambda x: ";".join(pd.Series(x, dtype = "str")))
    docsums["Projects"] = docsums.Projects.apply(json.dumps)
    docsums["Samples"] = docsums.Samples.apply(json.dumps)
    docsums.to_csv(f, mode = "a", header = not f.tell(), index = False)
