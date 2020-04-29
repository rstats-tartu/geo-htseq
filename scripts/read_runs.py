import os
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd
import time


def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i : i + n]


def beetroot(brokenxml):
    xml = "<?xml version='1.0' encoding='UTF-8' ?><root>" + brokenxml + "</root>"
    tree = ET.ElementTree(ET.fromstring(xml))
    return tree.getroot()


def spotify(accessions, retmax=20):
    # Get ids for accessions
    handle = Entrez.esearch(
        db="bioproject", term=" OR ".join(accessions), retmode="text", retmax=retmax
    )
    records = Entrez.read(handle)

    # Get project acc for ids
    handle = Entrez.esummary(
        db="bioproject", id=",".join(records["IdList"]), retmode="text", retmax=retmax
    )
    records = Entrez.read(handle)
    prjna = [
        ds["Project_Acc"] for ds in records["DocumentSummarySet"]["DocumentSummary"]
    ]

    # Get run ids
    handle = Entrez.esearch(
        db="sra", term=" OR ".join(prjna), retmode="text", retmax=retmax
    )
    records = Entrez.read(handle)

    # Get run metadata
    handle = Entrez.efetch(
        db="sra", id=",".join(records["IdList"]), rettype="docsum", retmax=retmax
    )
    records = Entrez.read(handle)

    # Parse run metadata
    df = pd.DataFrame()
    for record in records:
        runs = beetroot(record["Runs"])
        exps = beetroot(record["ExpXml"])
        spots = runs[0].attrib
        spots.update({"name": exps[2].attrib["name"]})
        spots.update({"platform": list(exps[6].attrib.values())[0]})
        df = df.append(spots, True)
    return df[["acc", "name", "total_bases", "total_spots", "platform"]]


if __name__ == "__main__":

    # Parse parameters
    Entrez.email = snakemake.params.get("email", None)
    if Entrez.email is None:
        raise ValueError("An email must be provided")

    Entrez.api_key = snakemake.params.get("api_key", None)
    if Entrez.api_key is None:
        print(
            "Personal API key from NCBI. If not set, only 3 queries per second are allowed. 10 queries per seconds otherwise with a valid API key."
        )

    retmax = snakemake.params.get("retmax", 100000)
    batch_size = snakemake.params.get("batch_size", 1000)
    sleep = snakemake.params.get("sleep", 0)

    # Parse input to chunks
    acc = pd.read_csv(snakemake.input[0], sep=",")["Accession"]
    chunked_acc = chunks(acc.to_list(), batch_size)

    # Fetch and parse summaries
    with open(snakemake.output[0], "a") as output_handle:
        for chunk in chunked_acc:
            time.sleep(sleep)
            spots = spotify(chunked_acc, retmax=retmax)
            spots.to_csv(
                output_handle, mode="a", header=not output_handle.tell(), index=False
            )
