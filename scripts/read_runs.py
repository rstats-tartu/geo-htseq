import os
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
from tqdm import tqdm
from requests import Request, Session

def empty_dataframe(acc, err, var):
    df = pd.DataFrame(columns=var)
    df.loc[0,"geo_accession"] = acc
    df.loc[0,"err"] = err
    return df

def spotify(acc, email, **kwargs):

    Entrez.email = email
    if "api_key" in kwargs:
        Entrez.api_key = kwargs.pop("api_key")
    Entrez.max_tries = kwargs.pop("max_tries", 3)
    Entrez.sleep_between_tries = kwargs.pop("sleep_between_tries", 15)
    retmax = kwargs.pop("retmax", 20)
    keep = "geo_accession,study_accession,run_accession,sample_name,sample_title,spots,bases,tax_id,organism,LIBRARY_STRATEGY,LIBRARY_SOURCE,LIBRARY_SELECTION,LIBRARY_LAYOUT,PLATFORM,MODEL,err".split(",")

    term = '({}[All Fields]) AND "gsm"[Filter]'.format(acc)
    handle = Entrez.esearch(db="gds", term=term, retmode="text", retmax = retmax)
    records=Entrez.read(handle)
    handle.close()
    gsm=records["IdList"]
    handle = Entrez.elink(dbfrom="gds", id=",".join(gsm), linkname="gds_sra", retmax=retmax)
    records = Entrez.read(handle)
    handle.close()

    if len(records[0]["LinkSetDb"]) == 0:
        handle = Entrez.esummary(db="gds", id=",".join(gsm), retmode="text", retmax=retmax)
        records = Entrez.read(handle)
        handle.close()
        gsm_acc = [i["Accession"] for i in records]
        try:
            srx = [i['ExtRelations'][0]['TargetObject'] for i in records]
        except IndexError as e:
            err = "SRA Experiment is not public: {}".format(e)
            return empty_dataframe(acc, err, keep)
        handle = Entrez.esearch(db="sra", term=" OR ".join(srx), retmode="text", retmax=retmax)
        records = Entrez.read(handle)
        handle.close()
        uids = records['IdList']
    else:
        uids = [link["Id"] for link in records[0]["LinkSetDb"][0]["Link"]]

    url_endpoint = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {"db": "sra", "id": ",".join(uids), "rettype": "full", "retmode": "xml", "email": Entrez.email}
    if Entrez.api_key:
        params.update({"api_key": Entrez.api_key})
    method = "GET"
    if len(uids) > 200:
        method = "POST"
    s = Session()
    req = Request(method, url_endpoint, params=params)
    prepped = s.prepare_request(req)
    resp = s.send(prepped)
    tree = ET.ElementTree(ET.fromstring(resp.text))
    root = tree.getroot()

    
    
    if root.findall(".//ERROR"):
        err = root.findall(".//ERROR")[0].text
        return empty_dataframe(acc, err, keep)
    
    platform = [{"PLATFORM": i.tag, "MODEL": i.find("INSTRUMENT_MODEL").text} for i in root.findall(".//EXPERIMENT_PACKAGE/EXPERIMENT/PLATFORM/")]
    design = []
    for e in root.iter("LIBRARY_DESCRIPTOR"):
        libdesc = {}
        libdesc.update({"LIBRARY_STRATEGY": e.find("LIBRARY_STRATEGY").text})
        libdesc.update({"LIBRARY_SOURCE": e.find("LIBRARY_SOURCE").text})
        libdesc.update({"LIBRARY_SELECTION": e.find("LIBRARY_SELECTION").text})
        libdesc.update({"LIBRARY_LAYOUT": e.find("LIBRARY_LAYOUT")[0].tag})
        design.append(libdesc)

    spots = [i.attrib for i in root.findall(".//EXPERIMENT_PACKAGE/RUN_SET/RUN/Pool/Member")]
    study = [{"study_accession": i.attrib.pop("accession")} for i in root.findall(".//EXPERIMENT_PACKAGE/EXPERIMENT/STUDY_REF")]
    df = pd.DataFrame([{**a, **b, **c, **d} for a,b,c,d in zip(study, platform, spots, design)])
    df.rename(columns = {"accession": "run_accession"}, inplace=True)
    df.insert(0, "geo_accession", acc)
    df.insert(0, "err", np.nan)
    return df[keep]


if __name__ == "__main__":

    # Parse parameters
    email = snakemake.params.get("email", None)
    if email is None:
        raise ValueError("An email must be provided")

    api_key = snakemake.params.get("api_key", None)
    if api_key is None:
        print(
            "Personal API key from NCBI. If not set, only 3 queries per second are allowed. 10 queries per seconds otherwise with a valid API key."
        )

    max_tries = snakemake.params.get("max_tries", 3)
    sleep_between_tries = snakemake.params.get("sleep_between_tries", 15)
    retmax = snakemake.params.get("retmax", 100000)

    # Parse input to chunks
    accessions = pd.read_csv(snakemake.input[0], sep=",")["Accession"]
    accessions = accessions.to_list()

    # Fetch and parse summaries
    with open(snakemake.output[0], "a") as output_handle:
        for i, acc in tqdm(enumerate(accessions), total=len(accessions)):
            print(acc)
            spots = spotify(acc, email=email, api_key=api_key, max_tries=max_tries, sleep_between_tries=sleep_between_tries, retmax=retmax)
            spots.to_csv(
                output_handle,
                sep=",",
                mode="a",
                header=not output_handle.tell(),
                index=False,
            )
