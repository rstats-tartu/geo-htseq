import pandas as pd
from urllib.parse import urlparse
import ftplib

input = snakemake.input[0]
output = snakemake.output[0]

ds = pd.read_csv(input, index_col="Accession")
paths = [urlparse(link).path for link in ds["FTPLink"]]

with ftplib.FTP("ftp.ncbi.nlm.nih.gov") as ftp:

    try:
        ftp.login("anonymous", "tapa741@gmail.com")
        with open(output, "a+") as f:
            for path in paths:
                print("Working on: ", path.split("/")[4])
                for dir in ["matrix", "suppl"]:
                    try:
                        ftp.cwd(path)
                        for item in ftp.nlst(dir):
                            f.write("{}\n".format(item))
                    except ftplib.error_temp as e:
                        print("FTP error:", e)
                        pass
    except ftplib.all_errors as e:
        print("FTP error:", e)
