import argparse
import pandas as pd
from urllib.parse import urlparse
import ftplib


def chunks(lst, size):
    """Yield successive n-sized chunks from lst."""
    n = len(lst) // int(size)
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def download_suppfilenames(input, output, email, dirs=["matrix", "suppl"], size=200):

    ds = pd.read_csv(input, index_col="Accession")
    ftplinks = chunks([urlparse(link).path for link in ds["FTPLink"]], size)
    suppfilenames = []
    for paths in ftplinks:
        with ftplib.FTP("ftp.ncbi.nlm.nih.gov") as ftp:
            try:
                ftp.login("anonymous", email)
                for path in paths:
                    print("Working on: ", path.split("/")[4])
                    if not isinstance(dirs, list):
                        dirs = [dirs]
                    for dir in dirs:
                        try:
                            ftp.cwd(path)
                            suppfilenames = suppfilenames + ftp.nlst(dir)
                        except ftplib.error_temp as e:
                            print("FTP error:", e)
                            pass
            except ftplib.all_errors as e:
                print("FTP error:", e)
    with open(output, "w") as f:
        for file in suppfilenames:
            f.write(f"{file}\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        metavar="FILE",
        help="a NCBI GEO document summary file with list of GEO accessions in column 'Accession' and FTP links in column 'FTPLink'",
        required=True,
    )
    parser.add_argument("--output", metavar="FILE", help="output file", required=True)
    parser.add_argument(
        "--email",
        metavar="EMAIL",
        help="email address for anonymous FTP",
        required=True,
    )
    parser.add_argument(
        "--dirs",
        metavar="DIR",
        help="FTP dirs to list",
    )
    parser.add_argument(
        "--size",
        metavar="INT",
        type=int,
        help="batch size",
    )
    args = parser.parse_args()

    download_suppfilenames(**args.__dict__)
