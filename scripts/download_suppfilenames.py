import argparse
import pandas as pd
from urllib.parse import urlparse
import ftplib


def download_suppfilenames(input, output, email):
    ds = pd.read_csv(input, index_col="Accession")
    paths = [urlparse(link).path for link in ds["FTPLink"]]
    with ftplib.FTP("ftp.ncbi.nlm.nih.gov") as ftp:

        try:
            ftp.login("anonymous", email)
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


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--list", metavar="FILE", help="input file with list of files to be downloaded"
    )
    parser.add_argument("--out", metavar="FILE", help="output file")
    parser.add_argument(
        "--email", metavar="EMAIL", help="email address for anonymous FTP"
    )
    args = parser.parse_args()

    download_suppfilenames(args.list, args.out, args.email)
