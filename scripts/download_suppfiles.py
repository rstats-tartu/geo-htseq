import ftplib
import os
import re
import argparse
from typing import Type


def chunks(lst, batchsize):
    """Yield successive n-sized chunks from lst."""
    n = max([len(lst) // int(batchsize), 1])
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


class download_suppfiles:
    def __init__(self, input=None, file=None, email=None, maxfilesize=1e9, batchsize=200, dir="."):
        self.input = input
        self.file = file
        self.email = email
        self.maxfilesize=maxfilesize
        self.batchsize = batchsize
        self.dir = dir
        self.p = re.compile("GSE\\d+")

        if self.email is None:
            raise TypeError("Please supply email address.")

    def from_list(self):
        if self.input is None:
            raise TypeError(
                "Please supply file with list of supplementary files to be downloaded."
            )
        with open(self.input, "r") as i:
            filenames = i.readlines()
        filenames = chunks(filenames, self.batchsize)
        for chunk in filenames:
            with ftplib.FTP("ftp.ncbi.nlm.nih.gov") as ftp:
                try:
                    ftp.login("anonymous", self.email)
                    for line in chunk:
                        try:
                            self.retr_ncbi(ftp, line)
                        except ftplib.all_errors as e:
                            print("FTP error:", e)
                            continue
                except ftplib.all_errors as e:
                    print("FTP error:", e)
                    continue


    def from_filename(self):
        if self.file is None:
            raise TypeError("Please supply supplementary file name to be downloaded.")
        with ftplib.FTP("ftp.ncbi.nlm.nih.gov") as ftp:
            try:
                ftp.login("anonymous", self.email)
                try:
                    line = self.file
                    self.retr_ncbi(ftp, line)
                except ftplib.all_errors as e:
                    print("FTP error:", e)
            except ftplib.all_errors as e:
                print("FTP error:", e)

    def retr_ncbi(self, ftp, line):
        path = os.path.join(self.dir, line.rstrip())
        if os.path.exists(path):
            restarg = {'rest': str(os.path.getsize(path))}
        else:
            restarg = {}
        filename = os.path.basename(path)
        id = self.p.search(filename).group(0)
        ftpdir = (
                        "/geo/series/"
                        + id[0:-3]
                        + "nnn/"
                        + id
                        + "/"
                        + os.path.dirname(line.rstrip())
                    )
        ftp.cwd(ftpdir)
        if ftp.size(filename) < self.maxfilesize:
            with open(path, "wb") as file:
                ftp.retrbinary("RETR " + filename, file.write, 1024, **restarg)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--file", metavar="FILENAME", help="supplementary file name to be downloaded"
    )
    group.add_argument(
        "--input", metavar="FILE", help="input file with list of files to be downloaded"
    )
    parser.add_argument(
        "--email",
        metavar="EMAIL",
        help="email address for anonymous FTP",
        required=True,
    )
    parser.add_argument(
        "--maxfilesize", metavar="INT", type=int, help="maximum supplementary file size to be downloaded", default=1e9
    )
    parser.add_argument(
        "--batchsize", metavar="INT", type=int, help="batch size", default=200
    )
    parser.add_argument("--dir", metavar="DIR", help="target directory", default=".")
    args = parser.parse_args()

    ds = download_suppfiles(**args.__dict__)
    if args.file:
        ds.from_filename()
    elif args.input:
        ds.from_list()
