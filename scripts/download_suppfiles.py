import ftplib
import os
import re


def download_suppfiles(input, output, email):
    p = re.compile("GSE\\d+")
    with ftplib.FTP("ftp.ncbi.nlm.nih.gov") as ftp:

        try:
            ftp.login("anonymous", email)
            with open(input) as i:
                for line in i:
                    path = "output/" + line.rstrip()
                    if not os.path.isfile(path):
                        filename = os.path.basename(path)
                        id = p.search(filename).group(0)
                        dir = (
                            "/geo/series/"
                            + id[0:-3]
                            + "nnn/"
                            + id
                            + "/"
                            + os.path.dirname(line.rstrip())
                        )
                        ftp.cwd(dir)
                        if ftp.size(filename) < 1e9:
                            with open(path, "wb") as file:
                                ftp.retrbinary("RETR " + filename, file.write, 1024)
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

    download_suppfiles(args.list, args.out, args.email)
