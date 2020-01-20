import ftplib
import os
import re

p = re.compile("GSE\\d+")

input = snakemake.input[0]
output = snakemake.output[0]

with ftplib.FTP("ftp.ncbi.nlm.nih.gov") as ftp:

    try:
        ftp.login("anonymous", "taavi.pall@ut.ee")
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
                    with open(path, "wb") as file:
                        ftp.retrbinary("RETR " + filename, file.write, 1024)
    except ftplib.all_errors as e:
        print("FTP error:", e)
