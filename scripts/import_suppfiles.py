import pandas as pd
import os
import re
import gzip
import tarfile

pv = re.compile("p.*val")
sp = re.compile("\s+")
xls = re.compile("xls")


def import_tabs(path):
    tabs = {}
    if xls.search(path):
        try:
            if path.endswith(".gz"):
                with gzip.open(path) as gz:
                    wb = pd.ExcelFile(gz)
            else:
                wb = pd.ExcelFile(path)
            sheets = wb.sheet_names
            for sheet in sheets:
                df = wb.parse(sheet, comment="#")
                if not df.empty:
                    tabs.update({os.path.basename(path) + "-sheet-" + sheet: df})
        except Exception as e:
            print("Error: ", e)
            pass
    else:
        try:
            r = pd.read_csv(
                path, sep=None, engine="python", iterator=True, skiprows=20, nrows=1000
            )
            sep = r._engine.data.dialect.delimiter
            df = pd.read_csv(path, sep=sep, comment="#")
            if all(["Unnamed" in i for i in list(df.columns)]):
                idx = find_header(df)
                if idx > 0:
                    df = pd.read_csv(path, sep=sep, comment="#", skiprows=idx)
            tabs.update({input: df})
            r.close()
        except Exception as e:
            print("Error: ", e)
            pass
    return tabs


def find_header(df, n=20):
    head = df.head(n)
    idx = 0
    for index, row in head.iterrows():
        if all([isinstance(i, str) for i in row]):
            idx = index
            break
    return idx


dir = "output/suppl/"
suppfiles = os.listdir(dir)

for input in suppfiles:
    print("Working on: ", input)
    path = dir + input
    tabs = import_tabs(path)
    for k, v in tabs.items():
        print("Table: ", k)
        print(v.head())
