import pandas as pd
import os
import re
import gzip
import tarfile
import io

xls = re.compile("xls")
keep = "|".join(
    ["\." + i + "(.gz)?$" for i in "tab xlsx diff tsv xls csv txt rtf".split(" ")]
)
keep = re.compile(keep)
gse = re.compile("GSE\d+_")
pv_str = "p.{0,4}val"
pv = re.compile(pv_str)
adj = re.compile("adj|fdr|corr")


def find_header(df, n=20):
    head = df.head(n)
    idx = 0
    for index, row in head.iterrows():
        if all([isinstance(i, str) for i in row]):
            idx = index
            break
    return idx


def read_csv(input, tar=None):
    input_name = input
    csv = input
    if isinstance(input, (tarfile.TarInfo)):
        input_name = gse.search(tar.name).group(0) + input.name.replace("/", "_")
        csv = io.StringIO(tar.extractfile(input).read().decode("unicode_escape"))
        input = tar.extractfile(input)
    r = pd.read_csv(
        csv, sep=None, engine="python", iterator=True, skiprows=20, nrows=1000
    )
    sep = r._engine.data.dialect.delimiter
    df = pd.read_csv(input, sep=sep, comment="#", encoding="unicode_escape")
    if all(["Unnamed" in i for i in list(df.columns)]):
        idx = find_header(df)
        if idx > 0:
            df = pd.read_csv(
                input, sep=sep, comment="#", skiprows=idx, encoding="unicode_escape"
            )
    return {os.path.basename(input_name): df}


def read_excel(input, tar=None):
    tabs = {}
    input_name = input
    if isinstance(input, (tarfile.TarInfo)):
        input_name = gse.search(tar.name).group(0) + input.name.replace("/", "_")
        input = tar.extractfile(input)
    if input_name.endswith(".gz"):
        with gzip.open(input) as gz:
            wb = pd.ExcelFile(gz)
    else:
        wb = pd.ExcelFile(input)
    sheets = wb.sheet_names
    for sheet in sheets:
        df = wb.parse(sheet, comment="#")
        if not df.empty:
            tabs.update({os.path.basename(input_name) + "-sheet-" + sheet: df})
    return tabs


def import_tar(path):
    out = {}
    with tarfile.open(path, "r:*") as tar:
        for member in tar:
            if member.isfile():
                if keep.search(member.name):
                    if not member.name.startswith("."):
                        try:
                            if xls.search(member.name):
                                out.update(read_excel(member, tar))
                            else:
                                out.update(read_csv(member, tar))
                        except Exception as e:
                            print("Error: ", e)
    return out


def import_flat(path):
    out = {}
    try:
        if xls.search(path):
            out.update(read_excel(path))
        else:
            out.update(read_csv(path))
    except Exception as e:
        print("Error: ", e)
    return out


def filter_pvalue_tables(input, pv=None, adj=None):
    return {
        k: v
        for k, v in input.items()
        if any([bool(pv.search(i) and not adj.search(i)) for i in v.columns])
    }


def summarise_pvalue_tables(df):
    expr = ["basemean", "value", "logcpm", "rpkm"]
    df.columns = map(str.lower, df.columns)
    pvalues = df.filter(regex=pv_str).copy()
    pvalues.columns = ["pvalue"]
    for e in expr:
        label = e
        if e is "value":
            e = "^value_\d"
        frames = df.filter(regex=e, axis=1).mean(axis=1, skipna=True)
        pvalues.loc[:, label] = frames
    return pvalues.dropna(how="all")


dir = "output/suppl/"
suppfiles = os.listdir(dir)
# path = "/Users/taavi/Downloads/GSE0_test.tar.gz"

out = {}
for input in suppfiles:
    print("Working on: ", input)
    path = dir + input
    if path.endswith("tar.gz"):
        dfs = filter_pvalue_tables(import_tar(path), pv, adj)
        dfsums = {k: summarise_pvalue_tables(v) for k, v in dfs.items()}
        out.update(dfsums)
    else:
        dfs = filter_pvalue_tables(import_flat(path), pv, adj)
        dfsums = {k: summarise_pvalue_tables(v) for k, v in dfs.items()}
        out.update(dfsums)

for k, v in out.items():
    print("Table: ", k)
    print(v)
