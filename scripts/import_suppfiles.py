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


def summarise_pvalue_tables(df, var=["basemean", "value", "logcpm", "rpkm"]):
    df.columns = map(str.lower, df.columns)
    pvalues = df.filter(regex=pv_str).copy()
    pvalues.columns = ["pvalue"]
    for v in var:
        label = v
        if v is "value":
            v = "^value_\d"
            label = "fpkm"
        frames = df.filter(regex=v, axis=1)
        if not frames.empty:
            frames = frames.mean(axis=1, skipna=True)
            pvalues.loc[:, label] = frames
    return pvalues.dropna(subset=["pvalue"])


# https://stackoverflow.com/a/32681075/1657871
def rle(inarray):
    """ run length encoding. Partial credit to R rle function. 
            Multi datatype arrays catered for including non Numpy
            returns: tuple (runlengths, startpositions, values) """
    ia = np.asarray(inarray)  # force numpy
    n = len(ia)
    if n == 0:
        return (None, None, None)
    else:
        y = np.array(ia[1:] != ia[:-1])  # pairwise unequal (string safe)
        i = np.append(np.where(y), n - 1)  # must include last element posi
        z = np.diff(np.append(-1, i))  # run lengths
        p = np.cumsum(np.append(0, z))[:-1]  # positions
        return (z, p, ia[i])


def get_hist_class(counts, breaks, fdr):
    qc = binom.ppf(1 - 1 / breaks * fdr, sum(counts), 1 / breaks)
    counts_over_qc = counts > qc
    for i in counts_over_qc:
        if all(~counts_over_qc):
            Class = "uniform"
        elif not any(counts_over_qc[rle(counts_over_qc)[0][0] + 2 :]):
            Class = "anti-conservative"
        elif not any(np.flip(counts_over_qc)[rle(np.flip(counts_over_qc))[0][0] + 2 :]):
            Class = "conservative"
        else:
            Class = "other"
    return Class



def summarise_pvalues(
    df,
    breaks=30,
    fdr=0.05,
    var={"basemean": 10, "fpkm": 0.5, "logcpm": -0.3, "rpkm": 0.5},
):
    bins = np.linspace(0, 1, breaks)
    center = (bins[:-1] + bins[1:]) / 2
    # Filter pvalues
    pf = pd.DataFrame()
    for k, v in var.items():
        if k in df.columns:
            pf = df.loc[df[k] >= v, ["pvalue"]]
            break
    # Make histogram
    hists = [np.histogram(i["pvalue"], bins=bins) for i in [df, pf] if not i.empty]
    counts = [counts for (counts, bins) in hists]
    # Assign class to histograms
    Class = [get_hist_class(i, breaks, fdr) for i in counts]
    # Calculate pi0


def write_to_csv(input, outpath):
    for k, v in input.items():
        v.to_csv(outpath + k + ".csv", sep=",", index=False)


dir = "output/suppl/"
suppfiles = os.listdir(dir)
# path = "/Users/taavi/Downloads/GSE0_test.tar.gz"


for input in suppfiles:
    print("Working on: ", input)
    path = dir + input
    if path.endswith("tar.gz"):
        frames = import_tar(path)
    else:
        frames = import_flat(path)
    frames = filter_pvalue_tables(frames, pv, adj)
    frames = {k: summarise_pvalue_tables(v) for k, v in frames.items()}
    write_to_csv(frames, "output/tmp/imported/")

# for k, v in out.items():
#     print("Table: ", k)
#     print(v)
