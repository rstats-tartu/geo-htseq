import os
import sys
import re
import gzip
import tarfile
import io
import scipy
import collections
import argparse
import pandas as pd
import numpy as np
from scipy.stats import binom
from pandas.api.types import is_string_dtype
from pathlib import Path
import numbers


xls = re.compile("xls")
keep = "|".join(
    ["\." + i + "(.gz)?$" for i in "tab xlsx diff tsv xls csv txt rtf".split(" ")]
)
keep = re.compile(keep)
gse = re.compile("GSE\d+_")
pv_str = "p[^a-zA-Z]{0,4}val"
pv = re.compile(pv_str)
adj = re.compile("adj|fdr|corr|thresh")
ws = re.compile(" ")
mtabs = re.compile("\w+\t{2,}\w+")
tab = re.compile("\t")
fields = ["Type", "Class", "Conversion", "pi0", "FDR_pval", "hist", "note"]
PValSum = collections.namedtuple("PValSum", fields, defaults=[np.nan] * 7)


def raw_pvalues(i):
    return bool(pv.search(i.lower()) and not adj.search(i.lower()))


def find_header(df, n=20):
    head = df.head(n)
    idx = 0
    for col in head:
        s = head[col]
        match = s.str.contains(pv_str, na=False)
        if any(match):
            idx = s.index[match].tolist()[0] + 1
            break
    if idx == 0:
        for index, row in head.iterrows():
            if all([isinstance(i, str) for i in row if i is not np.nan]):
                idx = index + 1
                break
    return idx


def csv_helper(input, input_name, csv, verbose=0):
    # Get comments and set rows to skip
    r = pd.read_csv(csv, sep=None, engine="python", iterator=True, nrows=1000)
    comment = None
    sep = r._engine.data.dialect.delimiter
    columns = r._engine.columns
    if isinstance(input, (tarfile.ExFileObject)):
        with csv as h:
            first_line = h.readline()
    elif re.search("gz$", input):
        with gzip.open(input, "rb") as h:
            first_line = h.readline().decode("utf-8").rstrip()
    else:
        with open(input, "r") as h:
            first_line = h.readline().rstrip()
    more_tabs_than_sep = len(tab.findall(first_line)) > len(re.findall(sep, first_line))
    if re.search("^#", first_line) or more_tabs_than_sep:
        comment = "#"
        # Get delimiter
        r = pd.read_csv(
            csv, sep=None, engine="python", iterator=True, skiprows=20, nrows=1000
        )
        sep = r._engine.data.dialect.delimiter
        columns = r._engine.columns
    if ws.search(sep):
        sep = "\s+"
    if mtabs.search(first_line):
        sep = "\t+"
    # Import file
    df = pd.read_csv(input, sep=sep, comment=comment, encoding="unicode_escape")
    # Check and fix column names
    # Case of extra level of delimiters in column names
    if len(df.columns) > len(columns):
        df = pd.read_csv(
            input,
            header=None,
            skiprows=[0],
            sep=sep,
            comment=comment,
            encoding="unicode_escape",
        ).drop([0])
        df.columns = columns
    unnamed = ["Unnamed" in i for i in df.columns]
    # Case of empty rows before header
    if all(unnamed):
        idx = find_header(df)
        if idx > 0:
            df = pd.read_csv(
                input, sep=sep, comment=comment, skiprows=idx, encoding="unicode_escape"
            )
    # Case of anonymous row names
    if unnamed[-1] & sum(unnamed) == 1:
        if any([pv.search(i) for i in df.columns]):
            df.columns = [df.columns[-1]] + list(df.columns[:-1])
    if verbose > 1:
        print("df after import:\n", df)
    return {os.path.basename(input_name): df}


def excel_helper(input, input_name, verbose=0):
    tabs = {}
    if input_name.endswith(".gz"):
        with gzip.open(input) as gz:
            wb = pd.ExcelFile(gz)
    else:
        wb = pd.ExcelFile(input)
    sheets = wb.sheet_names
    sheets = [i for i in sheets if "README" not in i]
    for sheet in sheets:
        df = wb.parse(sheet, comment="#")
        if verbose > 1:
            print("df after import:\n", df)
        if not df.empty:
            pu = sum(["Unnamed" in i for i in list(df.columns)]) / len(df.columns)
            if pu >= 2 / 3:
                idx = find_header(df)
                if idx > 0:
                    df = wb.parse(sheet, skiprows=idx)
            tabs.update({os.path.basename(input_name) + "-sheet-" + sheet: df})
    return tabs


def read_csv(input, tar=None):
    if isinstance(input, (tarfile.TarInfo)):
        input_name = os.path.basename(input.name)
        with tar.extractfile(input) as h:
            csv = io.StringIO(h.read().decode("unicode_escape"))
        with tar.extractfile(input) as h:
            out = csv_helper(h, input_name, csv)
    else:
        input_name = input
        csv = input
        out = csv_helper(input, input_name, csv)
    return out


def read_excel(input, tar=None):
    if isinstance(input, (tarfile.TarInfo)):
        input_name = os.path.basename(input.name)
        with tar.extractfile(input) as h:
            out = excel_helper(h, input_name)
    else:
        input_name = input
        out = excel_helper(input, input_name)
    return out


def import_flat(input, tar=None):
    out = {}
    try:
        if xls.search(input.name if tar else input):
            out.update(read_excel(input, tar=tar))
        else:
            d = read_csv(input, tar=tar)
            is_empty = [v.empty for v in d.values()][0]
            if is_empty:
                raise Exception("empty table")
            else:
                out.update(d)
    except Exception as e:
        if tar:
            key = os.path.basename(input.name)
        else:
            key = os.path.basename(input)
        out.update(note(key, e))
    return out


def import_tar(input):
    out = {}
    with tarfile.open(input, "r:*") as tar:
        for member in tar:
            if member.isfile():
                if keep.search(member.name):
                    out.update(import_flat(member, tar))
    return out


def filter_pvalue_tables(input, pv=None, adj=None):
    return {k: v for k, v in input.items() if any([raw_pvalues(i) for i in v.columns if not isinstance(i, numbers.Number)])}


def fix_column_dtype(df):
    for col in df.columns:
        s = df[col]
        if is_string_dtype(s):
            if "," in s[:5].astype(str).str.cat(sep=" "):
                df[col] = s.apply(lambda x: str(x).replace(",", "."))
            df[col] = pd.to_numeric(s, errors="coerce")
    return df


def summarise_pvalue_tables(
    df, var=["basemean", "value", "fpkm", "logcpm", "rpkm", "aveexpr"]
):
    df = df.filter(regex='^\D')
    df.columns = map(str.lower, df.columns)
    pval_cols = [i for i in df.columns if raw_pvalues(i)]
    pvalues = df[pval_cols].copy()
    # Check if there is ANOTHER(!!#?) level of ":" delimiters in p value column(s)
    extra_delim = ":"
    split_col = [i for i in pvalues.columns if extra_delim in i]
    if split_col:
        for index, col in enumerate(split_col):
            col_count = len(re.findall(extra_delim, col))
            obs_count = len(re.findall(extra_delim, str(pvalues.iloc[0, index])))
            if obs_count == 0:
                pass
            elif col_count == obs_count:
                new_cols = col.split(extra_delim)
                split_pval_col = [i for i in new_cols if raw_pvalues(i)]
                cols_split = pvalues.iloc[:, index].str.split(extra_delim, expand=True)
                try:
                    cols_split.columns = new_cols
                    pvalues[split_pval_col] = cols_split[split_pval_col]
                    pvalues.drop(col, axis=1, inplace=True)
                except ValueError:
                    pass
        pval_cols = [i for i in pvalues.columns if raw_pvalues(i)]
    pvalues_check = fix_column_dtype(pvalues)
    for v in var:
        label = v
        if v is "value":
            v = "^value_\d"
            label = "fpkm"
        exprs = df.filter(regex=v, axis=1)
        if not exprs.empty:
            exprs_check = fix_column_dtype(exprs)
            exprs_sum = exprs_check.mean(axis=1, skipna=True)
            pvalues_check.loc[:, label] = exprs_sum
            break
    pv_stacked = (
        pvalues_check.melt(id_vars=list(set(pvalues_check.columns) - set(pval_cols)))
        .set_index("variable")
        .rename(columns={"value": "pvalue"})
    )
    return pv_stacked.dropna()


# https://stackoverflow.com/a/32681075/1657871
def rle(inarray):
    """run length encoding. Partial credit to R rle function.
    Multi datatype arrays catered for including non Numpy
    returns: tuple (runlengths, startpositions, values)"""
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


def get_hist_class(counts, fdr):
    bins = len(counts)
    qc = binom.ppf(1 - 1 / bins * fdr, sum(counts), 1 / bins)
    counts_over_qc = counts > qc
    ru = rle(counts_over_qc)
    over_qc = ru[1][ru[2]]
    rufl = rle(np.flip(counts_over_qc))
    over_qc_fl = rufl[1][rufl[2]]
    if all(~counts_over_qc):
        Class = "uniform"
    elif len(over_qc) == 1:
        over_qc_prop = ru[0][ru[2]] / bins
        if over_qc == 0 and over_qc_prop < 1 / 3:
            Class = "anti-conservative"
        elif over_qc_fl == 0:
            Class = "conservative"
        else:
            Class = "other"
    elif len(over_qc) == 2:
        if over_qc[0] == 0 and over_qc_fl[0] == 0:
            Class = "bimodal"
        else:
            Class = "other"
    else:
        Class = "other"
    return Class


# https://gdsctools.readthedocs.io/en/master/_modules/gdsctools/qvalue.html#QValue
def estimate_pi0(
    pv,
    lambdas=None,
    pi0=None,
    df=3,
    method="smoother",
    smooth_log_pi0=False,
    verbose=True,
):
    """Estimate pi0 based on the pvalues"""
    try:
        pv = np.array(pv)
    except:
        pv = pv.copy()
    assert pv.min() >= 0 and pv.max() <= 1, "p-values should be between 0 and 1"
    if lambdas is None:
        epsilon = 1e-8
        lambdas = np.arange(0, 0.9 + 1e-8, 0.05)
    if len(lambdas) > 1 and len(lambdas) < 4:
        raise ValueError(
            """if length of lambda greater than 1, you need at least 4 values"""
        )
    if len(lambdas) >= 1 and (min(lambdas) < 0 or max(lambdas) >= 1):
        raise ValueError("lambdas must be in the range[0, 1[")
    m = float(len(pv))

    pv = pv.ravel()  # flatten array
    if pi0 is not None:
        pass
    elif len(lambdas) == 1:
        pi0 = np.mean(pv >= lambdas[0]) / (1 - lambdas[0])
        pi0 = min(pi0, 1)
    else:
        # evaluate pi0 for different lambdas
        pi0 = [np.mean(pv >= this) / (1 - this) for this in lambdas]
        # in R
        # lambda = seq(0,0.09, 0.1)
        # pi0 = c(1.0000000, 0.9759067, 0.9674164, 0.9622673, 0.9573241,
        #         0.9573241 0.9558824, 0.9573241, 0.9544406, 0.9457901)
        # spi0 = smooth.spline(lambda, pi0, df=3, all.knots=F, spar=0)
        # predict(spi0, x=max(lambda))$y  --> 0.9457946
        # spi0 = smooth.spline(lambda, pi0, df=3, all.knots=F)
        # predict(spi0, x=max(lambda))$y  --> 0.9485383
        # In this function, using pi0 and lambdas, we get 0.9457946
        # this is not too bad, the difference on the v17 data set
        # is about 0.3 %
        if method == "smoother":
            if smooth_log_pi0:
                pi0 = np.log(pi0)
            # In R, the interpolation is done with smooth.spline
            # within qvalue. However this is done with default
            # parameters, and this is different from the Python
            # code. Note, however, that smooth.spline has a parameter
            # called spar. If set to 0, then we would get the same
            # as in scipy. It looks like scipy has no equivalent of
            # the smooth.spline function in R if spar is not 0
            tck = scipy.interpolate.splrep(lambdas, pi0, k=df)
            pi0 = scipy.interpolate.splev(lambdas[-1], tck)
            if smooth_log_pi0:
                pi0 = np.exp(pi0)
            pi0 = min(pi0, 1.0)
        elif method == "lfdr":
            """Estimate proportion of null p-values
            by average local FDR
            Belinda Phipson and Gordon Smyth
            23 May 2012. Last revised 30 July 2012."""
            n = len(pv)
            i = np.array(list(range(1, n + 1)))
            i.sort()
            i = i[::-1]
            pv.sort()
            pv = pv[::-1]
            q = [min(i, 1) for i in n / np.array(i) * np.array(pv)]
            n1 = n + 1
            pi0 = sum(np.array(i) * q) / n / n1 * 2
        elif method == "bootstrap":
            raise NotImplementedError
            """minpi0 = min(pi0)
            mse = rep(0, len(lambdas))
            pi0.boot = rep(0, len(lambdas))
            for i in range(1,100):
                p.boot = sample(p, size = m, replace = TRUE)
                for i in range(0,len(lambdas)):
                    pi0.boot[i] <- mean(p.boot > lambdas[i])/(1 - lambdas[i])
                mse = mse + (pi0.boot - minpi0)^2
            pi0 = min(pi0[mse == min(mse)])
            pi0 = min(pi0, 1)"""
        if pi0 > 1:
            if verbose:
                print("got pi0 > 1 (%.3f), setting it to 1" % pi0)
            pi0 = 1.0
    assert pi0 >= 0 and pi0 <= 1, "pi0 is not between 0 and 1: %f" % pi0
    return pi0


def conversion(x, y):
    classes = pd.DataFrame.from_dict(
        {
            "uniform": ["same good", "improve, effects", "worsen", "worsen", "worsen"],
            "anti-conservative": [
                "effects lost",
                "same good",
                "worsen",
                "worsen",
                "worsen",
            ],
            "conservative": [
                "improvement; no effects",
                "improvement; effects",
                "same bad",
                "no improvement",
                "no improvement",
            ],
            "other": [
                "improvement; no effects",
                "improvement; effects",
                "no improvement",
                "same bad",
                "no improvement",
            ],
            "bimodal": [
                "improvement; no effects",
                "improvement; effects",
                "no improvement",
                "no improvement",
                "same bad",
            ],
        },
        orient="index",
        columns=["uniform", "anti-conservative", "conservative", "other", "bimodal"],
    )
    return classes.loc[x, y]


def summarise_pvalues(
    df,
    bins=30,
    fdr=0.05,
    var={
        "basemean": 10,
        "fpkm": 0.5,
        "logcpm": np.log2(0.5),
        "rpkm": 0.5,
        "aveexpr": np.log2(10),
    },
    pi0_method="lfdr",
    verbose=True,
):
    breaks = np.linspace(0, 1, bins)
    center = (breaks[:-1] + breaks[1:]) / 2
    out = {}
    grouped = df.groupby(level=0)
    for name, group in grouped:
        # Test if pvalues are in 0 to 1 range
        if group.min()["pvalue"] < 0 or group.max()["pvalue"] > 1:
            out.update(note(name, "p-values not in 0 to 1 range"))
            continue
        # Filter pvalues
        pf = pd.DataFrame()
        for k, v in var.items():
            if k in group.columns:
                pf = group.loc[group[k] >= v, ["pvalue"]]
                filt = k
                break
        # Make histogram
        pv_sets = [i for i in [group, pf] if not i.empty]
        hists = [np.histogram(i["pvalue"], bins=breaks) for i in pv_sets]
        counts = [counts.tolist() for (counts, bins) in hists]
        # Test if p-values are truncated
        truncated = rle([i == 0 for i in counts[0]])[1][-1] > 0
        if truncated:
            out.update(note(name, "p-values truncated or right-skewed"))
            continue
        # Assign class to histograms
        Class = [get_hist_class(i, fdr) for i in counts]
        # Conversion
        Type = ["raw"]
        conv = np.nan
        if len(Class) == 2:
            conv = conversion(Class[0], Class[1])
            Type = ["raw", filt]
        # Calculate pi0
        pi0 = []
        for i, c in zip(pv_sets, Class):
            if c in ["uniform", "anti-conservative"]:
                pi0_est = estimate_pi0(i["pvalue"], method=pi0_method, verbose=verbose)
                pi0.append(pi0_est)
            else:
                pi0.append(np.nan)
        # Number of effects < FDR
        fdr_effects = [sum(i["pvalue"] < fdr) for i in pv_sets]
        out.update(
            {
                name: pd.DataFrame(
                    PValSum(Type, Class, conv, pi0, fdr_effects, counts)._asdict()
                )
            }
        )
    return (
        pd.concat(
            [df for df in out.values()], keys=[k for k in out.keys()], names=["Set"]
        )
        .reset_index(level=["Set"])
        .astype(dtype={"FDR_pval": "Int64"})
    )


def note(filename, message):
    return {
        filename: pd.DataFrame(PValSum(note=str(message).rstrip())._asdict(), index=[0])
    }


def parse_key(k, filename):
    key = re.sub(r"^.*-(sheet-.*)", r"\1", k) + " from " + filename
    if k == filename:
        key = k
    return key


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--file", nargs="+", help="path to input file to be parsed")
    group.add_argument(
        "--list",
        metavar="FILE",
        type=argparse.FileType("r"),
        help="file with paths to input files, one per line",
    )
    parser.add_argument("--out", metavar="FILE", help="output file")
    parser.add_argument(
        "--vars",
        metavar="KEY=VALUE",
        nargs="*",
        default=["basemean=10", "fpkm=0.5", "logcpm=-0.3", "rpkm=0.5", "aveexpr=3.32"],
        help="variables for expression level filtering. Input 'key=value' pairs without spaces around equation mark",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=30,
        help="number of histogram bins, integer, default is 30",
    )
    parser.add_argument(
        "--fdr",
        type=float,
        default=0.05,
        help="false discovery rate, float, default is 0.05",
    )
    parser.add_argument(
        "--pi0method",
        type=str,
        default="lfdr",
        help="method to calculate pi0, string, default is 'lfdr'",
    )
    parser.add_argument(
        "--verbose", "-v", help="increase output verbosity", action="count", default=0
    )
    parser.add_argument(
        "--blacklist",
        metavar="FILE",
        type=argparse.FileType("r"),
        help="file with filenames to skip importing, one per line",
    )
    args = parser.parse_args()
    var = dict(map(lambda s: s.split("="), args.vars))
    VAR = {k: float(v) for k, v in var.items()}
    VAR.update({"value": VAR["fpkm"]})
    BINS = args.bins + 1
    FDR = args.fdr

    if args.file:
        input = args.file
    elif args.list:
        input = []
        with args.list as f:
            for line in f:
                input.append(line.rstrip())

    blacklist = []
    if args.blacklist:
        with args.blacklist as f:
            for line in f:
                file = os.path.basename(line.rstrip())
                blacklist.append(file)

    # Keep only inputs that exist
    input = [i for i in input if os.path.isfile(i)]

    # Drop files in blacklist
    if blacklist:
        input = [i for i in input if os.path.basename(i) not in blacklist]

    if len(input) == 0:
        Path(args.out).touch()
        sys.exit()
    out = {}
    for path in input:
        if args.verbose > 0:
            print("working on", path)
        filename = os.path.basename(path)
        if tarfile.is_tarfile(path):
            frames = import_tar(path)
        else:
            frames = import_flat(path)
        out.update(
            {
                parse_key(k, filename): v
                for k, v in frames.items()
                if all(i in fields for i in v.columns)
            }
        )
        frames = filter_pvalue_tables(frames, pv, adj)
        if len(frames) == 0:
            out.update(note(filename, "no pvalues"))
            continue
        else:
            frames = {
                k: summarise_pvalue_tables(v, var=VAR.keys()) for k, v in frames.items()
            }
            pv_stats = {
                k: summarise_pvalues(
                    v,
                    bins=BINS,
                    fdr=FDR,
                    var={k: v for k, v in VAR.items() if "value" not in k},
                    pi0_method=args.pi0method,
                    verbose=args.verbose,
                )
                if not v.empty
                else pd.DataFrame(
                    PValSum(note="all p-values are NaN")._asdict(), index=[0]
                )
                for k, v in frames.items()
            }
            for k, v in pv_stats.items():
                out.update({parse_key(k, filename): v})

    result = pd.concat(
        [df for df in out.values()],
        keys=[k for k in out.keys()],
        names=["id"],
        sort=False,
    )
    with open(args.out, "w") as f:
        result.reset_index(level="id").to_csv(f, index=False)
