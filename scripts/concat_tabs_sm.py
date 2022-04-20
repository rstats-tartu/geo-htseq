import pandas as pd
import argparse


def safely_read_csv(path, **kwargs):
    try:
        return pd.read_csv(path, **kwargs)
    except pd.errors.EmptyDataError:
        pass


def concatenate_tables(input, sep=","):
    frames = [safely_read_csv(f, sep=sep) for f in input]
    return pd.concat(frames, sort=False)


out = snakemake.output[0]
tabs = snakemake.input

with open(out, "w") as out:
    concatenate_tables(tabs).to_csv(out, index=False)
