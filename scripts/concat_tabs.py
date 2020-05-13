import pandas as pd


def safely_read_csv(path, **kwargs):
    try:
        return pd.read_csv(path, **kwargs)
    except pd.errors.EmptyDataError:
        pass


def concatenate_tables(input, sep=","):
    frames = [safely_read_csv(f, sep=sep) for f in input]
    return pd.concat(frames, sort=False)


with open(snakemake.output[0], "w") as out:
    concatenate_tables(snakemake.input).to_csv(out, index=False)
