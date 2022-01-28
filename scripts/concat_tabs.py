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


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--tabs", nargs="+", help="paths to input files to be concatenated")
    parser.add_argument("--out", metavar="FILE", help="output file", required=True)
    args = parser.parse_args()

    with open(args.out, "w") as out:
        concatenate_tables(args.tabs).to_csv(out, index=False)
