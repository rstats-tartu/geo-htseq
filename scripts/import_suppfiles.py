import pandas as pd
import os
import re
import gzip

pv = re.compile("p.*val")
sp = re.compile("\s+")

dir = "output/suppl/"
suppfiles = os.listdir(dir)

for input in suppfiles:
    print("Working on: ", input)
    path = dir + input
    tabs = {}
    if "xls" in input:
        try:
            if ".gz" in input:
                with gzip.open(path) as gz:
                    wb = pd.ExcelFile(gz)
            else:
                wb = pd.ExcelFile(path)
            sheets = wb.sheet_names
            for sheet in sheets:
                parsed_sheet = wb.parse(sheet)
                if not parsed_sheet.empty:
                    tabs.update({input + "-sheet-" + sheet: parsed_sheet})
        except Exception as e:
            print(input, " file: ", e)
            pass
    else:
        try:
            r = pd.read_csv(
                path, sep=None, engine="python", iterator=True, skiprows=20, nrows=1000
            )
            sep = r._engine.data.dialect.delimiter
            r.close()
            tabs.update({input: pd.read_csv(path, sep=sep, comment="#")})
        except Exception as e:
            print("Error: ", e)
            pass
    for k, v in tabs.items():
        print("Table: ", k)
        print(v)
