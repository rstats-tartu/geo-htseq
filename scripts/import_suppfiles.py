import pandas as pd
import os
import re

p = re.compile("p.*val")

path = "output/suppl/"
suppfiles = os.listdir(path)

for input in suppfiles:
    print("Working on: ", input)
    tabs = {}
    if "xls" in input:
        try:
            wb = pd.ExcelFile(path + input)
            sheets = wb.sheet_names
            for sheet in sheets:
                tabs.update({input + "_" + sheet : wb.parse(sheet, engine="xlrd")})
        except Exception as e:
            print(input, " file: ", e)
            pass
    else:
        try:
            r = pd.read_csv(path + input, sep=None, engine='python', iterator=True, skiprows = 20, nrows = 1000)
            sep = r._engine.data.dialect.delimiter
            r.close()
            tabs.update({input : pd.read_csv(path + input, sep=sep, comment="#")})
        except Exception as e:
            print("Error: ", e)
            pass
    for k,v in tabs.items():
        print("Table: ", k)
        print(v)
