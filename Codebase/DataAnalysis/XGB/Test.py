import pandas as pd
import matplotlib.pyplot as plt

from os import listdir
from os.path import isfile, join

mypath = "../../Data"
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
for file in onlyfiles:
    df = pd.read_hdf(f"../../Data/{file}")
    if "data" in file:
        df["channel"] = file[:-7]
    print(df)

