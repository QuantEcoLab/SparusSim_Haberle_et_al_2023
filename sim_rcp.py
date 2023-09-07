import pathlib
import pandas as pd
import geopandas as gpd
from sim import sim
from tqdm import tqdm
import numpy as np

# Please change this line to the path of the folder containing the data
files = list(pathlib.Path('<PATH TO FOLDER CONTAINING LOCAL TIMESERIES DATA>').glob('*.csv'))

files.sort()
files = files

# Starts of the simulation periods in the format YYYY-MM-DD
# Change according your data
# Multiple starts can be used to simulate multiple periods
periodS_s = [
    "2091-4-1",
]

periodS = pd.to_datetime(periodS_s, format="%Y-%m-%d").values

for per in periodS:
    df_4_5 = {
        "TTM": [],
        "FCR": [],
        "W2Y": [],
        "FCR2Y": [],
        "MTEMP2Y": [],
        "NBRD2Y": []
    }

    df_8_5 = {
        "TTM": [],
        "FCR": [],
        "W2Y": [],
        "FCR2Y": [],
        "MTEMP2Y": [],
        "NBRD2Y": []
    }

    for i in tqdm(range(len(files))):
        df = pd.read_csv(files[i], delimiter=",", index_col=0)
        df.index = pd.to_datetime(df.date, format="%Y-%m-%d")

        df = df.loc[per: per+pd.DateOffset(years=3)-pd.DateOffset(days=1)]
        temp = df["rcp4_5_2.5m"]+273.15
        rez = sim(temp)
        df_4_5["TTM"].append(rez["TTM"])
        df_4_5["FCR"].append(rez["FCR"])
        df_4_5["W2Y"].append(rez["W2Y"])
        df_4_5["FCR2Y"].append(rez["FCR2Y"])
        df_4_5["MTEMP2Y"].append(rez["MTEMP2Y"])
        df_4_5["NBRD2Y"].append(rez["NBRD2Y"])

        temp = df["rcp8_5_2.5m"]+273.15
        rez = sim(temp)
        df_8_5["TTM"].append(rez["TTM"])
        df_8_5["FCR"].append(rez["FCR"])
        df_8_5["W2Y"].append(rez["W2Y"])
        df_8_5["FCR2Y"].append(rez["FCR2Y"])
        df_8_5["MTEMP2Y"].append(rez["MTEMP2Y"])
        df_8_5["NBRD2Y"].append(rez["NBRD2Y"])

    df_4_5 = pd.DataFrame(df_4_5)
    df_8_5 = pd.DataFrame(df_8_5)

    index_ = [f.stem for f in files]
    df_4_5.index = [i[3:] for i in index_]
    df_8_5.index = [i[3:] for i in index_]


    df_4_5.to_csv(f"out/TTM_FCR_4_5_{np.datetime_as_string(per, unit='Y')}.csv", index_label="id")
    
   
    df_8_5.to_csv(f"out/TTM_FCR_8_5_{np.datetime_as_string(per, unit='Y')}.csv", index_label="id")

    print(f"Done with {per}")


