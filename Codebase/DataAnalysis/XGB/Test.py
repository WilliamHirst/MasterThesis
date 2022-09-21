import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import xgboost as XGB
from Plot_stuff.RHM import ROOT_Histo_Maker
from Plot_stuff.plot_set import *



from os import listdir
from os.path import isfile, join

mypath = "/scratch/William"
signal = "ttbar"
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
df = pd.DataFrame()
channels = []
for i in range(len(onlyfiles)):
    if "data" not in onlyfiles[i]:
        channel = onlyfiles[i][:-7]
        df_i = pd.read_hdf(f"{mypath}/{onlyfiles[i]}")
        df_i["isSignal"] = signal == channel
        df_i["channel"] = channel
        df = df.append(df_i, ignore_index=True)
        channels.append(channel)

xgb = XGB.XGBClassifier(
            max_depth=3, 
            n_estimators=120,
            learning_rate=0.1,
            n_jobs=4,
            tree_method="hist",
            objective='binary:logistic',
            missing=-999.0,
            use_label_encoder=False,
            eval_metric="error") 




df_train, df_test = train_test_split(df, test_size=0.2)

y_train = df_train.isSignal
df_train_channel = df_train.channel
df_train_weights = df_train.wgt_SG
df_train = df_train.drop(columns=["isSignal", "channel", "wgt_SG"])


xgb = xgb.fit(df_train, y_train, sample_weight = df_train_weights)
mc_predict = []
mc_weights = []
for channel in channels:
    isChannel = df_train_channel == channel
    channel_pred = xgb.predict_proba(df_train[isChannel])[:,1]
    mc_predict.append(channel_pred)
    mc_weights.append(df_train_weights[isChannel])



ROOT_Histo_Maker(mc_predict, 
                 mc_weights,channels, 
                 bin_max = 1, 
                 bin_min = 0,
                 nr_bins = 20, 
                 y_max= 1e8, 
                 y_min = 0.5, 
                 variable_name = r"$XGB-Output$", 
                 saveAs = "ttbarSearch.pdf")

