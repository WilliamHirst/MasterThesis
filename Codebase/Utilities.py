import numpy as np
import pandas as pd




def separateByChannel(prediction, weights, df_channel, channels):
    mc_predict = []
    mc_weights = []

    for channel in channels:
        isChannel = df_channel == channel
        channel_pred = prediction[isChannel]
        mc_predict.append(channel_pred)
        mc_weights.append(weights[isChannel])
    return mc_predict, mc_weights


def removeNegWeights(weight):
    P = np.sum(weight)/np.sum(np.abs(weight))
    return P*np.abs(weight)

def scaleWeights(weights, y):
    weights[y == 0.0] *= np.sum(weights[y == 1 ]) / np.sum(weights[y == 0])
    return weights

def getXYW(df, signal = "isSignal", channel = "channel", wgt = "wgt_SG") :
    y = df[signal]
    channels = df[channel]
    weights = df[wgt]

    df = df.drop(columns=[signal, channel, wgt])

    return df, y, weights, channels







def timer(start_time=None):
    from datetime import datetime
    if not start_time:
        start_time = datetime.now()
        return start_time
    elif start_time:
        thour, temp_sec = divmod((datetime.now() - start_time).total_seconds(), 3600)
        tmin, tsec = divmod(temp_sec, 60)
        print('\nTime taken: %i hours %i minutes and %s seconds.' % (thour, tmin, round(tsec, 2)))

def loadDf(location, signal):
    from os import listdir
    from os.path import isfile, join
    onlyfiles = [f for f in listdir(location) if isfile(join(location, f))]
    df = pd.DataFrame()
    channels = []
    for i in range(len(onlyfiles)):
        if "data" not in onlyfiles[i]:
            channel = onlyfiles[i][:-7]
            df_i = pd.read_hdf(f"{location}/{onlyfiles[i]}")
            df_i["isSignal"] = signal == channel
            df_i["channel"] = channel
            df = df.append(df_i, ignore_index=True)
            channels.append(channel)
    return df, channels
