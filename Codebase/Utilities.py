from locale import D_T_FMT
import numpy as np
import pandas as pd



def loadDf(location, signal):
    from os import listdir
    from os.path import isfile, join
    onlyfiles = [f for f in listdir(location) if isfile(join(location, f))]
    df = pd.DataFrame()
    y = np.array([])
    channels = []
    for i in range(len(onlyfiles)):
        if "data" not in onlyfiles[i]:
            channel = onlyfiles[i][:-7]
            df_i = pd.read_hdf(f"{location}/{onlyfiles[i]}")
            df_i["isSignal"] = channel == signal
            df_i["channel"] = channel
            df = df.append(df_i, ignore_index=True)
            channels.append(channel)
        else:
            df_data = pd.read_hdf(f"{location}/{onlyfiles[i]}")
            df_data = df_data.drop(columns = ["wgt_SG"])
    y = df.isSignal
    df = df.drop(columns = ["isSignal"])
    return df, y, df_data, channels

def mergeToRoot(MC, MC_wgt, Data, Channels):
    df = {}
    for i in range(len(Channels)):
        df_i = pd.DataFrame()
        df_i["ML_Val"] = np.asarray(MC[i])
        df_i["wgt"] = np.asarray(MC_wgt[i])
        df[Channels[i]] = df_i

    df_i = pd.DataFrame()
    df_i["ML_Val"] = np.asarray(Data)
    df_i["wgt"] = np.ones(len(Data))
    df["Data"] = df_i
    print(df)

    import ROOT
    rdf = ROOT.RDF.MakeNumpyDataFrame(df)

    return rdf


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


def splitData(X, Y, split_v = 0.2, isEven = False, split_b = 0.2):
    from sklearn.utils import shuffle

    X,Y = X.to_numpy(), Y.to_numpy()

    X,Y = shuffle(X,Y, random_state=2)

    x_s = X[Y == 1]
    x_b = X[Y == 0]
    y_s = Y[Y == 1] 
    y_b = Y[Y == 0]


    rng = np.random.default_rng(seed=2)
    size_s = len(x_s)

    if isEven:
        size_b = size_s
    else:
        size_b = int(len(x_b)*split_b)

    
    indx = rng.choice(len(x_b), size_b, replace=False)
    
    re_indx = np.delete(np.arange(len(x_b)), indx)

    split_indx  = int((size_s+size_b)*split_v)
   
    X_train = np.concatenate((x_s[:,:-2], x_b[indx,:-2]), axis=0)
    W_train = np.concatenate((x_s[:,-2],  x_b[indx,-2]))
    C_train = np.concatenate((x_s[:,-1],  x_b[indx,-1]))
    Y_train = np.concatenate((y_s, y_b[indx]))

    X_train, Y_train, W_train, C_train = shuffle(X_train, Y_train, W_train, C_train, random_state=2)
    X_train, Y_train, W_train, C_train = X_train[split_indx:], Y_train[split_indx:], W_train[split_indx:], C_train[split_indx:]
    X_val, Y_val, W_val, C_val = X_train[0:split_indx], Y_train[0:split_indx], W_train[0:split_indx], C_train[0:split_indx]

    if not isEven:
        W_train[Y_train == 0.0] *= np.sum(W_train[Y_train == 1 ]) / np.sum(W_train[Y_train == 0])
        W_val[Y_val == 0.0] *= np.sum(W_val[Y_val == 1 ]) / np.sum(W_val[Y_val == 0])

    C_test = x_b[:,-1]
    W_test = x_b[:,-2]#[re_indx,-2] 
    X_test = x_b[:,:-2]#[re_indx,:-2]
    Y_test = y_b#[re_indx]
    
    W_train = removeNegWeights(W_train)

    return (X_train, Y_train, W_train, C_train), (X_val, Y_val, W_val, C_val), (X_test, Y_test, W_test, C_test)