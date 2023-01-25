from locale import D_T_FMT
import pickle as pkl
import numpy as np
import pandas as pd
from DataAnalysis.FeatureSelection import lowFeats 
import json





"""
MERGE AND TRANSFORM DATAFRAMES.
"""
def loadDf(location, signal = None, incHigh = True, notInc = []):
    from os import listdir
    from os.path import isfile, join
    onlyfiles = [f for f in listdir(location) if isfile(join(location, f))]
    df = pd.DataFrame()
    df_data = pd.DataFrame()
    y = np.array([])
    channels = []
    
    for i in range(len(onlyfiles)):

        cont = 0
        for chan in notInc:
            if chan in onlyfiles[i]:
                cont = 1            
        if cont: continue
        if "data" not in onlyfiles[i]:
            channel = onlyfiles[i][:-7]
            df_i = pd.read_hdf(f"{location}/{onlyfiles[i]}")
            df_i["channel"] = channel
            df = df.append(df_i, ignore_index=True,sort=False)
            channels.append(channel)
        else:

            data_i = pd.read_hdf(f"{location}/{onlyfiles[i]}")
            data_i = data_i.drop(columns = ["wgt_SG", "type"])
            df_data = df_data.append(data_i, ignore_index=True,sort=False) 
    if signal is None:
        y = df.type
    else:
        df = df[df.type != 1]
        y = df["channel"] == signal

    df = df.drop(columns = ["type"])
    if not incHigh:
        df = df.drop(columns = [feat for feat in df.keys() if feat not in lowFeats])
        df_data = df_data.drop(columns = [feat for feat in df.keys() if feat not in lowFeats]) 
    return df, y, df_data, channels

def mergeToRoot(MC, MC_wgt, Data, Channels, CutOff = None):
    import ROOT
    df = {}
    if CutOff is None:
        CutOff = 0
    for i in range(len(Channels)):
        df_i = pd.DataFrame()
        ML_Val = np.array(MC[i],dtype=np.float64)
        wgt = np.array(MC_wgt[i], dtype=np.float64)        
        df_i = {"ML_Val": ML_Val[ML_Val >= CutOff], "wgt": wgt[ML_Val >= CutOff] }
        df[Channels[i]] = ROOT.RDF.MakeNumpyDataFrame(df_i)

    df_i = pd.DataFrame()
    ML_Val = np.array(Data,dtype=np.float64)
    wgt = np.ones(len(Data),dtype=np.float64)
    df_i = {"ML_Val": ML_Val[ML_Val >= CutOff], "wgt": wgt[ML_Val >= CutOff] }
    df["data18"] = ROOT.RDF.MakeNumpyDataFrame(df_i)
    return df

def saveLoad(name, data = None, type = "Numpy"):
    if data is None:
        if type ==   "Numpy":
            with open(f"{name}", 'rb') as file:
                output = np.load(file,allow_pickle=True)
        elif type == "Pandas":
            output = pd.read_hdf(f"{name}")
        return output
    else:
        if type ==   "Numpy":
            with open(f"{name}", 'wb') as f:
                np.save(f, np.asarray(data,dtype=object).ravel())
        elif type == "Pandas":
            data.to_hdf(f"{name}","mini")
        return 
def loadSamples(ex_path = ""):
    X_train = saveLoad(f"{ex_path}samples/X_train_sample.hdf5", type = "Pandas")
    Y_train = saveLoad(f"{ex_path}samples/Y_train_sample.hdf5", type = "Pandas")
    W_train = saveLoad(f"{ex_path}samples/W_train_sample.hdf5", type = "Pandas")
    C_train = saveLoad(f"{ex_path}samples/C_train_sample.hdf5", type = "Pandas")

    X_val = saveLoad(f"{ex_path}samples/X_val_sample.hdf5", type = "Pandas")
    Y_val = saveLoad(f"{ex_path}samples/Y_val_sample.hdf5", type = "Pandas")
    W_val = saveLoad(f"{ex_path}samples/W_val_sample.hdf5", type = "Pandas")
    C_val = saveLoad(f"{ex_path}samples/C_val_sample.hdf5", type = "Pandas")
    Tr = (X_train, Y_train, W_train, C_train)
    Va = (X_val, Y_val, W_val, C_val)
    return Tr, Va
def separateByChannel(prediction, weights, df_channel, channels):
    mc_predict = []
    mc_weights = []

    for channel in channels:
        isChannel = df_channel == channel
        mc_predict.append(prediction[isChannel].ravel())
        mc_weights.append(weights[isChannel].ravel())
    return mc_predict, mc_weights



"""
DATA-HANDLING FUNCTIONS 
"""
def splitAndPrepData(X,  
                     Y, 
                     split_v = 0.2, 
                     scaleWeight = True, 
                     scale = False, 
                     PCA = False, 
                     n_components = None, 
                     ret_scaleFactor = False,
                     addParam = False,
                     removeNeg = False):
    from sklearn.model_selection import train_test_split
   
    X_train, X_val, Y_train, Y_val = train_test_split(X, Y, test_size = split_v, random_state=42)

    W_train = X_train.wgt_SG.array
    W_val = X_val.wgt_SG.array

    C_train = X_train.channel
    C_val = X_val.channel
   
    if scaleWeight:
        W_train, _ = scaleWeights(W_train, Y_train)
        W_val, scaleFactor = scaleWeights(W_val, Y_val, ret_scaleFactor)

    if removeNeg: 
        W_train = removeNegWeights(W_train)
        W_val = removeNegWeights(W_val)

    W_train = pd.DataFrame(W_train)
    W_val = pd.DataFrame(W_val)

    X_train = X_train.drop(columns = ["channel", "wgt_SG"])
    X_val = X_val.drop(columns = ["channel", "wgt_SG"])

    X_train.index = np.arange(len(X_train))
    X_val.index = np.arange(len(X_val))
    Y_train.index = np.arange(len(Y_train))
    Y_val.index = np.arange(len(Y_val))

    print("Anna er kul")
    
    if scale:
        X_train, X_val = scaleData(X_train, X_val, scaler = "Standard")

    if PCA:
        X_train, X_val = PCAData(X_train, X_val, n_components = n_components)


    Tr = (X_train.astype('float32'), Y_train.astype('float32'), W_train.astype('float32'), C_train)
    if ret_scaleFactor:
        Va = (X_val.astype('float32'), Y_val.astype('float32'), W_val.astype('float32'), C_val, scaleFactor)
    else:
        Va = (X_val.astype('float32'), Y_val.astype('float32'), W_val.astype('float32'), C_val)

    return Tr, Va

def removeNegWeights(weight):
    P = np.sum(weight)/np.sum(np.abs(weight))
    return P*np.abs(weight)

def scaleWeights(weights, y, ret_scaleFactor = False):
    scaleFactor = np.sum(weights[y == 1]) / np.sum(weights[y == 0])
    weights[y == 0.0] *= scaleFactor
    if ret_scaleFactor:
        return weights, scaleFactor
    return weights, None


def scaleData(X_train, X_val = None, scaler = "Standard"):
    if scaler == "Standard":
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler()
    elif scaler == "MinMax":
        from sklearn.preprocessing import MinMaxScaler
        scaler = MinMaxScaler()
    scaler.fit(X_train)
    X_train[X_train.keys()] = scaler.transform(X_train[X_train.keys()])
    if X_val is None:
        return X_train

    X_val[X_val.keys()] = scaler.transform(X_val[X_val.keys()])
    return X_train, X_val

def PCAData(X_train, X_val = None, n_components = None):
    from sklearn.decomposition import PCA
    if n_components is None:
        pca = PCA()
    elif n_components < 1 and n_components > 0:
        pca = PCA(n_components = n_components, svd_solver = 'full' )
    else:
        pca = PCA(n_components)

    bf = nFeats(X_train)
    pca = pca.fit(X_train)

    X_train = pca.transform(X_train)
    af = nFeats(X_train)
    print(f"Removed {bf-af} features. {af} features left.")

    if X_val is None:
        return X_train
    
    X_val = pca.transform(X_val)
    return X_train, X_val

def AddParameters(df,Y,data):
    C = df["channel"] 
    columns_s = C[Y.to_numpy() == 1] 
    unique_c = columns_s.unique()
    dists = {}
    b_indx = Y.to_numpy() == 0
    df["param1"] = np.nan
    df["param2"] = np.nan

    for c in unique_c:
        txt = c.split("p0")
        m1 = txt[0][21:24]
        m2 = txt[2] 
        indx = (C == c).to_numpy()
        df.loc[indx,"param1"] = m1
        df.loc[indx,"param2"] = m2
        dists[f"{np.sum(indx)}"] = [m1,m2]

    prob = np.asarray([int(d) for d in dists.keys()])
    vals = np.asarray(list(dists.values()))
    keys = np.asarray(list(dists.keys()))

    bkg_params =  vals[np.random.choice(np.arange(len(keys)),np.sum(b_indx), p=prob/np.sum(prob) )]

    df.loc[b_indx,"param1"] = bkg_params[:,0]
    df.loc[b_indx,"param2"] = bkg_params[:,1]

    data_params =  vals[np.random.choice(np.arange(len(keys)),len(data), p=prob/np.sum(prob) )]
    data["param1"] = data_params[:,0]
    data["param2"] = data_params[:,1]
    return df, data


"""
EXTRA FUNCTIONS
"""
def timer(start_time=None):
    from datetime import datetime
    if not start_time:
        start_time = datetime.now()
        return start_time
    elif start_time:
        thour, temp_sec = divmod((datetime.now() - start_time).total_seconds(), 3600)
        tmin, tsec = divmod(temp_sec, 60)
        print('\nTime taken: %i hours %i minutes and %s seconds.' % (thour, tmin, round(tsec, 2)))

def nFeats(data):
    try:
        nF = len(data.keys())
    except:
        print(data)
        nF = len(data[0])
    return nF


def Calc_Sig(y_MC, y_label, sample_weight, y_Data = None,sf = None, best_threshold = None):
    if best_threshold is None:
        from sklearn.metrics import roc_curve

        fpr, tpr, thresholds = roc_curve(y_label, y_MC, sample_weight = sample_weight, pos_label=1)

        gmeans = np.sqrt(np.array(tpr) * (1-np.array(fpr)/np.max(np.array(fpr))))
        ix = np.argmax(gmeans)
        best_threshold = thresholds[ix]

    bkg_indx = y_label == 0

    sample_weight.loc[np.ravel(np.asarray(y_label==0))] /=  sf
  
    if y_Data is None:
        thresholds = np.concatenate((np.linspace(0.9,0.98,100),np.linspace(0.98,0.9999,100) ))
        max = 0
        for thresh in thresholds:
            nrB = np.sum(sample_weight[np.ravel(y_MC>thresh) * np.ravel(bkg_indx)].to_numpy() )
            nrS = np.sum(sample_weight[np.ravel(y_MC>thresh) * np.ravel(y_label == 1) ].to_numpy() )
            sig  = np.sqrt(2*((nrS + nrB)*np.log(1+nrS/nrB)-nrS))
            if sig >max:
                max = sig
                max_thresh = thresh 
        sig = max
    else:
        nrB = int(np.sum(sample_weight[np.ravel(y_MC>best_threshold)]))
        nrS = len(y_Data[np.ravel(y_Data>best_threshold)]) - nrB
        nrS *= nrS>0
        sig  = np.sqrt(2*((nrS + nrB)*np.log(1+nrS/nrB)-nrS))

    print(f"The significance: {sig}.")
    return sig

def saveToJson(score, m1, m2, metric, method):
    if metric == "Auc":
        file = 'AUC'
    else:
        file = 'SIG'
    with open(f'results/{file}.json', 'r') as openfile:
        # Reading from json file
        json_object = json.load(openfile)
    
    json_object[method][f"{m1}_{m2}"] = {'score': score, 'm1': m1, 'm2': m2}

    json_object = json.dumps(json_object, indent=4)
 
    # Writing to sample.json
    with open(f'results/{file}.json', "w") as outfile:
        outfile.write(json_object)

def EmptyJson( metric, method):
    if metric == "Auc":
        file = 'AUC'
    else:
        file = 'SIG'

    with open(f'results/{file}.json', 'r') as openfile:
        # Reading from json file
        json_object = json.load(openfile)
    
    json_object[method] = {}

    json_object = json.dumps(json_object, indent=4)
 
    # Writing to sample.json
    with open(f'results/{file}.json', "w") as outfile:
        outfile.write(json_object)




    

