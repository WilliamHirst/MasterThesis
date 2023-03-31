from locale import D_T_FMT
import pickle as pkl
import numpy as np
import pandas as pd
from DataAnalysis.FeatureSelection import lowFeats 
import json



"""
MERGE AND TRANSFORM DATAFRAMES.
"""
def loadDf(location, signal = None, incHigh = True, notInc = [], IncludeRange = None):
    from os import listdir
    from os.path import isfile, join
    onlyfiles = [f for f in listdir(location) if isfile(join(location, f))]
    df = pd.DataFrame()
    df_data = pd.DataFrame()
    y = np.array([])
    channels = []
    for i in range(len(onlyfiles)):
        if checkIfInclude(notInc, IncludeRange, onlyfiles[i]): continue
            
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



def mergeToRoot(MC, MC_wgt, Channels, Data=None, CutOff = None):
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
    if Data is not None:
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
                     removeNeg = False):
    """_summary_

    Args:
        X (_type_): Feature data set
        Y (_type_): Labels for data set
        split_v (float, optional): Ratio to split data in training and background. 
                                   Defaults to 0.2.
        scaleWeight (bool, optional): Bool . Defaults to True.
        scale (bool, optional): _description_. Defaults to False.
        PCA (bool, optional): _description_. Defaults to False.
        n_components (_type_, optional): _description_. Defaults to None.
        ret_scaleFactor (bool, optional): _description_. Defaults to False.
        removeNeg (bool, optional): _description_. Defaults to False.

    Returns:
        _type_: _description_
    """

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

def AddParameters(df,Y,data = None):
    C = df["channel"] 
    columns_s = C[Y.to_numpy() == 1] 
    unique_c = columns_s.unique()
    dists = {}
    b_indx = Y.to_numpy() == 0
    df["param1"] = np.nan
    df["param2"] = np.nan

    for c in unique_c:
        elem = c.split("WZ")
        m1 = elem[1][0:3]
        if "p0p0" in c:
            elem = c.split("p0p0")
        else:
            elem = c.split("p0")
        m2 = elem[1]  
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

    if data is None:
        return df
        
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


def Calc_Sig(y_MC, y_label, sample_weight, y_Data = None,sf = None, best_threshold = None, max_sig = None, returnNr = False):
    bkg_indx = y_label == 0

    sample_weight.loc[np.ravel(np.asarray(y_label==0))] /=  sf

    if max_sig is None:
        max_sig = 0.9999

    m_s = 0
    m_b = 0
    thresholds = np.concatenate((np.linspace(0.9,0.98,100),np.linspace(0.98,max_sig,100) ))
    max = 0
    for thresh in thresholds:
        nrB = np.sum(sample_weight[np.ravel(y_MC>thresh) * np.ravel(bkg_indx)].to_numpy() )
        nrS = np.sum(sample_weight[np.ravel(y_MC>thresh) * np.ravel(y_label == 1) ].to_numpy() )
        sig  = np.sqrt(2*((nrS + nrB)*np.log(1+nrS/nrB)-nrS))
        if sig >max:
            max = sig
            best_threshold = thresh
            m_s = nrS
            m_b = nrB
    sig = max
    if y_Data is not None:
        nrB = int(np.sum(sample_weight[np.ravel(y_MC>best_threshold)]))
        nrS = len(y_Data[np.ravel(y_Data>best_threshold)]) - nrB
        nrS *= nrS>0
        sig  = np.sqrt(2*((nrS + nrB)*np.log(1+nrS/nrB)-nrS))

    print(m_b)
    print(m_s)
    print(best_threshold)

    print(f"The significance: {sig}.")
    if returnNr:
        return sig, nrB, nrS
        
def saveToTxt(m1, m2, nbkg, nsig, method):
    from os.path import exists
    path = f'../results/{method}Sig.txt'
    if not exists(path):
        f = open(path, 'w')
        f.write(f"m1    m2    nbkg    nsig")
    else:
        f = open(path, 'a')
    f.write(f"{m1}    {m2}    {nbkg}    {nsig}")
    
    f.close()

        
    

def saveToJson(score, m1, m2, metric, method):
    if metric == "Auc":
        file = 'AUC'
    elif metric == "NrEvents":
        file = metric
    else:
        file = 'SIG'
        
    with open(f'../results/{file}.json', 'r') as openfile:
        # Reading from json file
        json_object = json.load(openfile)
    
    json_object[method][f"{m1}_{m2}"] = {'score': score, 'm1': m1, 'm2': m2}

    json_object = json.dumps(json_object, indent=4)
 
    # Writing to sample.json
    with open(f'../results/{file}.json', "w") as outfile:
        outfile.write(json_object)

def EmptyJson( metric, method):
    if metric == "Auc":
        file = 'AUC'
    else:
        file = 'SIG'

    with open(f'../results/{file}.json', 'r') as openfile:
        # Reading from json file
        json_object = json.load(openfile)
    
    json_object[method] = {}

    json_object = json.dumps(json_object, indent=4)
 
    # Writing to sample.json
    with open(f'../results/{file}.json', "w") as outfile:
        outfile.write(json_object)

def checkIfInclude(notInc, IncludeRange, channel): 
    for chan in notInc:
        if chan in channel:
            return 1

    # Checks if signal has mass inside inclusive mass range.
    if IncludeRange is not None and "MGPy8EGA14N23" in channel:
        m1Min, m1Max, m2Min, m2Max = IncludeRange
        m1, m2 = getMass(channel)
        if int(m1) >= m1Min and int(m1) <= m1Max and int(m2) >= m2Min and int(m2) <= m2Max:
            return 0
        return 1
    return 0

def getMass(string):
    elem = string.split("WZ")
    m1 = elem[1][0:3]
    if "p0p0" in string:
        elem = string.split("p0p0")
    else:
        elem = string.split("p0")
    m2 = elem[1]  
    return m1, m2
    



    

