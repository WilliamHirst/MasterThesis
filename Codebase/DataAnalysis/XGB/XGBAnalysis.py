from time import time
import xgboost as XGB
from sklearn import decomposition


import sys
sys.path.insert(1, "../")
from Plot_stuff.plot_set import *
from Plot_stuff.ROCM import *
from Plot_stuff.FI import *
from Plot_stuff.ROOTPlot import *

sys.path.insert(1, "../../")
from Utilities import *
from Plot_stuff.HM import *

myPath = "/storage/William_Sakarias/William_Data"

name = "XGBNoWeights"
signal = "SUSY"
train = True
IncludeRange = [400, 800, 0, 400] 
notInc=["ttbarHNLfull","LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75", "500p0", "550p0", "600p0", "450p0p0350p0", "450p0p0250p0","450p0p0150p0", "450p0p050p0", "400p0p00p0","400p0p0100p0","400p0p0200p0", "400p0p0300p0","650p050p0", "650p0150p0", "650p0250p0", "650p0350p0", "700p0100p0", "700p00p0", "700p0200p0", "700p0300p0", "700p0400p0"]

print(f"Starting test: Model = {name} -- Signal = {signal}")

df, y, df_data, channels = loadDf(myPath, notInc=notInc, IncludeRange=IncludeRange)
print("Preparing data....")
train, val = splitAndPrepData(df, y, scale = True, ret_scaleFactor=True, removeNeg=True)#, PCA=True, n_components=1-1e-3)
print("Done.")

X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val, scaleFactor = val
nrFeature = nFeats(X_train)

sum_wpos = W_train[Y_train == 1.0].sum()
sum_wneg = W_train[Y_train == 0.0].sum()

W_train[:] = 1 
W_val[:] = 1 

W_train[Y_train == 0.0] *= 1/np.sum(W_train[Y_train == 1.0])
W_val[Y_val == 0.0] *= 1/np.sum(W_val[Y_val == 1.0])

xgb = XGB.XGBClassifier(
            tree_method="hist",
            objective='binary:logistic',
            eval_metric = 'auc',
            use_label_encoder=False,
            #scale_pos_weight = np.float32(sum_wneg/sum_wpos)[0],
            ) 
            
xgb = xgb.fit(X_train, Y_train, early_stopping_rounds = 10, eval_set= [(X_val, Y_val)], sample_weight = W_train, sample_weight_eval_set = [W_val])

W = df.wgt_SG
C = df.channel
Y = y
df = df.drop(columns = ["channel", "wgt_SG"])
df, df_data = scaleData(df,df_data)

HM(xgb, df, Y, W, C, data = None, name = f"{name}Grid", metric="Sig", save = True, mlType = 'XGB')
    
    
    

    

# mc_predict, mc_weights = separateByChannel(prediction, weights, C, C.unique())

# saveLoad("results/predict_sorted_test.npy", mc_predict)
# saveLoad("results/weights_sorted_test.npy", mc_weights)
# saveLoad("results/predict_data_test.npy", model(df_data))
# saveLoad("results/channels_test.npy", C)