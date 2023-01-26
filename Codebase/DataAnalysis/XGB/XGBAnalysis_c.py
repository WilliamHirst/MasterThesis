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

name = "XGB"
signal = "SUSY"
train = True


print(f"Starting test: Model = {name} -- Signal = {signal}")

df, y, df_data, channels = loadDf(myPath, notInc=["ttbarHNLfull","LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75"])

print("Preparing data....")
train, val = splitAndPrepData(df, y, scale = True, ret_scaleFactor=True, removeNeg=True)#, PCA=True, n_components=1-1e-3)
print("Done.")

X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val, scaleFactor = val
nrFeature = nFeats(X_train)





sum_wpos = W_train[Y_train == 1.0].sum()
sum_wneg = W_train[Y_train == 0.0].sum()

xgb = XGB.XGBClassifier(
            max_depth=4, 
            n_estimators=100,
            learning_rate=0.1,
            n_jobs=4,
            subsample = 0.8,
            tree_method="hist",
            objective='binary:logistic',
            eval_metric = 'auc',
            scale_pos_weight = np.float32(sum_wneg/sum_wpos)[0],
            missing=-999.0,
            use_label_encoder=False,
            ) 
            
xgb = xgb.fit(X_train, Y_train, sample_weight = W_train, eval_set= [(X_val, Y_val)], sample_weight_eval_set = [W_val], early_stopping_rounds = 10)

W = df.wgt_SG
C = df.channel
Y = y
df = df.drop(columns = ["channel", "wgt_SG"])
df, df_data = scaleData(df,df_data)

HM(xgb, df, Y, W, C, data = None, name = f"SUSY/{name}Grid", metric="Sig", save = False, mlType = 'XGB')
    
    
    

    

# mc_predict, mc_weights = separateByChannel(prediction, weights, C, C.unique())

# saveLoad("results/predict_sorted_test.npy", mc_predict)
# saveLoad("results/weights_sorted_test.npy", mc_weights)
# saveLoad("results/predict_data_test.npy", model(df_data))
# saveLoad("results/channels_test.npy", C)