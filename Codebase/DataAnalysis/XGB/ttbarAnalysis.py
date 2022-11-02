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


myPath = "/storage/William_Sakarias/William_Data"

signal = "ttbar"

df, y, df_data, channels = loadDf(myPath, signal)

print("Preparing data....")
train, val = splitAndPrepData(df, y, scale = True)
print("Done.")

X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val = val

"""pca = decomposition.PCA(n_components=15)
pca.fit(X_train)
X_train = pca.transform(X_train)
X_val = pca.transform(X_val)"""


sum_wpos = W_train[Y_train == 1.0].sum()
sum_wneg = W_train[Y_train == 0.0].sum()
print(np.float32(sum_wneg/sum_wpos))



xgb = XGB.XGBClassifier(
            max_depth=3, 
            n_estimators=50,
            learning_rate=0.1,
            n_jobs=4,
            tree_method="hist",
            objective='binary:logistic',
            eval_metric = 'auc',
            scale_pos_weight = np.float32(sum_wneg/sum_wpos)[0],
            missing=-999.0,
            use_label_encoder=False,
            ) 
            
time = timer()
print("Training....")
xgb = xgb.fit(X_train, Y_train, sample_weight = W_train,eval_set= [(X_val, Y_val)], sample_weight_eval_set = [W_val])
print("Done")
timer(time)

plotRoc(Y_train, 
        xgb.predict_proba(X_train)[:,1], 
        W_train,
        "Training-data", 
        plot = True,
        return_score = True, 
        name = f"../../../thesis/Figures/MLResults/XGB/{signal}SearchROCTrain.pdf")

plotRoc(Y_val, 
        xgb.predict_proba(X_val)[:,1], 
        W_val,
        "Validation-data", 
        plot = True,
        return_score = True, 
        name = f"../../../thesis/Figures/MLResults/XGB/{signal}SearchROCVal.pdf")

"""
channel = df.channel
wgt = df.wgt_SG
df = df.drop(columns=["wgt_SG","channel"])


predict_sorted, weights_sorted =  separateByChannel(xgb.predict_proba(df)[:,1], wgt, channel, channels)
predict_data = xgb.predict_proba(df_data)[:,1]


print("Plotting...")
PlotRootHisto(predict_sorted, 
              weights_sorted, 
              predict_data, 
              channels, 
              title = f"XGB/{signal}SearchDist", 
              xlabel = "XGB-Output", 
              bins = 30)
PlotRootHisto(predict_sorted, 
              weights_sorted, 
              predict_data, 
              channels, 
              title = f"XGB/{signal}SearchDist_C7", 
              xlabel = "XGB-Output", 
              bins = 30,
              CutOff = 0.7)

plotFI(xgb, df.keys(), signal)
print("Finshed plots.")
"""



