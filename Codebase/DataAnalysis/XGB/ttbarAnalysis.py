from time import time
from sklearn.model_selection import train_test_split
import xgboost as XGB
from Plot_stuff.plot_set import *
from Plot_stuff.ROCM import *
from Plot_stuff.FI import *
from Plot_stuff.ROOTPlot import *

import sys
sys.path.insert(1, "../../")
from Utilities import *


myPath = "/storage/William_Sakarias/William_Data"

signal = "ttbar"


xgb = XGB.XGBClassifier(
            max_depth=4, 
            n_estimators=120,
            learning_rate=0.1,
            n_jobs=4,
            tree_method="hist",
            objective='binary:logistic',
            missing=-999.0,
            use_label_encoder=False,
            eval_metric="error") 

df, y, df_data, channels = loadDf(myPath, signal)
train, val, test = splitData(df, y)

X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val = val
X_test, Y_test, W_test, C_test = test

time = timer()
print("Training....")
xgb = xgb.fit(X_train, Y_train, eval_set= [(X_val, Y_val)], sample_weight = W_train)
print("Done")
timer(time)

plotRoc(Y_train, 
        xgb.predict_proba(X_train)[:,1], 
        W_train,
        "Training-data", 
        return_score = True, 
        name = f"../../../thesis/Figures/MLResults/XGB/{signal}SearchROCTrain.pdf")

plotRoc(Y_val, 
        xgb.predict_proba(X_val)[:,1], 
        W_val,
        "Validation-data", 
        return_score = True, 
        name = f"../../../thesis/Figures/MLResults/XGB/{signal}SearchROCVal.pdf")


channel = df.channel
wgt = df.wgt_SG
df = df.drop(columns=["wgt_SG","channel"])

predict_sorted, weights_sorted =  separateByChannel(xgb.predict_proba(df)[:,1], wgt, channel, channels)
predict_data = xgb.predict_proba(df_data)[:,1]

PlotRootHisto(predict_sorted, 
              weights_sorted, 
              predict_data, 
              channels, 
              title = f"{signal}SearchDist", 
              xlabel = "XGB-Output", 
              bins = 30)
PlotRootHisto(predict_sorted, 
              weights_sorted, 
              predict_data, 
              channels, 
              title = f"{signal}SearchDist_C7", 
              xlabel = "XGB-Output", 
              bins = 30,
              CutOff = 0.7)

plotFI(xgb, df.keys(), signal)




