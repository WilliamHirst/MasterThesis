from time import time
from sklearn.model_selection import train_test_split
import xgboost as XGB
from Plot_stuff.RHM import ROOT_Histo_Maker
from Plot_stuff.plot_set import *
from Plot_stuff.ROCM import *
import sys
sys.path.insert(1, "../../")
from Utilities import *


myPath = "/storage/William_Sakarias/William_Data"

signal = "topOther"

df, channels = loadDf(myPath, signal)

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




df_train, df_test = train_test_split(df, test_size=0.2)

y_train = df_train.isSignal
df_train_channel = df_train.channel
df_train_weights = df_train.wgt_SG

df_train = df_train.drop(columns=["isSignal", "channel", "wgt_SG"])

df_train_weights = removeNegWeights(df_train_weights)
df_train_weights = scaleWeights(df_train_weights, y_train)

time = timer()
print("Training....")
xgb = xgb.fit(df_train, y_train, sample_weight = df_train_weights.array)
print("Done")
timer(time)
prediction = xgb.predict_proba(df_train)[:,1]

predict_sorted, weights_sorted =  separateByChannel(prediction, df_train_weights, df_train_channel, channels)


ROOT_Histo_Maker(predict_sorted, 
                 weights_sorted,
                 channels, 
                 bin_max = 1, 
                 bin_min = 0,
                 nr_bins = 20, 
                 y_max= 1e5, 
                 y_min = 1e-2, 
                 variable_name = r"$XGB-Output$", 
                 saveAs = f"../../Figures/MLResults/XGB/{signal}SearchDist.pdf")


plotRoc(y_train, 
        prediction, 
        df_train_weights, 
        "Training-data", 
        return_score = True, 
        name = f"../../Figures/MLResults/XGB/{signal}ttbarSearchROC.pdf")