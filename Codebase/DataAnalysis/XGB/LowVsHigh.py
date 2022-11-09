import xgboost as XGB

import sys
sys.path.insert(1, "../")
from Plot_stuff.plot_set import *
from Plot_stuff.ROCM import *
from Plot_stuff.FI import *
from Plot_stuff.ROOTPlot import *

sys.path.insert(1, "../../")
from Utilities import *


myPath = "/storage/William_Sakarias/William_Data"

siglist = ["LRSMWR2400NR50",
          "WeHNL5040Glt01ddlepfiltch1",
          "WeHNL5060Glt01ddlepfiltch1",
          "WeHNL5070Glt01ddlepfiltch1",
          "WmuHNL5040Glt01ddlepfiltch1",
          "LRSMWR4500NR400",
          "WmuHNL5070Glt01ddlepfiltch1"]


signal = "ttbar"
df_h, y_h, df_data_h, channels_h = loadDf(myPath, incHigh = True, signal = signal)
df_l, y_l, df_data_l, channels_l = loadDf(myPath, incHigh = False, signal = signal)


print("Preparing data....")
train_h, val_h = splitAndPrepData(df_h, y_h, split_v = 0.2, scaleWeight = True)
X_train_h, Y_train_h, W_train_h, C_train_h = train_h
X_val_h, Y_val_h, W_val_h, C_val_h = val_h

train_l, val_l = splitAndPrepData(df_l, y_l, split_v = 0.2, scaleWeight = True)
X_train_l, Y_train_l, W_train_l, C_train_l = train_l
X_val_l, Y_val_l, W_val_l, C_val_l = val_l
print("Done.")



xgbShallow = XGB.XGBClassifier(
            max_depth=3, 
            n_estimators=30,
            learning_rate=0.1,
            n_jobs=4,
            tree_method="hist",
            objective='binary:logistic',
            missing=-999.0,
            use_label_encoder=False,
            eval_metric="error") 

xgbDeep = XGB.XGBClassifier(
            max_depth=3, 
            n_estimators=75,
            learning_rate=0.1,
            n_jobs=4,
            tree_method="hist",
            objective='binary:logistic',
            missing=-999.0,
            use_label_encoder=False,
            eval_metric="error") 

highScoresDeep = []
highScoresShallow = []
N = 10

for i in range(1,N):
    index = int(len(X_train_h)*i/N-1)
    time = timer()
    print("Training....")
    xgbDeepHist = xgbDeep.fit(X_train_l[:index], Y_train_l[:index], 
                  sample_weight = W_train_l[:index],
                  eval_set = [(X_val_l, Y_val_l)], 
                  sample_weight_eval_set = [W_val_l])
    xgbShallowHist = xgbShallow.fit(X_train_h[:index], Y_train_h[:index], 
                  sample_weight = W_train_h[:index],
                  eval_set = [(X_val_h, Y_val_h)], 
                  sample_weight_eval_set = [W_val_h])
    print("Done")
    timer(time)
    scoreDeep = plotRoc(Y_val_l, 
            xgbDeepHist.predict_proba(X_val_l)[:,1], 
            W_val_l,
            "", 
            return_score = True, 
            name = f"",
            plot = False,
            )
    scoreShallow = plotRoc(Y_val_h, 
            xgbShallowHist.predict_proba(X_val_h)[:,1], 
            W_val_h,
            "", 
            return_score = True, 
            name = f"",
            plot = False,
            )
    highScoresDeep.append(scoreDeep)
    highScoresShallow.append(scoreShallow)
    print(f"Completed: {i/N*100 :.2f} %")
import matplotlib.pyplot as plt
plt.plot(highScoresDeep, label = "Deep")
plt.plot(highScoresShallow, label = "Shallow")
plt.legend()
plt.show()
plt.savefig("lowHig.pdf")
plt.show()
    







