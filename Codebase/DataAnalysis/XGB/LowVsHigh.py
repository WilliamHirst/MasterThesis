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
df, y, df_data, channels = loadDf(myPath, incHigh = True, signal = signal)


print("Preparing data....")
train, val = splitAndPrepData(df, y, split_v = 0.2, scaleWeight = True)
X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val = val
print("Done.")



xgb = XGB.XGBClassifier(
            max_depth=4, 
            n_estimators=100,
            learning_rate=0.1,
            n_jobs=4,
            tree_method="hist",
            objective='binary:logistic',
            missing=-999.0,
            use_label_encoder=False,
            eval_metric="error") 

highScores = []


for i in range(1,6):
    index = int(len(X_train)*i/5-1)
    time = timer()
    print("Training....")
    xgb = xgb.fit(X_train[:index], Y_train[:index], 
                  sample_weight = W_train[:index],
                  eval_set = [(X_val, Y_val)], 
                  sample_weight_eval_set = [W_val])
    print("Done")
    timer(time)
    score = plotRoc(Y_val, 
            xgb.predict_proba(X_val)[:,1], 
            W_val,
            "", 
            return_score = True, 
            name = f"",
            plot = False,
            )
    highScores.append(score)
import matplotlib.pyplot as plt
plt.plot(score)
plt.savefig("lowHig.pdf")
plt.show()
    







