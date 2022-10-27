import xgboost as XGB

import sys
sys.path.insert(1, "../")
from Plot_stuff.plot_set import *
from Plot_stuff.ROCM import *
from Plot_stuff.FI import *
from Plot_stuff.ROOTPlot import *

sys.path.insert(1, "../../")
from Utilities import *
from sklearn.model_selection import train_test_split


myPath = "/storage/William_Sakarias/William_Data"

siglist = ["LRSMWR2400NR50",
          "WeHNL5040Glt01ddlepfiltch1",
          "WeHNL5060Glt01ddlepfiltch1",
          "WeHNL5070Glt01ddlepfiltch1",
          "WmuHNL5040Glt01ddlepfiltch1",
          "LRSMWR4500NR400",
          "WmuHNL5070Glt01ddlepfiltch1"]


signal = "ttbar"

df, y, df_data, channels = loadDf(myPath, incHigh = True)

X_train, X_test, y_train, y_test = train_test_split(df, y, test_size=0.2, random_state=42)

print(np.sum(X_test[X_test.channel == "Zeejets"].wgt_SG)/np.sum(X_test.wgt_SG),
      np.sum(df[df.channel == "Zeejets"].wgt_SG)/np.sum(df.wgt_SG)   
)
print(np.sum(X_test[X_test.channel == "ttbar"].wgt_SG)/np.sum(X_test.wgt_SG),
      np.sum(df[df.channel == "ttbar"].wgt_SG)/np.sum(df.wgt_SG)   
)
print(np.sum(X_test[X_test.channel == "Wjets"].wgt_SG)/np.sum(X_test.wgt_SG),
      np.sum(df[df.channel == "Wjets"].wgt_SG)/np.sum(df.wgt_SG)
)
print(np.sum(X_test[X_test.channel == "higgs"].wgt_SG)/np.sum(X_test.wgt_SG),
      np.sum(df[df.channel == "higgs"].wgt_SG)/np.sum(df.wgt_SG)
)


exit()

y = df["channel"] == signal

xgb = XGB.XGBClassifier(
            max_depth=4, 
            n_estimators=150,
            learning_rate=0.1,
            n_jobs=4,
            tree_method="hist",
            objective='binary:logistic',
            missing=-999.0,
            use_label_encoder=False,
            eval_metric="error") 



print("Preparing data....")
train, val, test = splitData(df, y)
print("Done.")

highScores = []



X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val = val
X_test, Y_test, W_test, C_test = test

print(len(X_train[:,:,0]))
print(len(X_train[:,:,1]))
print(len(X_train[:,:,2]))
print(len(X_train[:,:,3]))
print(len(X_train[:,:,4]))
exit()
for i in range(5):
    time = timer()
    print("Training....")
    xgb = xgb.fit(X_train[:,:,i], Y_train[:,:,i], 
                  sample_weight = W_train[:,:,i],
                  eval_set = [(X_val[:,:,i], Y_val[:,:,i])], 
                  sample_weight_eval_set = [W_val[:,:,i]])
    print("Done")
    timer(time)
    score = plotRoc(Y_val, 
            xgb.predict_proba(X_val[:,:,i])[:,1], 
            W_val,
            "", 
            return_score = True, 
            name = f"",
            plot = False,
            )
    







