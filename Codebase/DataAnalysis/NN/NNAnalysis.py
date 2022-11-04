from time import time
import xgboost as XGB
import sys
sys.path.insert(1, "../")

from Plot_stuff.plot_set import *
from Plot_stuff.ROCM import *
from Plot_stuff.FI import *
import tensorflow as tf


sys.path.insert(1, "../../")
from Utilities import loadDf, splitAndPrepData, timer, separateByChannel, saveLoad

from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession

config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)


myPath = "/storage/William_Sakarias/William_Data"

signal = "ttbar"

name = "test"

hypermodel = tf.keras.models.load_model(f"models/model_{name}.h5")

print(hypermodel.summary())


df, y, df_data, channels = loadDf(myPath,notInc=["LRS", "filtch"])

print("Preparing data....")
train, val = splitAndPrepData(df, y, scale = True)
print("Done.")

X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val = val

channel = df.channel
wgt = np.asarray(df.wgt_SG)
df = df.drop(columns=["wgt_SG","channel"])

print(X_train)
time = timer()
print("Training....")
with tf.device("/GPU:0"):
        history = hypermodel.fit(X_train, 
                                 Y_train,
                                 sample_weight = W_train, 
                                 epochs=20, 
                                 batch_size=8096, 
                                 validation_data=(X_val, Y_val, W_val),
                                 verbose = 1)
        pred_Train = hypermodel.predict(X_train)
        pred_Val = hypermodel.predict(X_val)
        pred_MC = hypermodel.predict(df)
        pred_Data = hypermodel.predict(df_data).ravel()

print("Done")
timer(time)

plotRoc(Y_train, 
        pred_Train, 
        W_train,
        "Training-data", 
        plot = True,
        return_score = True, 
        name = f"../../../thesis/Figures/MLResults/NN/{signal}SearchROCTrain.pdf")

plotRoc(Y_val, 
        pred_Val, 
        W_val,
        "Validation-data", 
        plot = True,
        return_score = True, 
        name = f"../../../thesis/Figures/MLResults/NN/{signal}SearchROCVal.pdf")
print("Plotted ROC-curves.")



predict_sorted, weights_sorted =  separateByChannel(pred_MC , wgt, channel, channels)
weights_sorted = np.asarray(weights_sorted, dtype = object)

print(predict_sorted,"\n")
print(weights_sorted,"\n")
print(pred_Data,"\n")
print(channels,"\n")

saveLoad("predict_sorted_test.npy", predict_sorted)
saveLoad("weights_sorted_test.npy", weights_sorted)
saveLoad("predict_data_test.npy", pred_Data)
saveLoad("channels_test.npy", channels)







