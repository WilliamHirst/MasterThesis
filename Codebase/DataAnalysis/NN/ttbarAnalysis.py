from time import time
import xgboost as XGB
import sys
sys.path.insert(1, "../")

from Plot_stuff.plot_set import *
from Plot_stuff.ROCM import *
from Plot_stuff.FI import *
import tensorflow as tf


sys.path.insert(1, "../../")
from Utilities import loadDf, saveLoad, splitAndPrepData, separateByChannel, timer

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


df, y, df_data, channels = loadDf(myPath, signal)

print("Preparing data....")
train, val = splitAndPrepData(df, y, scale = True)
print("Done.")

X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val = val

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
print("Done")
timer(time)

plotRoc(Y_train, 
        hypermodel.predict(X_train), 
        W_train,
        "Training-data", 
        plot = True,
        return_score = True, 
        name = f"../../../thesis/Figures/MLResults/NN/{signal}SearchROCTrain.pdf")

plotRoc(Y_val, 
        hypermodel.predict(X_val), 
        W_val,
        "Validation-data", 
        plot = True,
        return_score = True, 
        name = f"../../../thesis/Figures/MLResults/NN/{signal}SearchROCVal.pdf")
print("Plotted ROC-curves.")
"""
channel = df.channel
wgt = np.asarray(df.wgt_SG)
df = df.drop(columns=["wgt_SG","channel"])

predict_sorted, weights_sorted =  separateByChannel(hypermodel.predict(df), wgt, channel, channels)
predict_data = hypermodel.predict(df_data).ravel()
weights_sorted = np.asarray(weights_sorted, dtype = object)


saveLoad("predict_sorted_test.npy", predict_sorted)
saveLoad("weights_sorted_test.npy", weights_sorted)
saveLoad("predict_data_test.npy", predict_data)
saveLoad("channels_test.npy", channels)
"""






