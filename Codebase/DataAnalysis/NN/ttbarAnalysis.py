from time import time
import xgboost as XGB
import sys
sys.path.insert(1, "../")

from Plot_stuff.plot_set import *
from Plot_stuff.ROCM import *
from Plot_stuff.FI import *
import tensorflow as tf


sys.path.insert(1, "../../")
from Utilities import loadDf, saveLoad, splitData, separateByChannel, timer


myPath = "/storage/William_Sakarias/William_Data"

signal = "ttbar"

name = "test"

hypermodel = tf.keras.models.load_model(f"models/model_{name}.h5")


df, y, df_data, channels = loadDf(myPath)
train, val, test = splitData(df, y)

X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val = val
X_test, Y_test, W_test, C_test = test

time = timer()
print("Training....")
history = hypermodel.fit(X_train, 
                         Y_train,
                         sample_weight = W_train, 
                         epochs=200, 
                         batch_size=int(len(X_train)/30), 
                         validation_data=(X_val, Y_val, W_val))
print("Done")
timer(time)

plotRoc(Y_train, 
        hypermodel.predict(X_train), 
        W_train,
        "Training-data", 
        return_score = True, 
        name = f"../../../thesis/Figures/MLResults/NN/{signal}SearchROCTrain.pdf")

plotRoc(Y_val, 
        hypermodel.predict(X_val), 
        W_val,
        "Validation-data", 
        return_score = True, 
        name = f"../../../thesis/Figures/MLResults/NN/{signal}SearchROCVal.pdf")

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






