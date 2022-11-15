import tensorflow as tf
from tensorflow.keras import optimizers

from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession

from keras.layers import Activation
from keras import backend as K
from keras.utils.generic_utils import get_custom_objects



config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)

import sys
sys.path.insert(1, "../")

from Plot_stuff.plot_set import *
from Plot_stuff.ROCM import *

sys.path.insert(1, "../../")
from Utilities import *

def MaxChannel(input):
        for i in range(0,tf.size(input).eval(),3):
                max = input[i]
                indx = i
                for j in range(i+1,i+3):
                        if input(j) > max:
                                input[indx] = 0
                                max = input[j]
                                indx = j
                        else:
                                input[j] = 0
     
get_custom_objects().update({'MaxChannel': Activation(MaxChannel)})



myPath = "/storage/William_Sakarias/William_Data"

signal = "ttbarHNLMaxChannel"

df, y, df_data, channels = loadDf(myPath, notInc=["LRS", "filtch"])

print("Preparing data....")
train, val = splitAndPrepData(df, y, scale = True, PCA = False)
print("Done.")

X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val = val
nrFeature = nFeats(X_train)


model = tf.keras.Sequential()
model.add(tf.keras.layers.InputLayer(input_shape=(nrFeature,)))
model.add(tf.keras.layers.Dense(600,activation = MaxChannel))
model.add(tf.keras.layers.Dense(600,activation = MaxChannel))
model.add(tf.keras.layers.Dense(600,activation = MaxChannel))
model.add(tf.keras.layers.Dense(1, activation="sigmoid"))
optimizer = optimizers.Adam(learning_rate=1e-3)
model.compile(loss="binary_crossentropy", optimizer=optimizer, metrics=["AUC"])


with tf.device("/GPU:0"):
    callback = tf.keras.callbacks.EarlyStopping(monitor='val_auc', 
                                                patience=10, 
                                                restore_best_weights = True,
                                                verbose = 1)
    history = model.fit(X_train, 
                        Y_train,
                        sample_weight = W_train, 
                        epochs=15, 
                        batch_size=8096, 
                        #callbacks = [callback],
                        validation_data=(X_val, Y_val, W_val),
                        verbose = 1)
    pred_Train = model.predict(X_train)
    pred_Val = model.predict(X_val)

plotRoc(Y_train, 
        pred_Train, 
        W_train,
        "Training-data", 
        plot = True,
        return_score = True, 
        name = f"DNN/{signal}SearchROCTrain.pdf")

plotRoc(Y_val, 
        pred_Val, 
        W_val,
        "Validation-data", 
        plot = True,
        return_score = True, 
        name = f"DNN/{signal}SearchROCVal.pdf")

