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

def max_out(inputs, num_units = 200, axis=None):
    shape = inputs.get_shape().as_list()
    if shape[0] is None:
        shape[0] = -1
    if axis is None:  # Assume that channel is the last dimension
        axis = -1
    num_channels = shape[axis]
    if num_channels % num_units:
        raise ValueError('number of features({}) is not '
                         'a multiple of num_units({})'.format(num_channels, num_units))
    shape[axis] = num_units
    shape += [num_channels // num_units]
    outputs = tf.reduce_max(tf.reshape(inputs, shape), -1, keepdims=False)
    return outputs

def channel_out(inputs, num_units = 200, axis=None):
    shape = inputs.get_shape().as_list()
    if shape[0] is None:
        shape[0] = -1
    if axis is None:  # Assume that channel is the last dimension
        axis = -1
    num_channels = shape[axis]
    if num_channels % num_units:
        raise ValueError('number of features({}) is not '
                         'a multiple of num_units({})'.format(num_channels, num_units))
    shape[axis] = num_units
    shape += [num_channels // num_units]
    grouped = tf.reshape(inputs, shape)
    top_vals = tf.reduce_max(tf.reshape(inputs, shape), -1, keepdims=True)
    isMax = tf.reshape(tf.greater_equal(grouped, top_vals), [shape[0], num_channels])
    output = tf.multiply(tf.cast(isMax,tf.float32), inputs)

    return output



myPath = "/storage/William_Sakarias/William_Data"

signal = "ttbarHNLMaxOut"

print(f"Starting test: {signal}")

df, y, df_data, channels = loadDf(myPath, notInc=["LRS", "filtch"])

print("Preparing data....")
train, val = splitAndPrepData(df, y, scale = True, PCA = False)
print("Done.")

X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val = val
nrFeature = nFeats(X_train)

print("Compiling Model")
model = tf.keras.Sequential()
model.add(tf.keras.layers.InputLayer(input_shape=(nrFeature,)))
model.add(tf.keras.layers.Dense(600,activation = max_out))
model.add(tf.keras.layers.Dense(600,activation = max_out))
model.add(tf.keras.layers.Dense(600,activation = max_out))
model.add(tf.keras.layers.Dense(1, activation="sigmoid"))
optimizer = optimizers.Adam(learning_rate=1e-3)
model.compile(loss="binary_crossentropy", optimizer=optimizer, metrics=["AUC"])
print("Done compiling.")

with tf.device("/GPU:0"):
    callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', 
                                                patience=5, 
                                                restore_best_weights = True,
                                                verbose = 1)
    history = model.fit(X_train, 
                        Y_train,
                        sample_weight = W_train, 
                        epochs=15, 
                        batch_size=8096, 
                        callbacks = [callback],
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

state = True
while state == True:
    answ = input("Do you want to save model? (y/n) ")
    if answ == "y":
        name = "test_ExtraNodes"
        model.save(f"models/model_{name}.h5")
        model.save()
        model.save_weights(f'models/model_{name}_weights.h5')
        state = False
        print("Model saved")
    elif answ == "n":
        state = False
        print("Model not saved")

