import sys
import tensorflow as tf
tf.random.set_seed(42)
from tensorflow.keras import optimizers,regularizers


from tensorflow.compat.v1 import ConfigProto, InteractiveSession

from layers import MaxOut, ChannelOut, StochChannelOut

sys.path.insert(1, "../")
from Plot_stuff.plot_set import *
from Plot_stuff.ROCM import *
from Plot_stuff.HM import *

sys.path.insert(1, "../../")
from Utilities import *


config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)



myPath = "/storage/William_Sakarias/William_Data"

signal = "SUSY"

print(f"Starting test: {signal}")

df, y, df_data, channels = loadDf(myPath, notInc=["ttbarHNLfull","LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75"])

#df, df_data = AddParameters(df, y, df_data)


print("Preparing data....")
train, val = splitAndPrepData(df, y, scale = True, ret_scaleFactor=True)#, PCA=True, n_components=1-1e-2)
print("Done.")

# train, val = loadSamples()

X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val, scaleFactor = val
nrFeature = nFeats(X_train)


print("Compiling Model")
model = tf.keras.Sequential()
model.add(tf.keras.layers.InputLayer(input_shape=(nrFeature,)))
model.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
model.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
model.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
model.add(tf.keras.layers.Dense(1, activation="sigmoid"))

optimizer = optimizers.Adam(learning_rate=1e-3)
model.compile(loss="binary_crossentropy", optimizer=optimizer, weighted_metrics="AUC")
print(model.summary())
print("Done compiling.")

with tf.device("/GPU:0"):
    callback = tf.keras.callbacks.EarlyStopping(monitor='val_auc', 
                                                patience=10, 
                                                restore_best_weights = True,
                                                verbose = 1,
                                                mode = "max")
    history = model.fit(X_train, 
                        Y_train,
                        sample_weight = W_train, 
                        epochs=1, 
                        batch_size=8192, 
                        callbacks = [callback], #, CC],
                        validation_data=(X_val, Y_val, W_val),
                        verbose = 1)
    pred_Train = model.predict(X_train, batch_size=8192)
    pred_Val = model.predict(X_val, batch_size=8192)
    #Calc_Sig(Y_val, pred_Val, W_val/scaleFactor)
    
    #print(f"Optimal Validation AUC: {np.max(history.history['val_auc']):.5}")

    plotRoc(Y_val, 
            pred_Val, 
            W_val,
            "",
            plot = False)
    HM(model, X_val, Y_val, W_val, C_val, name = "../../../thesis/Figures/MLResults/NN/SUSY/NNGrid.pdf")
