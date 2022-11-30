import tensorflow as tf
from tensorflow.keras import optimizers
from tensorflow.keras.constraints import max_norm


from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession

from CustomeObjects import max_out, channel_out, global_max_out




config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)

import sys
sys.path.insert(1, "../")

from Plot_stuff.plot_set import *
from Plot_stuff.ROCM import *

sys.path.insert(1, "../../")
from Utilities import *




myPath = "/storage/William_Sakarias/William_Data"

signal = "ttbarHNLMaxChannel"

print(f"Starting test: {signal}")

df, y, df_data, channels = loadDf(myPath, notInc=["LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75"])

print("Preparing data....")
train, val = splitAndPrepData(df, y, scale = True, ret_scaleFactor=True, PCA=True, n_components=1-1e-2)
print("Done.")

X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val, scaleFactor = val
nrFeature = nFeats(X_train)

print("Compiling Model")
model = tf.keras.Sequential()
model.add(tf.keras.layers.InputLayer(input_shape=(nrFeature,)))
#model.add(tf.keras.layers.Dropout(0.2))
model.add(tf.keras.layers.Dense(600,activation = max_out))
#model.add(tf.keras.layers.Dropout(0.2))
model.add(tf.keras.layers.Dense(600,activation = max_out))
#model.add(tf.keras.layers.Dropout(0.2))
model.add(tf.keras.layers.Dense(600,activation = max_out))
#model.add(tf.keras.layers.Dropout(0.2))
model.add(tf.keras.layers.Dense(1, activation="sigmoid"))

optimizer = optimizers.Adam(learning_rate=1e-3)
model.compile(loss="binary_crossentropy", optimizer=optimizer, weighted_metrics="AUC")
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
                        epochs=100, 
                        batch_size=8096, 
                        callbacks = [callback], #, CC],
                        validation_data=(X_val, Y_val, W_val),
                        verbose = 1)
    pred_Train = model.predict(X_train, batch_size=8096)
    pred_Val = model.predict(X_val, batch_size=8096)
    Calc_Sig(Y_val, pred_Val, W_val/scaleFactor)
    
    print(f"Optimal Validation AUC: {np.max(history.history['val_auc']):.5}")

    plotRoc(Y_val, 
            pred_Val, 
            W_val,
            "",
            plot = False)




exit()
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
        model.save_weights(f'models/model_{name}_weights.h5')
        state = False
        print("Model saved")
    elif answ == "n":
        state = False
        print("Model not saved")



# tf.keras.layers.LeakyReLU(alpha=0.01)
# tf.keras.layers.Dropout(alpha=0.01)