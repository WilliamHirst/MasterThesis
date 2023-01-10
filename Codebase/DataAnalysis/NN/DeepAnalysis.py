import sys
import tensorflow as tf
tf.random.set_seed(42)
from tensorflow.keras import optimizers,regularizers


from tensorflow.compat.v1 import ConfigProto, InteractiveSession

from layers import MaxOut, ChannelOut, StochChannelOut

sys.path.insert(1, "../")
from Plot_stuff.plot_set import *
from Plot_stuff.ROCM import *

sys.path.insert(1, "../../")
from Utilities import *


config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)



myPath = "/storage/William_Sakarias/William_Data"

signal = "ttbarHNLMaxChannel"

print(f"Starting test: {signal}")

df, y, df_data, channels = loadDf(myPath, notInc=["ttbarHNLfull","LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75"])

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
model.add(tf.keras.layers.Dropout(0.2))
model.add(MaxOut(units=600, num_inputs=nrFeature, num_groups=200))
model.add(tf.keras.layers.Dropout(0.2))
model.add(MaxOut(units=600, num_inputs=200, num_groups=200))
model.add(tf.keras.layers.Dropout(0.2))
model.add(MaxOut(units=600, num_inputs=200, num_groups=200))
model.add(tf.keras.layers.Dropout(0.2))
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
                        epochs=100, 
                        batch_size=8192, 
                        callbacks = [callback], #, CC],
                        validation_data=(X_val, Y_val, W_val),
                        verbose = 1)
    pred_Train = model.predict(X_train, batch_size=8192)
    pred_Val = model.predict(X_val, batch_size=8192)
    Calc_Sig(Y_val, pred_Val, W_val/scaleFactor)
    
    print(f"Optimal Validation AUC: {np.max(history.history['val_auc']):.5}")

    plotRoc(Y_val, 
            pred_Val, 
            W_val,
            "",
            plot = False)


# saveLoad("samples/X_train_sample.hdf5", X_train[:100000], type = "Pandas")
# saveLoad("samples/Y_train_sample.hdf5", Y_train[:100000], type = "Pandas")
# saveLoad("samples/W_train_sample.hdf5", W_train[:100000], type = "Pandas")
# saveLoad("samples/C_train_sample.hdf5", C_train[:100000], type = "Pandas")

# saveLoad("samples/X_val_sample.hdf5", X_val[:100000], type = "Pandas")
# saveLoad("samples/Y_val_sample.hdf5", Y_val[:100000], type = "Pandas")
# saveLoad("samples/W_val_sample.hdf5", W_val[:100000], type = "Pandas")
# saveLoad("samples/C_val_sample.hdf5", C_val[:100000], type = "Pandas")



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