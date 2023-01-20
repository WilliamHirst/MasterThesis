import sys
import tensorflow as tf
tf.random.set_seed(42)
from tensorflow.keras import optimizers,regularizers
from tensorflow.compat.v1 import ConfigProto, InteractiveSession

config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)
from layers import MaxOut, ChannelOut, StochChannelOut

sys.path.insert(1, "../")
from Plot_stuff.plot_set import *
from Plot_stuff.ROCM import *
from Plot_stuff.HM import *

sys.path.insert(1, "../../")
from Utilities import *

myPath = "/storage/William_Sakarias/William_Data"

signal = "SUSY"

print(f"Starting test: {signal}")

df, y, df_data, channels = loadDf(myPath, notInc=["ttbarHNLfull","LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75"])

print("Preparing data....")
train, val = splitAndPrepData(df, y, scale = True, ret_scaleFactor=True)#, PCA=True, n_components=1-1e-2)
print("Done.")

X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val, scaleFactor = val
nrFeature = nFeats(X_train)


print("Compiling Model")
model = tf.keras.Sequential()
model.add(tf.keras.layers.InputLayer(input_shape=(nrFeature,)))
model.add(tf.keras.layers.Dropout(0.15))
model.add(MaxOut(units=600, num_inputs=nrFeature, num_groups=200))
model.add(tf.keras.layers.Dropout(0.15))
model.add(MaxOut(units=600, num_inputs=200, num_groups=200))
model.add(tf.keras.layers.Dropout(0.15))
model.add(MaxOut(units=600, num_inputs=200, num_groups=200))
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
                        callbacks = [callback],
                        validation_data=(X_val, Y_val, W_val),
                        verbose = 1)
    pred_Train = model.predict(X_train, batch_size=8192)
    pred_Val = model.predict(X_val, batch_size=8192)

    
    W = df.wgt_SG
    C = df.channel
    Y = y
    df = df.drop(columns = ["channel", "wgt_SG"])
    df, df_data = scaleData(df,df_data)
    
    HM(model, df, Y, W, C, data = None, name = "../../../thesis/Figures/MLResults/NN/SUSY/MaxOutGrid", metric="Sig")

# mc_predict, mc_weights = separateByChannel(prediction, weights, C, C.unique())

# saveLoad("results/predict_sorted_test.npy", mc_predict)
# saveLoad("results/weights_sorted_test.npy", mc_weights)
# saveLoad("results/predict_data_test.npy", model(df_data))
# saveLoad("results/channels_test.npy", C)