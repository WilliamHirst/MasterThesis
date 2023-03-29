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
from Plot_stuff.THP import *


sys.path.insert(1, "../../")
from Utilities import *

myPath = "/storage/William_Sakarias/William_Data"

name = "NN_Interpolation"
signal = "SUSY"

notInc = ["ttbarHNLfull","LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75"]
IncludeRange = [450, 650, 150, 350] 
TestMasses = {"500": ["200", "250","300"], "550": ["200", "300"], "600": ["200", "250","300"]} 

print(f"Starting test: Model = {name} -- Signal = {signal}")


df, y, df_data, channels = loadDf(myPath, notInc = notInc, IncludeRange = IncludeRange)

"""
Create one signal dataset
"""
# indx_m = (df.channel == "MGPy8EGA14N23LOC1N2WZ550p0250p03L2L7").to_numpy() 
# indx_1 = (y==0).to_numpy() + indx_m

# df_1 = df[indx_1] 
# y_1 = y[indx_1]

# print(df_1.channel.unique())
"""
Remove signals which will be used to test.
"""
indx_b = (y==0).to_numpy()
indx = indx_b.copy()
for chan in df[y == 1].channel.unique():
    m1, m2 = getMass(chan)
    if m1 in list(TestMasses.keys()):
        if m2 in TestMasses[m1]:
            print(m1,m2)
            continue
    indx += (df.channel == chan).to_numpy()

df_R = df[indx]
y_R = y[indx]

print(df_R.channel.unique())



print("Preparing data....")
train_1, val_1 = splitAndPrepData(df_R, y_R, scale = True, ret_scaleFactor=True)
print("Done.")

X_train_1, Y_train_1, W_train_1, C_train_1 = train_1
X_val_1, Y_val_1, W_val_1, C_val_1, scaleFactor_1 = val_1
nrFeature_1 = nFeats(X_train_1)


"""
NN
"""
model_1 = tf.keras.Sequential()
model_1.add(tf.keras.layers.InputLayer(input_shape=(nrFeature_1,)))
model_1.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
model_1.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
model_1.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
model_1.add(tf.keras.layers.Dense(1, activation="sigmoid"))

optimizer = optimizers.Adam(learning_rate=1e-3)
model_1.compile(loss="binary_crossentropy", optimizer=optimizer, weighted_metrics="AUC")
print("Done compiling.")

with tf.device("/GPU:0"):
    
    callback_1 = tf.keras.callbacks.EarlyStopping(monitor='val_auc', 
                                                patience=10, 
                                                restore_best_weights = True,
                                                verbose = 1,
                                                mode = "max")
    history = model_1.fit(X_train_1, 
                          Y_train_1,
                          sample_weight = W_train_1, 
                          epochs=20, 
                          batch_size=8192, 
                          callbacks = [callback_1],
                          validation_data=(X_val_1, Y_val_1, W_val_1),
                          verbose = 1)

    THP(history, name, signal)
    exit()
    df_test = df
    y_test = y 

    W = df_test.wgt_SG
    C = df_test.channel
    df_test = df_test.drop(columns = ["channel", "wgt_SG"])


    df_test = scaleData(df_test)


    name = "NN_OneMass_Overfitting15_Interpolation"
    HM(model_1, df_test, y_test, W, C, data = None, name = f"Interpolation/{name}Grid", metric="Sig", save = True)
    THP(history, name, signal)

    