import sys
import tensorflow as tf
tf.random.set_seed(42)
from tensorflow.keras import optimizers
from tensorflow.compat.v1 import ConfigProto, InteractiveSession


config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)

sys.path.insert(1, "../")
from Plot_stuff.plot_set import *
from Plot_stuff.ROCM import *
from Plot_stuff.HM import *
from Plot_stuff.THP import *

sys.path.insert(1, "../../")
from Utilities import *

myPath = "/storage/William_Sakarias/William_Data"

name = "PNN_FS_MLM"
signal = "SUSY"
#notInc=["ttbarHNLfull","LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75", "p01p0"]
notInc=["ttbarHNLfull","LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75", "p01p0", "WZ100p0p0","WZ150p0p050p0", "WZ150p0p00p0p", "WZ200p0p00p0", "WZ200p0p050p0"]



print(f"Starting test: Model = {name} -- Signal = {signal}")

df, y, df_data, channels = loadDf(myPath, notInc=notInc)
indx = (y ==0).to_numpy() + (df.channel=="MGPy8EGA14N23LOC1N2WZ250p050p03L2L7").to_numpy() + (df.channel=="MGPy8EGA14N23LOC1N2WZ350p0p00p0p03L2L7").to_numpy() + (df.channel=="MGPy8EGA14N23LOC1N2WZ300p0p0200p0p03L2L7").to_numpy()+(df.channel=="MGPy8EGA14N23LOC1N2WZ500p0p0100p0p03L2L7").to_numpy()

df = df[indx]
y = y[indx]
print(df[y==1].channel.unique())
df["param1"] = 250 
df["param2"] = 50 

W = df.wgt_SG
C = df.channel
Y = y
df = df.drop(columns = ["channel", "wgt_SG"])
df  = scaleData(df)
nrFeature = nFeats(df)


print("Compiling Model")
model = tf.keras.Sequential()
model.add(tf.keras.layers.InputLayer(input_shape=(nrFeature,)))
model.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
model.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
model.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
model.add(tf.keras.layers.Dense(1,   activation="sigmoid"))

model.load_weights(f"models/model_{name}.h5")

optimizer = optimizers.Adam(learning_rate=1e-3)
model.compile(loss="binary_crossentropy", optimizer=optimizer, weighted_metrics="AUC")
print(model.summary())
print("Done compiling.")

with tf.device("/GPU:0"):
    df = df
    W = W
    C = C
    C_u = C.unique()
    prediction = model.predict(df, batch_size=8192)
    mc_predict, mc_weights = separateByChannel(prediction, W, C, C_u)

    prediction = model.predict(df)
    weights = W


    
    
    

    

mc_predict, mc_weights = separateByChannel(prediction, W, C, C.unique())

saveLoad("results/predict_sorted_PNN.npy", mc_predict)
saveLoad("results/weights_sorted_PNN.npy", mc_weights)
saveLoad("results/predict_data_PNN.npy", model(df_data))
saveLoad("results/channels_PNN.npy", C)