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

category = "PNN"
name = "PNN_FS_MLM"
signal = "SUSY"

notInc=["ttbarHNLfull","LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75", "p01p0", "WZ100p0p0","WZ150p0p050p0", "WZ150p0p00p0p", "WZ200p0p00p0", "WZ200p0p050p0"]

print(f"Starting test: Model = {name} -- Signal = {signal}")

df, y, df_data, channels = loadDf(myPath, notInc=notInc)
df = AddParameters(df, y)


# df["param1"] = 3000 # 50
# df["param2"] = 2000 # 250


W = df.wgt_SG
C = df.channel
Y = y
df = df.drop(columns = ["channel", "wgt_SG"])
df  = scaleData(df)
print(df)
""" 
Setting the parameters manually to one signal comb.
"""
index = C == "MGPy8EGA14N23LOC1N2WZ250p050p03L2L7"
df_param1 = np.array(df[index]["param1"])
df_param2 = np.array(df[index]["param2"])
df["param1"] = df_param1[0]
df["param2"] = df_param2[0]
nrFeature = nFeats(df)
print(df)


print("Compiling Model")
model = tf.keras.Sequential()
model.add(tf.keras.layers.InputLayer(input_shape=(nrFeature,)))
model.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
model.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
model.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
model.add(tf.keras.layers.Dense(1,   activation="sigmoid"))

model.load_weights(f"models/{category}/model_{name}.h5")

optimizer = optimizers.Adam(learning_rate=1e-3)
model.compile(loss="binary_crossentropy", optimizer=optimizer, weighted_metrics="AUC")
print(model.summary())
print("Done compiling.")

with tf.device("/GPU:0"):
    C_u = C.unique()
    prediction = model.predict(df, batch_size=8192)
    mc_predict, mc_weights = separateByChannel(prediction, W, C, C_u)
    
saveLoad("results/PNNDistTest/predict_sorted_PNNV2.npy", mc_predict)
saveLoad("results/PNNDistTest/weights_sorted_PNNV2.npy", mc_weights)
saveLoad("results/PNNDistTest/channels_PNNV2.npy", C_u)