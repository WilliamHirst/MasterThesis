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

name = "PNNPCA_FS_MLM"
signal = "SUSY"
train = False
#notInc=["ttbarHNLfull","LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75", "p01p0"]
notInc=["ttbarHNLfull","LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75", "p01p0", "WZ100p0p0","WZ150p0p050p0", "WZ150p0p00p0p", "WZ200p0p00p0", "WZ200p0p050p0"]



print(f"Starting test: Model = {name} -- Signal = {signal}")

df, y, df_data, channels = loadDf(myPath, notInc=notInc)
df, df_data = AddParameters(df, y,df_data)

if train:
    print("Preparing data....")
    train, val = splitAndPrepData(df, y, scale = True, ret_scaleFactor=True, PCA=True, n_components=1-1e-3)
    print("Done.")

    X_train, Y_train, W_train, C_train = train
    X_val, Y_val, W_val, C_val, scaleFactor = val
    nrFeature = nFeats(X_train)
else:
    W = df.wgt_SG
    C = df.channel
    Y = y
    df = df.drop(columns = ["channel", "wgt_SG"])
    df, df_data = scaleData(df,df_data)
    df = PCAData(df, n_components=1-1e-3)
    nrFeature = nFeats(df)


print("Compiling Model")
model = tf.keras.Sequential()
model.add(tf.keras.layers.InputLayer(input_shape=(nrFeature,)))
model.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
model.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
model.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
model.add(tf.keras.layers.Dense(1,   activation="sigmoid"))

if not train:
    model.load_weights(f"models/model_{name}.h5")

optimizer = optimizers.Adam(learning_rate=1e-3)
model.compile(loss="binary_crossentropy", optimizer=optimizer, weighted_metrics="AUC")
print(model.summary())
print("Done compiling.")

with tf.device("/GPU:0"):
    if train:
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

        model.save_weights(f"models/model_{name}.h5")
        #THP(history=history, model = name ,signal = signal )


    else: 
        HM(model, df, Y, W, C, data = None, name = f"FS/{name}Grid", metric="Sig", save = True)
    
    
    

    

# mc_predict, mc_weights = separateByChannel(prediction, weights, C, C.unique())

# saveLoad("results/predict_sorted_test.npy", mc_predict)
# saveLoad("results/weights_sorted_test.npy", mc_weights)
# saveLoad("results/predict_data_test.npy", model(df_data))
# saveLoad("results/channels_test.npy", C)