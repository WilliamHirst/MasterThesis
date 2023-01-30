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

name = "NN_oneMass"
signal = "SUSY"
train = False


print(f"Starting test: Model = {name} -- Signal = {signal}")

df, y, df_data, channels = loadDf(myPath, notInc=["ttbarHNLfull","LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75"])
dfPNN, df_dataPNN = AddParameters(df, y,df_data)

"""
Create one signal dataset
"""
indx = (y==0).to_numpy() + (df.channel == "MGPy8EGA14N23LOC1N2WZ750p0p0250p0p03L2L7").to_numpy()
df_1 = df[indx] 
y_1 = y[(y==0).to_numpy() + (df.channel == "MGPy8EGA14N23LOC1N2WZ750p0p0250p0p03L2L7").to_numpy()]

"""
Remove 4 signals which will be used to test.
"""
indx =  (df.channel!="MGPy8EGA14N23LOC1N2WZ700p0p0250p0p03L2L7").to_numpy() + (df.channel!="MGPy8EGA14N23LOC1N2WZ800p0p0250p0p03L2L7").to_numpy() + (df.channel!="MGPy8EGA14N23LOC1N2WZ750p0p0200p0p03L2L7").to_numpy()+(df.channel!="MGPy8EGA14N23LOC1N2WZ750p0p0300p0p03L2L7").to_numpy()
dfPNN = dfPNN[ indx] 
y = y[ indx] 

print("Preparing data....")
trainPNN, valPNN = splitAndPrepData(dfPNN, y, scale = True, ret_scaleFactor=True)#, PCA=True, n_components=1-1e-3)
train_1, val_1 = splitAndPrepData(df_1, y_1, scale = True, ret_scaleFactor=True)#, PCA=True, n_components=1-1e-3)
print("Done.")


X_trainPNN, Y_trainPNN, W_trainPNN, C_trainPNN = trainPNN
X_valPNN, Y_valPNN, W_valPNN, C_valPNN, scaleFactorPNN = valPNN
nrFeaturePNN = nFeats(X_trainPNN)



X_train_1, Y_train_1, W_train_1, C_train_1 = train_1
X_val_1, Y_val_1, W_val_1, C_val_1, scaleFactor_1 = val_1
nrFeature_1 = nFeats(X_train_1)



"""
PNN
"""
print("Compiling Model")
modelPNN = tf.keras.Sequential()
modelPNN.add(tf.keras.layers.InputLayer(input_shape=(nrFeaturePNN,)))
modelPNN.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
modelPNN.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
modelPNN.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
modelPNN.add(tf.keras.layers.Dense(1, activation="sigmoid"))

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
print(model_1.summary())
print("Done compiling.")

with tf.device("/GPU:0"):
    if train:
        callback = tf.keras.callbacks.EarlyStopping(monitor='val_auc', 
                                                    patience=10, 
                                                    restore_best_weights = True,
                                                    verbose = 1,
                                                    mode = "max")
        history = model_1.fit(X_trainPNN, 
                            Y_trainPNN,
                            sample_weight = W_trainPNN, 
                            epochs=100, 
                            batch_size=8192, 
                            callbacks = [callback],
                            validation_data=(X_valPNN, Y_valPNN, W_valPNN),
                            verbose = 1)

        history = modelPNN.fit(X_trainPNN, 
                            Y_trainPNN,
                            sample_weight = W_trainPNN, 
                            epochs=100, 
                            batch_size=8192, 
                            callbacks = [callback],
                            validation_data=(X_valPNN, Y_valPNN, W_valPNN),
                            verbose = 1)

        # pred_Train = model.predict(X_train, batch_size=8192)
        # pred_Val = model.predict(X_val, batch_size=8192)
    
    indx =  (df.channel=="MGPy8EGA14N23LOC1N2WZ700p0p0250p0p03L2L7").to_numpy() + (df.channel=="MGPy8EGA14N23LOC1N2WZ800p0p0250p0p03L2L7").to_numpy() + (df.channel=="MGPy8EGA14N23LOC1N2WZ750p0p0200p0p03L2L7").to_numpy()+(df.channel=="MGPy8EGA14N23LOC1N2WZ750p0p0300p0p03L2L7").to_numpy()

    df_test = df[indx + (y==0).to_numpy()]
    y_test = y[indx + (y==0).to_numpy()]

    W = df_test.wgt_SG
    C = df_test.channel
    df_testPNN = AddParameters(df_test, y)

    df_test = df_test.drop(columns = ["channel", "wgt_SG"])
    df_testPNN = df_testPNN.drop(columns = ["channel", "wgt_SG"])


    df_test = scaleData(df_test)
    df_testPNN = scaleData(df_testPNN)


    name = "NN_oneMass"
    HM(model_1, df_test, y_test, W, C, data = None, name = f"SUSY/{name}Grid", metric="Sig", save = True)


    name = "PNN_oneMass"
    HM(modelPNN, df_testPNN, y_test, W, C, data = None, name = f"SUSY/{name}Grid", metric="Sig", save = True)
    
    
    

    

# mc_predict, mc_weights = separateByChannel(prediction, weights, C, C.unique())

# saveLoad("results/predict_sorted_test.npy", mc_predict)
# saveLoad("results/weights_sorted_test.npy", mc_weights)
# saveLoad("results/predict_data_test.npy", model(df_data))
# saveLoad("results/channels_test.npy", C)