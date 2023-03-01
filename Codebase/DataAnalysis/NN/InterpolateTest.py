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

name = "PNN_oneMass"
signal = "SUSY"
train = False
notInc = ["ttbarHNLfull","LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75"]
IncludeRange = [450, 650, 150, 350]

print(f"Starting test: Model = {name} -- Signal = {signal}")


df, y, df_data, channels = loadDf(myPath, notInc = notInc, IncludeRange = IncludeRange)
dfPNN, df_dataPNN = AddParameters(df, y,df_data)

"""
Create one signal dataset
"""
# indx_1 = (y==0).to_numpy() + (df.channel == "MGPy8EGA14N23LOC1N2WZ750p0p0250p0p03L2L7").to_numpy()
# df_1 = df[indx_1] 
# y_1 = y[indx_1]

"""
Remove 4 signals which will be used to test.
"""
indx =  (df.channel!="MGPy8EGA14N23LOC1N2WZ700p0p0250p0p03L2L7").to_numpy() + (df.channel!="MGPy8EGA14N23LOC1N2WZ800p0p0250p0p03L2L7").to_numpy() + (df.channel!="MGPy8EGA14N23LOC1N2WZ750p0p0200p0p03L2L7").to_numpy()+(df.channel!="MGPy8EGA14N23LOC1N2WZ750p0p0300p0p03L2L7").to_numpy()
dfPNN = dfPNN[indx] 
yPNN = y[indx] 

# dfMaxOut = df[indx]

print("Preparing data....")
trainPNN, valPNN = splitAndPrepData(dfPNN, yPNN, scale = True, ret_scaleFactor=True)#, PCA=True, n_components=1-1e-3)
# train_1, val_1 = splitAndPrepData(df_1, y_1, scale = True, ret_scaleFactor=True)#, PCA=True, n_components=1-1e-3)
# train_MaxOut, val_MaxOut = splitAndPrepData(dfMaxOut, yPNN, scale = True, ret_scaleFactor=True)#, PCA=True, n_components=1-1e-3)
print("Done.")


X_trainPNN, Y_trainPNN, W_trainPNN, C_trainPNN = trainPNN
X_valPNN, Y_valPNN, W_valPNN, C_valPNN, scaleFactorPNN = valPNN
nrFeaturePNN = nFeats(X_trainPNN)



# X_train_1, Y_train_1, W_train_1, C_train_1 = train_1
# X_val_1, Y_val_1, W_val_1, C_val_1, scaleFactor_1 = val_1
# nrFeature_1 = nFeats(X_train_1)

# X_trainMaxOut, Y_trainMaxOut, W_trainMaxOut, C_trainMaxOut = train_MaxOut
# X_valMaxOut, Y_valMaxOut, W_valMaxOut, C_valMaxOut, scaleFactorMaxOut = val_MaxOut
# nrFeatureMaxOut = nFeats(X_trainMaxOut)

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
MaxOut
"""
# print("Compiling Model")
# modelMaxOut = tf.keras.Sequential()
# modelMaxOut.add(tf.keras.layers.InputLayer(input_shape=(nrFeatureMaxOut,)))
# modelMaxOut.add(tf.keras.layers.Dropout(0.15))
# modelMaxOut.add(MaxOut(units=600, num_inputs=nrFeatureMaxOut, num_groups=200))
# modelMaxOut.add(tf.keras.layers.Dropout(0.15))
# modelMaxOut.add(MaxOut(units=600, num_inputs=200, num_groups=200))
# modelMaxOut.add(tf.keras.layers.Dropout(0.15))
# modelMaxOut.add(MaxOut(units=600, num_inputs=200, num_groups=200))
# modelMaxOut.add(tf.keras.layers.Dense(1, activation="sigmoid"))

# """
# NN
# """
# model_1 = tf.keras.Sequential()
# model_1.add(tf.keras.layers.InputLayer(input_shape=(nrFeature_1,)))
# model_1.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
# model_1.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
# model_1.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
# model_1.add(tf.keras.layers.Dense(1, activation="sigmoid"))

optimizer = optimizers.Adam(learning_rate=1e-3)
# model_1.compile(loss="binary_crossentropy", optimizer=optimizer, weighted_metrics="AUC")
modelPNN.compile(loss="binary_crossentropy", optimizer=optimizer, weighted_metrics="AUC")
# modelMaxOut.compile(loss="binary_crossentropy", optimizer=optimizer, weighted_metrics="AUC")
print("Done compiling.")

with tf.device("/GPU:0"):
    callback = tf.keras.callbacks.EarlyStopping(monitor='val_auc', 
                                                patience=10, 
                                                restore_best_weights = True,
                                                verbose = 1,
                                                mode = "max")
    # history = model_1.fit(X_train_1, 
    #                       Y_train_1,
    #                       sample_weight = W_train_1, 
    #                       epochs=100, 
    #                       batch_size=8192, 
    #                       callbacks = [callback],
    #                       validation_data=(X_val_1, Y_val_1, W_val_1),
    #                       verbose = 1)

    history = modelPNN.fit(X_trainPNN, 
                           Y_trainPNN,
                           sample_weight = W_trainPNN, 
                           epochs=100, 
                           batch_size=8192, 
                           callbacks = [callback],
                           validation_data=(X_valPNN, Y_valPNN, W_valPNN),
                           verbose = 1)

    # history = modelMaxOut.fit(X_trainMaxOut, 
    #                           Y_trainMaxOut,
    #                           sample_weight = W_trainMaxOut, 
    #                           epochs=100, 
    #                           batch_size=8192, 
    #                           callbacks = [callback],
    #                           validation_data=(X_valMaxOut, Y_valMaxOut, W_valMaxOut),
    #                           verbose = 1)
    
    
    indx =  (df.channel=="MGPy8EGA14N23LOC1N2WZ700p0p0250p0p03L2L7").to_numpy() + (df.channel=="MGPy8EGA14N23LOC1N2WZ800p0p0250p0p03L2L7").to_numpy() + (df.channel=="MGPy8EGA14N23LOC1N2WZ750p0p0200p0p03L2L7").to_numpy()+(df.channel=="MGPy8EGA14N23LOC1N2WZ750p0p0300p0p03L2L7").to_numpy()+ (df.channel == "MGPy8EGA14N23LOC1N2WZ750p0p0250p0p03L2L7").to_numpy()

    df_test = df[indx + (y==0).to_numpy()]
    y_test = y[indx + (y==0).to_numpy()]

    W = df_test.wgt_SG
    C = df_test.channel
    df_testPNN = AddParameters(df_test, y_test)

    # df_test = df_test.drop(columns = ["channel", "wgt_SG"])
    df_testPNN = df_testPNN.drop(columns = ["channel", "wgt_SG"])


    # df_test = scaleData(df_test)
    df_testPNN = scaleData(df_testPNN)

    # name = "MaxOut_oneMass"
    # HM(modelMaxOut, df_test, y_test, W, C, data = None, name = f"SUSY/Interpolation/{name}Grid", metric="Sig", save = False)

    # name = "NN_oneMass"
    # HM(model_1, df_test, y_test, W, C, data = None, name = f"SUSY/Interpolation/{name}Grid", metric="Sig", save = True)


    name = "PNN_oneMass"
    HM(modelPNN, df_testPNN, y_test, W, C, data = None, name = f"SUSY/Interpolation/{name}Grid", metric="Sig", save = False)
    