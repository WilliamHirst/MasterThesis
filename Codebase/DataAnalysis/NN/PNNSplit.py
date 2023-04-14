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


signal = "SUSY"
train = False
#notInc=["ttbarHNLfull","LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75", "p01p0"]
notInc=["ttbarHNLfull","LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75", "p01p0", "WZ100p0p0","WZ150p0p050p0", "WZ150p0p00p0p", "WZ200p0p00p0", "WZ200p0p050p0"] #FS_MLM

# split = [X2_min, X2_max, X1_min, X1_max] 
split1 = [[200, 350, 0, 250] ,[200, 400, 0, 250]] #Yellow
split2 = [[400, 600, 200, 400], [350, 650, 150, 400]] #Green
split3 = [[400, 600, 0, 150] ,[400, 650, 0, 250]] #Purple
split4 = [[650, 800, 200, 400],[550, 800, 150, 400]] #Cyan
split5 = [[650, 800, 0, 150],[550, 800, 0, 250]] #Red
split6 = [[200, 800, 0, 150],[200, 800, 0, 150]] #Cyan2
split7 = [[650, 800, 0, 150],[200, 600, 100, 350]] #Green2
split8 = [[650, 800, 0, 150],[300, 800, 150, 400]] #Red2
# splits = [split1,split2,split3,split4,split5]
splits = [split6,split7,split8]
models = ["PNNSplit1", "PNNSplit2", "PNNSplit3", "PNNSplit4", "PNNSplit5"]
models = ["PNNSplit6", "PNNSplit7", "PNNSplit8"]

for split, name in zip(splits,models):
    print(f"Starting test: Model = {name} -- Signal = {signal}")
    print(split)
    df, y, df_data, channels = loadDf(myPath, notInc=notInc,IncludeRange=split[1])
    df, df_data = AddParameters(df, y,df_data) # Only for PNN
    print(df.channel.unique())

    print("Preparing data....")
    train, val = splitAndPrepData(df, y, scale = True, ret_scaleFactor=True, PCA=True, n_components=1-1e-3) #Only for PCA
    print("Done.")

    X_train, Y_train, W_train, C_train = train
    X_val, Y_val, W_val, C_val, scaleFactor = val
    nrFeature = nFeats(X_train)


    # for channel in df.channel.unique():
    #     if checkIfInclude([],split[0],channel) and "MGPy8EGA14N23" in channel:
    #         y = y[df.channel != channel]
    #         df = df[df.channel != channel]
    # print(df.channel.unique())
    W = df.wgt_SG
    C = df.channel
    Y = y
    df = df.drop(columns = ["channel", "wgt_SG"])
    df, df_data = scaleData(df,df_data)
    df, df_data = PCAData(df, df_data, n_components=1-1e-3)


    print("Compiling Model")
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.InputLayer(input_shape=(nrFeature,)))
    model.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
    model.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
    model.add(tf.keras.layers.Dense(600, activation=tf.keras.layers.LeakyReLU(alpha=0.01)))
    model.add(tf.keras.layers.Dense(1,   activation="sigmoid"))

    # if not train:
    #     model.load_weights(f"models/PNNSplit/model_{name}.h5")

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

        model.save_weights(f"models/PNNSplit/model_{name}.h5")


        HM(model, df, Y, W, C, data = df_data, name = f"FS/{name}Grid", metric="Sig", save = False, saveTxt=True)