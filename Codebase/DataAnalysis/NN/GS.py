import tensorflow as tf
from tensorflow.keras import optimizers
#Added section to remove error. Test later if can remove.
from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession

config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)

import keras_tuner as kt
import sys
sys.path.insert(1, "../../")
from Utilities import *



myPath = "/storage/William_Sakarias/William_Data"

procUsed = int(100/10)

signal = "ttbarHNL"

df, y, df_data, channels = loadDf(myPath, notInc=["LRS", "filtch"])

print("Preparing data....")
train, val = splitAndPrepData(df, y, scale = True, PCA = False)#, n_components = 1. - 1e-3)
print("Done.")

X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val = val

nrFeature = nFeats(X_train)
class MyHyperModel(kt.HyperModel):
    def build(self, hp):
        model = tf.keras.Sequential()
        model.add(tf.keras.layers.InputLayer(input_shape=(nrFeature,)))
        model.add(tf.keras.layers.Dropout( hp.Choice("drop_out", values=[0.0, .1,.2]),))
        for i in range(hp.Int('num_layers', 2, 5)):
            act = hp.Choice('act_' + str(i), 
                            ['relu',
                             'leaky_relu',
                             'sigmoid', 
                             'tanh'])

            if act == "leaky_relu":
                act = tf.keras.layers.LeakyReLU(alpha=0.01)

            model.add(tf.keras.layers.Dense(units=hp.Int('units_' + str(i),
                                                min_value=10,
                                                max_value=30,
                                                step=5),
                               activation = act))

        model.add(tf.keras.layers.Dense(1, activation="sigmoid"))
        hp_learning_rate = hp.Choice("learning_rate", values=[1e-2, 5e-3, 1e-3, 5e-4])
        optimizer = optimizers.Adam(learning_rate=hp_learning_rate)

        model.compile(loss="binary_crossentropy", optimizer=optimizer, metrics=["AUC"])
        return model

start_time = timer(None)




tuner = kt.Hyperband(
    MyHyperModel(),
    objective=kt.Objective('val_auc', direction='max'),
    max_epochs=50,
    factor=3,
    directory="GridSearches",
    project_name="NN",
    overwrite=True,
)
callback = tf.keras.callbacks.EarlyStopping(monitor='AUC', patience=3)
with tf.device("/GPU:0"):
    tuner.search(X_train[::procUsed], 
                    Y_train[::procUsed], 
                    sample_weight = W_train[::procUsed], 
                    batch_size =8096, 
                    epochs=50, 
                    validation_data=(X_val[::procUsed],
                                    Y_val[::procUsed],
                                    W_val[::procUsed]),
                    #callbacks = [callback]
                )
timer(start_time)
best_hps = tuner.get_best_hyperparameters(num_trials=1)[0]


print(f"The hyperparameter search is complete.")
print(f"The optimal number of layers: {best_hps.get('num_layers')}")
for i in range(best_hps.get('num_layers')):
    print(f"Hidden layer {i} has optimal nodes {best_hps.get('units_' + str(i))} and activation {best_hps.get('act_' + str(i))}")
print(f"The optimal learning rate for the optimize is {best_hps.get('learning_rate')}")

state = True
while state == True:
    answ = input("Do you want to save model? (y/n) ")
    if answ == "y":
        name = "test_ExtraLayers2"
        tuner.hypermodel.build(best_hps).save(f"models/model_{name}.h5")
        state = False
        print("Model saved")
    elif answ == "n":
        state = False
        print("Model not saved")