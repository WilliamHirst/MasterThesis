import tensorflow as tf
from tensorflow.keras import optimizers

import keras_tuner as kt
import sys
sys.path.insert(1, "../../")
from Utilities import *



myPath = "/storage/William_Sakarias/William_Data"

signal = "ttbar"

procUsed = int(100/20)

df, y, df_data, channels = loadDf(myPath, signal)
train, val, test = splitData(df, y)

X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val = val
X_test, Y_test, W_test, C_test = test

nrFeature = len(X_train[0])
nData = len(X_train[0,:])
class MyHyperModel(kt.HyperModel):
    def build(self, hp):
        model = tf.keras.Sequential(
            [
                tf.keras.layers.Dense(
                    units=hp.Int("num_of_neurons0", min_value=10, max_value=30, step=10),
                    activation=tf.keras.layers.LeakyReLU(alpha=0.01),
                    input_shape=(nrFeature,),
                ),
                tf.keras.layers.Dense(
                    units=hp.Int("num_of_neurons1", min_value=10, max_value=30, step=10),
                    activation=tf.keras.layers.LeakyReLU(alpha=0.01),
                ),
                tf.keras.layers.Dense(1, activation="sigmoid"),
            ]
        )
        hp_learning_rate = hp.Choice("learning_rate", values=[1e-2, 5e-3, 1e-3, 5e-4])
        optimizer = optimizers.Adam(learning_rate=hp_learning_rate)

        model.compile(loss="binary_crossentropy", optimizer=optimizer, metrics=["accuracy"])
        return model
    

start_time = timer(None)


tuner = kt.Hyperband(
    MyHyperModel(),
    objective="val_accuracy",
    max_epochs=50,
    factor=3,
    directory="GridSearches",
    project_name="NN",
    overwrite=True,
)
callback = tf.keras.callbacks.EarlyStopping(monitor='val_accuracy', patience=3)
with tf.device("/GPU:0"):
    tuner.search(X_train[::procUsed], 
                    Y_train[::procUsed], 
                    sample_weight = W_train[::procUsed], 
                    batch_size =int(nData/40), 
                    epochs=50, 
                    validation_data=(X_val[::procUsed],
                                    Y_val[::procUsed],
                                    W_val[::procUsed]),
                    callbacks = [callback]
                )
timer(start_time)
best_hps = tuner.get_best_hyperparameters(num_trials=1)[0]

print(
    f"""
    The hyperparameter search is complete. \
    The optimal number nodes in start, \
    first, second layer is {best_hps.get('num_of_neurons0')}, \
    {best_hps.get('num_of_neurons1')} and \
    the optimal learning rate for the optimizer\
    is {best_hps.get('learning_rate')}."""
)

state = True
while state == True:
    answ = input("Do you want to save model? (y/n) ")
    if answ == "y":
        name = "test"#input("name: ")
        tuner.hypermodel.build(best_hps).save(f"models/model_{name}.h5")
        state = False
        print("Model saved")
    elif answ == "n":
        state = False
        print("Model not saved")