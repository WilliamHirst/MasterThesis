import tensorflow as tf
from tensorflow.keras import optimizers

import keras_tuner as kt
import sys
sys.path.insert(1, "../../")
from Utilities import *



myPath = "/storage/William_Sakarias/William_Data"

signal = "ttbar"



df, y, df_data, channels = loadDf(myPath, signal)
train, val, test = splitData(df, y)

X_train, Y_train, W_train, C_train = train
X_val, Y_val, W_val, C_val = val
X_test, Y_test, W_test, C_test = test

nrFeature = len(X_train[0])
nData = len(X_train[0,:])*0.8
class MyHyperModel(kt.HyperModel):
    def build(self, hp):
        model = tf.keras.Sequential(
            [
                tf.keras.layers.Dense(
                    units=hp.Int("num_of_neurons0", min_value=5, max_value=40, step=5),
                    activation=tf.keras.layers.LeakyReLU(alpha=0.01),
                    input_shape=(nrFeature,),
                ),
                tf.keras.layers.Dense(
                    units=hp.Int("num_of_neurons1", min_value=5, max_value=40, step=5),
                    activation=tf.keras.layers.LeakyReLU(alpha=0.01),
                ),
                tf.keras.layers.Dense(1, activation="sigmoid"),
            ]
        )
        hp_learning_rate = hp.Choice("learning_rate", values=[9e-2, 9.5e-2, 1e-3, 1.5e-3])
        optimizer = optimizers.Adam(learning_rate=hp_learning_rate)

        model.compile(loss="binary_crossentropy", optimizer=optimizer, metrics=["accuracy"])
        return model
    def fit(self, hp, model, *args, **kwargs):
        return model.fit(
            *args,
            batch_size=int(nData/40),
            **kwargs,
        )

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
with tf.device("/GPU:0"):
    tuner.search(X_train.astype('float32'), Y_train.astype('float32'), epochs=100, validation_data=(X_val.astype('float32'),Y_val.astype('float32')))
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
        name = input("name: ")
        tuner.hypermodel.build(best_hps).save(f"tf_models/model_{name}.h5")
        state = False
        print("Model saved")
    elif answ == "n":
        state = False
        print("Model not saved")