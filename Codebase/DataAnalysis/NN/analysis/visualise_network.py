"""Contains the code used to plot visualise LWTA NN architecture
and activation."""
import matplotlib.pyplot as plt
import tensorflow as tf

import network_plot_tools as npt
import plot_utils

import sys

sys.path.insert(1, "../../../")
from Utilities import *

def plot_channelout_architecture(network: tf.keras.Model,
                                 inputs,
                                 targets,
                                 ax=None,
                                 **plot_kwargs):
    """Plots a channel out network on specified ax, with pathways indicating
    the active nodes for every datapoint in the input.

    Args:
        network (tf.keras.Model): A tf neural network model with
                                  ChannelOut hidden layers.
        input (list): Iterable of the feature inputs that generate the
                      layer activations.
        ax (ax, optional): plt ax on which to plot. Defaults to None.
    """
    if ax is None:
        fig, ax = plt.subplots()
        ax.set_xticks([])
        ax.set_yticks([])

    npt.plot_nodes(network.layers, ax=ax)
    for i in range(len(inputs)):
        all_activations = npt.get_all_activations(network,
                                                  inputs[i].reshape(1, -1))
        isactive = [(activations != 0.).reshape([activations.shape[-1]])
                    for activations in all_activations]
        # isactive.append(np.where(
        #     all_activations[-1] == np.max(all_activations[-1]),
        #     True, False
        # ).reshape(all_activations[-1].shape[-1]))
        plot_kwargs = dict(
            color=plot_utils.colors[int(targets[i])+1],
            lw=1,
            alpha=0.005
        )
        ax = npt.plot_pathways(network.layers, isactive, ax=ax, **plot_kwargs)
        
        ax = npt.plot_value_line(network.layers, isactive, all_activations, ax=ax, **plot_kwargs)
    
    ax = npt.plotAxis(network.layers, isactive, all_activations, ax=ax)
    
    return fig,ax


if __name__ == "__main__":
    import tensorflow as tf
    from tensorflow.keras import optimizers
    import numpy as np
    import pandas as pd

    sys.path.insert(1, "../")
    from tensorno.layers import MaxOut


    myPath = "/storage/William_Sakarias/William_Data"

    signal = "ttbarHNLMaxChannel"

    print(f"Starting test: {signal}")

    df, y, df_data, channels = loadDf(myPath, notInc=["LRS", "filtch", "LepMLm15","LepMLp15","LepMLm75"])

    print("Preparing data....")
    train, val = splitAndPrepData(df, y, scale = True, ret_scaleFactor=True)
    print("Done.")

    X_train, Y_train, W_train, C_train = train
    X_val, Y_val, W_val, C_val, scaleFactor = val
    nrFeature = nFeats(X_train)
    # X_val = np.random.rand(500,32)
    # nrFeature = nFeats(X_val)

    model = tf.keras.Sequential()
    model.add(tf.keras.layers.InputLayer(input_shape=(nrFeature,)))
    model.add(MaxOut(units=10, num_inputs=nrFeature, num_groups=5))
    model.add(MaxOut(units=10, num_inputs=5, num_groups=5))
    model.add(MaxOut(units=10, num_inputs=5, num_groups=5))
    model.add(tf.keras.layers.Dense(1, activation="sigmoid"))

    optimizer = optimizers.Adam(learning_rate=1e-3)
    model.compile(loss="binary_crossentropy", optimizer=optimizer, weighted_metrics="AUC")
    print("Done compiling.")

    index_1 = Y_val[Y_val == 1].index[:500]
    index_2 = Y_val[Y_val == 0].index[:500]
    index = index_1.append(index_2)

    X_viz = X_val.loc[index].sample(frac = 1).reset_index(drop=True)
    Y_viz = Y_val.loc[index].sample(frac = 1).reset_index(drop=True)

    print(Y_viz)

    fig, ax = plot_channelout_architecture(model,
                                 X_viz.values,
                                 Y_viz.values,
                                 )
    ax.plot()
    fig.savefig("BeforeTraining.pdf")
    #ax.show()
    
    with tf.device("/GPU:0"):
        callback = tf.keras.callbacks.EarlyStopping(monitor='val_auc', 
                                                    patience=3, 
                                                    restore_best_weights = True,
                                                    verbose = 1,
                                                    mode = "max")
        history = model.fit(X_train, 
                            Y_train,
                            sample_weight = W_train, 
                            epochs=100, 
                            batch_size=8096, 
                            callbacks = [callback], #, CC],
                            validation_data=(X_val, Y_val, W_val),
                            verbose = 1)
    
    fig, ax = plot_channelout_architecture(model,
                                 X_viz.values,
                                 Y_viz.values,
                                 )
    ax.plot()
    fig.savefig("AfterTraining.pdf")
