"""Contains useful functions for plotting neural network architecture."""
import tensorflow as tf
import matplotlib.pyplot as plt
import numpy as np

import context
import plot_utils
from tensorno.layers import MaxOut, ChannelOut


def get_layer_names(network: tf.keras.Model) -> list:
    """Get the string names of the layers of a tf Model instance.

    Args:
        network (tf.keras.Model): Model from which to get layer names.

    Returns:
        list: list of the string names of the layers in network.
    """
    return list(map(lambda x: x.name, network.layers))


def get_layer_activations(network: tf.keras.Model,
                          layer,
                          input):
    """Get the output from a given layer in a tf Model given input.

    Args:
        network (tf.keras.Model): Model from which to extract activation.
        layer (str or tf.keras.layers.Layer): Either the name of the layer or
                                              the Layer instance from which to
                                              get the outputs.
        input (np.ndarray or tf.Tensor): The initial feature input to generate
                                         the activations.

    Returns:
        np.ndarray or tf.Tensor: Output of the given layer.
    """
    if isinstance(layer, str):
        intermediate_model = tf.keras.Model(network.input,
                                            network.get_layer(layer).output,
                                            name="intermediate")
    elif isinstance(layer, tf.keras.layers.Layer):
        intermediate_model = tf.keras.Model(network.input,
                                            layer.output,
                                            name="intermediate")
    else:
        raise ValueError("Given layer is not a string or Layer instance"
                         f" (is type {type(layer)})")
    intermediate_model(input)

    layer_i = intermediate_model.layers[-1]

    return np.asarray(layer_i.counter)
    
    # if isinstance(input, np.ndarray):
    #     return intermediate_model(input).numpy()
    # else:
    #     return intermediate_model(input)


def get_all_activations(network: tf.keras.Model, input) -> list:
    """Gets the output from every layer in a tf Model instance given input.

    Args:
        network (tf.keras.Model): Model from which to extract activation.
        input (np.ndarray or tf.Tensor): The initial feature input to generate
                                         the activations.

    Returns:
        list: list of np.ndarrays or tf.Tensors with the output from
              every layer.
    """
    activations = list()
    for layer in network.layers:
        if isinstance(layer, (MaxOut, ChannelOut)):
            activations.append(get_layer_activations(network, layer, input))

    return activations


def group_nodes(x: float, nodes: int, groups: int, ax):
    """Plots rectangles around each group of competitors in a LWTA layer.

    Args:
        x (float): x-position of the layer on ax.
        nodes (int): Number of nodes in the layer.
        groups (int): Number of groups in the layer.
        ax (ax): plt ax object on which to plot.

    Returns:
        tuple: ax object on which the function plotted.
    """
    competitors = nodes // groups
    for group in range(groups):
        ax.add_patch(plt.Rectangle((x-0.1,
                                    group*competitors - nodes/2 - 0.3),
                                   0.2, competitors-0.4,
                                   fc='none',
                                   ec=plot_utils.colors[-1],
                                   lw=1.2,
                                   clip_on=False))

    return ax


def plot_nodes(layers: list, ax=None):
    """Plots all the nodes in every layer in layers, grouping LWTA groups.

    Args:
        layers (list): list of tf.keras.layers.Layer instances.
        ax (ax, optional): plt ax on which to plot. Defaults to None.

    Returns:
        tuple: ax object on which the function plotted.
    """
    if ax is None:
        _, ax = plt.subplots()
        ax.set_facecolor("white")
        ax.set_xticks([])
        ax.set_yticks([])

    for x, layer in enumerate(layers):
        if isinstance(layer, (MaxOut, ChannelOut)):
            nodes = layer.units
            for node in np.arange(-nodes/2, nodes/2):
                ax.scatter([x], [node], color=plot_utils.colors[0])

            if isinstance(layer, (MaxOut, ChannelOut)):
                group_nodes(x, nodes, layer.num_groups, ax=ax)

    return ax


def plot_active_nodes(layers: list, isactive: list, ax=None):
    """Plots all the nodes in every layer in layers, grouping LWTA
    and indicating whether the nodes are active or not.

    Args:
        layers (list): list of tf.keras.layers.Layer instances.
        isactive (list): list of boolean np.ndarrays indicating whether node
                         is active or not by layer.
        ax (ax, optional): plt ax on which to plot. Defaults to None.

    Returns:
        tuple: ax object on which the function plotted.
    """
    if ax is None:
        _, ax = plt.subplots()
        ax.set_facecolor("white")
        ax.set_xticks([])
        ax.set_yticks([])

    for x, (layer, active) in enumerate(layers, isactive):
        colors = np.where(active, "r", "b")
        nodes = layer.units
        for node, color in zip(np.arange(-nodes/2, nodes/2), colors):
            ax.scatter([x], [node], c=color)

        if isinstance(layer, (MaxOut, ChannelOut)):
            group_nodes(x, nodes, layer.num_groups, ax=ax)

    return ax


def plot_pathways(layers: list, isactive: list, ax=None, **plot_kwargs):
    """Plots lines going through all active nodes of all layers.

    Args:
        layers (list): list of tf.keras.layers.Layer instances.
        isactive (list): list of boolean np.ndarrays indicating whether node
                         is active or not by layer.
        ax (ax, optional): plt ax on which to plot. Defaults to None.

    Returns:
        tuple: ax object on which the function plotted.
    """
    if ax is None:
        _, ax = plt.subplots()
        ax.set_facecolor("white")
        ax.set_xticks([])
        ax.set_yticks([])

    line_kwargs = dict(
        color=plot_utils.colors[1],
        lw=1,
        alpha=0.005
    )
    line_kwargs.update(plot_kwargs)

    active_nodes_old = [np.nan]
    for x, (layer, active) in enumerate(zip(layers, isactive)):
        nodes = layer.units
        active_nodes_new = np.where(active,
                                    np.arange(-nodes/2, nodes/2),
                                    np.nan)
        for old_node in active_nodes_old:
            for new_node in active_nodes_new:
                ax.plot([x-1, x], [old_node, new_node], **line_kwargs)
        active_nodes_old = active_nodes_new

    return ax

def plot_value_line(layers: list, isactive: list, pred: float, ax=None, **plot_kwargs):
    """Plots lines going through all active nodes of all layers.

    Args:
        layers (list): list of tf.keras.layers.Layer instances.
        isactive (list): list of boolean np.ndarrays indicating whether node
                         is active or not by layer.
        ax (ax, optional): plt ax on which to plot. Defaults to None.

    Returns:
        tuple: ax object on which the function plotted.
    """
    if ax is None:
        _, ax = plt.subplots()
        ax.set_facecolor("white")
        ax.set_xticks([])
        ax.set_yticks([])

    
    line_kwargs = dict(
        color=plot_utils.colors[1],
        lw=1,
        alpha=0.1
    )
    line_kwargs.update(plot_kwargs)


    x = len(layers)
    nodes = layers[-1].units
    last_nodes = np.where(isactive[-1],
                            np.arange(-nodes/2, nodes/2),
                            np.nan)
    for p in pred[-1][0]:
            for node in last_nodes: 
                ax.plot([x-1, x], [node, p], **line_kwargs)

    ax.arrow(x = x , dx = 0, y = last_nodes[0],  dy = last_nodes[-1]-last_nodes[0], color=plot_utils.colors[1], alpha = 1)
    return ax