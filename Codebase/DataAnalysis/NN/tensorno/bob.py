"""Contains builders used for creating maxout and channelout neural networks
in tensorflow"""
import tensorflow as tf
import keras_tuner

if __name__ == "__main__":
    from layers import MaxOut, ChannelOut
    from utils import get_custom_initializers
else:
    from tensorno.layers import MaxOut, ChannelOut
    from tensorno.utils import get_custom_initializers


def build_LWTA_regressor(
    num_layers: int,
    units: tuple,
    num_groups: tuple,
    Layer: str,
    num_features: int = 2,
    **compile_kwargs
) -> tf.keras.Sequential:
    """Builds a LWTA FFNN for regression with the specified architecture.

    Args:
        num_layers (int): Number of hidden layers.
        units (tuple): Iterable with the number of nodes for each hidden layer.
        num_groups (tuple): Iterable with the number of competing groups
                            for each hidden layer.
        Layer (str): Specifies whether to use MaxOut or ChannelOut layers.
                     Either "max_out" or "channel_out".
        num_features (int, optional): Number of input features. Defaults to 2.

    Returns:
        tf.keras.Sequential: Compiled Sequential model as regressor.
    """
    # Get the appropriate Layer class
    if Layer == "max_out":
        Layer = MaxOut
    elif Layer == "channel_out":
        Layer = ChannelOut
    else:
        raise ValueError(f"Given Layer ({Layer}) is not supported. "
                         "must be `max_out` or `channel_out`.")

    model = tf.keras.Sequential()
    model.add(tf.keras.layers.InputLayer(
        input_shape=(num_features,),
        name="input"
    ))

    # Add hidden layers.
    num_inputs = num_features
    for layer in range(num_layers):
        model.add(Layer(
            units=units[layer],
            num_inputs=num_inputs,
            num_groups=num_groups[layer],
            **get_custom_initializers(
                num_features if layer == 1 else num_groups[layer-1]
            ),
            name=f"{Layer.__name__.lower()}_{layer+1}"
        ))
        # MaxOut and ChannelOut have different number of outputs.
        if isinstance(Layer, MaxOut):
            num_inputs = num_groups[layer]
        else:
            num_inputs = units[layer]

    # Add output layer.
    model.add(tf.keras.layers.Dense(
        units=1,
        activation="linear",
        **get_custom_initializers(num_inputs),
        name="output"
    ))

    kwargs = dict(
        optimizer="adam",
        loss="mse"
    )
    kwargs.update(compile_kwargs)
    model.compile(**kwargs)

    return model


def build_LWTA_classifier(
    num_layers: int,
    units: tuple,
    num_groups: tuple,
    Layer: str,
    num_features: int = 2,
    num_categories: int = 1,
    **compile_kwargs
) -> tf.keras.Sequential:
    """Builds a LWTA FFNN for classification with the specified architecture.

    Args:
        num_layers (int): Number of hidden layers.
        units (tuple): Iterable with the number of nodes for each hidden layer.
        num_groups (tuple): Iterable with the number of competing groups
                            for each hidden layer.
        Layer (str): Specifies whether to use MaxOut or ChannelOut layers.
                     Either "max_out" or "channel_out".
        num_features (int, optional): Number of input features. Defaults to 2.
        num_categories (int, optional): Number of output categories.
                                        Defaults to 2.

    Returns:
        tf.keras.Sequential: Compiled Sequential model as classifier.
    """
    # Get the appropriate Layer class
    if Layer == "max_out":
        Layer = MaxOut
    elif Layer == "channel_out":
        Layer = ChannelOut
    else:
        raise ValueError(f"Given Layer ({Layer}) is not supported. "
                         "must be `max_out` or `channel_out`.")

    model = tf.keras.Sequential()
    model.add(tf.keras.layers.InputLayer(
        input_shape=(num_features,),
        name="input"
    ))

    # Add hidden layers.
    num_inputs = num_features
    for layer in range(num_layers):
        if num_groups[layer] > units[layer]:
            num_groups[layer] = units[layer]
        model.add(Layer(
            units=units[layer],
            num_inputs=num_inputs,
            num_groups=num_groups[layer],
            **get_custom_initializers(
                num_features if layer == 1 else num_groups[layer-1]
            ),
            name=f"{Layer.__name__.lower()}_{layer+1}"
        ))
        # MaxOut and ChannelOut have different number of outputs.
        if Layer is MaxOut:
            num_inputs = num_groups[layer]
        else:
            num_inputs = units[layer]

    # Add output layer.
    model.add(tf.keras.layers.Dense(
        units=num_categories,
        activation="sigmoid" if num_categories == 1 else "softmax",
        **get_custom_initializers(num_inputs),
        name="output"
    ))

    kwargs = dict(
        optimizer="adam",
        loss=("binary" if num_categories == 1 else "categorical")
        + "_crossentropy",
        metrics=["accuracy"]
    )
    kwargs.update(compile_kwargs)
    model.compile(**kwargs)

    return model


def build_LWTA_architecture(
    hp: keras_tuner.HyperParameters,
    Layer: str,
    layer_choices: list = [2, 3, 4, 5],
    node_choices: list = [4, 8, 16, 32],
    group_choices: list = [1, 2, 4, 8],
    isregressor: bool = True,
    num_features: int = None,
    num_categories: int = None
) -> tf.keras.Sequential:
    """Builder for interfacing with keras_tuner to tune model architecture.

    Args:
        hp (keras_tuner.HyperParameters): Stores and generates hyperparameters
                                          for the model.
        Layer (str): Specifies whether to use MaxOut or ChannelOut layers.
                     Either "max_out" or "channel_out".
        isregressor (bool, optional): Whether the model is used for regression
                                      or classification. Defaults to True.

    Returns:
        tf.keras.Sequential: Sequential model with a choice for architecture.
    """
    num_layers = hp.Choice("num_layers", layer_choices)
    units = [hp.Choice(f"num_nodes{i}", node_choices)
             for i in range(1, num_layers+1)]
    num_groups = [hp.Choice(f"num_groups{i}", group_choices)
                  for i in range(1, num_layers+1)]

    if isregressor:
        return build_LWTA_regressor(num_layers, units, num_groups, Layer,
                                    num_features)
    else:
        return build_LWTA_classifier(num_layers, units, num_groups, Layer,
                                     num_features, num_categories)


def LWTA_architecture_builder(
    Layer: str,
    layer_choices: list = [2, 3, 4, 5],
    node_choices: list = [4, 8, 16, 32],
    group_choices: list = [1, 2, 4, 8],
    isregressor: bool = True,
    num_features: int = None,
    num_categories: int = None
):
    """Returns appropriate architecture builder for interfacing with
    keras_tuner for tuning the architecture of a LWTA model.

    Args:
        Layer (str): Specifies whether to use MaxOut or ChannelOut layers.
                     Either "max_out" or "channel_out".
        isregressor (bool, optional): Whether the model is used for regression
                                      or classification. Defaults to True.

    Returns:
        function: Builder that keras_tuner can use for hyperparameter search.
    """
    return lambda hp: build_LWTA_architecture(hp,
                                              Layer=Layer,
                                              layer_choices=layer_choices,
                                              node_choices=node_choices,
                                              group_choices=group_choices,
                                              isregressor=isregressor,
                                              num_features=num_features,
                                              num_categories=num_categories)
