"""Contains various utility functions used in the tensorno library."""

import tensorflow as tf
from math import sqrt


def get_custom_initializers(num_inputs: int, bias: float = 0.) -> dict:
    """Returns a kwargs dict that can be passed to tf.keras.Dense-type layers
    with kernel and bias initializers. The kernel is initialized with a normal
    distribution weighted with the number of inputs. the bias is initialized to
    the same constant value.

    Args:
        num_inputs (int): Number of inputs to the layer.
        bias (float, optional): Bias terms' initial value. Defaults to 0..

    Returns:
        dict: kwargs to be passed to contructor of tf.keras.layers.Dense.
    """
    result = dict(
        kernel_initializer=tf.keras.initializers.RandomNormal(
            stddev=sqrt(2./num_inputs)
        ),
        bias_initializer=tf.keras.initializers.Constant(bias)
    )

    return result
