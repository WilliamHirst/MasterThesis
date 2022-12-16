"""Contains the MaxOut and ChannelOut Layer classes
interfaced with tensorflow"""

import tensorflow as tf

# from typing import Optional
from tensorflow.keras.layers import Layer
from tensorflow.python.ops import gen_math_ops, nn_ops


class MaxOut(Layer):
    """``MaxOut``."""

    def __init__(
        self,
        units: int,
        num_inputs: int,
        num_groups: int = 2,
        kernel_initializer="glorot_uniform",
        bias_initializer="zeros",
        regularizer = None,
        **kwargs
    ):
        # Using Layer initialization
        super().__init__(**kwargs)

        # Initializing units and groups
        if num_groups > units:
            num_groups = units
        elif units % num_groups != 0:
            raise ValueError(f"Number of units ({units}) "
                             "is not divisible by number of groups "
                             f"({num_groups})")
        self.units = int(units) if not isinstance(units, int) else units
        if self.units < 0:
            raise ValueError(
                "Received an invalid value for `units`, expected "
                f"a positive integer. Received: units={units}"
            )
        self.num_groups = int(num_groups) if not isinstance(num_groups, int) \
            else num_groups
        if self.num_groups < 0:
            raise ValueError(
                "Received an invalid value for `num_groups`, expected "
                f"a positive integer. Received: num_groups={num_groups}"
            )

        # Initializing weights and biases
        self.kernel_initializer = tf.keras.initializers.get(kernel_initializer)
        self.bias_initializer = tf.keras.initializers.get(bias_initializer)

        num_competitors = self.units // self.num_groups
        self.new_shape = [num_inputs, self.num_groups, num_competitors]

        self.counter = tf.zeros([num_inputs, self.units], dtype=tf.float32)

        self.kernel = self.add_weight(
            "kernel",
            shape=[num_inputs, self.units],
            initializer=self.kernel_initializer,
            dtype=tf.float32,
            regularizer = regularizer,
            trainable=True,
        )
        self.bias = self.add_weight(
            "bias",
            shape=[self.units, ],
            initializer=self.bias_initializer,
            dtype=tf.float32,
            trainable=True,
        )

    def call(self, inputs: tf.Tensor) -> tf.Tensor:
        # Passing input through weight kernel and adding bias terms
        inputs = gen_math_ops.MatMul(a=inputs, b=self.kernel)
        inputs = nn_ops.bias_add(inputs, self.bias)

        num_inputs = inputs.shape[0]
        if num_inputs is None:
            num_inputs = -1
        num_competitors = self.units // self.num_groups
        new_shape = [num_inputs, self.num_groups, num_competitors]

        # Reshaping outputs such that they are grouped correctly
        inputs = tf.reshape(inputs, new_shape)
        # Finding maximum activation in each group
        outputs = tf.math.reduce_max(inputs, axis=-1,keepdims=True)

        counter = tf.where(tf.equal(inputs, outputs), outputs, 0.)
        # Reshaping outputs to original input shape
        self.counter = tf.reshape(counter, [num_inputs, self.units])

        return tf.reshape(outputs,[num_inputs, self.num_groups])

    def get_config(self):
        config = super().get_config()
        config.update(
            {
                "units": self.units,
                "num_groups": self.num_groups
            }
        )
        return config


class ChannelOut(Layer):
    """``ChannelOut``."""

    def __init__(
        self,
        units: int,
        num_inputs: int,
        num_groups: int = 2,
        kernel_initializer="glorot_uniform",
        bias_initializer="zeros",
        regularizer = None,
        **kwargs
    ):
        super().__init__(**kwargs)

        # Initializing units and groups
        if num_groups > units:
            num_groups = units
        elif units % num_groups != 0:
            raise ValueError(f"Number of units ({units}) "
                             "is not divisible by number of groups "
                             f"({num_groups})")
        self.units = int(units) if not isinstance(units, int) else units
        if self.units < 0:
            raise ValueError(
                "Received an invalid value for `units`, expected "
                f"a positive integer. Received: units={units}"
            )
        self.num_groups = int(num_groups) if not isinstance(num_groups, int) \
            else num_groups
        if self.num_groups < 0:
            raise ValueError(
                "Received an invalid value for `num_groups`, expected "
                f"a positive integer. Received: num_groups={num_groups}"
            )

        # Initializing weights and biases
        self.kernel_initializer = tf.keras.initializers.get(kernel_initializer)
        self.bias_initializer = tf.keras.initializers.get(bias_initializer)

        self.counter = tf.zeros([num_inputs, self.units], dtype=tf.float32)

        self.kernel = self.add_weight(
            "kernel",
            shape=[num_inputs, self.units],
            initializer=self.kernel_initializer,
            dtype=tf.float32,
            regularizer = regularizer,
            trainable=True,
        )
        self.bias = self.add_weight(
            "bias",
            shape=[self.units, ],
            initializer=self.bias_initializer,
            dtype=tf.float32,
            trainable=True,
        )

    def call(self,
             inputs: tf.Tensor,
             mask: tf.Tensor = None
             ) -> tf.Tensor:
        # Pass input through weight kernel and adding bias terms.
        inputs = gen_math_ops.MatMul(a=inputs, b=self.kernel)
        inputs = nn_ops.bias_add(inputs, self.bias)

        num_inputs = inputs.shape[0]
        if num_inputs is None:
            num_inputs = -1

        # Reshaping inputs such that they are grouped correctly
        num_competitors = self.units // self.num_groups
        new_shape = [num_inputs, self.num_groups, num_competitors]
        inputs = tf.reshape(inputs, new_shape)

        # Finding maximum activations and setting losers to 0.
        outputs = tf.math.reduce_max(inputs, axis=-1, keepdims=True)
        outputs = tf.where(tf.equal(inputs, outputs), outputs, 0.)
        # Reshaping outputs to original input shape
        outputs = tf.reshape(outputs, [num_inputs, self.units])

        self.counter = outputs
        return outputs

    def get_config(self):
        config = super().get_config()
        config.update(
            {
                "units": self.units,
                "num_groups": self.num_groups
            }
        )
        return config

class StochChannelOut(Layer):
    """``ChannelOut``."""

    def __init__(
        self,
        units: int,
        num_inputs: int,
        num_groups: int = 2,
        kernel_initializer="glorot_uniform",
        bias_initializer="zeros",
        regularizer = None,
        **kwargs
    ):
        super().__init__(**kwargs)

        # Initializing units and groups
        if num_groups > units:
            num_groups = units
        elif units % num_groups != 0:
            raise ValueError(f"Number of units ({units}) "
                             "is not divisible by number of groups "
                             f"({num_groups})")
        self.units = int(units) if not isinstance(units, int) else units
        if self.units < 0:
            raise ValueError(
                "Received an invalid value for `units`, expected "
                f"a positive integer. Received: units={units}"
            )
        self.num_groups = int(num_groups) if not isinstance(num_groups, int) \
            else num_groups
        if self.num_groups < 0:
            raise ValueError(
                "Received an invalid value for `num_groups`, expected "
                f"a positive integer. Received: num_groups={num_groups}"
            )

        self.index = tf.range(start = 0, limit = self.units,  dtype=tf.int32)
        self.zeros = tf.zeros_like(self.index)

        # Initializing weights and biases
        self.kernel_initializer = tf.keras.initializers.get(kernel_initializer)
        self.bias_initializer = tf.keras.initializers.get(bias_initializer)


        self.kernel = self.add_weight(
            "kernel",
            shape=[num_inputs, self.units],
            initializer=self.kernel_initializer,
            dtype=tf.float32,
            trainable=True,
            regularizer = regularizer
        )
        self.bias = self.add_weight(
            "bias",
            shape=[self.units, ],
            initializer=self.bias_initializer,
            dtype=tf.float32,
            trainable=True,
        )

    def call(self,
             inputs: tf.Tensor,
             mask: tf.Tensor = None
             ) -> tf.Tensor:
        # Pass input through weight kernel and adding bias terms.
        inputs = gen_math_ops.MatMul(a=inputs, b=self.kernel)
        inputs = nn_ops.bias_add(inputs, self.bias)
        #tf.print(inputs)
        num_inputs = inputs.shape[0]

        if num_inputs is None:
            num_inputs = -1
        
        shuffle_index = tf.random.shuffle(self.index)
        unshuffle_index = tf.tensor_scatter_nd_update(tensor = self.zeros , indices = tf.reshape(shuffle_index, [inputs.shape[1],1]), updates = self.index)
        inputs_s = tf.gather(inputs, shuffle_index, axis = 1)
  

        # Reshaping inputs such that they are grouped correctly
        num_competitors = self.units // self.num_groups
        new_shape = [num_inputs, self.num_groups, num_competitors]
        inputs_s = tf.reshape(inputs_s, new_shape)

        # Finding maximum activations and setting losers to 0.
        outputs = tf.math.reduce_max(inputs_s, axis=-1, keepdims=True)

        outputs = tf.where(tf.equal(inputs_s, outputs), 1.0, 0.)
        # Reshaping outputs to original input shape
        outputs = tf.reshape(outputs, [num_inputs, self.units])

        outputs = tf.gather(outputs, unshuffle_index, axis = 1) 
        outputs = tf.multiply(inputs, outputs)
        self.counter = outputs
        return outputs

    def get_config(self):
        config = super().get_config()
        config.update(
            {
                "units": self.units,
                "num_groups": self.num_groups
            }
        )
        return config
