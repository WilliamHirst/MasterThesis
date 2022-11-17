import tensorflow as tf
import keras as keras
import numpy as np

from sklearn.metrics import roc_curve, roc_auc_score, auc




def max_out(inputs, num_units = 200, axis=None):
    shape = inputs.get_shape().as_list()
    if shape[0] is None:
        shape[0] = -1
    if axis is None:  # Assume that channel is the last dimension
        axis = -1
    num_channels = shape[axis]
    if num_channels % num_units:
        raise ValueError('number of features({}) is not '
                         'a multiple of num_units({})'.format(num_channels, num_units))
    shape[axis] = num_units
    shape += [num_channels // num_units]
    outputs = tf.reduce_max(tf.reshape(inputs, shape), -1, keepdims=False)
    return outputs

def channel_out(inputs, num_units = 200, axis=None):
    shape = inputs.get_shape().as_list()
    if shape[0] is None:
        shape[0] = -1
    if axis is None:  # Assume that channel is the last dimension
        axis = -1
    num_channels = shape[axis]
    if num_channels % num_units:
        raise ValueError('number of features({}) is not '
                         'a multiple of num_units({})'.format(num_channels, num_units))
    shape[axis] = num_units
    shape += [num_channels // num_units]
    grouped = tf.reshape(inputs, shape)
    top_vals = tf.reduce_max(tf.reshape(inputs, shape), -1, keepdims=True)
    isMax = tf.reshape(tf.greater_equal(grouped, top_vals), [shape[0], num_channels])
    output = tf.multiply(tf.cast(isMax,tf.float32), inputs)

    return output

class Cust_Callback(keras.callbacks.Callback):
    def on_train_begin(self, logs={}):
        self._data = []

    def on_epoch_end(self, batch, logs={}):
        X_val, y_val, sample_weight = self.validation_data[0], self.validation_data[1],  self.validation_data[2]

        y_pred = np.asarray(self.model.predict(X_val))

        fpr, tpr, thresholds = roc_curve(y_val,y_pred, sample_weight = sample_weight, pos_label=1)


        gmeans = np.sqrt(np.array(tpr) * (1-np.array(fpr)/np.max(np.array(fpr))))
        ix = np.argmax(gmeans)
        best_threshold = thresholds[ix]

        nrB = np.sum(sample_weight[y_pred < best_threshold])
        print(nrB)
        nrS = np.sum(sample_weight[y_pred > best_threshold])
        print(nrS)
        sig = nrS/np.sqrt(nrB)

        self._data.append({
            'val_sig': sig,
        })
        print(f"The significance: {sig} ")

        return

    def get_data(self):
        return self._data

class Cust_Metric(tf.keras.metrics.Metric):
    def __init__(self, num_classes):
        super().__init__()
        self.num_classes = num_classes

    def AUC(self, y_true, y_pred, sample_weight):
        fpr, tpr, _ = roc_curve(y_true,y_pred, sample_weight = sample_weight, pos_label=1)

        roc_auc = auc(fpr,tpr)

        return roc_auc

    def Sig(self, y_true, y_pred, sample_weight):
        fpr, tpr, _ = roc_curve(y_true,y_pred, sample_weight = sample_weight, pos_label=1)

        roc_auc = auc(fpr,tpr)

        return roc_auc

        


