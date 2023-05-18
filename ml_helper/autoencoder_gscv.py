import numpy as np

from sklearn.base import BaseEstimator, RegressorMixin, TransformerMixin
from sklearn.metrics import mean_squared_error

import tensorflow as tf
from tensorflow.keras.layers import Input, Dense
from tensorflow.keras.models import Model
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras import backend as K
from keras import regularizers

def weighted_mse(feature_weights):
    """
    Computes a custom weighted mean squared error (MSE) loss function.

    Args:
        feature_weights (numpy.ndarray): Weights calculated for each feature during feature extraction.

    Returns:
        function: Custom loss function that calculates the weighted MSE.
    """
    feature_weights = tf.constant(feature_weights, dtype=tf.float32)
    def loss(X_true, X_pred):
        squared_difference = tf.square(X_true - X_pred)
        weighted_squared_difference = squared_difference * feature_weights
        return tf.math.reduce_mean(weighted_squared_difference)
    return loss

def create_layer_sizes(num_features, num_layers, min_bottleneck_nodes):
    """
    Creates layer sizes for the autoencoder's encoding and decoding layers.

    Args:
        num_features (int): Number of input features.
        num_layers (int): Number of hidden layers in the autoencoder.
        min_bottleneck_nodes (int): Minimum number of nodes in the bottleneck layer.

    Returns:
        tuple: Tuple containing lists of encoding and decoding layer sizes.
    """
    # Below if-else statements may not be necessary in adopting this code for future integration into other models.
    # It was designed to test feature subsetting in voting classifier.

    def nearest_multiple_of_5(x):
        return (x + 4) // 5 * 5

    if num_layers == 1:
        bottleneck_nodes = min(min_bottleneck_nodes, int(num_features / 5))
        bottleneck_nodes = nearest_multiple_of_5(bottleneck_nodes)
        encoding_layers = [bottleneck_nodes]
        decoding_layers = []
    else:
        encoding_layers = []
        last_layer_size = num_features
        for i in range(1, num_layers):
            next_layer_size = max(min_bottleneck_nodes, int(num_features / (5 ** i)))
            next_layer_size = nearest_multiple_of_5(next_layer_size)

            encoding_layers.append(next_layer_size)
            last_layer_size = next_layer_size
            if len(encoding_layers) == num_layers - 1 or next_layer_size == min_bottleneck_nodes:
                break

        if encoding_layers[-1] > min_bottleneck_nodes:
            encoding_layers.append(min_bottleneck_nodes)
        else:
            encoding_layers[-1] = min_bottleneck_nodes

        decoding_layers = encoding_layers[-2::-1]

    return encoding_layers, decoding_layers

class Autoencoder(BaseEstimator, RegressorMixin, TransformerMixin):
    """
    Autoencoder model for feature extraction and data reconstruction.

    Args:
        encode_decode_pairs (dict): Dictionary specifying encoding and decoding layer sizes.
        epochs (int): Number of training epochs.
        batch_size (int): Batch size for training.
        activation_functions (dict): Dictionary specifying activation functions for encoder, decoder, and final output layers.
        feature_weights (numpy.ndarray, optional): Weights for each feature to prioritize "important" features during training. Defaults to None.
    """
    def __init__(self, encode_decode_pairs={'encoding_layers':[32], 'decoding_layers':[32]}, epochs=100, batch_size=32,
                 activation_functions={'encoder_activation': 'sigmoid', 'decoder_activation': 'sigmoid', 'final_activation': 'linear'}, feature_weights=None):
        self.encode_decode_pairs = encode_decode_pairs
        self.activation_functions = activation_functions
        self.epochs = epochs
        self.batch_size = batch_size
        self.feature_weights = feature_weights
        self.autoencoder = None
        self.encoder = None

    def _create_autoencoder(self, X):
        """
        Creates the autoencoder model.

        Args:
            X (numpy.ndarray): Input data.

        Returns:
            tuple: Tuple containing the autoencoder model and the encoder model.
        """
        input_dim = X.shape[1]
        encoding_layers = self.encode_decode_pairs['encoding_layers']
        decoding_layers = self.encode_decode_pairs['decoding_layers']
        encoder_activation = self.activation_functions['encoder_activation']
        decoder_activation = self.activation_functions['decoder_activation']
        final_activation = self.activation_functions['final_activation']

        # Encoder
        input_layer = Input(shape=(input_dim,))
        x = input_layer
        for units in encoding_layers:
            x = Dense(units, activation=encoder_activation)(x)
        encoded = x

        # Decoder
        for units in decoding_layers:
            x = Dense(units, activation=decoder_activation)(x)
        decoded = Dense(input_dim, activation=final_activation)(x)

        autoencoder = Model(input_layer, decoded)
        encoder = Model(input_layer, encoded)

        if self.feature_weights is not None:
            loss_function = weighted_mse(self.feature_weights)
        else:
            loss_function = 'mse'

        autoencoder.compile(optimizer='adam', loss=loss_function)

        return(autoencoder, encoder)

    def fit(self, X, y=None):
        """
        Fits the autoencoder to the input data.

        Args:
            X (numpy.ndarray): Input data.

        Returns:
            self
        """
        early_stopping_callback = EarlyStopping(monitor='loss', patience=10, restore_best_weights=True)

        self.autoencoder, self.encoder = self._create_autoencoder(X)
        self.autoencoder.fit(X, X, epochs=self.epochs, batch_size=self.batch_size, verbose=0, callbacks=[early_stopping_callback])
        return(self)

    def predict(self, X):
        """
        Reconstructs the input data using the autoencoder.

        Args:
            X (numpy.ndarray): Input data.

        Returns:
            numpy.ndarray: Reconstructed data.
        """
        if self.autoencoder is None:
            raise RuntimeError("The autoencoder has not been fitted yet.")
        return(self.autoencoder.predict(X))

    def score(self, X, y=None):
        """
        Computes the mean squared error (MSE) score for the input data.
        Need to transform the score returned by the autoencoder to adhere to GridSearchCV.
        GSCV takes the highest score as the "best". For MSE, lowest represent better scores.

        Args:
            X (numpy.ndarray): Input data.

        Returns:
            float: Negative mean squared error score.
        """
        X_pred = self.autoencoder.predict(X)
        if self.feature_weights is not None:
            loss_function = weighted_mse(self.feature_weights)
            reconst_error = loss_function(K.constant(X), K.constant(X_pred)).numpy()
        else:
            reconst_error = np.mean(((X - X_pred) ** 2))
        return(-reconst_error)

    def transform(self, X):
        """
        Encodes the input data using the encoder component of the autoencoder.

        Args:
            X (numpy.ndarray): Input data.

        Returns:
            numpy.ndarray: Encoded data.
        """
        if self.encoder is None:
            raise RuntimeError("The autoencoder has not been fitted yet.")
        return(self.encoder.predict(X))

    def get_params(self, deep=None):
        """
        Returns the parameters of the autoencoder.

        Args:
            deep (bool, optional): Whether to get parameters recursively. Defaults to None.

        Returns:
            dict: Autoencoder parameters.
        """
        return({'encode_decode_pairs': self.encode_decode_pairs,
                'activation_functions': self.activation_functions,
                'epochs': self.epochs,
                'batch_size': self.batch_size,
                'feature_weights': self.feature_weights})

    def summary(self):
        """
        Prints the summary of the autoencoder model.

        Raises:
            RuntimeError: If the autoencoder has not been fitted yet.
        """
        if self.autoencoder is not None:
            self.autoencoder.summary()
        else:
            raise(RuntimeError("The autoencoder has not been fitted yet."))