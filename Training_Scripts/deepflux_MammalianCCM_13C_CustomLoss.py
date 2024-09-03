import tensorflow as tf
tf.random.set_seed(0)

import gc
import numpy as np
import keras
from keras.models import model_from_json
from keras.callbacks import ModelCheckpoint
from keras.layers import BatchNormalization
from sklearn.model_selection import train_test_split

def create_ANN(input_shape, output_shape):
    model = keras.Sequential(name="model_ANN")
                                                                                                                                                 
    model.add(keras.layers.Input(shape=input_shape))

    model.add(keras.layers.Flatten())

    ## Nodes configuration 2
    model.add(keras.layers.Dense(1024, activation = 'relu', name="Dense_1"))
    model.add(keras.layers.Dense(512, activation = 'relu', name="Dense_2"))
    model.add(keras.layers.Dense(256, activation = 'relu', name="Dense_3"))
    model.add(keras.layers.Dense(128, activation = 'relu', name="Dense_4"))
    model.add(keras.layers.Dense(64, activation = 'relu', name="Dense_5"))

    model.add(keras.layers.Dense(output_shape, name="output"))

    print(model.summary())
    return model

label = np.transpose(np.loadtxt("labeling.dat", delimiter=',', skiprows=0))
flux = np.loadtxt("fluxes.dat", skiprows=0)

net_flux = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
exchange_flux = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54]

X_train, X_test, y_train, y_test = train_test_split(label, flux, test_size=0.2, random_state=0)
X_test, X_val, y_test, y_val = train_test_split(X_test, y_test, test_size=0.5, random_state=0)

ANN_regression = create_ANN(4800, 55)

def custom_loss(y_true, y_pred):
    # Load Kernels and get full fluxes
    kernelNet = tf.transpose(np.loadtxt("../TrainedModels/ANN_Labels_to_Fluxes/Kernels/KernelNet_CCM.txt"))
    kernelXch = tf.transpose(np.loadtxt("../TrainedModels/ANN_Labels_to_Fluxes/Kernels/KernelXch_CCM.txt"))
    kernelNet = tf.cast(kernelNet, dtype=tf.float32)
    kernelXch = tf.cast(kernelXch, dtype=tf.float32)

    net_flux = tf.constant([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
    exchange_flux = tf.constant([21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54])
    
    fullNet_pred = tf.matmul(tf.gather(y_pred,net_flux,axis=1),kernelNet)
    fullXch_pred = tf.matmul(tf.gather(y_pred,exchange_flux,axis=1),kernelXch)
    fullFluxes_pred = tf.concat((fullNet_pred,fullXch_pred),axis=1)
    
    fullNet_true = tf.matmul(tf.gather(y_true,net_flux,axis=1),kernelNet)
    fullXch_true = tf.matmul(tf.gather(y_true,exchange_flux,axis=1),kernelXch)
    fullFluxes_true = tf.concat((fullNet_true,fullXch_true),axis=1)
    
    squaredError = tf.square(fullFluxes_true - fullFluxes_pred)

    weights = np.loadtxt("CCM_Weights.txt")
    weightedSquaredError = weights*squaredError
    
    loss = tf.reduce_mean(weightedSquaredError)
    
    return loss

ANN_regression.compile(
    optimizer='adam', loss=custom_loss, metrics=['mse'])

# Create callback
filepath = 'CCM_ANN.h5'
checkpoint = ModelCheckpoint(filepath=filepath,
                             monitor='val_loss',
                             verbose=1,
                             save_best_only=True,
                             save_weights_only=True,
                             mode='min')
callbacks = [checkpoint]

# Save the ANN architecture to json
ANN_regression_json = ANN_regression.to_json()
with open("CCM_ANN.json", "w") as json_file: # Change for each configuration
    json_file.write(ANN_regression_json)

# Evaluate the model
loss, acc = ANN_regression.evaluate(X_test, y_test, verbose=2)
print("Untrained model, accuracy: {:5.2f}%".format(100 * acc))

loss, acc = ANN_regression.evaluate(X_test, y_test, verbose=2)
print("Restored model, accuracy: {:5.2f}%".format(100 * acc))

# Fit the ANN to the training set
history = ANN_regression.fit(
  x=X_train, y=y_train, shuffle=True, validation_data=(X_val, y_val),
  epochs=250, callbacks=callbacks, batch_size=32, verbose=2
)