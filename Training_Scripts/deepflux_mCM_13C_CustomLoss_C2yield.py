# %%

import tensorflow as tf
tf.random.set_seed(0)

import gc
import numpy as np
import keras
from keras.models import model_from_json
from keras.callbacks import ModelCheckpoint
from keras.layers import BatchNormalization
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
# %%


def create_ANN(input_shape, output_shape):
    model = keras.Sequential(name="model_ANN")
                                                                                                                                                 
    model.add(keras.layers.Input(shape=input_shape))

    model.add(keras.layers.Flatten())
  
    ## Nodes configuration 
    model.add(keras.layers.Dense(884, activation = 'relu', name="Dense_1")) 
    model.add(keras.layers.Dense(512, activation = 'relu', name="Dense_2")) 
    model.add(keras.layers.Dense(214, activation = 'relu', name="Dense_3")) 

    model.add(keras.layers.Dense(output_shape, name="output"))

    print(model.summary())
    return model


# %%
label = np.loadtxt( r"fulllabeling_11162025_clean_noNegSerBiosynth_v2_noErr.dat", delimiter=',', dtype=np.float32, skiprows=0)
flux  = np.loadtxt(r"FreeFluxes_11162025_clean_noNegSerBiosynth_v2_noErr.dat", delimiter=None,skiprows=0, dtype=np.float32)

max_flux = np.max(flux, axis = 0) # max flux along each column
#max_flux = np.loadtxt("max_free_flux.dat",skiprows=0, dtype=np.float32)
flux = flux/max_flux

np.savetxt("max_free_flux.dat",max_flux,fmt="%.6f")


net_flux = list(range(39)) 
exchange_flux = list(range(39,114))


X_train, X_test, y_train, y_test = train_test_split(label, flux, test_size=0.2, random_state=0)
X_test, X_val, y_test, y_val = train_test_split(X_test, y_test, test_size=0.5, random_state=0)

C2_yield_train = y_train[:,30].reshape(-1,1)
C2_yield_test = y_test[:,30].reshape(-1,1)
C2_yield_val = y_val[:,30].reshape(-1,1)

X_train = np.concatenate((X_train,C2_yield_train),axis=1)
X_test = np.concatenate((X_test,C2_yield_test),axis=1)
X_val = np.concatenate((X_val,C2_yield_val),axis=1)


# %%
ANN_regression = create_ANN(10609, 114)

kernelNet = np.loadtxt("kernel_net.txt").astype(np.float32).T
kernelXch = np.loadtxt("kernel_xch.txt").astype(np.float32).T

kernelNet = tf.constant(kernelNet)
kernelXch = tf.constant(kernelXch)

net_flux = tf.constant(list(range(39)))
exchange_flux = tf.constant(list(range(39, 114)))


def custom_loss_free(y_true, y_pred):
    squaredError = tf.square(y_true - y_pred)
    weights = np.ones(114, dtype=np.float32) #equal weights to start
    weights[2] = 15   # cs, g6pdh, ppc, PGA_serd
    weights[32] = 5
    weights[35] = 5
    weights[38] = 5
    weights[36] = 5 # eda
    weights[34] = 15 # oaadc
    weights[33] = 15 # me
    weights[30] = 5 #EX_AC
    weights[42] = 5 #xch_FBA
    weights[43] = 5 #xch_TPI
    weights[56] = 5 #xch_ICDH
    weights[61] = 5 #xch_MDH
    weights[67] = 5 # xch_PGA_Gly
    weights[0] = 0.05 # EX_CO2



    print(weights)
    weightedSquaredError = weights*squaredError
    
    loss = tf.reduce_mean(weightedSquaredError)
    
    return loss


def custom_loss(y_true, y_pred):
    
    fullNet_pred = tf.matmul(tf.gather(y_pred,net_flux,axis=1),kernelNet)
    fullXch_pred = tf.matmul(tf.gather(y_pred,exchange_flux,axis=1),kernelXch)
    fullFluxes_pred = tf.concat((fullNet_pred,fullXch_pred),axis=1)
    
    fullNet_true = tf.matmul(tf.gather(y_true,net_flux,axis=1),kernelNet)
    fullXch_true = tf.matmul(tf.gather(y_true,exchange_flux,axis=1),kernelXch)
    fullFluxes_true = tf.concat((fullNet_true,fullXch_true),axis=1)
    
    squaredError = tf.square(fullFluxes_true - fullFluxes_pred)
    weightedSquaredError = weights*squaredError
    
    loss = tf.reduce_mean(weightedSquaredError)
    
    return loss

#adjust optimizer learning rate
optimizer = keras.optimizers.Adam(learning_rate=1e-3)

ANN_regression.compile(
    optimizer=optimizer, loss=custom_loss_free, metrics=['mse'])

# Create callback
filepath = 'mCM_ANN_C2yield.h5'
checkpoint = ModelCheckpoint(filepath=filepath,
                             monitor='val_loss',
                             verbose=1,
                             save_best_only=True,
                             save_weights_only=True,
                             mode='min')
callbacks = [checkpoint]

# Save the ANN architecture to json
ANN_regression_json = ANN_regression.to_json()
with open("mCM_ANN_C2yield.json", "w") as json_file: 
    json_file.write(ANN_regression_json)

# Evaluate the model pre-training
loss, acc = ANN_regression.evaluate(X_test, y_test, verbose=2)
print("Untrained model, accuracy: {:5.2f}%".format(100 * acc))



# Fit the ANN to the training set
history = ANN_regression.fit(
  x=X_train, y=y_train, validation_data=(X_val, y_val),
  #x=X_train, y=y_train,
  epochs=50, callbacks=callbacks, batch_size=32, verbose=2  
) 



# Evaluate the model post-training
loss, acc = ANN_regression.evaluate(X_test, y_test, verbose=2)
print("Restored model, accuracy: {:5.2f}%".format(100 * acc))



