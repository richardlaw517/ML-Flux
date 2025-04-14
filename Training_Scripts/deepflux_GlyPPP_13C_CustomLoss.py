# Fix random seed number for reproducibility
from numpy.random import seed
seed(0)
import tensorflow as tf
tf.random.set_seed(1)

import numpy as np
import keras
from keras.models import model_from_json
from keras.callbacks import ModelCheckpoint
from keras.layers import BatchNormalization
from sklearn.model_selection import train_test_split
import gc

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

label = np.loadtxt("../../MFEA/training/GlyPPP_13C/20250219/GlyPPP_labeling_Seed0_Constraints1of3_20pct_20250219.dat")
label = np.concatenate((label,np.loadtxt("../../MFEA/training/GlyPPP_13C/20250219/GlyPPP_labeling_Seed1_Constraints2of3_30pct_20250219.dat")))
label = np.concatenate((label,np.loadtxt("../../MFEA/training/GlyPPP_13C/20250219/GlyPPP_labeling_Seed2_Constraints3of3_50pct_20250219.dat")))
flux = np.loadtxt("../../MFEA/training/GlyPPP_13C/20250219/GlyPPP_fluxes_Seed0_Constraints1of3_20pct_20250219.dat")
flux = np.concatenate((flux,np.loadtxt("../../MFEA/training/GlyPPP_13C/20250219/GlyPPP_fluxes_Seed1_Constraints2of3_30pct_20250219.dat")))
flux = np.concatenate((flux,np.loadtxt("../../MFEA/training/GlyPPP_13C/20250219/GlyPPP_fluxes_Seed2_Constraints3of3_50pct_20250219.dat")))

net_flux = [0,1,2,3,4,5,6,7,8,9,10,11,12]
exchange_flux = [13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]

# First, separate each subset
subset1_X = label[:200000]
subset1_y = flux[:200000]
subset2_X = label[200000:500000]
subset2_y = flux[200000:500000]
subset3_X = label[500000:]
subset3_y = flux[500000:]

# Delete unused variable
del label
gc.collect()

# Then, perform train-test split for each subset separately
subset1_X_train, subset1_X_test, subset1_y_train, subset1_y_test = train_test_split(subset1_X, subset1_y, test_size=0.2, random_state=0)
subset2_X_train, subset2_X_test, subset2_y_train, subset2_y_test = train_test_split(subset2_X, subset2_y, test_size=0.2, random_state=0)
subset3_X_train, subset3_X_test, subset3_y_train, subset3_y_test = train_test_split(subset3_X, subset3_y, test_size=0.2, random_state=0)
subset1_X_test, subset1_X_val, subset1_y_test, subset1_y_val = train_test_split(subset1_X_test, subset1_y_test, test_size=0.5, random_state=0)
subset2_X_test, subset2_X_val, subset2_y_test, subset2_y_val = train_test_split(subset2_X_test, subset2_y_test, test_size=0.5, random_state=0)
subset3_X_test, subset3_X_val, subset3_y_test, subset3_y_val = train_test_split(subset3_X_test, subset3_y_test, test_size=0.5, random_state=0)

del subset1_X
del subset1_y
del subset2_X
del subset2_y
del subset3_X
del subset3_y
gc.collect()

# Concatenate the subsets back together
X_train = np.concatenate((subset1_X_train, subset2_X_train, subset3_X_train), axis=0)
X_test = np.concatenate((subset1_X_test, subset2_X_test, subset3_X_test), axis=0)
X_val = np.concatenate((subset1_X_val, subset2_X_val, subset3_X_val), axis=0)
y_train = np.concatenate((subset1_y_train, subset2_y_train, subset3_y_train), axis=0)
y_test = np.concatenate((subset1_y_test, subset2_y_test, subset3_y_test), axis=0)
y_val = np.concatenate((subset1_y_val, subset2_y_val, subset3_y_val), axis=0)

del subset1_X_test
del subset1_X_train
del subset1_X_val
del subset2_X_test
del subset2_X_train
del subset2_X_val
del subset3_X_test
del subset3_X_train
del subset3_X_val
gc.collect


# X_train, X_test, y_train, y_test = train_test_split(label, flux, test_size=0.2, random_state=0)
# X_test, X_val, y_test, y_val = train_test_split(X_test, y_test, test_size=0.5, random_state=0)

# np.savetxt("flux_test_GlyPPP_1M_20250219.dat",y_test,fmt="%.6f")
# np.savetxt("label_test_GlyPPP_1M_20250219.dat",X_test,fmt="%.6f")
# np.savetxt("label_train_GlyPPP_1M_20250219.dat",X_train,fmt="%.6f")
# killer = 190/0

## Data transformation configuration 4
# y_train[:,net_flux] = np.piecewise(y_train[:,net_flux],[y_train[:,net_flux]<3.89048,y_train[:,net_flux]>=3.89048],[lambda y_train: 1/(1+np.exp(-y_train)),lambda y_train: np.log10(y_train)/np.log10(4)])
# y_train[:,exchange_flux] = np.piecewise(y_train[:,exchange_flux],[y_train[:,exchange_flux]<0,(0<=y_train[:,exchange_flux])&(y_train[:,exchange_flux]<1E-4),1E-4<=y_train[:,exchange_flux]],[-5,lambda y_train: y_train*1E4-5,lambda y_train: np.log10(y_train)])
# y_val[:,net_flux] = np.piecewise(y_val[:,net_flux],[y_val[:,net_flux]<3.89048,y_val[:,net_flux]>=3.89048],[lambda y_val: 1/(1+np.exp(-y_val)),lambda y_val: np.log10(y_val)/np.log10(4)])
# y_val[:,exchange_flux] = np.piecewise(y_val[:,exchange_flux],[y_val[:,exchange_flux]<0,(0<=y_val[:,exchange_flux])&(y_val[:,exchange_flux]<1E-4),1E-4<=y_val[:,exchange_flux]],[-5,lambda y_val: y_val*1E4-5,lambda y_val: np.log10(y_val)])
# y_test[:,net_flux] = np.piecewise(y_test[:,net_flux],[y_test[:,net_flux]<3.89048,y_test[:,net_flux]>=3.89048],[lambda y_test: 1/(1+np.exp(-y_test)),lambda y_test: np.log10(y_test)/np.log10(4)])
# y_test[:,exchange_flux] = np.piecewise(y_test[:,exchange_flux],[y_test[:,exchange_flux]<0,(0<=y_test[:,exchange_flux])&(y_test[:,exchange_flux]<1E-4),1E-4<=y_test[:,exchange_flux]],[-5,lambda y_test: y_test*1E4-5,lambda y_test: np.log10(y_test)])


ANN_regression = create_ANN(3456, 33) #num of columns in single matrix and number of free fluxes

def custom_loss(y_true, y_pred):
    # Load Kernels and get full fluxes
    kernelNet = tf.transpose(np.loadtxt("../model/KernelNet_GlyPPP.txt"))
    kernelXch = tf.transpose(np.loadtxt("../model/KernelXch_GlyPPP.txt"))
    kernelNet = tf.cast(kernelNet, dtype=tf.float32)
    kernelXch = tf.cast(kernelXch, dtype=tf.float32)

    net_flux = [0,1,2,3,4,5,6,7,8,9,10,11,12]
    exchange_flux = [13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]

    fullNet_pred = tf.matmul(tf.gather(y_pred,net_flux,axis=1),kernelNet)
    fullXch_pred = tf.matmul(tf.gather(y_pred,exchange_flux,axis=1),kernelXch)
    fullFluxes_pred = tf.concat((fullNet_pred,fullXch_pred),axis=1)
    
    fullNet_true = tf.matmul(tf.gather(y_true,net_flux,axis=1),kernelNet)
    fullXch_true = tf.matmul(tf.gather(y_true,exchange_flux,axis=1),kernelXch)
    fullFluxes_true = tf.concat((fullNet_true,fullXch_true),axis=1)
    
    # squaredError = tf.square(fullFluxes_true - fullFluxes_pred)

    # weights = 1/(abs(fullFluxes_true)+0.01)
    # weightedSquaredError = weights*squaredError
    
    # loss = tf.reduce_mean(weightedSquaredError)

    squaredError = tf.square(fullFluxes_true - fullFluxes_pred)
    
    weights = np.loadtxt("GlyPPP_13C/GlyPPP_Weights_20250219.txt")
    weightedSquaredError = weights*squaredError
    

    loss = tf.reduce_mean(weightedSquaredError)
    
    return loss

def custom_loss_dG(y_true, y_pred):
    # Load Kernels and get full fluxes
    kernelNet = tf.transpose(np.loadtxt("../model/KernelNet_GlyPPP.txt"))
    kernelXch = tf.transpose(np.loadtxt("../model/KernelXch_GlyPPP.txt"))
    kernelNet = tf.cast(kernelNet, dtype=tf.float32)
    kernelXch = tf.cast(kernelXch, dtype=tf.float32)

    net_flux = [0,1,2,3,4,5,6,7,8,9,10,11,12]
    exchange_flux = [13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]

    fullNet_pred = tf.matmul(tf.gather(y_pred,net_flux,axis=1),kernelNet)
    fullXch_pred = tf.matmul(tf.gather(y_pred,exchange_flux,axis=1),kernelXch)
    fullFluxes_pred = tf.concat((fullNet_pred,fullXch_pred),axis=1)
    
    fullNet_true = tf.matmul(tf.gather(y_true,net_flux,axis=1),kernelNet)
    fullXch_true = tf.matmul(tf.gather(y_true,exchange_flux,axis=1),kernelXch)
    fullFluxes_true = tf.concat((fullNet_true,fullXch_true),axis=1)

    squaredError = tf.square(fullFluxes_true - fullFluxes_pred)
    
    weights = np.loadtxt("GlyPPP_13C/GlyPPP_Weights_20250219.txt")
    weightedSquaredError = weights*squaredError

    r = 8.314462618
    t = 310

    pred_for_fluxes = tf.abs(fullXch_pred + tf.maximum(0.0, fullNet_pred))
    pred_rev_fluxes = tf.abs(fullXch_pred - tf.minimum(0.0, fullNet_pred))
    dG_Pred = (r * t * tf.math.log(pred_rev_fluxes / pred_for_fluxes)) / 1000.0

    true_for_fluxes = fullXch_true + tf.maximum(0.0, fullNet_true)
    true_rev_fluxes = fullXch_true - tf.minimum(0.0, fullNet_true)
    dG_True = (r * t * tf.math.log(true_rev_fluxes / true_for_fluxes)) / 1000.0

    squaredError_dG = tf.square(dG_True - dG_Pred)
    
    weights_dG = np.loadtxt("GlyPPP_13C/GlyPPP_Weights_20250219.txt")
    weightedSquaredError_dG = weights_dG #*squaredError_dG

    # loss = tf.reduce_mean(weightedSquaredError) + tf.reduce_mean(weightedSquaredError_dG)
    loss = tf.reduce_mean(tf.cast(weightedSquaredError, tf.float32)) + tf.reduce_mean(tf.cast(weightedSquaredError_dG, tf.float32))


    return loss

# Compile the ANN model
ANN_regression.compile(
    optimizer='adam', loss=custom_loss_dG, metrics=['mape']
)

# Create callback
filepath = '../../MFEA/GlyPPP_1M_20250219_CustomLoss_dG.h5'
checkpoint = ModelCheckpoint(filepath=filepath,
                             monitor='val_loss',
                             verbose=1,
                             save_best_only=True,
                             save_weights_only=True,
                             mode='min')
callbacks = [checkpoint]

# Save the ANN architecture to json
ANN_regression_json = ANN_regression.to_json()
with open("../model/GlyPPP_1M_20250219_CustomLoss_dG.json", "w") as json_file: # Change for each configuration
    json_file.write(ANN_regression_json)

# Evaluate the model
loss, acc = ANN_regression.evaluate(X_test, y_test, verbose=2)
print("Untrained model, accuracy: {:5.2f}%".format(100 * acc))

# ANN_regression.load_weights('../../MFEA/GlyPPP_1M_20240612.h5')
    
loss, acc = ANN_regression.evaluate(X_test, y_test, verbose=2)
print("Restored model, accuracy: {:5.2f}%".format(100 * acc))

# Fit the ANN to the training set
history = ANN_regression.fit(
  x=X_train, y=y_train, shuffle=True, validation_data=(X_val, y_val),
  epochs=250, callbacks=callbacks, batch_size=32, verbose=2
)
