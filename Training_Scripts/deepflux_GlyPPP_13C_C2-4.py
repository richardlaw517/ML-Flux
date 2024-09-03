# Fix random seed number for reproducibility
from numpy.random import seed
seed(0)
import tensorflow as tf
tf.random.set_seed(0)

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
    model.add(keras.layers.Dense(4851, activation = 'relu', name="Dense_1"))
    model.add(keras.layers.Dense(70, activation = 'relu', name="Dense_2"))
    model.add(keras.layers.Dense(70, activation = 'relu', name="Dense_3"))
    model.add(keras.layers.Dense(70, activation = 'relu', name="Dense_4"))
    model.add(keras.layers.Dense(693, activation = 'relu', name="Dense_5"))

    model.add(keras.layers.Dense(output_shape, name="output"))

    print(model.summary())
    return model

label = np.loadtxt("GlyPPP_13C/labeling_1of10_GlyPPP.dat")
label = np.concatenate((label,np.loadtxt("labeling_2of10_GlyPPP.dat")))
label = np.concatenate((label,np.loadtxt("labeling_3of10_GlyPPP.dat")))
label = np.concatenate((label,np.loadtxt("labeling_4of10_GlyPPP.dat")))
label = np.concatenate((label,np.loadtxt("labeling_5of10_GlyPPP.dat")))
label = np.concatenate((label,np.loadtxt("labeling_6of10_GlyPPP.dat")))
label = np.concatenate((label,np.loadtxt("labeling_7of10_GlyPPP.dat")))
label = np.concatenate((label,np.loadtxt("labeling_8of10_GlyPPP.dat")))
label = np.concatenate((label,np.loadtxt("labeling_9of10_GlyPPP.dat")))
label = np.concatenate((label,np.loadtxt("labeling_10of10_GlyPPP.dat")))
flux = np.loadtxt("fluxes_GlyPPP_1M.dat")

net_flux = [0,1,2,3,4,5,6,7,8,9,10,11,12]
exchange_flux = [13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]

X_train, X_test, y_train, y_test = train_test_split(label, flux, test_size=0.2, random_state=0)
X_test, X_val, y_test, y_val = train_test_split(X_test, y_test, test_size=0.5, random_state=0)

#np.savetxt("flux_test_GlyPPP_1M.dat",y_test,fmt="%.6f")
#np.savetxt("label_test_GlyPPP_1M.dat",X_test,fmt="%.6f")

## Data transformation configuration 4
y_train[:,net_flux] = np.piecewise(y_train[:,net_flux],[y_train[:,net_flux]<3.89048,y_train[:,net_flux]>=3.89048],[lambda y_train: 1/(1+np.exp(-y_train)),lambda y_train: np.log10(y_train)/np.log10(4)])
y_train[:,exchange_flux] = np.piecewise(y_train[:,exchange_flux],[y_train[:,exchange_flux]<0,(0<=y_train[:,exchange_flux])&(y_train[:,exchange_flux]<1E-4),1E-4<=y_train[:,exchange_flux]],[-5,lambda y_train: y_train*1E4-5,lambda y_train: np.log10(y_train)])
y_val[:,net_flux] = np.piecewise(y_val[:,net_flux],[y_val[:,net_flux]<3.89048,y_val[:,net_flux]>=3.89048],[lambda y_val: 1/(1+np.exp(-y_val)),lambda y_val: np.log10(y_val)/np.log10(4)])
y_val[:,exchange_flux] = np.piecewise(y_val[:,exchange_flux],[y_val[:,exchange_flux]<0,(0<=y_val[:,exchange_flux])&(y_val[:,exchange_flux]<1E-4),1E-4<=y_val[:,exchange_flux]],[-5,lambda y_val: y_val*1E4-5,lambda y_val: np.log10(y_val)])
y_test[:,net_flux] = np.piecewise(y_test[:,net_flux],[y_test[:,net_flux]<3.89048,y_test[:,net_flux]>=3.89048],[lambda y_test: 1/(1+np.exp(-y_test)),lambda y_test: np.log10(y_test)/np.log10(4)])
y_test[:,exchange_flux] = np.piecewise(y_test[:,exchange_flux],[y_test[:,exchange_flux]<0,(0<=y_test[:,exchange_flux])&(y_test[:,exchange_flux]<1E-4),1E-4<=y_test[:,exchange_flux]],[-5,lambda y_test: y_test*1E4-5,lambda y_test: np.log10(y_test)])


ANN_regression = create_ANN(3456, 33) #num of columns in single matrix and number of free fluxes

# Compile the ANN model
ANN_regression.compile(
    optimizer='adam', loss='mae', metrics=['mae']
)

# Create callback
filepath = 'GlyPPP_1M_Seed1.h5'
checkpoint = ModelCheckpoint(filepath=filepath,
                             monitor='val_loss',
                             verbose=1,
                             save_best_only=True,
                             save_weights_only=True,
                             mode='min')
callbacks = [checkpoint]

# Evaluate the model
loss, acc = ANN_regression.evaluate(X_test, y_test, verbose=2)
print("Untrained model, accuracy: {:5.2f}%".format(100 * acc))

ANN_regression.load_weights('GlyPPP_1M_Seed1.h5')
    
loss, acc = ANN_regression.evaluate(X_test, y_test, verbose=2)
print("Restored model, accuracy: {:5.2f}%".format(100 * acc))

# Fit the ANN to the training set
history = ANN_regression.fit(
  x=X_train, y=y_train, shuffle=True, validation_data=(X_val, y_val),
  epochs=500, callbacks=callbacks, batch_size=32, verbose=2
)
