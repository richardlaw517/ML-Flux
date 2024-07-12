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
    model.add(keras.layers.Dense(9020, activation = 'relu', name="Dense_1"))
    model.add(keras.layers.Dense(95, activation = 'relu', name="Dense_2"))
    model.add(keras.layers.Dense(95, activation = 'relu', name="Dense_3"))
    model.add(keras.layers.Dense(95, activation = 'relu', name="Dense_4"))
    model.add(keras.layers.Dense(1705, activation = 'relu', name="Dense_5"))

    model.add(keras.layers.Dense(output_shape, name="output"))

    print(model.summary())
    return model

label = np.loadtxt("MammalianCCM_13C/labeling_1of10_CCM_1M.dat")
label = np.concatenate((label,np.loadtxt("MammalianCCM_13C/labeling_2of10_CCM_1M.dat")))
label = np.concatenate((label,np.loadtxt("MammalianCCM_13C/labeling_3of10_CCM_1M.dat")))
label = np.concatenate((label,np.loadtxt("MammalianCCM_13C/labeling_4of10_CCM_1M.dat")))
label = np.concatenate((label,np.loadtxt("MammalianCCM_13C/labeling_5of10_CCM_1M.dat")))
label = np.concatenate((label,np.loadtxt("MammalianCCM_13C/labeling_6of10_CCM_1M.dat")))
label = np.concatenate((label,np.loadtxt("MammalianCCM_13C/labeling_7of10_CCM_1M.dat")))
label = np.concatenate((label,np.loadtxt("MammalianCCM_13C/labeling_8of10_CCM_1M.dat")))
label = np.concatenate((label,np.loadtxt("MammalianCCM_13C/labeling_9of10_CCM_1M.dat")))
label = np.concatenate((label,np.loadtxt("MammalianCCM_13C/labeling_10of10_CCM_1M.dat")))
flux = np.loadtxt("MammalianCCM_13C/fluxes_CCM_1M.dat")

net_flux = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
exchange_flux = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54]

X_train, X_test, y_train, y_test = train_test_split(label, flux, test_size=0.2, random_state=0)
X_test, X_val, y_test, y_val = train_test_split(X_test, y_test, test_size=0.5, random_state=0)

#np.savetxt("flux_test_CCM_1M.dat",y_test,fmt="%.6f")
#np.savetxt("label_test_CCM_1M.dat",X_test,fmt="%.6f")
#np.savetxt("label_train_CCM_1M.dat",X_train,fmt="%.6f")

## Data transformation configuration 4
y_train[:,net_flux] = np.piecewise(y_train[:,net_flux],[y_train[:,net_flux]<3.89048,y_train[:,net_flux]>=3.89048],[lambda y_train: 1/(1+np.exp(-y_train)),lambda y_train: np.log10(y_train)/np.log10(4)])
y_train[:,exchange_flux] = np.piecewise(y_train[:,exchange_flux],[y_train[:,exchange_flux]<0,(0<=y_train[:,exchange_flux])&(y_train[:,exchange_flux]<1E-4),1E-4<=y_train[:,exchange_flux]],[-5,lambda y_train: y_train*1E4-5,lambda y_train: np.log10(y_train)])
y_val[:,net_flux] = np.piecewise(y_val[:,net_flux],[y_val[:,net_flux]<3.89048,y_val[:,net_flux]>=3.89048],[lambda y_val: 1/(1+np.exp(-y_val)),lambda y_val: np.log10(y_val)/np.log10(4)])
y_val[:,exchange_flux] = np.piecewise(y_val[:,exchange_flux],[y_val[:,exchange_flux]<0,(0<=y_val[:,exchange_flux])&(y_val[:,exchange_flux]<1E-4),1E-4<=y_val[:,exchange_flux]],[-5,lambda y_val: y_val*1E4-5,lambda y_val: np.log10(y_val)])
y_test[:,net_flux] = np.piecewise(y_test[:,net_flux],[y_test[:,net_flux]<3.89048,y_test[:,net_flux]>=3.89048],[lambda y_test: 1/(1+np.exp(-y_test)),lambda y_test: np.log10(y_test)/np.log10(4)])
y_test[:,exchange_flux] = np.piecewise(y_test[:,exchange_flux],[y_test[:,exchange_flux]<0,(0<=y_test[:,exchange_flux])&(y_test[:,exchange_flux]<1E-4),1E-4<=y_test[:,exchange_flux]],[-5,lambda y_test: y_test*1E4-5,lambda y_test: np.log10(y_test)])

ANN_regression = create_ANN(4800, 55)

# Compile the ANN model
ANN_regression.compile(
    optimizer='adam', loss='mae', metrics=['mae']
)

# Create callback
filepath = 'CCM_1M_Seed3.h5'
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

ANN_regression.load_weights('CCM_1M_Seed3.h5')
    
loss, acc = ANN_regression.evaluate(X_test, y_test, verbose=2)
print("Restored model, accuracy: {:5.2f}%".format(100 * acc))

# Fit the ANN to the training set
history = ANN_regression.fit(
  x=X_train, y=y_train, shuffle=True, validation_data=(X_val, y_val),
  epochs=500, callbacks=callbacks, batch_size=32, verbose=2
)
