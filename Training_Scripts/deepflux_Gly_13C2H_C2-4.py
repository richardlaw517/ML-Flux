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
    model.add(keras.layers.Dense(432, activation = 'relu', name="Dense_1"))
    model.add(keras.layers.Dense(59, activation = 'relu', name="Dense_2"))
    model.add(keras.layers.Dense(59, activation = 'relu', name="Dense_3"))
    model.add(keras.layers.Dense(59, activation = 'relu', name="Dense_4"))
    model.add(keras.layers.Dense(432, activation = 'relu', name="Dense_5"))

    model.add(keras.layers.Dense(output_shape, name="output"))

    print(model.summary())
    return model

label = np.loadtxt("Glycolysis/Gly_Labels.dat")
flux = np.loadtxt("Glycolysis/Gly_Fluxes.dat")
net_flux = [0,1,2,3,4,5,6,7]
exchange_flux = [8,9,10,11,12,13,14,15,16,17]

X_train, X_test, y_train, y_test = train_test_split(label, flux, test_size=0.2, random_state=0)
X_test, X_val, y_test, y_val = train_test_split(X_test, y_test, test_size=0.5, random_state=0)

np.savetxt("Glycolysis/flux_test_Gly.dat",y_test,fmt="%.3f")
np.savetxt("Glycolysis/label_test_Gly.dat",X_test,fmt="%.3f")
np.savetxt("Glycolysis/flux_train_Gly.dat",y_train,fmt="%.3f")
np.savetxt("Glycolysis/label_train_Gly.dat",X_train,fmt="%.3f")
np.savetxt("Glycolysis/flux_val_Gly.dat",y_val,fmt="%.3f")
np.savetxt("Glycolysis/label_val_Gly.dat",X_val,fmt="%.3f")

## Data transformation configuration 4
y_train[:,net_flux] = np.piecewise(y_train[:,net_flux],[y_train[:,net_flux]<3.89048,y_train[:,net_flux]>=3.89048],[lambda y_train: 1/(1+np.exp(-y_train)),lambda y_train: np.log10(y_train)/np.log10(4)])
y_train[:,exchange_flux] = np.piecewise(y_train[:,exchange_flux],[y_train[:,exchange_flux]<0,(0<=y_train[:,exchange_flux])&(y_train[:,exchange_flux]<1E-4),1E-4<=y_train[:,exchange_flux]],[-5,lambda y_train: y_train*1E4-5,lambda y_train: np.log10(y_train)])
y_val[:,net_flux] = np.piecewise(y_val[:,net_flux],[y_val[:,net_flux]<3.89048,y_val[:,net_flux]>=3.89048],[lambda y_val: 1/(1+np.exp(-y_val)),lambda y_val: np.log10(y_val)/np.log10(4)])
y_val[:,exchange_flux] = np.piecewise(y_val[:,exchange_flux],[y_val[:,exchange_flux]<0,(0<=y_val[:,exchange_flux])&(y_val[:,exchange_flux]<1E-4),1E-4<=y_val[:,exchange_flux]],[-5,lambda y_val: y_val*1E4-5,lambda y_val: np.log10(y_val)])
y_test[:,net_flux] = np.piecewise(y_test[:,net_flux],[y_test[:,net_flux]<3.89048,y_test[:,net_flux]>=3.89048],[lambda y_test: 1/(1+np.exp(-y_test)),lambda y_test: np.log10(y_test)/np.log10(4)])
y_test[:,exchange_flux] = np.piecewise(y_test[:,exchange_flux],[y_test[:,exchange_flux]<0,(0<=y_test[:,exchange_flux])&(y_test[:,exchange_flux]<1E-4),1E-4<=y_test[:,exchange_flux]],[-5,lambda y_test: y_test*1E4-5,lambda y_test: np.log10(y_test)])


ANN_regression = create_ANN(192, 18) #num of columns in single matrix and number of free fluxes

# Compile the ANN model
ANN_regression.compile(
    optimizer='adam', loss='mae', metrics=None
)

# Create callback
filepath = 'Gly_13C2H_Seed0.h5'
checkpoint = ModelCheckpoint(filepath=filepath,
                             monitor='val_loss',
                             verbose=1,
                             save_best_only=True,
                             save_weights_only=True,
                             mode='min')
callbacks = [checkpoint]

# Save the ANN architecture to json
ANN_regression_json = ANN_regression.to_json()
with open("Gly_13C2H_Seed0.json", "w") as json_file:
    json_file.write(ANN_regression_json)
    
# Fit the ANN to the training set
history = ANN_regression.fit(
  x=X_train, y=y_train, shuffle=True, validation_data=(X_val, y_val),
  epochs=2000, callbacks=callbacks, batch_size=32, verbose=2
)

# Load the ANN architecture from json
pred_model = model_from_json(ANN_regression_json)

# Load weights from the best model into ANN model
pred_model.load_weights("Gly_13C2H_Seed0.h5")

# Compile the loaded ANN model
pred_model.compile(optimizer='adam', metrics=['mae'])

## Evaluate the performace
score = pred_model.evaluate(X_test, y_test, verbose=0)
print("%s: %.5f" % (pred_model.metrics_names[1], score[1]))

# Predict the test data
y_pred = pred_model.predict(X_test)

## Data inverse-transformation Configuration 4
y_pred[:,net_flux] = np.piecewise(y_pred[:,net_flux],[y_pred[:,net_flux]<0.97997,y_pred[:,net_flux]>=0.97997],[lambda y_pred: np.log(y_pred/(1-y_pred)),lambda y_pred: 4**y_pred])
y_pred[:,exchange_flux] = np.piecewise(y_pred[:,exchange_flux],[y_pred[:,exchange_flux]>-4,(-5<y_pred[:,exchange_flux])&(y_pred[:,exchange_flux]<=-4),y_pred[:,exchange_flux]<=-5],[lambda y_pred: 10**y_pred,lambda y_pred: (y_pred+5)/10000,0])
y_test[:,net_flux] = np.piecewise(y_test[:,net_flux],[y_test[:,net_flux]<0.97997,y_test[:,net_flux]>=0.97997],[lambda y_test: np.log(y_test/(1-y_test)),lambda y_test: 4**y_test])
y_test[:,exchange_flux] = np.piecewise(y_test[:,exchange_flux],[y_test[:,exchange_flux]>-4,(-5<y_test[:,exchange_flux])&(y_test[:,exchange_flux]<=-4),y_test[:,exchange_flux]<=-5],[lambda y_test: 10**y_test,lambda y_test: (y_test+5)/10000,0])

y_pred[:,exchange_flux] = np.maximum(0, y_pred[:,exchange_flux]) # Ensure exchange fluxes >= 0
print("mae: ",np.mean(np.abs(y_test - y_pred)))

np.savetxt("Gly_13C2H/flux_pred_Gly_13C2H.dat", y_pred, fmt="%.6f")
