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
    model.add(keras.layers.Dense(216, activation = 'relu', name="Dense_1"))
    model.add(keras.layers.Dense(15, activation = 'relu', name="Dense_2"))
    model.add(keras.layers.Dense(15, activation = 'relu', name="Dense_3"))
    model.add(keras.layers.Dense(15, activation = 'relu', name="Dense_4"))
    model.add(keras.layers.Dense(36, activation = 'relu', name="Dense_5"))

    model.add(keras.layers.Dense(output_shape, name="output"))

    print(model.summary())
    return model

label = np.loadtxt("Simple/labeling_Simple.dat")
flux = np.loadtxt("Simple/fluxes_Simple.dat")
net_flux = [0,1]
exchange_flux = [2,3,4,5]

X_train, X_test, y_train, y_test = train_test_split(label, flux, test_size=0.2, random_state=0)
X_test, X_val, y_test, y_val = train_test_split(X_test, y_test, test_size=0.5, random_state=0)

np.savetxt("Simple/flux_test_Simple.dat", y_test, fmt="%.4f")
np.savetxt("Simple/flux_train_Simple.dat", y_train, fmt="%.4f")
np.savetxt("Simple/flux_val_Simple.dat", y_val, fmt="%.4f")
np.savetxt("Simple/labeling_test_Simple.dat", X_test, fmt="%.4f")
np.savetxt("Simple/labeling_train_Simple.dat", X_train, fmt="%.4f")
np.savetxt("Simple/labeling_val_Simple.dat", X_val, fmt="%.4f")

## Data transformation configuration 4
y_train[:,net_flux] = np.piecewise(y_train[:,net_flux],[y_train[:,net_flux]<3.89048,y_train[:,net_flux]>=3.89048],[lambda y_train: 1/(1+np.exp(-y_train)),lambda y_train: np.log10(y_train)/np.log10(4)])
y_train[:,exchange_flux] = np.piecewise(y_train[:,exchange_flux],[y_train[:,exchange_flux]<0,(0<=y_train[:,exchange_flux])&(y_train[:,exchange_flux]<1E-4),1E-4<=y_train[:,exchange_flux]],[-5,lambda y_train: y_train*1E4-5,lambda y_train: np.log10(y_train)])
y_val[:,net_flux] = np.piecewise(y_val[:,net_flux],[y_val[:,net_flux]<3.89048,y_val[:,net_flux]>=3.89048],[lambda y_val: 1/(1+np.exp(-y_val)),lambda y_val: np.log10(y_val)/np.log10(4)])
y_val[:,exchange_flux] = np.piecewise(y_val[:,exchange_flux],[y_val[:,exchange_flux]<0,(0<=y_val[:,exchange_flux])&(y_val[:,exchange_flux]<1E-4),1E-4<=y_val[:,exchange_flux]],[-5,lambda y_val: y_val*1E4-5,lambda y_val: np.log10(y_val)])
y_test[:,net_flux] = np.piecewise(y_test[:,net_flux],[y_test[:,net_flux]<3.89048,y_test[:,net_flux]>=3.89048],[lambda y_test: 1/(1+np.exp(-y_test)),lambda y_test: np.log10(y_test)/np.log10(4)])
y_test[:,exchange_flux] = np.piecewise(y_test[:,exchange_flux],[y_test[:,exchange_flux]<0,(0<=y_test[:,exchange_flux])&(y_test[:,exchange_flux]<1E-4),1E-4<=y_test[:,exchange_flux]],[-5,lambda y_test: y_test*1E4-5,lambda y_test: np.log10(y_test)])

ANN_regression = create_ANN(120, 6) # Number of inputs (MID vector size), Number of outputs (Fluxes)

# Compile the ANN model
ANN_regression.compile(
    optimizer='adam', loss='mae', metrics=None
)

# Create callback
filepath = 'Simple_ANN.h5' # Change for each configuration
checkpoint = ModelCheckpoint(filepath=filepath,
                             monitor='val_loss',
                             verbose=1,
                             save_best_only=True,
                             save_weights_only=True,
                             mode='min')
callbacks = [checkpoint]

# Save the ANN architecture to json
ANN_regression_json = ANN_regression.to_json()
with open("Simple_ANN_C2-4.json", "w") as json_file: # Change for each configuration
    json_file.write(ANN_regression_json)
    
# Fit the ANN to the training set
history = ANN_regression.fit(
  x=X_train, y=y_train, shuffle=True, validation_data=(X_val, y_val),
  epochs=2000, callbacks=callbacks, batch_size=1024, verbose=2
)

# Load the ANN architecture from json
pred_model = model_from_json(ANN_regression_json)

# Load weights from the best model into ANN model
pred_model.load_weights("Simple_ANN.h5") # Change for each configuration

# Compile the loaded ANN model
pred_model.compile(optimizer='adam', metrics=['mae'])

## Evaluate the performace
#score = pred_model.evaluate(X_test, y_test, verbose=0)
#print("%s: %.5f" % (pred_model.metrics_names[1], score[1]))

# Predict the test data
y_pred = pred_model.predict(X_test)
print(y_pred)

## Data inverse-transformation Configuration 4
y_pred[:,net_flux] = np.piecewise(y_pred[:,net_flux],[y_pred[:,net_flux]<0.97997,y_pred[:,net_flux]>=0.97997],[lambda y_pred: np.log(y_pred/(1-y_pred)),lambda y_pred: 4**y_pred])
y_pred[:,exchange_flux] = np.piecewise(y_pred[:,exchange_flux],[y_pred[:,exchange_flux]>-4,(-5<y_pred[:,exchange_flux])&(y_pred[:,exchange_flux]<=-4),y_pred[:,exchange_flux]<=-5],[lambda y_pred: 10**y_pred,lambda y_pred: (y_pred+5)/10000,0])
y_test[:,net_flux] = np.piecewise(y_test[:,net_flux],[y_test[:,net_flux]<0.97997,y_test[:,net_flux]>=0.97997],[lambda y_test: np.log(y_test/(1-y_test)),lambda y_test: 4**y_test])
y_test[:,exchange_flux] = np.piecewise(y_test[:,exchange_flux],[y_test[:,exchange_flux]>-4,(-5<y_test[:,exchange_flux])&(y_test[:,exchange_flux]<=-4),y_test[:,exchange_flux]<=-5],[lambda y_test: 10**y_test,lambda y_test: (y_test+5)/10000,0])

y_pred[:,exchange_flux] = np.maximum(0, y_pred[:,exchange_flux]) # Ensure exchange fluxes >= 0
print("mae: ",np.mean(np.abs(y_test - y_pred)))

np.savetxt("Simple/flux_pred_Simple.dat", y_pred, fmt="%.6f")