from imputeLabelsFromScratched import imputeLabelsFromScratched
from inpaintLabelsFromScratched import inpaintLabelsFromScratched
from predictFluxesFromLabels import predictFluxesFromLabels

import numpy as np
from keras.models import model_from_json

model_in = 'Simple' #CHANGE THIS
label_in = np.loadtxt('') #CHANGE THIS


# KNN
#fullLabelSet = imputeLabelsFromScratched(scratchedLabelSet=label_in,modelType=model_in)

# Inpainting
fullLabelSet = inpaintLabelsFromScratched(scratchedLabelSet=label_in,modelType=model_in)

freeFluxes, freeList, fullFluxes, fullList = predictFluxesFromLabels(fullLabelSet=fullLabelSet,modelType=model_in)

