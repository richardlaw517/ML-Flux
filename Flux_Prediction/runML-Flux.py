from imputeLabelsFromScratched import imputeLabelsFromScratched
from inpaintLabelsFromScratched import inpaintLabelsFromScratched
from predictFluxesFromLabels import predictFluxesFromLabels
from inpaintLabelsFromScratched_mCM_ANN import inpaintLabelsFromScratched_mCM_ANN

import numpy as np
from keras.models import model_from_json

model_in = 'mCM' #CHANGE THIS
label_in = np.loadtxt('') #CHANGE THIS
# C2_yield = np.loadtxt('') #CHANGE THIS IF USING C2 YIELD MODEL FOR mCM, ALSO CHANGE MODEL IN predictFluxesFromLabels



# KNN
#fullLabelSet = imputeLabelsFromScratched(scratchedLabelSet=label_in,modelType=model_in)

# Inpainting
if model_in == 'mCM': 
    fullLabelSet = inpaintLabelsFromScratched_mCM_ANN(
    scratchedLabelSet=label_in,
    modelType=model_in
)
else: 

    fullLabelSet = inpaintLabelsFromScratched(scratchedLabelSet=label_in,modelType=model_in)

# fullLabelSet = np.concatenate(fullLabelSet,C2_yield),axis=1)

freeFluxes, freeList, fullFluxes, fullList = predictFluxesFromLabels(fullLabelSet=fullLabelSet,modelType=model_in)

