## Runs a script to predict fluxes from a masked label set
# 0. Upload the relevent files and set options
# 1. Predicts fluxes using the full labeling set

from predictFluxesFromLabels import predictFluxesFromLabels
import time
import numpy as np
from keras.models import model_from_json


# 0. Upload files with labeling data and fluxes
fullLabelSet = np.loadtxt("")
modelType = "Gly13C2H" # (Simple, UpperGly, Gly13C2H, GlyPPP, or MammalianCCM)

# 1. Predicts fluxes using the predicted full labeling set
startTime = time.time()
freeFluxes, freeList, fullFluxes, fullList = predictFluxesFromLabels(fullLabelSet=fullLabelSet,modelType=modelType)
endTime = time.time()
print("The time of execution of flux prediction is :", (endTime-startTime) * 10**3, "ms")
np.savetxt(str(modelType)+"_FullPredictedFluxes_"+time.strftime("%Y%m%d-%H%M%S")+".txt", fullFluxes)