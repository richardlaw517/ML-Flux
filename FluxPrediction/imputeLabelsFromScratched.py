## function imputeLabelsFromScrated gives a full set of all labels for a chosen model given an incomplete labeling set for the model
## Inputs (in order):
    # scratchedLabelSet = vector of labeling information with missing values filled with NaN
    # modelType = string that classifies which model is being used to predict fluxes ('Simple','UpperGly','Gly13C2H','GlyPPP', or 'MammalianCCM')
## Outputs (in order):
    # fullLabelSet: vector of all labeling information that is fully filled out by the user or imputation
## Note: This is currently written to only accept a single scratched vector at a time (i.e. 1 experiment), unlike predictFluxesFromLabels

import numpy as np

def imputeLabelsFromScratched(scratchedLabelSet,modelType):
    # Prepare training labels based on which model is used
    if modelType == 'Simple':
        label_train = np.loadtxt('')
    elif modelType == 'UpperGly':
        label_train = np.loadtxt('')
    elif modelType == 'Gly13C2H':
        label_train = np.loadtxt('')
    elif modelType == 'GlyPPP':
        label_train = np.loadtxt('')
    elif modelType == 'MammalianCCM':
        label_train = np.loadtxt('')
    else: raise ValueError("Unexpected model name. Accepted models are 'Simple','UpperGly','Gly13C2H','GlyPPP', or 'MammalianCCM'")

    n_label = len(scratchedLabelSet)

    missing_list = np.argwhere(np.isnan(scratchedLabelSet))
    missing_list = missing_list.reshape(len(missing_list))
    known_list = list(range(n_label))
    [known_list.remove(x) for x in missing_list]

    X = label_train[:,known_list]
    y = label_train[:,missing_list]

    from sklearn.neighbors import KNeighborsRegressor
    neigh = KNeighborsRegressor(n_neighbors=5, weights='distance')
    neigh.fit(X, y)

    pred = neigh.predict(scratchedLabelSet[known_list].reshape(1, -1))
    neigh_dist, neigh_idx = neigh.kneighbors(scratchedLabelSet[known_list].reshape(1, -1))

    label_all = np.zeros(n_label)

    label_all[known_list] = scratchedLabelSet[known_list] 
    label_all[missing_list] = pred
    fullLabelSet = label_all.reshape(1,n_label)

    
    return fullLabelSet
