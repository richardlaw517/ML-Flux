## function predictFluxes gives the predicted set of all fluxes for a chosen model given a full labeling set for the model
## Inputs (in order):
    # fullLabelSet = array of all labeling information that is fully filled out by the user or imputation
    # modelType = string that classifies which model is being used to predict fluxes ('Simple','UpperGly','Gly13C2H','GlyPPP', or 'MammalianCCM')
## Outputs (in order):
    # freeFluxes = array of the predicted free/independent fluxes, with net fluxes first followed by exchange fluxes
    # freeList = string vector listing the names of the free  fluxes in order of appearance of freeFluxes
    # fullFluxes = array of all predicted net fluxes , with net fluxes first followed by exchange fluxes
    # fullList = string vector listing the names of all fluxes in order of appearance of fullFluxes
## Note: This is currently written to be able to accept multiple label inputs at a time (i.e. multiple experiments)

import numpy as np
from keras.models import model_from_json

def predictFluxesFromLabels(fullLabelSet,modelType):
    # Prepare function calls based on which model is used
    if modelType == 'Simple':
        modelArchitectureFile = 'Trained_Models/ANN_Labels_to_Fluxes/Simple_ANN.json'
        modelWeightsFiles = 'Trained_Models/ANN_Labels_to_Fluxes/Simple_ANN.h5'
        net_flux = [0,1]
        exchange_flux = [2,3,4,5]
        freeList = ['Ex_E','Ex_F','v2','v3','v4','v5']
        fullList = ['v1','v2','v3','v4','v5','Ex_E','Ex_F','v1_x','v2_x','v3_x','v4_x','v_5_x','Ex_E_x','Ex_F_x']
        kernelNet = np.loadtxt("Trained_Models/ANN_Labels_to_Fluxes/Kernels/KernelNet_Simple.txt")
        kernelXch = np.loadtxt("Trained_Models/ANN_Labels_to_Fluxes/Kernels/KernelXch_Simple.txt")
    elif modelType == 'UpperGly':
        modelArchitectureFile = 'Trained_Models/ANN_Labels_to_Fluxes/UpperGly_ANN.json'
        modelWeightsFiles = 'Trained_Models/ANN_Labels_to_Fluxes/UpperGly_ANN.h5'
        net_flux = [0,1,2]
        exchange_flux = [3,4,5,6]
        freeList = ['G6P_EX','F6P_EX','DHAP_EX','PGI','PFK','FBA','TPI']
        fullList = ['G6P_EX','F6P_EX','DHAP_EX','PGI','PFK','FBA','TPI']
        kernelNet = np.loadtxt("Trained_Models/ANN_Labels_to_Fluxes/Kernels/KernelNet_UpperGly.txt")
        kernelXch = np.loadtxt("Trained_Models/ANN_Labels_to_Fluxes/Kernels/KernelXch_UpperGly.txt")
    elif modelType == 'Gly13C2H':
        modelArchitectureFile = 'Trained_Models/ANN_Labels_to_Fluxes/Glycolysis_ANN.json'
        modelWeightsFiles = 'Trained_Models/ANN_Labels_to_Fluxes/Glycolysis_ANN.h5'
        net_flux = [0,1,2,3,4,5,6,7]
        exchange_flux = [8,9,10,11,12,13,14,15,16,17]
        freeList = ['H_EX','NADH_EX','G6P_EX','F6P_EX','DHAP_EX','PG3_EX','PEP_EX','PYR_EX','PGI','PFK','FBA','TPI','GAPD','PGK','PGM1','PGM2','ENO','PYK']
        fullList = ['GLC_IN','H_IN','NADH_IN','PGI','PFK','FBA','TPI','GAPD','PGK','PGM1','PGM2','ENO','PYK','H_EX','NADH_EX','G6P_EX','F6P_EX','DHAP_EX','PG3_EX','PEP_EX',
                    'GLC_IN','H_IN','NADH_IN','PGI','PFK','FBA','TPI','GAPD','PGK','PGM1','PGM2','ENO','PYK','H_EX','NADH_EX','G6P_EX','F6P_EX','DHAP_EX','PG3_EX','PEP_EX']
        kernelNet = np.loadtxt("Trained_Models/ANN_Labels_to_Fluxes/Kernels/KernelNet_Glycolysis.txt")
        kernelXch = np.loadtxt("Trained_Models/ANN_Labels_to_Fluxes/Kernels/KernelXch_Glycolysis.txt")
    elif modelType == 'GlyPPP':
        modelArchitectureFile = 'Trained_Models/ANN_Labels_to_Fluxes/GlyPPP_ANN.json'
        modelWeightsFiles = 'Trained_Models/ANN_Labels_to_Fluxes/GlyPPP_ANN.h5'
        net_flux = [0,1,2,3,4,5,6,7,8,9,10,11,12]
        exchange_flux = [13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]
        freeList = ['TAL','SBPASE','CO2_EX','H_EX','NADPH_EX','NADH_EX','G6P_EX','F6P_EX','DHAP_EX','R5P_EX','E4P_EX','PG3_EX','PEP_EX'
                    'PGI','PFK','FBA','TPI','GAPD','PGK','PGM1','PGM2','ENO','G6PDH','GND','RPI','RPE','TKT2','TKT1','TAL','SBA','SBPASE','PGI_leak','RPI_leak']
        fullList = ['GLC_IN','CO2_IN','H_IN','NADPH_IN','NADH_IN','PGI','PFK','FBA','TPI','GAPD','PGK','PGM1','PGM2','ENO','G6PDH','GND','RPI','RPE','TKT2','TKT1','TAL','SBA','SBPASE','CO2_EX','H_EX','NADPH_EX','NADH_EX','G6P_EX','F6P_EX','DHAP_EX','R5P_EX',
                    'E4P_EX','PG3_EX','PEP_EX','PGI_leak','RPI_leak',
                    'GLC_IN','CO2_IN','H_IN','NADPH_IN','NADH_IN','PGI','PFK','FBA','TPI','GAPD','PGK','PGM1','PGM2','ENO','G6PDH','GND','RPI','RPE','TKT2','TKT1','TAL','SBA','SBPASE','CO2_EX','H_EX','NADPH_EX','NADH_EX','G6P_EX','F6P_EX','DHAP_EX','R5P_EX','E4P_EX','PG3_EX','PEP_EX','PGI_leak','RPI_leak']
        kernelNet = np.loadtxt("Trained_Models/ANN_Labels_to_Fluxes/Kernels/KernelNet_GlyPPP.txt")
        kernelXch = np.loadtxt("Trained_Models/ANN_Labels_to_Fluxes/Kernels/KernelXch_GlyPPP.txt")
    elif modelType == 'MammalianCCM':
        modelArchitectureFile = 'Trained_Models/ANN_Labels_to_Fluxes/CCM_ANN.json'
        modelWeightsFiles = 'Trained_Models/ANN_Labels_to_Fluxes/CCM_ANN.h5'
        net_flux = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
        exchange_flux = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54]
        freeList = ['EX_CO2','EX_DHAP','EX_PGA','EX_PYR','EX_LAC','EX_G6P','pc','tal','SBPase','EX_R5P','EX_OAA','IN_AC','fum','mdh','OGA_Glu','IN_Glu','OGA_Gln','EX_Gln','IN_OAA','EX_OGA_Glu','EX_AC_cyt',
                    'hk','pgi','pfk','fba','tpi','gapd','pgk','eno','pyk','ldh','ppck','me','pc','g6pdh','gnd','rpi','rpe','tkt2','tkt1','tal','SBA','SBPase','pdh','cs','acitl','icdh','akgdh','sucoas','sucd','fum','mdh','PYR_Ala','OGA_Glu','OGA_Gln']
        fullList = ['IN_GLC','IN_CO2','EX_CO2','hk','pgi','pfk','fba','tpi','gapd','pgk','eno','pyk','ldh','EX_DHAP','EX_PGA','EX_PYR','EX_LAC','EX_G6P','ppck','me','pc','g6pdh','gnd','rpi','rpe','tkt2','tkt1','tal','SBA','SBPase','EX_R5P','EX_OAA','pdh','IN_AC','cs','acitl','icdh','akgdh','sucoas','sucd','fum','mdh','PYR_Ala','EX_PYR_Ala','OGA_Glu','IN_Gln','IN_Glu','OGA_Gln','EX_Gln','IN_OAA','EX_OGA_Glu','EX_AC_cyt',
                    'IN_GLC','IN_CO2','EX_CO2','hk','pgi','pfk','fba','tpi','gapd','pgk','eno','pyk','ldh','EX_DHAP','EX_PGA','EX_PYR','EX_LAC','EX_G6P','ppck','me','pc','g6pdh','gnd','rpi','rpe','tkt2','tkt1','tal','SBA','SBPase','EX_R5P','EX_OAA','pdh','IN_AC','cs','acitl','icdh','akgdh','sucoas','sucd','fum','mdh','PYR_Ala','EX_PYR_Ala','OGA_Glu','IN_Gln','IN_Glu','OGA_Gln','EX_Gln','IN_OAA','EX_OGA_Glu','EX_AC_cyt']
        kernelNet = np.loadtxt("Trained_Models/ANN_Labels_to_Fluxes/Kernels/KernelNet_CCM.txt")
        kernelXch = np.loadtxt("Trained_Models/ANN_Labels_to_Fluxes/Kernels/KernelXch_CCM.txt")
    else: raise ValueError("Unexpected model name. Accepted models are 'Simple','UpperGly','Gly13C2H','GlyPPP', or 'MammalianCCM'")
    
    # Load model architecture
    json_file = open(modelArchitectureFile, 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    loadedModel = model_from_json(loaded_model_json)

    # Load weights into model
    loadedModel.load_weights(modelWeightsFiles) 

    # Compile model and optimizer
    loadedModel.compile(optimizer='adam', loss='mae', metrics=[None])

    # Predict fluxes and transform output into real flux values
    if modelType == 'MammalianCCM':
        freeFluxes = loadedModel.predict(fullLabelSet)
    else:
        freeFluxes = loadedModel.predict(fullLabelSet)
        freeFluxes[:,net_flux] = np.piecewise(freeFluxes[:,net_flux],[freeFluxes[:,net_flux]<0,freeFluxes[:,net_flux]<0.97997,freeFluxes[:,net_flux]>=0.97997],[-4.5,lambda freeFluxes: np.log(freeFluxes/(1-freeFluxes)),lambda freeFluxes: 4**freeFluxes])
        freeFluxes[:,exchange_flux] = np.piecewise(freeFluxes[:,exchange_flux],[freeFluxes[:,exchange_flux]>-4,(-5<freeFluxes[:,exchange_flux])&(freeFluxes[:,exchange_flux]<=-4),freeFluxes[:,exchange_flux]<=-5],[lambda freeFluxes: 10**freeFluxes,lambda freeFluxes: (freeFluxes+5)/10000,0])
    
    freeFluxes[:,exchange_flux] = np.maximum(0, freeFluxes[:,exchange_flux]) # Ensures all exchange fluxes are positive

    # Generate the full set of net and exchange fluxes
    fullNet = np.matmul(kernelNet,np.transpose(freeFluxes[:,net_flux]))
    fullXch = np.matmul(kernelXch,np.transpose(freeFluxes[:,exchange_flux]))
    fullFluxes = np.transpose(np.concatenate((fullNet,fullXch)))

    return freeFluxes, freeList, fullFluxes, fullList