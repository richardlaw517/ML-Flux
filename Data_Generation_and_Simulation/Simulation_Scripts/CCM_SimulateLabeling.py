#-------------------------------------------------------------------------------
# Name:        mfapy CCM_SimulateLabeling.py
#              Code used to simulate mass isotopomer distributions adapted from code in MFApy
#              The mfapy python package is required to run this script.
#
# Citation:    https://doi.org/10.1016/j.mec.2021.e00177
#-------------------------------------------------------------------------------
import mfapy
import csv
import numpy as np
#from mfapyio import *

if __name__ == '__main__':

    # Function to read CSV file
    def read_csv(file_path):
        data = []
        with open(file_path, 'r', newline='') as file:
            reader = csv.reader(file)
            for row in reader:
                data.append(row)
        return data

    # Function to write CSV file
    def write_csv(file_path, data):
        with open(file_path, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(data)

    # Load metabolic model from txt file
    reactions, reversible, metabolites, target_fragments = mfapy.mfapyio.load_metabolic_model("CCM_model.txt", format = "text")
    # Construction of MetabolicModel instance
    model = mfapy.metabolicmodel.MetabolicModel(reactions, reversible, metabolites, target_fragments)
    # Configurations
    model.set_configuration(callbacklevel = 0)

    # Load all fluxes to simulate (will be replaced with a function to generate the fluxes)
    fluxdata = np.loadtxt('FluxDataforMFApy.dat')
    fluxdata.reshape(157,-1)
    
    # Define tracers
    tracers_glc = [[1,1,1,1,1,1],
                [1,0,0,0,0,0],
                [0,1,0,0,0,0],
                [1,1,0,0,0,0],
                [0,0,0,0,0,1],
                [1,0,0,0,0,1],
                [0,0,1,0,0,0],
                [0,0,0,0,1,0],
                [0,0,0,1,0,0],
                [0,0,0,0,1,1]]
    tracers_gln = [[0,0,0,0,0],
                [1,1,1,1,1]]
    tracers = []
    for row_glc in tracers_glc:
        for row_gln in tracers_gln:
            tracers.append([row_glc, row_gln])
    tracers[1] = [[0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1]]

    Data_all = []
    for fluxNum in range(0,fluxdata.shape[1]):
        print("Starting simulation "+str(fluxNum))
        # Define which set of fluxes to iterate over
        flux2sim = fluxdata[:,fluxNum]
        csv_data = read_csv("CCM_status_1.csv")
        for i,row in enumerate(csv_data[1:]):
            row[3] = flux2sim[i]
        write_csv('modelfluxinput.csv',csv_data)
        
        # Add fluxes into  model
        flux = model.load_states('modelfluxinput.csv', format = 'csv')
        # flux = model.load_states('CCM_status_NCB_ML.csv', format = 'csv')
        # Generation of instances of CarbonSource class from model
        cs = model.generate_carbon_source_templete()

        Data_1FluxSet = []
        # Simulate over all tracer configurations for 1 flux set
        for tracerConfig in range(len(tracers)):
            glc_tracer = {'#' + ''.join(str(i) for i in tracers[tracerConfig][0]): 1.0}
            gln_tracer = {'#' + ''.join(str(i) for i in tracers[tracerConfig][1]): 1.0}
            cs.set_each_isotopomer('GLCIN', glc_tracer, correction = "no")
            cs.set_each_isotopomer('GlnIN', gln_tracer, correction = "no")

            # Generation of MDV instances from metabolic state and carbon source. (Simultation)
            mdv = model.generate_mdv(flux, cs)

            # Append results to matrix
            fragments = list(mdv.mdv.keys())[:-1]  
            for fragment in fragments:
                for isotope_number in range(8):  # Iterate over numbers 0 to 7
                    if isotope_number in mdv.mdv[fragment]:
                        ratio = mdv.mdv[fragment][isotope_number]['ratio']
                    else:
                        ratio = 0
                    Data_1FluxSet.append(ratio)
        
        
        # Define constants
        numMets = 30
        numExp = 20
        lmid = 8

        arrayData = np.array(Data_1FluxSet)
        # Reshape and rearrange the column
        reshaped_column = arrayData.reshape(numExp, numMets,lmid)
        rearranged_column = reshaped_column.transpose(1, 0, 2).reshape(-1)

        result = []
        for i in range(0, len(rearranged_column), lmid):
            subarr = rearranged_column[i:i+lmid]
            rounded_subarr = np.around(subarr, decimals=3)
            sum_subarr = np.sum(rounded_subarr)
            if sum_subarr != 1:
                diff = 1 - sum_subarr
                rounded_subarr[0] += diff
                rounded_subarr[0] = round(rounded_subarr[0], 3)  # Round again after adjustment
            result.extend(rounded_subarr)
        
        Data_all.append(result)
    
    transData = list(zip(*Data_all))
    print("Start writing")
    with open('CCM_SimulatedFluxes.csv', 'w', newline='') as f:
        import csv
        writer = csv.writer(f)
        writer.writerows(transData)
    print("Done")